#!/usr/bin/env python
# SPDX-License-Identifier: LGPL-3.0-or-later
# SPDX-FileCopyrightText: Copyright CERN for the benefit of the SHiP Collaboration

import json
import os

import geometry_config
import ROOT
import shipRoot_conf
import shipunit as u
from heavyFlavourScaling import (
    check_run_type_override,
    derive_cross_sections,
    format_summary,
)

mcEngine = "TGeant4"
simEngine = "Pythia8"
checkOverlap = True
outputDir = "."
dy = 6.0  # 10.
ds = 8  # 9 # 5=TP muon shield, 6=magnetized hadron, 7=short magnet design, 9=optimised with T4 as constraint, 8=requires config file

# example for primary interaction, nobias: python $FAIRSHIP/muonShieldOptimization/run_fixedTarget.py -n 10000 -e 10 -f -r 10
#                                                               10000 events, energy cut 10GeV, run nr 10, override existing output folder
# example for charm decays, python $FAIRSHIP/muonShieldOptimization/run_fixedTarget.py -C -M -n 10000 -e 10  -r 60 -b 50 -f
#                                                               10000 events, charm decays, energy cut 10GeV, run nr 60, override existing output folder
#                                                               increase di-muon BRs for resonances < 1.1GeV by a factor 50

# ----------------------------- Yandex production ------------------------------
import argparse
import logging
import shutil

logging.info("")
logger = logging.getLogger(os.path.splitext(os.path.basename(os.sys.argv[0]))[0])
logger.setLevel(logging.INFO)


def get_work_dir(run_number, tag: str | None = None) -> str:
    import socket

    host = socket.gethostname()
    job_base_name = os.path.splitext(os.path.basename(os.sys.argv[0]))[0]
    if tag:
        out_dir = f"{host}_{job_base_name}_{run_number}_{tag}"
    else:
        out_dir = f"{host}_{job_base_name}_{run_number}"
    return out_dir


logger.info("SHiP proton-on-taget simulator (C) Thomas Ruf, 2017")

ap = argparse.ArgumentParser(description='Run SHiP "pot" simulation')
ap.add_argument("-d", "--debug", action="store_true")
ap.add_argument("-f", "--force", action="store_true", help="force overwriting output directory")
ap.add_argument("-r", "--run-number", type=int, dest="runnr", default=1)
ap.add_argument(
    "--reproducible",
    action="store_true",
    help="Reduce nondeterministic log output for reproducibility/testing",
)
ap.add_argument(
    "-e", "--ecut", type=float, help="energy cut", default=0.5
)  # GeV   with 1 : ~1sec / event, with 2: 0.4sec / event, 10: 0.13sec
ap.add_argument("-n", "--num-events", type=int, help="number of events to generate", dest="nev", default=100)
ap.add_argument(
    "-G",
    "--G4only",
    action=argparse.BooleanOptionalAction,
    default=False,
    help="Whether or not to use Geant4 directly, no Pythia8 (--no-G4only or --G4only). Default set to False.",
)
ap.add_argument(
    "-P",
    "--pythiaDecay",
    action=argparse.BooleanOptionalAction,
    default=False,
    help="Whether or not to use Pythia8 for decays (--no-PythiaDecay or --PythiaDecay). Default set to False.",
)
ap.add_argument("-t", "--tau-only", action=argparse.BooleanOptionalAction, dest="tauOnly", default=False)
ap.add_argument("-J", "--Jpsi-mainly", action=argparse.BooleanOptionalAction, dest="JpsiMainly", default=False)
ap.add_argument("-b", "--boostDiMuon", type=float, default=1.0, help="boost Di-muon branching ratios")
ap.add_argument("-X", "--boostFactor", type=float, default=1.0, help="boost Di-muon prod cross sections")
ap.add_argument(
    "--kaon-pion-splits",
    type=int,
    default=0,
    help="splitting factor for kaons and pions, in order to boost the number of muons stemming from their decays",
)
ap.add_argument(
    "--multiple-kpi-splits", action="store_true", help="split kaons and pions multiple times along the track path"
)

ap.add_argument("-C", "--charm", action=argparse.BooleanOptionalAction, default=False, help="generate charm decays")
ap.add_argument("-B", "--beauty", action=argparse.BooleanOptionalAction, default=False, help="generate beauty decays")
ap.add_argument(
    "-M",
    "--storeOnlyMuons",
    action=argparse.BooleanOptionalAction,
    default=False,
    help="store only muons, ignore neutrinos",
)
ap.add_argument("-N", "--skipNeutrinos", action=argparse.BooleanOptionalAction, default=False, help="skip neutrinos")
ap.add_argument(
    "-D",
    "--4darkPhoton",
    action=argparse.BooleanOptionalAction,
    dest="FourDP",
    default=False,
    help="enable ntuple production",
)
# for charm production
# A run is either charm or beauty, so only one ratio override is ever meaningful.
cross_section = ap.add_mutually_exclusive_group()
cross_section.add_argument(
    "-cc", "--chicc", type=float, default=None, help="ccbar over mbias cross section (overrides target-derived value)"
)
cross_section.add_argument(
    "-bb", "--chibb", type=float, default=None, help="bbbar over mbias cross section (overrides target-derived value)"
)
ap.add_argument(
    "--target-composition",
    default="W",
    choices=["W", "Mo"],
    help="Target composition. Default is Tungsten (W); Molybdenum (Mo) is the other preset.",
)
ap.add_argument(
    "-A",
    type=float,
    default=None,
    help=(
        "Target mass number; overrides --target-composition preset. "
        "Used to scale chicc/chibb as (A/A_Mo)^(heavyflavour_Ascale-mbias_Ascale) "
        "(default exponent: 0.29)."
    ),
)
ap.add_argument("-p", "--pot", default=4e13, help="number of protons on target per spill to normalize on")
ap.add_argument("-S", "--nStart", type=int, help="first event of input file to start", dest="nStart", default=0)
DEFAULT_CHARM_INPUT = (
    ROOT.gSystem.Getenv("EOSSHIP")
    + "/eos/experiment/ship/data/Charm/"
    "Cascade-parp16-MSTP82-1-MSEL4-76Mpot_1.root"
)

ap.add_argument(
    "-I",
    "--InputFile",
    type=str,
    dest="charmInputFile",
    default=None,
    help=(
        "External heavy-flavour input ROOT file. "
        "For -C/--charm, the standard Charm/Cascade file is used if omitted. "
        "For -B/--beauty, this option is REQUIRED and must point to a beauty "
        "input file containing B hadrons."
    ),
)
ap.add_argument("-o", "--output", type=str, help="output directory", dest="work_dir", default=None)
ap.add_argument(
    "-rs", "--seed", type=int, help="random seed; default value is 0, see TRrandom::SetSeed documentation", default=0
)
ap.add_argument(
    "--DecayVolumeMedium",
    help="Set Decay Volume Medium. Choices are helium (default) or vacuums.",
    default="helium",
    choices=["helium", "vacuums"],
)
ap.add_argument(
    "--shieldName",
    help="Name of the shield in the database.",
    default="TRY_2026",
    choices=["TRY_2025", "TRY_2026"],
)
ap.add_argument(
    "--AddMuonShield",
    help="Whether or not to add the muon shield. Default set to False.",
    default=False,
    action=argparse.BooleanOptionalAction,
)
ap.add_argument(
    "--AddMuonShieldField",
    help="Whether or not to add the muon shield magnetic field. Default set to False.",
    default=False,
    action=argparse.BooleanOptionalAction,
)
ap.add_argument(
    "--AddHadronAbsorberOnly",
    help="Whether to only add the hadron absorber part of the muon shield. Default set to True.",
    default=True,
    action=argparse.BooleanOptionalAction,
)

ap.add_argument(
    "--z-offset", type=float, dest="z_offset", default=-84.0, help="z-offset for the FixedTargetGenerator [mm]"
)
ap.add_argument(
    "--x-offset", type=float, dest="x_offset", default=0.0, help="x-offset for the FixedTargetGenerator [mm]"
)
ap.add_argument(
    "--y-offset", type=float, dest="y_offset", default=0.0, help="y-offset for the FixedTargetGenerator [mm]"
)
ap.add_argument(
    "--beam-smear", type=float, dest="beam_smear", default=16.0, help="beam smearing for the FixedTargetGenerator [mm]"
)
ap.add_argument(
    "--beam-paint",
    type=float,
    dest="beam_paint",
    default=50.0,
    help="beam painting radius for the FixedTargetGenerator [mm]",
)
ap.add_argument(
    "--TARGET_YAML",
    dest="TARGET_YAML",
    help="File for target configuration",
    default=os.path.expandvars("$FAIRSHIP/geometry/target_config.yaml"),
)

ap.add_argument(
    "--AddCylindricalSensPlane",
    action="store_true",
    help="Whether or not to add cylindrical sensitive plane around the target. False by default.",
)
ap.add_argument(
    "--AddPostTargetSensPlane",
    action="store_true",
    help="Whether or not to add sensitive plane after the target. False by default.",
)

# --- charged B/Bc -> mu HNL production study -------------------------------
hnl_signal_mode = ap.add_mutually_exclusive_group()
hnl_signal_mode.add_argument(
    "--hnl-parent",
    choices=["b", "bc"],
    default=None,
    help=(
        "Force the charged-parent decay and retain BOTH charge conjugates: "
        "'b' = B+/- (|PDG|=521), 'bc' = Bc+/- (|PDG|=541). "
        "Requires -B/--beauty and -P/--pythiaDecay. "
        "Bc belongs to the beauty/MSEL=5 cascade."
    ),
)
# Legacy charge-specific B switches are kept for backward compatibility.
hnl_signal_mode.add_argument(
    "--bplus-hnl",
    action="store_true",
    help="Legacy: retain only B+ -> mu+ HNL.",
)
hnl_signal_mode.add_argument(
    "--bminus-hnl",
    action="store_true",
    help="Legacy: retain only B- -> mu- HNL.",
)
ap.add_argument("--hnl-mass", type=float, default=1.0, help="HNL mass [GeV]")
ap.add_argument("--hnl-pdg", type=int, default=9900015, help="HNL PDG code")
ap.add_argument(
    "--hnl-ctau-mm",
    type=float,
    default=1.0e-3,
    help=(
        "Technical Pythia proper decay length [mm]. For this production-only study the HNL is decayed "
        "to nu_mu anti-nu_mu so it is retained in MCTrack but is not transported as an unknown Geant4 particle."
    ),
)

args = ap.parse_args()
if args.debug:
    logger.setLevel(logging.DEBUG)

if args.kaon_pion_splits < 0:
    ap.error("--kaon-pion-splits must be >= 0")
if args.multiple_kpi_splits and args.kaon_pion_splits == 0:
    ap.error("--multiple-kpi-splits requires --kaon-pion-splits > 0")
# Canonical HNL mode used throughout this script.
if args.hnl_parent is not None:
    hnl_mode = args.hnl_parent
elif args.bplus_hnl:
    hnl_mode = "bplus"
elif args.bminus_hnl:
    hnl_mode = "bminus"
else:
    hnl_mode = None

PARENT_ABS_PDG = {
    "b": 521,
    "bc": 541,
    "bplus": 521,
    "bminus": 521,
}
PARENT_MASS_GEV = {
    "b": 5.27934,
    "bc": 6.27447,
    "bplus": 5.27934,
    "bminus": 5.27934,
}

if hnl_mode is not None:
    if not args.beauty:
        ap.error(
            f"HNL mode {hnl_mode!r} requires -B/--beauty. "
            "Bc is part of the beauty/MSEL=5 cascade, not the charm/MSEL=4 cascade."
        )
    if not args.pythiaDecay:
        ap.error(
            f"HNL mode {hnl_mode!r} requires -P/--pythiaDecay; "
            "EvtGen would otherwise own the charged B/Bc decay."
        )
    if args.hnl_mass <= 0.0:
        ap.error("--hnl-mass must be positive")

    threshold = PARENT_MASS_GEV[hnl_mode] - 0.105658
    if args.hnl_mass >= threshold:
        ap.error(
            f"--hnl-mass={args.hnl_mass:g} GeV is above the two-body threshold "
            f"for {hnl_mode} -> mu HNL ({threshold:.3f} GeV)."
        )

    if args.hnl_ctau_mm <= 0.0:
        ap.error("--hnl-ctau-mm must be positive")


if args.G4only:
    args.charm = False
    args.beauty = False
    withEvtGen = False
    args.pythiaDecay = False
elif args.pythiaDecay:
    withEvtGen = False
    logger.info("use Pythia8 as primary decayer")
else:
    withEvtGen = True
    logger.info("use EvtGen as primary decayer")
# withEvtGen = args.withEvtGen
if args.charm and args.beauty:
    logger.warning("charm and beauty decays are set! Beauty gets priority")
    args.charm = False

# Heavy-flavour external-input validation.
if args.charm and args.charmInputFile is None:
    args.charmInputFile = DEFAULT_CHARM_INPUT

if args.beauty and args.charmInputFile is None:
    ap.error(
        "-B/--beauty requires -I/--InputFile pointing to a BEAUTY input "
        "ROOT file containing B hadrons. The historical default input is a "
        "Charm/Cascade file and cannot produce B+/- -> mu+/- HNL."
    )

charmInputFile = args.charmInputFile

if args.beauty:
    normalized_input = str(charmInputFile).replace("\\", "/")
    if "/Charm/" in normalized_input and "Cascade" in normalized_input:
        ap.error(
            "Beauty mode received a file that looks like the standard "
            f"Charm/Cascade input: {charmInputFile}. "
            "Use -I with your beauty/MSEL=5 production ROOT input instead."
        )

    print("=" * 72)
    print("Beauty external input:")
    print(f"  {charmInputFile}")

    if hnl_mode is not None:
        required_parent = PARENT_ABS_PDG[hnl_mode]
        print(f"Required charged parent: |PDG|={required_parent}")

        check_file = ROOT.TFile.Open(charmInputFile, "READ")
        if not check_file or check_file.IsZombie():
            ap.error(f"Could not open beauty input file: {charmInputFile}")

        check_tree = check_file.Get("pythia6")
        if check_tree:
            n_parent = check_tree.Draw(
                "id", f"abs(id)=={required_parent}", "goff"
            )
            if n_parent <= 0:
                check_file.Close()
                ap.error(
                    f"Beauty input contains no |PDG|={required_parent} entries. "
                    "Generate/use an MSEL=5 cascade containing the requested parent."
                )
            print(f"Found {n_parent} input rows with |PDG|={required_parent}")
        else:
            print(
                "WARNING: no 'pythia6' tree found; parent-species presence "
                "could not be checked here."
            )
        check_file.Close()

    print("=" * 72)

if args.work_dir is None:
    if args.charm:
        args.work_dir = get_work_dir(args.runnr, "charm")
    if args.beauty:
        beauty_tag = "beauty"
        if hnl_mode == "b":
            beauty_tag = "beauty_BpmHNL"
        elif hnl_mode == "bc":
            beauty_tag = "beauty_BcpmHNL"
        elif hnl_mode == "bplus":
            beauty_tag = "beauty_BplusHNL"
        elif hnl_mode == "bminus":
            beauty_tag = "beauty_BminusHNL"
        args.work_dir = get_work_dir(args.runnr, beauty_tag)
    else:
        args.work_dir = get_work_dir(args.runnr)

logger.debug("work_dir: %s" % args.work_dir)
logger.debug("command line arguments: %s", args)
if os.path.exists(args.work_dir):
    logger.warning("output directory '%s' already exists." % args.work_dir)
    if args.force:
        logger.warning("...cleaning")
        for root, dirs, files in os.walk(args.work_dir):
            for f in files:
                os.unlink(os.path.join(root, f))
            for d in dirs:
                shutil.rmtree(os.path.join(root, d))
    else:
        logger.warning("...use '-f' option to overwrite it")
else:
    os.makedirs(args.work_dir)

os.chdir(args.work_dir)
# -------------------------------------------------------------------
# PYTHIA8 requires Random:seed to be in range [0, 900000000]
# When seed=0, ROOT generates a seed from system time which can exceed this limit
seed = args.seed
if seed == 0:
    ROOT.gRandom.SetSeed(0)  # Generate time-based seed
    seed = ROOT.gRandom.GetSeed()
# Clamp to PYTHIA8's maximum allowed seed value
if seed > 900000000:
    seed = seed % 900000000
ROOT.gRandom.SetSeed(seed)
shipRoot_conf.configure()  # load basic libraries, prepare atexit for python
if args.reproducible and not args.debug:
    ROOT.gErrorIgnoreLevel = ROOT.kWarning
ship_geo_kwargs = {
    "Yheight": dy,
    "DecayVolumeMedium": args.DecayVolumeMedium,
    "shieldName": args.shieldName,
    "TARGET_YAML": args.TARGET_YAML,
}
ship_geo = geometry_config.create_config(**ship_geo_kwargs)

txt = "pythia8_Geant4_"
if withEvtGen:
    txt = "pythia8_evtgen_Geant4_"

if hnl_mode == "b":
    txt += "BpmHNL_"
elif hnl_mode == "bc":
    txt += "BcpmHNL_"
elif hnl_mode == "bplus":
    txt += "BplusHNL_"
elif hnl_mode == "bminus":
    txt += "BminusHNL_"

outFile = f"{outputDir}/{txt}{args.runnr}_{args.ecut}.root"
parFile = f"{outputDir}/ship.params.{txt}{args.runnr}_{args.ecut}.root"

# -----Timer--------------------------------------------------------
timer = ROOT.TStopwatch()
timer.Start()

# -----Create simulation run----------------------------------------
run = ROOT.FairRunSim()
run.SetName(mcEngine)  # Transport engine
if hasattr(run, "SetRunId"):
    run.SetRunId(args.runnr)
sink = ROOT.FairRootFileSink(outFile)
run.SetSink(sink)
ROOT.SetOwnership(sink, False)  # C++ FairRun takes ownership
if args.boostFactor > 1:
    # Turn off UseGeneralProcess to access GammaToMuons directly when cross-sections need to be changed
    os.environ["SET_GENERAL_PROCESS_TO_FALSE"] = "1"
if args.kaon_pion_splits > 0:
    os.environ["KAON_PION_SPLITS"] = str(args.kaon_pion_splits)
run.SetUserConfig("g4Config.C")  # user configuration file default g4Config.C
rtdb = run.GetRuntimeDb()

# -----Materials----------------------------------------------
run.SetMaterials("media.geo")
# -----Create geometry----------------------------------------------
cave = ROOT.ShipCave("CAVE")
cave.SetGeometryFileName("caveWithAir.geo")

run.AddModule(cave)
ROOT.SetOwnership(cave, False)  # C++ FairRunSim takes ownership

TargetStation = ROOT.ShipTargetStation(
    name="TargetStation",
    tl=ship_geo.target.length,
    tz=ship_geo.target.z,
    nS=ship_geo.target.nS,
    HeT=ship_geo.target.HeT,
)
TargetStation.SetLayerPosMat(
    d=ship_geo.target.xy,
    L=ship_geo.target.slices_length,
    G=ship_geo.target.slices_gap,
    M=ship_geo.target.slices_material,
)
target_version = getattr(ship_geo.target, "version", 1)
TargetStation.SetDesign(target_version)
if target_version >= 2:
    TargetStation.SetLastDiskDiameter(ship_geo.target.xy2)
TargetStation.SetShieldingReferenceLength(ship_geo.target.length_fixed)
run.AddModule(TargetStation)
ROOT.SetOwnership(TargetStation, False)  # C++ FairRunSim takes ownership


if args.AddPostTargetSensPlane:
    sensPlanePostT = ROOT.exitHadronAbsorber()
    sensPlanePostT.SetEnergyCut(args.ecut * u.GeV)
    sensPlanePostT.SetVetoPointName("PlanePostT")
    # by default, if the z-position is not set, the positioning is behind the hadron abosorber and the tracks are stopped when they hit the sens plane
    # if the z-position is set and has a reasonable value (below 1E8), then the tracks are not stopped and continue to the last plane after the hadron absorber
    sensPlanePostT.SetZposition(
        ship_geo.target.length + 7.6 * u.cm + 300 * u.mm
    )  # target length + vessel shift + shielding length
    sensPlanePostT.SetUseCaveCoordinates()  # position set from the cave to avoid extrusions since the plane is larger than the target vacuum box

    if args.storeOnlyMuons:
        sensPlanePostT.SetOnlyMuons()
    if args.skipNeutrinos:
        sensPlanePostT.SkipNeutrinos()
    if args.FourDP:
        sensPlanePostT.SetOpt4DP()
    run.AddModule(sensPlanePostT)
    ROOT.SetOwnership(sensPlanePostT, False)  # C++ FairRunSim takes ownership


if args.AddMuonShield or args.AddHadronAbsorberOnly:
    n_params = 15
    if not args.AddMuonShieldField:
        for i in range(ship_geo.muShield.nMagnets):
            ship_geo.muShield.params[i * n_params + 14] = 0  # set B field to 0
    if args.AddHadronAbsorberOnly:
        ship_geo.muShield.params = ship_geo.muShield.params[:15]  # set dXIn to 0

    MuonShield = ROOT.ShipMuonShield(
        in_params=list(ship_geo.muShield.params),
        z=ship_geo.muShield.z,
        WithConstShieldField=True,
        SC_key=ship_geo.SC_mag,
    )
    # MuonShield.SetSupports(False) # otherwise overlap with sensitive Plane
    run.AddModule(MuonShield)  # needs to be added because of magn hadron shield.
    ROOT.SetOwnership(MuonShield, False)  # C++ FairRunSim takes ownership


sensPlaneHA = ROOT.exitHadronAbsorber()
sensPlaneHA.SetNSplits(args.kaon_pion_splits)  # type: ignore[missing-attribute]
if args.multiple_kpi_splits:
    sensPlaneHA.SetSplitMultipleTimes()  # type: ignore[missing-attribute]
sensPlaneHA.SetEnergyCut(args.ecut * u.GeV)
sensPlaneHA.SetVetoPointName("PlaneHA")

sensPlaneT = None
if args.AddCylindricalSensPlane:  # add additional sensitive plane around target
    sensPlaneT = ROOT.exitHadronAbsorber()
    sensPlaneT.SetEnergyCut(args.ecut * u.GeV)
    sensPlaneT.SetVetoPointName("PlaneT")
    sensPlaneT.SetCylindricalPlane()
    # by default, if the z-position is not set, the positioning is behind the hadron abosorber and the tracks are stopped when they hit the sens plane
    # if the z-position is set and has a reasonable value (below 1E8), then the tracks are not stopped and continue to the last plane after the hadron absorber
    sensPlaneT.SetZposition(ship_geo.target.length)

if args.storeOnlyMuons:
    sensPlaneHA.SetOnlyMuons()
    if sensPlaneT is not None:
        sensPlaneT.SetOnlyMuons()
if args.skipNeutrinos:
    sensPlaneHA.SkipNeutrinos()
    if sensPlaneT is not None:
        sensPlaneT.SkipNeutrinos()
if args.FourDP:  # in case a ntuple should be filled with pi0,etas,omega
    sensPlaneHA.SetOpt4DP()
    if sensPlaneT is not None:
        sensPlaneT.SetOpt4DP()

run.AddModule(sensPlaneHA)
ROOT.SetOwnership(sensPlaneHA, False)  # C++ FairRunSim takes ownership

if args.AddCylindricalSensPlane:
    run.AddModule(sensPlaneT)
    ROOT.SetOwnership(sensPlaneT, False)  # C++ FairRunSim takes ownership

# -----Create PrimaryGenerator--------------------------------------
primGen = ROOT.FairPrimaryGenerator()
P8gen = ROOT.FixedTargetGenerator()
P8gen.SetZoffset(args.z_offset * u.mm)
P8gen.SetXoffset(args.x_offset * u.mm)
P8gen.SetYoffset(args.y_offset * u.mm)
P8gen.SetSmearBeam(args.beam_smear * u.mm)
P8gen.SetPaintRadius(args.beam_paint * u.mm)
# Use geometry constants instead of fragile TGeo navigation
P8gen.SetTargetCoordinates(ship_geo.target.z0, ship_geo.target.z0 + ship_geo.target.length)
P8gen.SetMom(400.0 * u.GeV)
P8gen.SetEnergyCut(args.ecut * u.GeV)
P8gen.SetDebug(args.debug)
P8gen.SetHeartBeat(100000)
if args.G4only:
    P8gen.SetG4only()
if args.JpsiMainly:
    P8gen.SetJpsiMainly()
if args.tauOnly:
    P8gen.SetTauOnly()
if withEvtGen:
    P8gen.WithEvtGen()
if args.boostDiMuon > 1:
    P8gen.SetBoost(
        args.boostDiMuon
    )  # will increase BR for rare eta,omega,rho ... mesons decaying to 2 muons in Pythia8
    # and later copied to Geant4
P8gen.SetSeed(seed)
# for charm/beauty
#        print ' for experts: p pot= number of protons on target per spill to normalize on'
#        print '            : c chicc= ccbar over mbias cross section'
if args.charm or args.beauty:
    check_run_type_override(args.beauty, args.chicc, args.chibb)
    cs = derive_cross_sections(args.target_composition, args.A, args.chicc, args.chibb)
    P8gen.SetChicc(cs.chicc)
    P8gen.SetChibb(cs.chibb)
    print(format_summary(cs, None if args.A is not None else args.target_composition))
    print("--- process heavy flavours ---")
    P8gen.InitForCharmOrBeauty(charmInputFile, args.nev, args.pot, args.nStart)
primGen.AddGenerator(P8gen)
ROOT.SetOwnership(P8gen, False)  # C++ FairPrimaryGenerator takes ownership
#
run.SetGenerator(primGen)
ROOT.SetOwnership(primGen, False)  # C++ FairRunSim takes ownership

# -----Initialize simulation run------------------------------------
run.Init()

# Configure charged B/Bc -> mu HNL after FixedTargetGenerator has created
# its Pythia8 object. PYTHIA configures particle/antiparticle decays through
# the positive particle-data entry; the negative parent is charge conjugated.
if hnl_mode is not None:
    p8 = P8gen.GetPythia()
    hnl = args.hnl_pdg
    parent_abs = PARENT_ABS_PDG[hnl_mode]

    commands = [
        (
            f"{hnl}:new = N2 N2 2 0 0 "
            f"{args.hnl_mass:.12g} 0.0 0.0 0.0 "
            f"{args.hnl_ctau_mm:.12g} 0 1 0 1 0"
        ),
        f"{hnl}:isResonance = false",
        f"{hnl}:mayDecay = on",
        f"{hnl}:oneChannel = 1 1.0 0 14 -14",
        # positive charged parent -> mu+ (-13) + HNL;
        # the negative charged parent is automatic charge conjugation.
        f"{parent_abs}:oneChannel = 1 1.0 0 -13 {hnl}",
    ]

    for command in commands:
        ok = p8.readString(command)
        if not ok:
            raise RuntimeError(f"PYTHIA rejected command: {command}")

    pdg = ROOT.TDatabasePDG.Instance()
    if not pdg.GetParticle(hnl):
        pdg.AddParticle(
            "N2", "N2", args.hnl_mass,
            True, 0.0, 0.0, "HNL", hnl
        )

    if hnl_mode == "b":
        selected_decay = "B+/- -> mu+/- HNL (both charges)"
    elif hnl_mode == "bc":
        selected_decay = "Bc+/- -> mu+/- HNL (both charges)"
    elif hnl_mode == "bplus":
        selected_decay = "B+ -> mu+ HNL"
    else:
        selected_decay = "B- -> mu- HNL"

    print("=" * 72)
    print(f"Requested charged-parent signal: {selected_decay}")
    print(
        f"PYTHIA configured through +{parent_abs}; "
        "the negative parent uses the automatic charge-conjugate decay."
    )
    print(
        f"HNL mass={args.hnl_mass:g} GeV, "
        f"PDG={hnl}, technical ctau={args.hnl_ctau_mm:g} mm"
    )
    print("=" * 72)
    p8.particleData.list(parent_abs)
    p8.particleData.list(hnl)

gMC = ROOT.TVirtualMC.GetMC()
fStack = gMC.GetStack()
fStack.SetMinPoints(1)
fStack.SetEnergyCut(-1.0)
if args.kaon_pion_splits > 0:
    fStack.SetSplitting()
#
import AddDiMuonDecayChannelsToG4

AddDiMuonDecayChannelsToG4.Initialize(P8gen.GetPythia())

# boost gamma2muon conversion
if args.boostFactor > 1:
    ROOT.gROOT.ProcessLine('#include "Geant4/G4ProcessTable.hh"')
    ROOT.gROOT.ProcessLine('#include "Geant4/G4AnnihiToMuPair.hh"')
    ROOT.gROOT.ProcessLine('#include "Geant4/G4GammaConversionToMuons.hh"')
    gProcessTable = ROOT.G4ProcessTable.GetProcessTable()
    procAnnihil = gProcessTable.FindProcess(ROOT.G4String("AnnihiToMuPair"), ROOT.G4String("e+"))
    procGMuPair = gProcessTable.FindProcess(ROOT.G4String("GammaToMuPair"), ROOT.G4String("gamma"))
    procAnnihil.SetCrossSecFactor(args.boostFactor)
    procGMuPair.SetCrossSecFactor(args.boostFactor)

# -----Start run----------------------------------------------------
run.Run(args.nev)

# -----Finish-------------------------------------------------------
timer.Stop()
rtime = timer.RealTime()
ctime = timer.CpuTime()
print(" ")
print("Macro finished successfully.")
print(f"Output file is {outFile}")
if not args.reproducible:
    print(f"Real time {rtime} s, CPU time {ctime} s")
# ---post processing--- remove empty events --- save histograms
tmpFile = outFile + "tmp"
if ROOT.gROOT.GetListOfFiles().GetEntries() > 0:
    fin = ROOT.gROOT.GetListOfFiles()[0]
else:
    fin = ROOT.TFile.Open(outFile)
fHeader = fin.Get("FileHeader")
if fHeader:
    fHeader.SetRunId(args.runnr)
else:
    print("WARNING: FileHeader not found in simulation output; skipped FileHeader RunID update")
if args.charm or args.beauty:
    # normalization for charm
    poteq = P8gen.GetPotForCharm()
    info = "POT equivalent = %7.3G" % (poteq)
else:
    info = f"POT = {args.nev}"

conditions = " with ecut=" + str(args.ecut)
if args.JpsiMainly:
    conditions += " J"
if args.tauOnly:
    conditions += " T"
if withEvtGen:
    conditions += " V"
if args.boostDiMuon > 1:
    conditions += " diMu" + str(args.boostDiMuon)
if args.boostFactor > 1:
    conditions += " X" + str(args.boostFactor)

info += conditions
if fHeader:
    fHeader.SetTitle(info)
    print(f"Data generated {fHeader.GetTitle()}")
else:
    print(f"Data generated {info}")

nt = fin.Get("4DP")
if nt:
    nt = fin["4DP"]
    tf = ROOT.TFile("FourDP.root", "recreate")
    tnt = nt.CloneTree(0)
    for i in range(nt.GetEntries()):
        rc = nt.GetEvent(i)
        rc = tnt.Fill(nt.id, nt.px, nt.py, nt.pz, nt.x, nt.y, nt.z)
    tnt.Write()
    tf.Close()

t = fin["cbmsim"]
fout = ROOT.TFile(tmpFile, "recreate")
sTree = t.CloneTree(0)

nEvents = 0
nHNLRaw = 0
nSignalPlusRaw = 0
nSignalMinusRaw = 0


def event_has_signal_chain(tracks, mother_pdg, muon_pdg, hnl_pdg):
    """Require parent -> mu + HNL with both daughters from the same direct mother."""
    hnl_mothers = set()

    for i_tr in range(len(tracks)):
        tr = tracks[i_tr]
        if tr.GetPdgCode() != hnl_pdg:
            continue
        mother_id = tr.GetMotherId()
        if mother_id >= 0:
            hnl_mothers.add(mother_id)

    for i_tr in range(len(tracks)):
        mu = tracks[i_tr]
        if mu.GetPdgCode() != muon_pdg:
            continue

        mother_id = mu.GetMotherId()
        if mother_id < 0 or mother_id >= len(tracks):
            continue
        if tracks[mother_id].GetPdgCode() != mother_pdg:
            continue
        if mother_id in hnl_mothers:
            return True

    return False


for n in range(t.GetEntries()):
    rc = t.GetEvent(n)

    has_positive_signal = False
    has_negative_signal = False

    if hasattr(t, "MCTrack"):
        tracks = t.MCTrack

        for i_tr in range(len(tracks)):
            if tracks[i_tr].GetPdgCode() == args.hnl_pdg:
                nHNLRaw += 1

        if hnl_mode is not None:
            parent_abs = PARENT_ABS_PDG[hnl_mode]
            has_positive_signal = event_has_signal_chain(
                tracks, parent_abs, -13, args.hnl_pdg
            )
            has_negative_signal = event_has_signal_chain(
                tracks, -parent_abs, 13, args.hnl_pdg
            )

            if has_positive_signal:
                nSignalPlusRaw += 1
            if has_negative_signal:
                nSignalMinusRaw += 1

    if hnl_mode == "bplus":
        keep_event = has_positive_signal
    elif hnl_mode == "bminus":
        keep_event = has_negative_signal
    elif hnl_mode in ("b", "bc"):
        keep_event = has_positive_signal or has_negative_signal
    else:
        keep_event = (
            (len(t.PlaneHAPoint) > 0)
            or (
                args.AddCylindricalSensPlane
                and len(t.PlaneTPoint) > 0
            )
            or (
                args.AddPostTargetSensPlane
                and len(t.PlanePostTPoint) > 0
            )
        )

    if keep_event:
        rc = sTree.Fill()
        nEvents += 1

if hnl_mode is not None:
    parent_abs = PARENT_ABS_PDG[hnl_mode]
    parent_label = "B" if parent_abs == 521 else "Bc"
    print("=" * 72)
    print("RAW cbmsim diagnostic before signal filtering")
    print(f"  HNL MCTracks (PDG {args.hnl_pdg}) : {nHNLRaw}")
    print(f"  {parent_label}+ -> mu+ HNL events             : {nSignalPlusRaw}")
    print(f"  {parent_label}- -> mu- HNL events             : {nSignalMinusRaw}")
    print(f"  Retained signal events              : {nEvents}")
    print("=" * 72)

fout.cd()
for k in fin.GetListOfKeys():
    x = fin.Get(k.GetName())
    className = x.Class().GetName()
    if className.find("TTree") < 0 and className.find("TNtuple") < 0:
        xcopy = x.Clone()
        rc = xcopy.Write()
sTree.AutoSave()
if fHeader:
    ff = fHeader.Clone(fout.GetName())
    fout.cd()
    ff.Write("FileHeader", ROOT.TObject.kSingleKey)
sTree.Write()
fout.Close()

rc1 = os.system("rm  " + outFile)
rc2 = os.system("mv " + tmpFile + " " + outFile)
print("removed out file, moved tmpFile to out file", rc1, rc2)

if rc1 == 0 and rc2 == 0:
    print("INFO: Adding file summary")
    fsr = vars(args)
    with ROOT.TFile.Open(outFile, "UPDATE") as _of:
        _of.WriteObject(ROOT.TString(json.dumps(fsr)), "FileSummary")
        if hnl_mode is not None:
            parent_abs = PARENT_ABS_PDG[hnl_mode]
            if hnl_mode == "b":
                process_name = "B+/- -> mu+/- HNL"
            elif hnl_mode == "bc":
                process_name = "Bc+/- -> mu+/- HNL"
            elif hnl_mode == "bplus":
                process_name = "B+ -> mu+ HNL"
            else:
                process_name = "B- -> mu- HNL"

            hnl_meta = {
                "process": process_name,
                "hnl_mode": hnl_mode,
                "parent_abs_pdg": parent_abs,
                "both_charges": hnl_mode in ("b", "bc"),
                "hnl_pdg": args.hnl_pdg,
                "hnl_mass_GeV": args.hnl_mass,
                "hnl_ctau_mm": args.hnl_ctau_mm,
                "beam_momentum_GeV": 400.0,
                "target_composition": args.target_composition,
                "target_z0_cm": float(ship_geo.target.z0 / u.cm),
                "target_z_end_cm": float(
                    (ship_geo.target.z0 + ship_geo.target.length) / u.cm
                ),
                "target_transverse_size_cm": float(
                    ship_geo.target.xy / u.cm
                ),
                "note": (
                    f"PYTHIA is configured through +{parent_abs}; "
                    "the negative parent is the automatic charge-conjugate decay. "
                    "Modes 'b' and 'bc' retain both parent charges."
                ),
            }
            _of.WriteObject(
                ROOT.TString(json.dumps(hnl_meta)),
                "HNLFixedTargetConfig",
            )
else:
    print("WARNING: tempFile mv or rm not successful. No attempt at FileSummary writing")

fin.SetWritable(False)  # bpyass flush error

if hnl_mode == "b":
    print(f"Number of retained B+/- -> mu+/- HNL events: {nEvents}")
elif hnl_mode == "bc":
    print(f"Number of retained Bc+/- -> mu+/- HNL events: {nEvents}")
elif hnl_mode == "bplus":
    print(f"Number of retained B+ -> mu+ HNL events: {nEvents}")
elif hnl_mode == "bminus":
    print(f"Number of retained B- -> mu- HNL events: {nEvents}")
else:
    print(f"Number of events produced with activity after hadron absorber: {nEvents}")

if checkOverlap:
    sGeo = ROOT.gGeoManager
    sGeo.CheckOverlaps()
    sGeo.PrintOverlaps()
    run.CreateGeometryFile("%s/geofile_full.root" % (outputDir))
    import saveBasicParameters

    saveBasicParameters.execute("%s/geofile_full.root" % (outputDir), ship_geo)
