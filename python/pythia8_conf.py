# SPDX-License-Identifier: LGPL-3.0-or-later
# SPDX-FileCopyrightText: Copyright CERN for the benefit of the SHiP Collaboration

import contextlib
import copy
import os

import hnl
import readDecayTable
import ROOT
import rpvsusy
import shipunit as u
import yaml
from method_logger import MethodLogger
from pythia8_conf_utils import (
    add_channel,
    add_particles,
    add_tau_channel,
    addHNLtoROOT,
    compute_max_total_br,
    exit_if_zero_br,
    fill_missing_channels,
    get_br,
    getbr_rpvsusy,
    getmaxsumbrrpvsusy,
    gettotalbrrpvsusy,
    make_interpolators,
    make_particles_stable,
    print_scale_factor,
)


def configurerpvsusy(
    P8gen, mass, couplings, sfermionmass, benchmark, inclusive, deepCopy: bool = False, debug: bool = True
) -> None:
    # configure pythia8 for Ship usage
    _exit_stack = contextlib.ExitStack()
    if debug:
        pythia_log = _exit_stack.enter_context(open("pythia8_conf.txt", "w"))  # noqa: SIM115
        P8gen = MethodLogger(P8gen, sink=pythia_log)
    h = make_interpolators(os.path.expandvars(f"$FAIRSHIP/shipgen/branchingratiosrpvsusybench{benchmark}.dat"))
    P8gen.UseRandom3()
    P8gen.SetMom(400)  # beam momentum in GeV
    if deepCopy:
        P8gen.UseDeepCopy()
    ROOT.TDatabasePDG.Instance()
    # let strange particle decay in Geant4
    make_particles_stable(P8gen, above_lifetime=1)

    if inclusive is True:  # For backward compatibility with a boolean argument
        inclusive = "True"

    if inclusive == "True":
        setup_pythia_inclusive(P8gen)

    # generate RPV neutralino from inclusive charm hadrons
    if inclusive == "c":
        P8gen.SetParameters("HardQCD::hardccbar  = on")
        # add RPVSUSY
        rpvsusy_instance = rpvsusy.RPVSUSY(mass, couplings, sfermionmass, benchmark, debug=True)
        ctau = rpvsusy_instance.computeNLifetime(system="FairShip") * u.c_light * u.cm
        print("RPVSUSY ctau ", ctau)
        P8gen.SetParameters(f"9900015:new = N2 N2 2 0 0 {mass:.12} 0.0 0.0 0.0 {ctau / u.mm:.12}  0   1   0   1   0")
        P8gen.SetParameters("9900015:isResonance = false")
        P8gen.SetParameters("Next:numberCount    =  0")
        # Configuring decay modes...
        rpvsusy_instance.AddChannelsToPythia(P8gen)

        # Finish HNL setup...
        P8gen.SetParameters("9900015:mayDecay = on")
        P8gen.SetHNLId(9900015)
        # also add to PDG
        gamma = u.hbarc / float(ctau)  # 197.3269631e-16 / float(ctau) # hbar*c = 197 MeV*fm = 197e-16 GeV*cm
        addHNLtoROOT(pid=9900015, m=mass, g=gamma)
        # 12 14 16 neutrinos replace with N2
        charmhistograms = ["d_mu", "ds_mu"]
        # no tau decay here to consider
        maxsumBR = getmaxsumbrrpvsusy(h, charmhistograms, mass, couplings)
        exit_if_zero_br(maxsumBR, inclusive, mass, particle="RPV neutralino")
        gettotalbrrpvsusy(h, charmhistograms, mass, couplings)

        # overwrite D_s+ decays
        P8gen.SetParameters(
            "431:new  D_s+  D_s-    1   3   0    1.96849"
            "    0.00000    0.00000    0.00000  1.49900e-01   0   1   0   1   0"
        )
        sumBR = 0.0
        br_ds_mu = getbr_rpvsusy(h, "ds_mu", mass, couplings[1])
        if br_ds_mu > 0.0:
            P8gen.SetParameters(f"431:addChannel      1  {br_ds_mu / maxsumBR:.12}    0      -13       9900015")
            sumBR += float(br_ds_mu / maxsumBR)
        if sumBR < 1.0 and sumBR > 0.0:
            P8gen.SetParameters(f"431:addChannel      1   {1.0 - sumBR:.12}    0       22      -11")

        # overwrite D+ decays
        P8gen.SetParameters(
            "411:new  D+ D-    1   3   0    1.86962    0.00000    0.00000    0.00000  3.11800e-01   0   1   0   1   0"
        )
        sumBR = 0.0
        br_d_mu = getbr_rpvsusy(h, "d_mu", mass, couplings[1])
        if br_d_mu > 0.0:
            P8gen.SetParameters(f"411:addChannel      1  {br_d_mu / maxsumBR:.12}    0      -13       9900015")
            sumBR += float(br_d_mu / maxsumBR)
        if sumBR < 1.0 and sumBR > 0.0:
            P8gen.SetParameters(f"411:addChannel      1   {1.0 - sumBR:.12}    0       22      -11")

        P8gen.List(9900015)

    if inclusive == "b":
        P8gen.SetParameters("HardQCD::hardbbbar  = on")
        # add RPVSUSY
        rpvsusy_instance = rpvsusy.RPVSUSY(mass, couplings, sfermionmass, benchmark, debug=True)
        ctau = rpvsusy_instance.computeNLifetime(system="FairShip") * u.c_light * u.cm
        P8gen.SetParameters(f"9900015:new = N2 N2 2 0 0 {mass:.12} 0.0 0.0 0.0 {ctau / u.mm:.12}  0   1   0   1   0")
        P8gen.SetParameters("9900015:isResonance = false")
        # Configuring decay modes...
        rpvsusy_instance.AddChannelsToPythia(P8gen)
        # Finish HNL setup...
        P8gen.SetParameters("9900015:mayDecay = on")
        P8gen.SetHNLId(9900015)
        # also add to PDG
        gamma = u.hbarc / float(ctau)  # 197.3269631e-16 / float(ctau) # hbar*c = 197 MeV*fm = 197e-16 GeV*cm
        addHNLtoROOT(pid=9900015, m=mass, g=gamma)
        # 12 14 16 neutrinos replace with N2
        beautyhistograms = ["b_mu", "b_tau", "b0_nu_mu", "b0_nu_tau"]
        maxsumBR = getmaxsumbrrpvsusy(h, beautyhistograms, mass, couplings)
        exit_if_zero_br(maxsumBR, inclusive, mass, particle="RPV neutralino")
        gettotalbrrpvsusy(h, beautyhistograms, mass, couplings)

        # overwrite B+ decays
        P8gen.SetParameters(
            "521:new  B+               B-    1   3   0    5.27925"
            "    0.00000    0.00000    0.00000  4.91100e-01   0   1   0   1   0"
        )
        sumBR = 0.0
        br_b_tau = getbr_rpvsusy(h, "b_tau", mass, couplings[1])
        if br_b_tau > 0.0:
            P8gen.SetParameters(f"521:addChannel      1  {br_b_tau / maxsumBR:.12}    0       9900015      -15")
            sumBR += float(br_b_tau / maxsumBR)
        if sumBR < 1.0 and sumBR > 0.0:
            P8gen.SetParameters(f"521:addChannel      1   {1.0 - sumBR:.12}    0       22      22")

        # overwrite B0 decays
        P8gen.SetParameters(
            "511:new  B0  Bbar0    1   0   0    5.27958"
            "    0.00000    0.00000    0.00000  4.58700e-01   0   1   0   1   0"
        )
        sumBR = 0.0
        br_b0_nu_tau = getbr_rpvsusy(h, "b0_nu_tau", mass, couplings[1])
        if br_b0_nu_tau > 0.0:
            P8gen.SetParameters(f"511:addChannel      1  {br_b0_nu_tau / maxsumBR:.12}   22       9900015      16")
            sumBR += float(br_b0_nu_tau / maxsumBR)
        if sumBR < 1.0 and sumBR > 0.0:
            P8gen.SetParameters(f"511:addChannel      1   {1.0 - sumBR:.12}    0       22      22")

        P8gen.List(9900015)

    _exit_stack.close()



def _normalize_beauty_process_selection(process_selection):
    """Normalise charged-B/Bc selection aliases."""
    aliases = {
        "B": "b",
        "b+": "bplus",
        "B+": "bplus",
        "b-": "bminus",
        "B-": "bminus",
        "Bc": "bc",
        "BC": "bc",
        "b_c": "bc",
    }
    return aliases.get(process_selection, process_selection)


def _resolve_beauty_selection(data, process_selection):
    """Resolve b, bc, bplus or bminus from hnl_production YAML."""
    mode = _normalize_beauty_process_selection(process_selection)

    if mode not in ("b", "bc", "bplus", "bminus"):
        raise ValueError(f"Unsupported beauty selection: {mode!r}")

    if mode not in data["selections"]:
        raise KeyError(
            f"HNL production YAML has no selections[{mode!r}]"
        )

    selection = data["selections"][mode]
    particles = list(selection.get("particles", []))

    if mode == "bplus" and particles != [521]:
        raise ValueError(
            "selections.bplus.particles must be [521]"
        )
    if mode == "bminus" and particles != [-521]:
        raise ValueError(
            "selections.bminus.particles must be [-521]"
        )
    if mode == "bc" and particles != [541]:
        raise ValueError(
            "selections.bc.particles must be [541]; PYTHIA supplies Bc- by charge conjugation"
        )

    return mode, selection


def _canonical_particle_ids_for_pythia(particles):
    """Map signed semantic IDs to PYTHIA particle-data IDs."""
    result = []
    for particle in particles:
        if isinstance(particle, int):
            result.append(abs(particle))
        else:
            result.append(particle)
    return result


def _channel_for_pythia(channel):
    """Convert signed B-/Bc- YAML semantics to PYTHIA's positive parent entry."""
    out = copy.deepcopy(channel)
    parent_id = int(out["id"])

    if parent_id >= 0:
        return out

    if parent_id not in (-521, -541):
        raise ValueError(
            "Negative-parent conversion is implemented only for B-/Bc-; "
            f"got {parent_id}"
        )

    out["id"] = abs(parent_id)

    if "idlepton" in out:
        lepton = int(out["idlepton"])
        if lepton != 13:
            raise ValueError(
                "Negative charged parent -> mu- HNL must use idlepton: 13 in YAML; "
                f"got {lepton}"
            )
        out["idlepton"] = -13

    if "idhadron" in out:
        hadron = int(out["idhadron"])
        if hadron != 0:
            out["idhadron"] = -hadron

    return out


def _load_hnl_production_yaml(fairship_root):
    """Load the HNL production YAML, optionally overridden by environment."""
    override = os.environ.get("HNL_PRODUCTION_YAML")
    if override:
        datafile = os.path.expandvars(override)
    else:
        datafile = fairship_root + "/python/hnl_production_test.yaml"

    if not os.path.exists(datafile):
        raise FileNotFoundError(
            f"HNL production YAML not found: {datafile}"
        )

    with open(datafile) as f:
        data = yaml.load(f, Loader=yaml.FullLoader)

    for key in ("particles", "selections", "channels"):
        if key not in data:
            raise ValueError(
                f"{datafile} is missing required top-level key {key!r}"
            )

    print(f"HNL production YAML: {datafile}")
    return datafile, data


def configure(
    P8gen, mass, production_couplings, decay_couplings, process_selection, deepCopy: bool = False, debug: bool = True
) -> None:
    """
    This function configures a HNLPythia8Generator instance for SHiP usage.
    """

    if process_selection is True:  # For backward compatibility
        process_selection = "inclusive"

    # Wrap the Pythia8 object into a class logging all of its method calls
    _exit_stack = contextlib.ExitStack()
    if debug:
        pythia_log = _exit_stack.enter_context(open("pythia8_conf.txt", "w"))  # noqa: SIM115
        P8gen = MethodLogger(P8gen, sink=pythia_log)

    fairship_root = os.environ["FAIRSHIP"]
    histograms = make_interpolators(fairship_root + "/shipgen/branchingratios.dat")
    P8gen.UseRandom3()  # TRandom1 or TRandom3 ?
    P8gen.SetMom(400)  # beam momentum in GeV
    if deepCopy:
        P8gen.UseDeepCopy()
    ROOT.TDatabasePDG.Instance()
    P8gen.SetParameters("Next:numberCount    =  0")
    # let strange particle decay in Geant4
    make_particles_stable(P8gen, above_lifetime=1)

    # Load particle & decay data
    # ==========================

    datafile, data = _load_hnl_production_yaml(fairship_root)
    all_channels = data["channels"]

    # Inclusive
    # =========

    if process_selection == "inclusive":
        setup_pythia_inclusive(P8gen)

    # Charm decays only (with secondary production from tau)
    # ======================================================

    if process_selection == "c":
        selection = data["selections"]["c"]
        for cmd in selection["parameters"]:
            P8gen.SetParameters(cmd)
        add_hnl(P8gen, mass, decay_couplings)

        # Add new charmed particles
        # -------------------------

        # Select all charmed particles (+ tau lepton)
        c_particles = selection["particles"]
        tau_id = 15  # tau- Monte-Carlo ID
        add_particles(P8gen, c_particles + [tau_id], data)

        # Add HNL production channels from charmed particles
        # --------------------------------------------------

        # Find charm and tau decays to HNLs
        c_channels = [ch for ch in all_channels if ch["id"] in c_particles]
        tau_channels = [ch for ch in all_channels if ch["id"] == tau_id]
        # Standard model process: tau+ production from D_s+ decay
        ds_id = 431  # D_s+ Monte-Carlo ID
        ds_tau_br = 0.0548  # SM branching ratio Br(D_s+ -> tau+ nu_tau) (source: PDG 2018)

        # Compute the branching ratio scaling factor, taking into account
        # secondary production from tau
        # Decay chains are encoded as follows:
        #     [(top level id A, [br A -> B, br B -> C, ...]), ...]

        # Most charm particles directly decay to HNLs
        primary_decays = [(ch["id"], [get_br(histograms, ch, mass, production_couplings)]) for ch in c_channels]
        # The D_s+ can indirectly produce a HNL by first producing a tau+
        secondary_decays = [
            (ds_id, [ds_tau_br, get_br(histograms, ch, mass, production_couplings)]) for ch in tau_channels
        ]
        all_decays = primary_decays + secondary_decays

        # Compute maximum total branching ratio (to rescale all BRs)
        max_total_br = compute_max_total_br(all_decays)
        exit_if_zero_br(max_total_br, process_selection, mass)
        print_scale_factor(1 / max_total_br)

        # Add charm decays
        for ch in c_channels:
            add_channel(P8gen, ch, histograms, mass, production_couplings, 1 / max_total_br)
        # Add tau production from D_s+
        # We can freely rescale Br(Ds -> tau) and Br(tau -> N X...) as long as
        # Br(Ds -> tau -> N X...) remains the same.
        # Here, we set Br(tau -> N) to unity to make event generation more efficient.
        # The implicit assumption here is that we will disregard the tau during the analysis.
        total_tau_br = sum(branching_ratios[1] for (_, branching_ratios) in secondary_decays)
        assert ds_tau_br * total_tau_br <= max_total_br + 1e-12
        P8gen.SetParameters(
            f"431:addChannel      1  {ds_tau_br * total_tau_br / max_total_br:.12}    0      -15       16"
        )
        # Add secondary HNL production from tau
        for ch in tau_channels:
            # Rescale branching ratios only if some are non-zero. Otherwise leave them at zero.
            add_tau_channel(P8gen, ch, histograms, mass, production_couplings, 1 / (total_tau_br or 1))

        # Add dummy channels in place of SM processes
        fill_missing_channels(P8gen, max_total_br, all_decays)

        # List channels to confirm that Pythia has been properly set up
        P8gen.List(9900015)

    # B/Bc/B+/B- decays only
    # ======================

    beauty_mode = _normalize_beauty_process_selection(process_selection)

    if beauty_mode in ["b", "bc", "bplus", "bminus"]:
        resolved_mode, selection = _resolve_beauty_selection(
            data,
            beauty_mode,
        )

        for cmd in selection["parameters"]:
            P8gen.SetParameters(cmd)

        add_hnl(P8gen, mass, decay_couplings)

        # The YAML selection uses signed PHYSICS IDs, while PYTHIA particle
        # definitions use the positive particle-data entry.
        semantic_particles = list(selection["particles"])
        pythia_particles = _canonical_particle_ids_for_pythia(
            semantic_particles
        )
        add_particles(P8gen, pythia_particles, data)

        # Select by signed YAML ID first: this distinguishes B+ from B-.
        semantic_channels = [
            ch
            for ch in all_channels
            if int(ch["id"]) in semantic_particles
        ]

        if not semantic_channels:
            raise ValueError(
                f"No decay channels found for {resolved_mode!r} "
                f"with particles={semantic_particles} in {datafile}"
            )

        # Convert B- semantics to valid +521 commands before using the
        # existing FairShip helper functions.
        pythia_channels = [
            _channel_for_pythia(ch)
            for ch in semantic_channels
        ]

        print(
            f"Resolved beauty selection {resolved_mode!r}: "
            f"semantic particles={semantic_particles}, "
            f"PYTHIA particles={pythia_particles}"
        )
        for semantic, canonical in zip(
            semantic_channels,
            pythia_channels,
        ):
            print(
                "  YAML channel "
                f"id={semantic['id']}, "
                f"idlepton={semantic.get('idlepton')} "
                "-> PYTHIA channel "
                f"id={canonical['id']}, "
                f"idlepton={canonical.get('idlepton')}"
            )

        # Optional forced two-body channels are useful for detector-occupancy
        # studies where the production kinematics are wanted but the physical
        # B/Bc -> mu HNL branching ratio is applied later as an event weight.
        forced_channels = [ch for ch in pythia_channels if "forced_br" in ch]

        if forced_channels:
            if len(forced_channels) != len(pythia_channels):
                raise ValueError(
                    "Do not mix forced_br and branching-ratio-table channels "
                    f"inside selection {resolved_mode!r}"
                )
            for ch in forced_channels:
                parent_id = abs(int(ch["id"]))
                br = float(ch["forced_br"])
                if not (0.0 < br <= 1.0):
                    raise ValueError(f"forced_br must be in (0,1], got {br}")
                if br != 1.0:
                    raise ValueError(
                        "Current forced two-body helper expects forced_br=1.0; "
                        "apply physical branching-ratio weights in analysis."
                    )
                idlepton = int(ch["idlepton"])
                P8gen.SetParameters(
                    f"{parent_id}:oneChannel = 1 1.0 0 {idlepton} 9900015"
                )
                print(
                    f"Forced occupancy channel: {parent_id} -> "
                    f"{idlepton} + 9900015 (BR set to 1 for generation)"
                )
        else:
            decays = [
                (
                    ch["id"],
                    [
                        get_br(
                            histograms,
                            ch,
                            mass,
                            production_couplings,
                        )
                    ],
                )
                for ch in pythia_channels
            ]

            max_total_br = compute_max_total_br(decays)
            exit_if_zero_br(max_total_br, resolved_mode, mass)
            print_scale_factor(1 / max_total_br)

            for ch in pythia_channels:
                add_channel(
                    P8gen,
                    ch,
                    histograms,
                    mass,
                    production_couplings,
                    1 / max_total_br,
                )

            fill_missing_channels(
                P8gen,
                max_total_br,
                decays,
            )

        for pid in sorted(set(pythia_particles)):
            P8gen.List(pid)
        P8gen.List(9900015)

    _exit_stack.close()


def add_hnl(P8gen, mass, decay_couplings) -> None:
    "Adds the HNL to Pythia and ROOT"
    hnl_instance = hnl.HNL(mass, decay_couplings, debug=True)
    ctau = hnl_instance.computeNLifetime(system="FairShip") * u.c_light * u.cm
    print(f"HNL ctau {ctau}")
    P8gen.SetParameters(f"9900015:new = N2 N2 2 0 0 {mass:.12} 0.0 0.0 0.0 {ctau / u.mm:.12}  0   1   0   1   0")
    P8gen.SetParameters("9900015:isResonance = false")
    # Configuring decay modes...
    readDecayTable.addHNLdecayChannels(
        P8gen, hnl_instance, conffile=os.path.expandvars("$FAIRSHIP/python/DecaySelection.conf"), verbose=False
    )
    # Finish HNL setup...
    P8gen.SetParameters("9900015:mayDecay = on")
    P8gen.SetHNLId(9900015)
    # also add to PDG
    gamma = u.hbarc / float(ctau)  # 197.3269631e-16 / float(ctau) # hbar*c = 197 MeV*fm = 197e-16 GeV*cm
    addHNLtoROOT(pid=9900015, m=mass, g=gamma)


def setup_pythia_inclusive(P8gen) -> None:
    P8gen.SetParameters("SoftQCD:inelastic = on")
    P8gen.SetParameters("PhotonCollision:gmgm2mumu = on")
    P8gen.SetParameters("PromptPhoton:all = on")
    P8gen.SetParameters("WeakBosonExchange:all = on")
