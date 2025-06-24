# example for accessing smeared hits and fitted tracks
import os
import sys
import ROOT
import ctypes
import rootUtils as ut
import shipunit as u
from ShipGeoConfig import ConfigRegistry
from rootpyPickler import Unpickler
from decorators import *
import shipRoot_conf
from argparse import ArgumentParser

ROOT.gROOT.SetBatch(True)

shipRoot_conf.configure()
PDG = ROOT.TDatabasePDG.Instance()
#chi2CutOff  = 4.
fiducialCut = False
#measCutFK = 25
#measCutPR = 22
#docaCut = 2.

parser = ArgumentParser()

parser.add_argument("-f", "--inputFile", dest="inputFile", help="Input file", required=True)
parser.add_argument("-n", "--nEvents",   dest="nEvents",   help="Number of events to analyze", required=False,  default=999999,type=int)
parser.add_argument("-g", "--geoFile",   dest="geoFile",   help="ROOT geofile", required=True)
parser.add_argument("--Debug",           dest="Debug", help="Switch on debugging", required=False, action="store_true")
options = parser.parse_args()

if not options.inputFile.find(',')<0 :
  sTree = ROOT.TChain("cbmsim")
  for x in options.inputFile.split(','):
    sTree.AddFile(x)
else:
  f = ROOT.TFile(options.inputFile)
  sTree = f.Get("cbmsim")

if not options.geoFile:
  options.geoFile = options.inputFile.replace('ship.','geofile_full.').replace('_rec.','.')
else:
  fgeo = ROOT.TFile(options.geoFile)
  
# new geofile, load Shipgeo dictionary written by run_simScript.py
upkl    = Unpickler(fgeo)
ShipGeo = upkl.load('ShipGeo')
dy = ShipGeo.Yheight/u.m
# -----Create geometry----------------------------------------------
#import shipDet_conf
#run = ROOT.FairRunSim()
#run.SetName("TGeant4")  # Transport engine
#run.SetSink(ROOT.FairRootFileSink(ROOT.TMemFile('output', 'recreate')))  # Dummy output file
#run.SetUserConfig("g4Config_basic.C") # geant4 transport not used, only needed for the mag field
#rtdb = run.GetRuntimeDb()
# -----Create geometry----------------------------------------------
#modules = shipDet_conf.configure(run,ShipGeo)

import geomGeant4
#if hasattr(ShipGeo.Bfield,"fieldMap"):
#  fieldMaker = geomGeant4.addVMCFields(ShipGeo, '', True, withVirtualMC = False)
#else:
#  print("no fieldmap given, geofile too old, not anymore support")
#  exit(-1)
sGeo = fgeo.Get("FAIRGeom")
#geoMat =  ROOT.genfit.TGeoMaterialInterface()
#ROOT.genfit.MaterialEffects.getInstance().init(geoMat)
#bfield = ROOT.genfit.FairShipFields()
#bfield.setField(fieldMaker.getGlobalField())
#fM = ROOT.genfit.FieldManager.getInstance()
#fM.init(bfield)

blob_radius = 50.0  # cm
blob_z = ShipGeo.UpstreamTagger.z
target_z = ShipGeo.target.z0
#vessel_r = min(sGeo.vetoStation.XMax, sGeo.vetoStation.YMax)


#print(blob_radius, blob_z, target_z, vessel_r)

volDict = {}
i=0
for x in ROOT.gGeoManager.GetListOfVolumes():
 volDict[i]=x.GetName()
 i+=1

# prepare veto decisions
import shipVeto
veto = shipVeto.Task(sTree)
vetoDets={}
log={}
h = {}

ut.bookHist(h,'x_coordinate_true_vertex_decay', 'x-coordinate decay vertex (true); x [cm]; Events', 2000, -1000.0, 1000.0)
ut.bookHist(h,'y_coordinate_true_vertex_decay', 'y-coordinate decay vertex (true); y [cm]; Events', 2000, -1000.0, 1000.0)
ut.bookHist(h,'z_coordinate_true_vertex_decay', 'z-coordinate decay vertex (true); z [cm]; Events', 2000, -5000.0, 5000.0)
ut.bookHist(h,'xz_coordinate_true_vertex_decay', 'xz-plane coordinates decay vertex (true); x [cm]; z [cm]; Events', 2000, -1000.0, 1000.0, 2000, -1000.0, 1000.0)
ut.bookHist(h,'yz_coordinate_true_vertex_decay', 'yz-plane coordinates decay vertex (true); y [cm]; z [cm]; Events', 2000, -1000.0, 1000.0, 2000, -1000.0, 1000.0)
ut.bookHist(h,'xy_coordinate_true_vertex_decay', 'xy-plane coordinates decay vertex (true); x [cm]; y [cm]; Events', 2000, -1000.0, 1000.0, 2000, -1000.0, 1000.0)

ut.bookHist(h,'angle_theta_true_daughter_1_wrt_HNL', '#theta daughter 1 wrt HNL (true); #theta [degree]; Events', 400, 0.0, 25.0)
ut.bookHist(h,'angle_theta_true_daughter_2_wrt_HNL', '#theta daughter 2 wrt HNL (true); #theta [degree]; Events', 400, 0.0, 25.0)

ut.bookHist(h,'px_true_HNL', 'p_{x} HNL (true); p_{x} [GeV]; Events', 400, -1.0, 1.0)
ut.bookHist(h,'py_true_HNL', 'p_{y} HNL (true); p_{y} [GeV]; Events', 400, -1.0, 1.0)
ut.bookHist(h,'pz_true_HNL', 'p_{z} HNL (true); p_{z} [GeV]; Events', 400, 0.0, 100.0)
ut.bookHist(h,'pt_true_HNL', 'p_{T} HNL (true); p_{T} [GeV]; Events', 400, 0.0, 3.0)
ut.bookHist(h,'p_true_HNL', 'p HNL (true); p [GeV]; Events', 400, 0.0, 110)
ut.bookHist(h,'E_true_HNL', 'E HNL (true); E [GeV]; Events', 400, 0.0, 110)
ut.bookHist(h,'theta_true_HNL', '#theta HNL (true); #theta [degree]; Events', 400, 0.0, 3.0)
ut.bookHist(h,'phi_true_HNL', '#phi HNL (true); #phi [degree]; Events', 400, -180.0, 180.0)

ut.bookHist(h,'px_true_daughter1', 'p_{x} daughter 1 (true); p_{x} [GeV]; Events', 400, -1.0, 1.0)
ut.bookHist(h,'py_true_daughter1', 'p_{y} daughter 1 (true); p_{y} [GeV]; Events', 400, -1.0, 1.0)
ut.bookHist(h,'pz_true_daughter1', 'p_{z} daughter 1 (true); p_{z} [GeV]; Events', 400, 0.0, 100.0)
ut.bookHist(h,'pt_true_daughter1', 'p_{T} daughter 1 (true); p_{T} [GeV]; Events', 400, 0.0, 3.0)
ut.bookHist(h,'p_true_daughter1', 'p daughter 1 (true); p [GeV]; Events', 400, 0.0, 110)
ut.bookHist(h,'E_true_daughter1', 'E daughter 1 (true); E [GeV]; Events', 400, 0.0, 110)
ut.bookHist(h,'theta_true_daughter1', '#theta daughter 1 (true); #theta [degree]; Events', 400, 0.0, 3.0)
ut.bookHist(h,'phi_true_daughter1', '#phi daughter 1 (true); #phi [degree]; Events', 400, -180.0, 180.0)

ut.bookHist(h,'px_true_daughter2', 'p_{x} daughter 2 (true); p_{x} [GeV]; Events', 400, -1.0, 1.0)
ut.bookHist(h,'py_true_daughter2', 'p_{y} daughter 2 (true); p_{y} [GeV]; Events', 400, -1.0, 1.0)
ut.bookHist(h,'pz_true_daughter2', 'p_{z} daughter 2 (true); p_{z} [GeV]; Events', 400, 0.0, 100.0)
ut.bookHist(h,'pt_true_daughter2', 'p_{T} daughter 2 (true); p_{T} [GeV]; Events', 400, 0.0, 3.0)
ut.bookHist(h,'p_true_daughter2', 'p daughter 2 (true); p [GeV]; Events', 400, 0.0, 110)
ut.bookHist(h,'E_true_daughter2', 'E daughter 2 (true); E [GeV]; Events', 400, 0.0, 110)
ut.bookHist(h,'theta_true_daughter2', '#theta daughter 2 (true); #theta [degree]; Events', 400, 0.0, 3.0)
ut.bookHist(h,'phi_true_daughter2', '#phi daughter 2 (true); #phi [degree]; Events', 400, -180.0, 180.0)

ut.bookHist(h,'cos_theta_star', 'cos(#theta^{*}) in HNL rest frame; cos(#theta^{*}); Events', 400, 0.9, 1.0)

def makePlots():
  ut.bookCanvas(h, key='cos_theta_star', title='cos(#theta^{*}) in HNL rest frame', nx=800, ny=600, cx=1, cy=1)
  h['cos_theta_star'].Draw()

  ut.bookCanvas(h, key='x_coordinate_true_vertex_decay', title='x-coordinate decay vertex (true)', nx=800, ny=600, cx=1, cy=1)
  h['x_coordinate_true_vertex_decay'].Draw()
  ut.bookCanvas(h, key='y_coordinate_true_vertex_decay', title='y-coordinate decay vertex (true)', nx=800, ny=600, cx=1, cy=1)
  h['y_coordinate_true_vertex_decay'].Draw()
  ut.bookCanvas(h, key='z_coordinate_true_vertex_decay', title='z-coordinate decay vertex (true)', nx=800, ny=600, cx=1, cy=1)
  h['z_coordinate_true_vertex_decay'].Draw()

  ut.bookCanvas(h, key='xz_coordinate_true_vertex_decay', title='xz-plane coordinates decay vertex (true)', nx=800, ny=600, cx=1, cy=1)
  h['xz_coordinate_true_vertex_decay'].Draw()
  ut.bookCanvas(h, key='yz_coordinate_true_vertex_decay', title='yz-plane coordinates decay vertex (true)', nx=800, ny=600, cx=1, cy=1)
  h['yz_coordinate_true_vertex_decay'].Draw()
  ut.bookCanvas(h, key='xy_coordinate_true_vertex_decay', title='xy-plane coordinates decay vertex (true)', nx=800, ny=600, cx=1, cy=1)
  h['xy_coordinate_true_vertex_decay'].Draw()

  ut.bookCanvas(h, key='angle_theta_true_daughter_1_wrt_HNL', title='#theta daughter wrt HNL (true)', nx=800, ny=600, cx=1, cy=1)
  h['angle_theta_true_daughter_1_wrt_HNL'].Draw()
  ut.bookCanvas(h, key='angle_theta_true_daughter_2_wrt_HNL', title='#theta daughter wrt HNL (true)', nx=800, ny=600, cx=1, cy=1)
  h['angle_theta_true_daughter_2_wrt_HNL'].Draw()
  
  ut.bookCanvas(h, key='px_true_HNL', title='p_{x} HNL (true)', nx=800, ny=600, cx=1, cy=1)
  h['px_true_HNL'].Draw()
  ut.bookCanvas(h, key='py_true_HNL', title='p_{y} HNL (true)', nx=800, ny=600, cx=1, cy=1)
  h['py_true_HNL'].Draw()
  ut.bookCanvas(h, key='pz_true_HNL', title='p_{z} HNL (true)', nx=800, ny=600, cx=1, cy=1)
  h['pz_true_HNL'].Draw()
  ut.bookCanvas(h, key='pt_true_HNL', title='p_{T} HNL (true)', nx=800, ny=600, cx=1, cy=1)
  h['pt_true_HNL'].Draw()
  ut.bookCanvas(h, key='p_true_HNL', title='p HNL (true)', nx=800, ny=600, cx=1, cy=1)
  h['p_true_HNL'].Draw()
  ut.bookCanvas(h, key='E_true_HNL', title='E HNL (true)', nx=800, ny=600, cx=1, cy=1)
  h['E_true_HNL'].Draw()
  ut.bookCanvas(h, key='theta_true_HNL', title='#theta HNL (true)', nx=800, ny=600, cx=1, cy=1)
  h['theta_true_HNL'].Draw()
  ut.bookCanvas(h, key='phi_true_HNL', title='#phi HNL (true)', nx=800, ny=600, cx=1, cy=1)
  h['phi_true_HNL'].Draw()

  ut.bookCanvas(h, key='px_true_daughter1', title='p_{x} daughter 1 (true)', nx=800, ny=600, cx=1, cy=1)
  h['px_true_daughter1'].Draw()
  ut.bookCanvas(h, key='py_true_daughter1', title='p_{y} daughter 1 (true)', nx=800, ny=600, cx=1, cy=1)
  h['py_true_daughter1'].Draw()
  ut.bookCanvas(h, key='pz_true_daughter1', title='p_{z} daughter 1 (true)', nx=800, ny=600, cx=1, cy=1)
  h['pz_true_daughter1'].Draw()
  ut.bookCanvas(h, key='pt_true_daughter1', title='p_{T} daughter 1 (true)', nx=800, ny=600, cx=1, cy=1)
  h['pt_true_daughter1'].Draw()
  ut.bookCanvas(h, key='p_true_daughter1', title='p daughter 1 (true)', nx=800, ny=600, cx=1, cy=1)
  h['p_true_daughter1'].Draw()
  ut.bookCanvas(h, key='E_true_daughter1', title='E daughter 1 (true)', nx=800, ny=600, cx=1, cy=1)
  h['E_true_daughter1'].Draw()
  ut.bookCanvas(h, key='theta_true_daughter1', title='#theta daughter 1 (true)', nx=800, ny=600, cx=1, cy=1)
  h['theta_true_daughter1'].Draw()
  ut.bookCanvas(h, key='phi_true_daughter1', title='#phi daughter 1 (true)', nx=800, ny=600, cx=1, cy=1)
  h['phi_true_daughter1'].Draw()

  ut.bookCanvas(h, key='px_true_daughter2', title='p_{x} daughter 2 (true)', nx=800, ny=600, cx=1, cy=1)
  h['px_true_daughter2'].Draw()
  ut.bookCanvas(h, key='py_true_daughter2', title='p_{y} daughter 2 (true)', nx=800, ny=600, cx=1, cy=1)
  h['py_true_daughter2'].Draw()
  ut.bookCanvas(h, key='pz_true_daughter2', title='p_{z} daughter 2 (true)', nx=800, ny=600, cx=1, cy=1)
  h['pz_true_daughter2'].Draw()
  ut.bookCanvas(h, key='pt_true_daughter2', title='p_{T} daughter 2 (true)', nx=800, ny=600, cx=1, cy=1)
  h['pt_true_daughter2'].Draw()
  ut.bookCanvas(h, key='p_true_daughter2', title='p daughter 2 (true)', nx=800, ny=600, cx=1, cy=1)
  h['p_true_daughter2'].Draw()
  ut.bookCanvas(h, key='E_true_daughter2', title='E daughter 2 (true)', nx=800, ny=600, cx=1, cy=1)
  h['E_true_daughter2'].Draw()
  ut.bookCanvas(h, key='theta_true_daughter2', title='#theta daughter 2 (true)', nx=800, ny=600, cx=1, cy=1)
  h['theta_true_daughter2'].Draw()
  ut.bookCanvas(h, key='phi_true_daughter2', title='#phi daughter 2 (true)', nx=800, ny=600, cx=1, cy=1)
  h['phi_true_daughter2'].Draw()
  
  print('finished making plots')

import TrackExtrapolateTool

def getAngularDistributionsLLPs():
  for i, HNLtrack in enumerate(sTree.MCTrack):
    if HNLtrack.GetPdgCode() != 999: continue
    
    # decay products
    daughters = {}
    for j, daughterTracks in enumerate(sTree.MCTrack):
      if daughterTracks.GetMotherId() == i:
        reached_ubt = any(hit.GetTrackID() == j for hit in sTree.UpstreamTaggerPoint)
        reached_sbt = any(hit.GetTrackID() == j for hit in sTree.vetoPoint)
        reached_sst = any(hit.GetTrackID() == j for hit in sTree.strawtubesPoint)
        daughters[j] = [reached_ubt, reached_sbt, reached_sst]

    if len(daughters) != 2: continue

    daughter1_PDGcode = 11
    daughter2_PDGcode =	211
    daughter1 = daughter2 = None
    for daughter_index, daughter_detector_hit in daughters.items():
      d = sTree.MCTrack.At(daughter_index)
      pdg = abs(d.GetPdgCode())
      hitUBT = daughter_detector_hit[0]
      hitSBT = daughter_detector_hit[1]
      hitSST = daughter_detector_hit[2]
      if pdg == daughter1_PDGcode:
        daughter1 = d
      elif pdg == daughter2_PDGcode:
        daughter2 = d

    HNL_true_p4 = ROOT.Math.PxPyPzEVector(HNLtrack.GetPx(), HNLtrack.GetPy(), HNLtrack.GetPz(), HNLtrack.GetEnergy())
    daughter1_true_p4 = ROOT.Math.PxPyPzEVector(daughter1.GetPx(), daughter1.GetPy(), daughter1.GetPz(), daughter1.GetEnergy())
    daughter2_true_p4 = ROOT.Math.PxPyPzEVector(daughter2.GetPx(), daughter2.GetPy(), daughter2.GetPz(), daughter2.GetEnergy())

    HNL_true_lab_p4 = ROOT.TLorentzVector(HNL_true_p4.Px(), HNL_true_p4.Py(), HNL_true_p4.Pz(), HNL_true_p4.E())
    daughter1_true_lab_p4 = ROOT.TLorentzVector(daughter1_true_p4.Px(), daughter1_true_p4.Py(), daughter1_true_p4.Pz(), daughter1_true_p4.E())

    beta = -HNL_true_lab_p4.BoostVector()
    electron_rest = ROOT.TLorentzVector(daughter1_true_lab_p4)
    electron_rest.Boost(beta)

    hnl_direction_lab = HNL_true_lab_p4.Vect().Unit()
    electron_direction_rest = daughter1_true_p4.Vect().Unit()
    cosThetaStar = electron_direction_rest.Dot(hnl_direction_lab)

    x_coordinate_true_vertex_decay, y_coordinate_true_vertex_decay, z_coordinate_true_vertex_decay = daughter1.GetStartX(), daughter1.GetStartY(), daughter1.GetStartZ()
  
    cos_angle_1 = (daughter1_true_p4.Px()*HNL_true_p4.Px() + daughter1_true_p4.Py()*HNL_true_p4.Py() + daughter1_true_p4.Pz()*HNL_true_p4.Pz()) / (daughter1_true_p4.P() * HNL_true_p4.P())
    angle_1 = ROOT.TMath.ACos(cos_angle_1)*ROOT.TMath.RadToDeg()

    cos_angle_2 = (daughter2_true_p4.Px()*HNL_true_p4.Px() + daughter2_true_p4.Py()*HNL_true_p4.Py() + daughter2_true_p4.Pz()*HNL_true_p4.Pz()) / (daughter2_true_p4.P() * HNL_true_p4.P())
    angle_2 = ROOT.TMath.ACos(cos_angle_2)*ROOT.TMath.RadToDeg()

    h['cos_theta_star'].Fill(cosThetaStar)
    
    h['x_coordinate_true_vertex_decay'].Fill(x_coordinate_true_vertex_decay)
    h['y_coordinate_true_vertex_decay'].Fill(y_coordinate_true_vertex_decay)
    h['z_coordinate_true_vertex_decay'].Fill(z_coordinate_true_vertex_decay)
    h['xz_coordinate_true_vertex_decay'].Fill(x_coordinate_true_vertex_decay, z_coordinate_true_vertex_decay)
    h['yz_coordinate_true_vertex_decay'].Fill(y_coordinate_true_vertex_decay, z_coordinate_true_vertex_decay)
    h['xy_coordinate_true_vertex_decay'].Fill(x_coordinate_true_vertex_decay, y_coordinate_true_vertex_decay)

    h['angle_theta_true_daughter_1_wrt_HNL'].Fill(angle_1)
    h['angle_theta_true_daughter_2_wrt_HNL'].Fill(angle_2)
    
    h['px_true_HNL'].Fill(HNL_true_p4.Px())
    h['py_true_HNL'].Fill(HNL_true_p4.Py())
    h['pz_true_HNL'].Fill(HNL_true_p4.Pz())
    h['pt_true_HNL'].Fill(HNL_true_p4.Pt())
    h['p_true_HNL'].Fill(HNL_true_p4.P())
    h['E_true_HNL'].Fill(HNL_true_p4.E())
    h['theta_true_HNL'].Fill(HNL_true_p4.Theta()*ROOT.TMath.RadToDeg())
    h['phi_true_HNL'].Fill(HNL_true_p4.Phi()*ROOT.TMath.RadToDeg())
    
    h['px_true_daughter1'].Fill(daughter1_true_p4.Px())
    h['py_true_daughter1'].Fill(daughter1_true_p4.Py())
    h['pz_true_daughter1'].Fill(daughter1_true_p4.Pz())
    h['pt_true_daughter1'].Fill(daughter1_true_p4.Pt())
    h['p_true_daughter1'].Fill(daughter1_true_p4.P())
    h['E_true_daughter1'].Fill(daughter1_true_p4.E())
    h['theta_true_daughter1'].Fill(daughter1_true_p4.Theta()*ROOT.TMath.RadToDeg())
    h['phi_true_daughter1'].Fill(daughter1_true_p4.Phi()*ROOT.TMath.RadToDeg())
    
    h['px_true_daughter2'].Fill(daughter2_true_p4.Px())
    h['py_true_daughter2'].Fill(daughter2_true_p4.Py())
    h['pz_true_daughter2'].Fill(daughter2_true_p4.Pz())
    h['pt_true_daughter2'].Fill(daughter2_true_p4.Pt())
    h['p_true_daughter2'].Fill(daughter2_true_p4.P())
    h['E_true_daughter2'].Fill(daughter2_true_p4.E())
    h['theta_true_daughter2'].Fill(daughter2_true_p4.Theta()*ROOT.TMath.RadToDeg())
    h['phi_true_daughter2'].Fill(daughter2_true_p4.Phi()*ROOT.TMath.RadToDeg())
    
# start event loop
def myEventLoop(n):
  rc = sTree.GetEntry(n)
  sTree.GetEntry(n)

  getAngularDistributionsLLPs()
  
sTree.GetEvent(0)
options.nEvents = min(sTree.GetEntries(),options.nEvents)

for n in range(options.nEvents):
  myEventLoop(n)
  sTree.FitTracks.Delete()
makePlots()

# output histograms
hfile = options.inputFile.split(',')[0].replace('_rec','_ana')
if "/eos" in hfile or not options.inputFile.find(',')<0:
# do not write to eos, write to local directory
  tmp = hfile.split('/')
  hfile = tmp[len(tmp)-1]
ROOT.gROOT.cd()
ut.writeHists(h,hfile)
