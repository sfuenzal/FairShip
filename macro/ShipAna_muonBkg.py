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

ut.bookHist(h,'Digi_SBTHits_YvsX', 'Y vs X digitalized SBT Hits', 400, -400., 400., 400, -400., 400.)
ut.bookHist(h,'Digi_SBTHits_Z', 'Z digitalized SBT Hits', 400, -3000., 3000.)

ut.bookHist(h,'digiSBT2MC', 'Digi SBT to MC', 400, 0., 4000.)

ut.bookHist(h,'UpstreamTaggerPoint_trackID', 'Upstream tagger point trackID', 400, 0, 400)
ut.bookHist(h,'UpstreamTaggerPoint_fY_versus_fX', 'Upstream tagger point Y versus X', 400, -1000, 1000, 400, 1000, 1000)
ut.bookHist(h,'UpstreamTaggerPoint_fZ', 'Upstream tagger point Z', 400, -1000, 1000)
ut.bookHist(h,'UpstreamTaggerPoint_fPy_versus_fPx', 'Upstream tagger point Py versus Px', 400, -10, 10, 400, -10, 10)
ut.bookHist(h,'UpstreamTaggerPoint_fPz', 'Upstream tagger point Pz', 400, 0, 400)
ut.bookHist(h,'UpstreamTaggerPoint_fPt', 'Upstream tagger point Pt', 400, 0, 100)
ut.bookHist(h,'UpstreamTaggerPoint_fTime', 'Upstream tagger point time', 400, 0, 1000)
ut.bookHist(h,'UpstreamTaggerPoint_fEloss', 'Upstream tagger point energy loss', 400, 0, 10)
ut.bookHist(h,'UpstreamTaggerPoint_fPdgCode', 'Upstream tagger point PDG code', 400, 0, 20)

ut.bookHist(h,'vetoPoint_trackID', 'Veto point trackID', 400, 0, 400)
ut.bookHist(h,'vetoPoint_fY_versus_fX', 'Veto point Y versus X', 400, -1000, 1000, 400, 1000, 1000)
ut.bookHist(h,'vetoPoint_fZ', 'Veto point Z', 400, -3000, 3000)
ut.bookHist(h,'vetoPoint_fPy_versus_fPx', 'Veto point Py versus Px', 400, -0.3, 0.3, 400, -0.3, 0.3)
ut.bookHist(h,'vetoPoint_fPz', 'Veto point Pz', 400, 0, 0.01)
ut.bookHist(h,'vetoPoint_fPt', 'Veto point Pt', 400, 0, 0.01)
ut.bookHist(h,'vetoPoint_fTime', 'Veto point time', 400, 0, 1000)
ut.bookHist(h,'vetoPoint_fEloss', 'Veto point energy loss', 400, 0, 0.01)
ut.bookHist(h,'vetoPoint_fPdgCode', 'Veto point PDG code', 400, 0, 20)

ut.bookHist(h,'strawtubesPoint_trackID', 'Strawtubes point trackID', 400, 0, 400)
ut.bookHist(h,'strawtubesPoint_fY_versus_fX', 'Strawtubes point Y versus X', 400, -200, 200, 400, -200, 200)
ut.bookHist(h,'strawtubesPoint_fZ', 'Strawtubes point Z', 400, -3000, 3000)
ut.bookHist(h,'strawtubesPoint_fPy_versus_fPx', 'Strawtubes point Py versus Px', 400, -0.06, 0.06, 400, -0.06, 0.06)
ut.bookHist(h,'strawtubesPoint_fPz', 'Strawtubes point Pz', 400, 0, 0.01)
ut.bookHist(h,'strawtubesPoint_fPt', 'Strawtubes point Pt', 400, 0, 0.01)
ut.bookHist(h,'strawtubesPoint_fTime', 'Strawtubes point time', 400, 0, 1000)
ut.bookHist(h,'strawtubesPoint_fEloss', 'Strawtubes point energy loss', 400, 0, 0.01)
ut.bookHist(h,'strawtubesPoint_fPdgCode', 'Strawtubes point PDG code', 400, 0, 20)

ut.bookHist(h,'MCTrack_MotherId', 'MC tracks mother ID', 400, 0, 400)
ut.bookHist(h,'MCTrack_fStartY_versus_fStartX', 'MC tracks start Y versus start X', 400, -1000, 1000, 400, 1000, 1000)
ut.bookHist(h,'MCTrack_fStartZ', 'MC tracks start Z', 400, -3000, 3000)
ut.bookHist(h,'MCTrack_fPy_versus_fPx', 'MC tracks Py versus Px', 400, -10, 10, 400, -10, 10)
ut.bookHist(h,'MCTrack_fPz', 'MC tracks Pz', 400, 0, 0.01)
ut.bookHist(h,'MCTrack_fPt', 'MC tracks Pt', 400, 0, 0.01)
ut.bookHist(h,'MCTrack_fM', 'MC tracks M', 400, 0, 20)
ut.bookHist(h,'MCTrack_fPdgCode', 'MC tracks PDG code', 400, 0, 20)

ut.bookHist(h,'UBT_Y_versus_X_muon_and_pion_hits', 'UBT Y versus X muon and pion hits; x [cm]; y [cm]; Muon and pion hits', 400, -1000, 1000, 400, 1000, 1000)
ut.bookHist(h,'UBT_Z_versus_X_muon_and_pion_hits', 'UBT Z versus X muon and pion hits; x [cm]; z [cm]; Muon and pion hits', 400, -1000, 1000, 400, 1000, 1000)
ut.bookHist(h,'UBT_Z_versus_Y_muon_and_pion_hits', 'UBT Z versus Y muon and pion hits; y [cm]; z [cm]; Muon and pion hits', 400, -1000, 1000, 400, 1000, 1000)

ut.bookHist(h,'UBT_Y_versus_X_electron_and_positron_hits', 'UBT Y versus X electron and positron hits; x [cm]; y [cm]; Electron and positron hits', 400, -1000, 1000, 400, 1000, 1000)
ut.bookHist(h,'UBT_Z_versus_X_electron_and_positron_hits', 'UBT Z versus X electron and positron hits; x [cm]; z [cm]; Electron and positron hits', 400, -1000, 1000, 400, 1000, 1000)
ut.bookHist(h,'UBT_Z_versus_Y_electron_and_positron_hits', 'UBT Z versus Y electron and positron hits; y [cm]; z [cm]; Electron and positron hits', 400, -1000, 1000, 400, 1000, 1000)
ut.bookHist(h,'UBT_electron_and_positron_momentum_versus_UBT_Z', 'UBT electron and positron momentum versus UBT Z; z [cm]; P [GeV]; Electron and positron hits', 400, -1000, 1000, 400, 0, 10)
ut.bookHist(h,'UBT_electron_and_positron_momentum', 'UBT electron and positron momentum; P [GeV]; Electron and positron hits', 400, 0, 10)

ut.bookHist(h,'UBT_Y_versus_X_photon_hits', 'UBT Y versus X photon hits; x [cm]; y [cm]; Photon hits', 400, -1000, 1000, 400, 1000, 1000)
ut.bookHist(h,'UBT_Z_versus_X_photon_hits', 'UBT Z versus X photon hits; x [cm]; z [cm]; Photon hits', 400, -1000, 1000, 400, 1000, 1000)
ut.bookHist(h,'UBT_Z_versus_Y_photon_hits', 'UBT Z versus Y photon hits; y [cm]; z [cm]; Photon hits', 400, -1000, 1000, 400, 1000, 1000)
ut.bookHist(h,'UBT_photon_energy_versus_UBT_Z', 'UBT photon energy versus UBT Z; z [cm]; Energy [GeV]; Photon hits', 400, -1000, 1000, 400, 0, 0.05)
ut.bookHist(h,'UBT_photon_energy', 'UBT photon energy; Energy [GeV]; Photon hits', 400, 0, 0.05)
ut.bookHist(h,'UBT_photon_momentum_versus_UBT_Z', 'UBT photon momentum versus UBT Z; z [cm]; P [GeV]; Photon hits', 400, -1000, 1000, 400, 0, 1000)
ut.bookHist(h,'UBT_photon_momentum', 'UBT photon momentum; P [GeV]; Photon hits', 400, 0, 10)

ut.bookHist(h,'SBT_Y_versus_X_photon_hits', 'SBT Y versus X photon hits; x [cm]; y [cm]; Photon hits', 400, -1000, 1000, 400, 1000, 1000)
ut.bookHist(h,'SBT_Z_versus_X_photon_hits', 'SBT Z versus X photon hits; x [cm]; z [cm]; Photon hits', 400, -1000, 1000, 400, 1000, 1000)
ut.bookHist(h,'SBT_Z_versus_Y_photon_hits', 'SBT Z versus Y photon hits; y [cm]; z [cm]; Photon hits', 400, -1000, 1000, 400, 1000, 1000)
ut.bookHist(h,'SBT_photon_energy_versus_SBT_Z', 'SBT photon energy versus SBT Z; z [cm]; Energy [GeV]; Photon hits', 400, -1000, 1000, 400, 0, 0.05)
ut.bookHist(h,'SBT_photon_energy', 'SBT photon energy; Energy [GeV]; Photon hits', 400, 0, 0.05)

ut.bookHist(h,'pT_muons_tracks_that_miss_SST', 'p_{T} muons tracks that miss SST; p_{T} [GeV]; Number of muon tracks', 400, 0, 0.1)

def makePlots():
  ut.bookCanvas(h, key='Digi_SBTHits_YvsX_canvas', title='Digi SBT Y vs X hit map', nx=800, ny=600, cx=1, cy=1)
  h['Digi_SBTHits_YvsX'].Draw('colz')
  ut.bookCanvas(h, key='Digi_SBTHits_Z_canvas', title='Digi SBT Z hit map', nx=800, ny=600, cx=1, cy=1)
  h['Digi_SBTHits_Z'].Draw()
  
  ut.bookCanvas(h, key='digiSBT2MC_canvas', title='Digi SBT to MC', nx=800, ny=600, cx=1, cy=1)
  h['digiSBT2MC'].Draw()

  ut.bookCanvas(h, key='UpstreamTaggerPoint_trackID_canvas', title='Upstream tagger point trackID', nx=800, ny=600, cx=1, cy=1)
  h['UpstreamTaggerPoint_trackID'].Draw()
  ut.bookCanvas(h, key='UpstreamTaggerPoint_fY_versus_fX_canvas', title='Upstream tagger point Y versus X', nx=800, ny=600, cx=1, cy=1)
  h['UpstreamTaggerPoint_fY_versus_fX'].Draw()
  ut.bookCanvas(h, key='UpstreamTaggerPoint_fZ_canvas', title='Upstream tagger point Z', nx=800, ny=600, cx=1, cy=1)
  h['UpstreamTaggerPoint_fZ'].Draw()
  ut.bookCanvas(h, key='UpstreamTaggerPoint_fPy_versus_fPx_canvas', title='Upstream tagger point Py versus Px', nx=800, ny=600, cx=1, cy=1)
  h['UpstreamTaggerPoint_fPy_versus_fPx'].Draw()
  ut.bookCanvas(h, key='UpstreamTaggerPoint_fPz_canvas', title='Upstream tagger point Pz', nx=800, ny=600, cx=1, cy=1)
  h['UpstreamTaggerPoint_fPz'].Draw()
  ut.bookCanvas(h, key='UpstreamTaggerPoint_fPt_canvas', title='Upstream tagger point Pt', nx=800, ny=600, cx=1, cy=1)
  h['UpstreamTaggerPoint_fPt'].Draw()
  ut.bookCanvas(h, key='UpstreamTaggerPoint_fTime_canvas', title='Upstream tagger point time', nx=800, ny=600, cx=1, cy=1)
  h['UpstreamTaggerPoint_fTime'].Draw()
  ut.bookCanvas(h, key='UpstreamTaggerPoint_fEloss_canvas', title='Upstream tagger point energy loss', nx=800, ny=600, cx=1, cy=1)
  h['UpstreamTaggerPoint_fEloss'].Draw()
  ut.bookCanvas(h, key='UpstreamTaggerPoint_fPdgCode_canvas', title='Upstream tagger point PDG code', nx=800, ny=600, cx=1, cy=1)
  h['UpstreamTaggerPoint_fPdgCode'].Draw()

  ut.bookCanvas(h, key='vetoPoint_trackID_canvas', title='Veto point trackID', nx=800, ny=600, cx=1, cy=1)
  h['vetoPoint_trackID'].Draw()
  ut.bookCanvas(h, key='vetoPoint_fY_versus_fX_canvas', title='Veto point Y versus X', nx=800, ny=600, cx=1, cy=1)
  h['vetoPoint_fY_versus_fX'].Draw()
  ut.bookCanvas(h, key='vetoPoint_fZ_canvas', title='Veto point Z', nx=800, ny=600, cx=1, cy=1)
  h['vetoPoint_fZ'].Draw()
  ut.bookCanvas(h, key='vetoPoint_fPy_versus_fPx_canvas', title='Veto point Py versus Px', nx=800, ny=600, cx=1, cy=1)
  h['vetoPoint_fPy_versus_fPx'].Draw()
  ut.bookCanvas(h, key='vetoPoint_fPz_canvas', title='Veto point Pz', nx=800, ny=600, cx=1, cy=1)
  h['vetoPoint_fPz'].Draw()
  ut.bookCanvas(h, key='vetoPoint_fPt_canvas', title='Veto point Pt', nx=800, ny=600, cx=1, cy=1)
  h['vetoPoint_fPt'].Draw()
  ut.bookCanvas(h, key='vetoPoint_fTime_canvas', title='Veto point time', nx=800, ny=600, cx=1, cy=1)
  h['vetoPoint_fTime'].Draw()
  ut.bookCanvas(h, key='vetoPoint_fEloss_canvas', title='Veto point energy loss', nx=800, ny=600, cx=1, cy=1)
  h['vetoPoint_fEloss'].Draw()
  ut.bookCanvas(h, key='vetoPoint_fPdgCode_canvas', title='Veto point PDG code', nx=800, ny=600, cx=1, cy=1)
  h['vetoPoint_fPdgCode'].Draw()
  
  ut.bookCanvas(h, key='strawtubesPoint_trackID_canvas', title='Strawtubes point trackID', nx=800, ny=600, cx=1, cy=1)
  h['strawtubesPoint_trackID'].Draw()
  ut.bookCanvas(h, key='strawtubesPoint_fY_versus_fX_canvas', title='Strawtubes point Y versus X', nx=800, ny=600, cx=1, cy=1)
  h['strawtubesPoint_fY_versus_fX'].Draw()
  ut.bookCanvas(h, key='strawtubesPoint_fZ_canvas', title='Strawtubes point Z', nx=800, ny=600, cx=1, cy=1)
  h['strawtubesPoint_fZ'].Draw()
  ut.bookCanvas(h, key='strawtubesPoint_fPy_versus_fPx_canvas', title='Strawtubes point Py versus Px', nx=800, ny=600, cx=1, cy=1)
  h['strawtubesPoint_fPy_versus_fPx'].Draw()
  ut.bookCanvas(h, key='strawtubesPoint_fPz_canvas', title='Strawtubes point Pz', nx=800, ny=600, cx=1, cy=1)
  h['strawtubesPoint_fPz'].Draw()
  ut.bookCanvas(h, key='strawtubesPoint_fPt_canvas', title='Strawtubes point Pt', nx=800, ny=600, cx=1, cy=1)
  h['strawtubesPoint_fPt'].Draw()
  ut.bookCanvas(h, key='strawtubesPoint_fTime_canvas', title='Strawtubes point time', nx=800, ny=600, cx=1, cy=1)
  h['strawtubesPoint_fTime'].Draw()
  ut.bookCanvas(h, key='strawtubesPoint_fEloss_canvas', title='Strawtubes point energy loss', nx=800, ny=600, cx=1, cy=1)
  h['strawtubesPoint_fEloss'].Draw()
  ut.bookCanvas(h, key='strawtubesPoint_fPdgCode_canvas', title='Strawtubes point PDG code', nx=800, ny=600, cx=1, cy=1)
  h['strawtubesPoint_fPdgCode'].Draw()

  ut.bookCanvas(h, key='MCTrack_MotherId_canvas', title='MC tracks Mother ID', nx=800, ny=600, cx=1, cy=1)
  h['MCTrack_MotherId'].Draw()
  ut.bookCanvas(h, key='MCTrack_fStartY_versus_fStartX_canvas', title='MC tracks start Y versus start X', nx=800, ny=600, cx=1, cy=1)
  h['MCTrack_fStartY_versus_fStartX'].Draw()
  ut.bookCanvas(h, key='MCTrack_fZ_canvas', title='MC tracks Z', nx=800, ny=600, cx=1, cy=1)
  h['MCTrack_fStartZ'].Draw()
  ut.bookCanvas(h, key='MCTrack_fPy_versus_fPx_canvas', title='MC tracks Py versus Px', nx=800, ny=600, cx=1, cy=1)
  h['MCTrack_fPy_versus_fPx'].Draw()
  ut.bookCanvas(h, key='MCTrack_fPz_canvas', title='MC tracks Pz', nx=800, ny=600, cx=1, cy=1)
  h['MCTrack_fPz'].Draw()
  ut.bookCanvas(h, key='MCTrack_fPt_canvas', title='MC tracks Pt', nx=800, ny=600, cx=1, cy=1)
  h['MCTrack_fPt'].Draw()
  ut.bookCanvas(h, key='MCTrack_fM_canvas', title='MC tracks M', nx=800, ny=600, cx=1, cy=1)
  h['MCTrack_fM'].Draw()
  ut.bookCanvas(h, key='MCTrack_fPdgCode_canvas', title='MC tracks PDG code', nx=800, ny=600, cx=1, cy=1)
  h['MCTrack_fPdgCode'].Draw()

  ut.bookCanvas(h, key='UBT_Y_versus_X_muon_and_pion_hits', title='UBT Y versus X muon and pion hits', nx=800, ny=600, cx=1, cy=1)
  h['UBT_Y_versus_X_muon_and_pion_hits'].Draw()
  ut.bookCanvas(h, key='UBT_Z_versus_X_muon_and_pion_hits', title='UBT Z versus X muon and pion hits', nx=800, ny=600, cx=1, cy=1)
  h['UBT_Y_versus_X_muon_and_pion_hits'].Draw()
  ut.bookCanvas(h, key='UBT_Z_versus_Y_muon_and_pion_hits', title='UBT Z versus Y muon and pion hits', nx=800, ny=600, cx=1, cy=1)
  h['UBT_Y_versus_X_muon_and_pion_hits'].Draw()

  ut.bookCanvas(h, key='UBT_Y_versus_X_electron_and_positron_hits', title='UBT Y versus X electron and positron hits', nx=800, ny=600, cx=1, cy=1)
  h['UBT_Y_versus_X_electron_and_positron_hits'].Draw()
  ut.bookCanvas(h, key='UBT_Z_versus_X_electron_and_positron_hits', title='UBT Z versus X electron and positron hits', nx=800, ny=600, cx=1, cy=1)
  h['UBT_Y_versus_X_electron_and_positron_hits'].Draw()
  ut.bookCanvas(h, key='UBT_Z_versus_Y_electron_and_positron_hits', title='UBT Z versus Y electron and positron hits', nx=800, ny=600, cx=1, cy=1)
  h['UBT_Y_versus_X_electron_and_positron_hits'].Draw()
  ut.bookCanvas(h, key='UBT_electron_and_positron_momentum_versus_UBT_Z', title='UBT electron and positron momentum versus UBT Z', nx=800, ny=600, cx=1, cy=1)
  h['UBT_electron_and_positron_momentum_versus_UBT_Z'].Draw()
  ut.bookCanvas(h, key='UBT_electron_and_positron_momentum', title='UBT electron and positron momentum', nx=800, ny=600, cx=1, cy=1)
  h['UBT_electron_and_positron_momentum'].Draw()
  
  ut.bookCanvas(h, key='UBT_Y_versus_X_photon_hits', title='UBT Y versus X photon hits', nx=800, ny=600, cx=1, cy=1)
  h['UBT_Y_versus_X_photon_hits'].Draw()
  ut.bookCanvas(h, key='UBT_Z_versus_X_photon_hits', title='UBT Z versus X photon hits', nx=800, ny=600, cx=1, cy=1)
  h['UBT_Y_versus_X_photon_hits'].Draw()
  ut.bookCanvas(h, key='UBT_Z_versus_Y_photon_hits', title='UBT Z versus Y photon hits', nx=800, ny=600, cx=1, cy=1)
  h['UBT_Y_versus_X_photon_hits'].Draw()
  ut.bookCanvas(h, key='UBT_photon_energy', title='UBT photon energy', nx=800, ny=600, cx=1, cy=1)
  h['UBT_photon_energy'].Draw()
  
  ut.bookCanvas(h, key='SBT_Y_versus_X_photon_hits', title='SBT Y versus X photon hits', nx=800, ny=600, cx=1, cy=1)
  h['SBT_Y_versus_X_photon_hits'].Draw()
  ut.bookCanvas(h, key='SBT_Z_versus_X_photon_hits', title='SBT Z versus X photon hits', nx=800, ny=600, cx=1, cy=1)
  h['SBT_Y_versus_X_photon_hits'].Draw()
  ut.bookCanvas(h, key='SBT_Z_versus_Y_photon_hits', title='SBT Z versus Y photon hits', nx=800, ny=600, cx=1, cy=1)
  h['SBT_Y_versus_X_photon_hits'].Draw()
  ut.bookCanvas(h, key='SBT_photon_energy_versus_SBT_Z', title='SBT photon energy versus SBT Z', nx=800, ny=600, cx=1, cy=1)
  h['SBT_photon_energy_versus_SBT_Z'].Draw()
  ut.bookCanvas(h, key='SBT_photon_energy', title='SBT photon energy', nx=800, ny=600, cx=1, cy=1)
  h['SBT_photon_energy'].Draw()

  ut.bookCanvas(h, key='PT_muons_miss_SST', title='p_{T} muons tracks that miss SST', nx=800, ny=600, cx=1, cy=1)
  h['PT_muons_miss_SST'].Draw()
  
  print('finished making plots')

import TrackExtrapolateTool

def myVertex(t1,t2,PosDir):
  # closest distance between two tracks
  # d = |pq . u x v|/|u x v|
  a = ROOT.TVector3(PosDir[t1][0](0), PosDir[t1][0](1), PosDir[t1][0](2))
  u = ROOT.TVector3(PosDir[t1][1](0), PosDir[t1][1](1), PosDir[t1][1](2))
  c = ROOT.TVector3(PosDir[t2][0](0), PosDir[t2][0](1), PosDir[t2][0](2))
  v = ROOT.TVector3(PosDir[t2][1](0), PosDir[t2][1](1), PosDir[t2][1](2))
  pq = a - c
  uCrossv = u.Cross(v)
  dist  = pq.Dot(uCrossv)/(uCrossv.Mag()+1E-8)
  # u.a - u.c + s*|u|**2 - u.v*t = 0
  # v.a - v.c + s*v.u - t*|v|**2 = 0
  E = u.Dot(a) - u.Dot(c)
  F = v.Dot(a) - v.Dot(c)
  A, B = u.Mag2(), -u.Dot(v)
  C, D = u.Dot(v), -v.Mag2()
  t = -(C*E-A*F)/(B*C-A*D)
  X = c.x()+v.x()*t
  Y = c.y()+v.y()*t
  Z = c.z()+v.z()*t

  return X, Y, Z, abs(dist)

# Position resolution: Trace back muon from SST to UBT (w/o) SBT requirement
def position_resolution_muon_from_SST_to_UBT():
  # Check how many muons miss SST
  track = -1
  PosDir = {}
  for i, mc in enumerate(sTree.MCTrack):
    if abs(mc.GetPdgCode()) != 13: continue
    p_start = ROOT.TVector3(mc.GetPx(), mc.GetPy(), mc.GetPz())
    pt = p_start.Perp()
    reached_ubt = any(hit.GetTrackID() == i for hit in sTree.UpstreamTaggerPoint)
    reached_sst = any(hit.GetTrackID() == i for hit in sTree.strawtubesPoint)
    if not reached_sst:
      h['pT_muons_tracks_that_miss_SST'].Fill(pt)

    if reached_ubt and reached_sst:
      track += 1
      if track >= len(sTree.FitTracks): continue
      FitTracksObj = sTree.FitTracks[track]
      fitStatus = FitTracksObj.getFitStatus()
      if not fitStatus.isFitConverged(): continue
      rep = FitTracksObj.getFittedState()
      PosDir[track] = [rep.getPos(), rep.getDir()]
      
  if len(PosDir) > 1:
    x, y, z, doca = myVertex(0, 1, PosDir)
    if doca < 2:
      print(x, y, z, doca)

      #print("vertex is at  < 5cm from vessel wall")
        
      #ip_pos = rep.getPos()
      #ip = ip_pos.Perp()
      #if ip > 250.0: continue
      
      #if hasattr(ShipGeo, "UpstreamTagger"):
      #rc, pos, mom = TrackExtrapolateTool.extrapolateToPlane(FitTracksObj, ShipGeo.UpstreamTagger.z)
      #print(FitTracksObj, rc, pos, mom)
      #if pos is None and mom is None: continue
      #print(pos.X(), vtx.X())
      
# Hit Map of muon crossing UBT (extended) with 1s spill, to obtain Max rate capability. 
# What is the e-/e+ flux at UBT Z position and what is their momentum.
# what is the P of the photons.
def UBT_hit_maps_and_particle_flux():
  for UpstreamTaggerPoint_it in sTree.UpstreamTaggerPoint:
    trackID = UpstreamTaggerPoint_it.GetTrackID()
    if trackID < 0 or trackID >= sTree.MCTrack.GetEntries(): continue
    mc = sTree.MCTrack[trackID]
    px, py, pz = mc.GetPx(), mc.GetPy(), mc.GetPz()
    mom = ROOT.TMath.Sqrt(px**2 + py**2 + pz**2)
    if abs(mc.GetPdgCode()) == 13 or abs(mc.GetPdgCode()) == 221:
      h['UBT_Y_versus_X_muon_and_pion_hits'].Fill(UpstreamTaggerPoint_it.GetX(), UpstreamTaggerPoint_it.GetY())
      h['UBT_Z_versus_X_muon_and_pion_hits'].Fill(UpstreamTaggerPoint_it.GetX(), UpstreamTaggerPoint_it.GetZ())
      h['UBT_Z_versus_Y_muon_and_pion_hits'].Fill(UpstreamTaggerPoint_it.GetY(), UpstreamTaggerPoint_it.GetZ())
    elif abs(mc.GetPdgCode()) == 11:
      h['UBT_Y_versus_X_electron_and_positron_hits'].Fill(UpstreamTaggerPoint_it.GetX(), UpstreamTaggerPoint_it.GetY())
      h['UBT_Z_versus_X_electron_and_positron_hits'].Fill(UpstreamTaggerPoint_it.GetX(), UpstreamTaggerPoint_it.GetZ())
      h['UBT_Z_versus_Y_electron_and_positron_hits'].Fill(UpstreamTaggerPoint_it.GetY(), UpstreamTaggerPoint_it.GetZ())
      h['UBT_electron_and_positron_momentum_versus_UBT_Z'].Fill(UpstreamTaggerPoint_it.GetZ(), mom)
      h['UBT_electron_and_positron_momentum'].Fill(mom)
    elif mc.GetPdgCode() == 22:
      h['UBT_photon_momentum_versus_UBT_Z'].Fill(UpstreamTaggerPoint_it.GetZ(), mom)
      h['UBT_photon_momentum'].Fill(mom)

# EM Debris: How many photons cross UBT (UpstreamTaggerPoint) and SBT (vetoPoint)
# Check E_gamma spectra at UBT z position
def EM_debris_UBT_SBT():
  for UpstreamTaggerPoint_it in sTree.UpstreamTaggerPoint:
    trackID = UpstreamTaggerPoint_it.GetTrackID()
    if trackID < 0 or trackID >= sTree.MCTrack.GetEntries(): continue
    mc = sTree.MCTrack[trackID]
    if mc.GetPdgCode() == 22:
      h['UBT_Y_versus_X_photon_hits'].Fill(UpstreamTaggerPoint_it.GetX(), UpstreamTaggerPoint_it.GetY())
      h['UBT_Z_versus_X_photon_hits'].Fill(UpstreamTaggerPoint_it.GetX(), UpstreamTaggerPoint_it.GetZ())
      h['UBT_Z_versus_Y_photon_hits'].Fill(UpstreamTaggerPoint_it.GetY(), UpstreamTaggerPoint_it.GetZ())
      h['UBT_photon_energy_versus_SBT_Z'].Fill(UpstreamTaggerPoint_it.GetZ(), mc.GetEnergy())
      h['UBT_photon_energy'].Fill(mc.GetEnergy())
    
  for vetoPoint_it in sTree.vetoPoint:
    trackID = vetoPoint_it.GetTrackID()
    if trackID < 0 or trackID >= sTree.MCTrack.GetEntries(): continue
    mc = sTree.MCTrack[trackID]
    if mc.GetPdgCode() == 22:
      h['SBT_Y_versus_X_photon_hits'].Fill(vetoPoint_it.GetX(), vetoPoint_it.GetY())
      h['SBT_Z_versus_X_photon_hits'].Fill(vetoPoint_it.GetX(), vetoPoint_it.GetZ())
      h['SBT_Z_versus_Y_photon_hits'].Fill(vetoPoint_it.GetY(), vetoPoint_it.GetZ())
      h['SBT_photon_energy_versus_SBT_Z'].Fill(vetoPoint_it.GetZ(), mc.GetEnergy())
      h['SBT_photon_energy'].Fill(mc.GetEnergy())
      
# Fill digi SBT Hits plots
def fill_digi_SBT_Hits_plots():
  for Digi_SBTHits_it in sTree.Digi_SBTHits:
    X_Digi_SBTHits = Digi_SBTHits_it.GetXYZ().X()
    Y_Digi_SBTHits = Digi_SBTHits_it.GetXYZ().Y()
    Z_Digi_SBTHits = Digi_SBTHits_it.GetXYZ().Z()

    h['Digi_SBTHits_YvsX'].Fill(X_Digi_SBTHits, Y_Digi_SBTHits)
    h['Digi_SBTHits_Z'].Fill(Z_Digi_SBTHits)

  for digiSBT2MC_it in sTree.digiSBT2MC:
    for digiSBT2MC_it_it in digiSBT2MC_it:
      h['digiSBT2MC'].Fill(digiSBT2MC_it_it)

# Fill UpstreamTaggerPoint (UBT) plots
def fill_UpstreamTaggerPoint_plots():
    for UpstreamTaggerPoint_it in sTree.UpstreamTaggerPoint:
      UpstreamTaggerPoint_trackID = UpstreamTaggerPoint_it.GetTrackID()
      UpstreamTaggerPoint_fX = UpstreamTaggerPoint_it.GetX()
      UpstreamTaggerPoint_fY = UpstreamTaggerPoint_it.GetY()
      UpstreamTaggerPoint_fZ = UpstreamTaggerPoint_it.GetZ()
      UpstreamTaggerPoint_fPx = UpstreamTaggerPoint_it.GetPx()
      UpstreamTaggerPoint_fPy = UpstreamTaggerPoint_it.GetPy()
      UpstreamTaggerPoint_fPz = UpstreamTaggerPoint_it.GetPz()
      UpstreamTaggerPoint_fPt = ROOT.TMath.Sqrt(UpstreamTaggerPoint_fPx**2 + UpstreamTaggerPoint_fPy**2)
      UpstreamTaggerPoint_fTime = UpstreamTaggerPoint_it.GetTime()
      UpstreamTaggerPoint_fEloss = UpstreamTaggerPoint_it.GetEnergyLoss()
      UpstreamTaggerPoint_fPdgCode = UpstreamTaggerPoint_it.PdgCode()
      h['UpstreamTaggerPoint_trackID'].Fill(UpstreamTaggerPoint_trackID)
      h['UpstreamTaggerPoint_fY_versus_fX'].Fill(UpstreamTaggerPoint_fX, UpstreamTaggerPoint_fY)
      h['UpstreamTaggerPoint_fZ'].Fill(UpstreamTaggerPoint_fZ)
      h['UpstreamTaggerPoint_fPy_versus_fPx'].Fill(UpstreamTaggerPoint_fPx, UpstreamTaggerPoint_fPy)
      h['UpstreamTaggerPoint_fPz'].Fill(UpstreamTaggerPoint_fPz)
      h['UpstreamTaggerPoint_fPt'].Fill(UpstreamTaggerPoint_fPt)
      h['UpstreamTaggerPoint_fTime'].Fill(UpstreamTaggerPoint_fTime)
      h['UpstreamTaggerPoint_fEloss'].Fill(UpstreamTaggerPoint_fEloss)
      h['UpstreamTaggerPoint_fPdgCode'].Fill(UpstreamTaggerPoint_fPdgCode)

# Fill vetoPoint (SBT) plots
def fill_vetoPoint_plots():
  for vetoPoint_it in sTree.vetoPoint:
    vetoPoint_trackID = vetoPoint_it.GetTrackID()
    vetoPoint_fX = vetoPoint_it.GetX()
    vetoPoint_fY = vetoPoint_it.GetY()
    vetoPoint_fZ = vetoPoint_it.GetZ()
    vetoPoint_fPx = vetoPoint_it.GetPx()
    vetoPoint_fPy = vetoPoint_it.GetPy()
    vetoPoint_fPz = vetoPoint_it.GetPz()
    vetoPoint_fPt = ROOT.TMath.Sqrt(vetoPoint_fPx**2 + vetoPoint_fPy**2)
    vetoPoint_fTime = vetoPoint_it.GetTime()
    vetoPoint_fEloss = vetoPoint_it.GetEnergyLoss()
    vetoPoint_fPdgCode = vetoPoint_it.PdgCode()
    h['vetoPoint_trackID'].Fill(vetoPoint_trackID)
    h['vetoPoint_fY_versus_fX'].Fill(vetoPoint_fX, vetoPoint_fY)
    h['vetoPoint_fZ'].Fill(vetoPoint_fZ)
    h['vetoPoint_fPy_versus_fPx'].Fill(vetoPoint_fPx, vetoPoint_fPy)
    h['vetoPoint_fPz'].Fill(vetoPoint_fPz)
    h['vetoPoint_fPt'].Fill(vetoPoint_fPt)
    h['vetoPoint_fTime'].Fill(vetoPoint_fTime)
    h['vetoPoint_fEloss'].Fill(vetoPoint_fEloss)
    h['vetoPoint_fPdgCode'].Fill(vetoPoint_fPdgCode)

# Fill strawtubesPoint (SST) plots
def fill_strawtubesPoint_plots():
    for strawtubesPoint_it in sTree.strawtubesPoint:
      strawtubesPoint_trackID = strawtubesPoint_it.GetTrackID()
      strawtubesPoint_fX = strawtubesPoint_it.GetX()
      strawtubesPoint_fY = strawtubesPoint_it.GetY()
      strawtubesPoint_fZ = strawtubesPoint_it.GetZ()
      strawtubesPoint_fPx = strawtubesPoint_it.GetPx()
      strawtubesPoint_fPy = strawtubesPoint_it.GetPy()
      strawtubesPoint_fPz = strawtubesPoint_it.GetPz()
      strawtubesPoint_fPt = ROOT.TMath.Sqrt(strawtubesPoint_fPx**2 + strawtubesPoint_fPy**2)
      strawtubesPoint_fTime = strawtubesPoint_it.GetTime()
      strawtubesPoint_fEloss = strawtubesPoint_it.GetEnergyLoss()
      strawtubesPoint_fPdgCode = strawtubesPoint_it.PdgCode()
      h['strawtubesPoint_trackID'].Fill(strawtubesPoint_trackID)
      h['strawtubesPoint_fY_versus_fX'].Fill(strawtubesPoint_fX, strawtubesPoint_fY)
      h['strawtubesPoint_fZ'].Fill(strawtubesPoint_fZ)
      h['strawtubesPoint_fPy_versus_fPx'].Fill(strawtubesPoint_fPx, strawtubesPoint_fPy)
      h['strawtubesPoint_fPz'].Fill(strawtubesPoint_fPz)
      h['strawtubesPoint_fPt'].Fill(strawtubesPoint_fPt)
      h['strawtubesPoint_fTime'].Fill(strawtubesPoint_fTime)
      h['strawtubesPoint_fEloss'].Fill(strawtubesPoint_fEloss)
      h['strawtubesPoint_fPdgCode'].Fill(strawtubesPoint_fPdgCode)

# Fill MCtrack plots
def fill_MCtrack_plots():
  for MCTrack_it in sTree.MCTrack:
    MCTrack_MotherId = MCTrack_it.GetMotherId()
    MCTrack_fStartX = MCTrack_it.GetStartX()
    MCTrack_fStartY = MCTrack_it.GetStartY()
    MCTrack_fStartZ = MCTrack_it.GetStartZ()
    MCTrack_fPx = MCTrack_it.GetPx()
    MCTrack_fPy = MCTrack_it.GetPy()
    MCTrack_fPz = MCTrack_it.GetPz()
    MCTrack_fPt = ROOT.TMath.Sqrt(MCTrack_fPx**2 + MCTrack_fPy**2)
    MCTrack_fM = MCTrack_it.GetMass()
    MCTrack_fPdgCode = MCTrack_it.GetPdgCode()
    h['MCTrack_MotherId'].Fill(MCTrack_MotherId)
    h['MCTrack_fStartY_versus_fStartX'].Fill(MCTrack_fStartX, MCTrack_fStartY)
    h['MCTrack_fStartZ'].Fill(MCTrack_fStartZ)
    h['MCTrack_fPy_versus_fPx'].Fill(MCTrack_fPx, MCTrack_fPy)
    h['MCTrack_fPz'].Fill(MCTrack_fPz)
    h['MCTrack_fPt'].Fill(MCTrack_fPt)
    h['MCTrack_fM'].Fill(MCTrack_fM)
    h['MCTrack_fPdgCode'].Fill(MCTrack_fPdgCode)
    
# start event loop
def myEventLoop(n):
  rc = sTree.GetEntry(n)
  sTree.GetEntry(n)

  fill_digi_SBT_Hits_plots()
  fill_UpstreamTaggerPoint_plots()
  fill_vetoPoint_plots()
  fill_strawtubesPoint_plots()
  fill_MCtrack_plots()
  UBT_hit_maps_and_particle_flux()
  EM_debris_UBT_SBT()
  position_resolution_muon_from_SST_to_UBT()
  
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
