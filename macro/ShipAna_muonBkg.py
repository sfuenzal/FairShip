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

chi2CutOff  = 4.
fiducialCut = False
measCutFK = 25
measCutPR = 22
docaCut = 2.

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
import shipDet_conf
run = ROOT.FairRunSim()
run.SetName("TGeant4")  # Transport engine
run.SetSink(ROOT.FairRootFileSink(ROOT.TMemFile('output', 'recreate')))  # Dummy output file
run.SetUserConfig("g4Config_basic.C") # geant4 transport not used, only needed for the mag field
rtdb = run.GetRuntimeDb()
# -----Create geometry----------------------------------------------
modules = shipDet_conf.configure(run,ShipGeo)

import geomGeant4
if hasattr(ShipGeo.Bfield,"fieldMap"):
  fieldMaker = geomGeant4.addVMCFields(ShipGeo, '', True, withVirtualMC = False)
else:
  print("no fieldmap given, geofile too old, not anymore support")
  exit(-1)
sGeo = fgeo.Get("FAIRGeom")
geoMat =  ROOT.genfit.TGeoMaterialInterface()
ROOT.genfit.MaterialEffects.getInstance().init(geoMat)
bfield = ROOT.genfit.FairShipFields()
bfield.setField(fieldMaker.getGlobalField())
fM = ROOT.genfit.FieldManager.getInstance()
fM.init(bfield)

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

ut.bookHist(h,'vetoPoint_fPx', 'Veto point Px', 400, -0.1, 0.1)
ut.bookHist(h,'vetoPoint_fPy', 'Veto point Py', 400, -0.1, 0.1)
ut.bookHist(h,'vetoPoint_fPz', 'Veto point Pz', 400, -0.1, 0.1)
ut.bookHist(h,'vetoPoint_fPt', 'Veto point Pt', 400, 0., 0.1)
ut.bookHist(h,'vetoPoint_fX', 'Veto point X', 400, -3000., 3000.)
ut.bookHist(h,'vetoPoint_fY', 'Veto point Y', 400, -3000., 3000.)
ut.bookHist(h,'vetoPoint_fZ', 'Veto point Z', 400, -3000., 3000.)

ut.bookHist(h,'strawtubesPoint_fPx', 'Straw tube point Px', 400, -0.1, 0.1)
ut.bookHist(h,'strawtubesPoint_fPy', 'Straw tube point Py', 400, -0.1, 0.1)
ut.bookHist(h,'strawtubesPoint_fPz', 'Straw tube point Pz', 400, -0.1, 0.1)
ut.bookHist(h,'strawtubesPoint_fPt', 'Straw tube point Pt', 400, 0., 0.1)
ut.bookHist(h,'strawtubesPoint_fX', 'Straw tube point X', 400, -3000., 3000.)
ut.bookHist(h,'strawtubesPoint_fY', 'Straw tube point Y', 400, -3000., 3000.)
ut.bookHist(h,'strawtubesPoint_fZ', 'Straw tube point Z', 400, -3000., 3000.)

ut.bookHist(h,'MCTrack_fPx', 'MC track Px', 400, -0.1, 0.1)
ut.bookHist(h,'MCTrack_fPy', 'MC track Py', 400, -0.1, 0.1)
ut.bookHist(h,'MCTrack_fPz', 'MC track Pz', 400, -0.1, 0.1)
ut.bookHist(h,'MCTrack_fPt', 'MC track Pt', 400, 0., 0.1)
ut.bookHist(h,'MCTrack_fStartX', 'MC track start X', 400, -3000., 3000.)
ut.bookHist(h,'MCTrack_fStartY', 'MC track start Y', 400, -3000., 3000.)
ut.bookHist(h,'MCTrack_fStartZ', 'MC track start Z', 400, -3000., 3000.)
ut.bookHist(h,'MCTrack_fM', 'MC track M', 400, 0., 10.)

def makePlots():
  ut.bookCanvas(h, key='Digi_SBTHits_YvsX_canvas', title='Digi SBT Y vs X hit map', nx=800, ny=600, cx=1, cy=1)
  h['Digi_SBTHits_YvsX'].Draw('colz')
  ut.bookCanvas(h, key='Digi_SBTHits_Z_canvas', title='Digi SBT Z hit map', nx=800, ny=600, cx=1, cy=1)
  h['Digi_SBTHits_Z'].Draw()
  
  ut.bookCanvas(h, key='digiSBT2MC_canvas', title='Digi SBT to MC', nx=800, ny=600, cx=1, cy=1)
  h['digiSBT2MC'].Draw()
  
  ut.bookCanvas(h, key='vetoPoint_fPx_canvas', title='Veto point Px', nx=800, ny=600, cx=1, cy=1)
  h['vetoPoint_fPx'].Draw()
  ut.bookCanvas(h, key='vetoPoint_fPy_canvas', title='Veto point Py', nx=800, ny=600, cx=1, cy=1)
  h['vetoPoint_fPy'].Draw()
  ut.bookCanvas(h, key='vetoPoint_fPz_canvas', title='Veto point Pz', nx=800, ny=600, cx=1, cy=1)
  h['vetoPoint_fPz'].Draw()
  ut.bookCanvas(h, key='vetoPoint_fPt_canvas', title='Veto point Pt', nx=800, ny=600, cx=1, cy=1)
  h['vetoPoint_fPt'].Draw()
  ut.bookCanvas(h, key='vetoPoint_fX_canvas', title='Veto point X', nx=800, ny=600, cx=1, cy=1)
  h['vetoPoint_fX'].Draw()
  ut.bookCanvas(h, key='vetoPoint_fY_canvas', title='Veto point Y', nx=800, ny=600, cx=1, cy=1)
  h['vetoPoint_fY'].Draw()
  ut.bookCanvas(h, key='vetoPoint_fZ_canvas', title='Veto point Z', nx=800, ny=600, cx=1, cy=1)
  h['vetoPoint_fZ'].Draw()
  
  ut.bookCanvas(h, key='strawtubesPoint_fPx_canvas', title='Straw tube point Px', nx=800, ny=600, cx=1, cy=1)
  h['strawtubesPoint_fPx'].Draw()
  ut.bookCanvas(h, key='strawtubesPoint_fPy_canvas', title='Straw tube point Py', nx=800, ny=600, cx=1, cy=1)
  h['strawtubesPoint_fPy'].Draw()
  ut.bookCanvas(h, key='strawtubesPoint_fPz_canvas', title='Straw tube point Pz', nx=800, ny=600, cx=1, cy=1)
  h['strawtubesPoint_fPz'].Draw()
  ut.bookCanvas(h, key='strawtubesPoint_fPt_canvas', title='Straw tube point Pt', nx=800, ny=600, cx=1, cy=1)
  h['strawtubesPoint_fPt'].Draw()
  ut.bookCanvas(h, key='strawtubesPoint_fX_canvas', title='Straw tube point X', nx=800, ny=600, cx=1, cy=1)
  h['strawtubesPoint_fX'].Draw()
  ut.bookCanvas(h, key='strawtubesPoint_fY_canvas', title='Straw tube point Y', nx=800, ny=600, cx=1, cy=1)
  h['strawtubesPoint_fY'].Draw()
  ut.bookCanvas(h, key='strawtubesPoint_fZ_canvas', title='Straw tube point Z', nx=800, ny=600, cx=1, cy=1)
  h['strawtubesPoint_fZ'].Draw()

  ut.bookCanvas(h, key='MCTrack_fPx_canvas', title='MC track Px', nx=800, ny=600, cx=1, cy=1)
  h['MCTrack_fPx'].Draw()
  ut.bookCanvas(h, key='MCTrack_fPy_canvas', title='MC track Py', nx=800, ny=600, cx=1, cy=1)
  h['MCTrack_fPy'].Draw()
  ut.bookCanvas(h, key='MCTrack_fPz_canvas', title='MC track Pz', nx=800, ny=600, cx=1, cy=1)
  h['MCTrack_fPz'].Draw()
  ut.bookCanvas(h, key='MCTrack_fPt_canvas', title='MC track Pt', nx=800, ny=600, cx=1, cy=1)
  h['MCTrack_fPt'].Draw()
  ut.bookCanvas(h, key='MCTrack_fStartX_canvas', title='MC track start X', nx=800, ny=600, cx=1, cy=1)
  h['MCTrack_fStartX'].Draw()
  ut.bookCanvas(h, key='MCTrack_fStartY_canvas', title='MC track start Y', nx=800, ny=600, cx=1, cy=1)
  h['MCTrack_fStartY'].Draw()
  ut.bookCanvas(h, key='MCTrack_fStartZ_canvas', title='MC track start Z', nx=800, ny=600, cx=1, cy=1)
  h['MCTrack_fStartZ'].Draw()
  ut.bookCanvas(h, key='MCTrack_fM_canvas', title='MC track M', nx=800, ny=600, cx=1, cy=1)
  h['MCTrack_fStartZ'].Draw()
  
  print('finished making plots')

# start event loop
def myEventLoop(n):
  rc = sTree.GetEntry(n)
  sTree.GetEntry(n)
  for Digi_SBTHits_it in sTree.Digi_SBTHits:
    X_Digi_SBTHits = Digi_SBTHits_it.GetXYZ().X()
    Y_Digi_SBTHits = Digi_SBTHits_it.GetXYZ().Y()
    Z_Digi_SBTHits = Digi_SBTHits_it.GetXYZ().Z()

    h['Digi_SBTHits_YvsX'].Fill(X_Digi_SBTHits, Y_Digi_SBTHits)
    h['Digi_SBTHits_Z'].Fill(Z_Digi_SBTHits)

  for digiSBT2MC_it in sTree.digiSBT2MC:
    for digiSBT2MC_it_it in digiSBT2MC_it:
      h['digiSBT2MC'].Fill(digiSBT2MC_it_it)

  for vetoPoint_it in sTree.vetoPoint:
    vetoPoint_fPx = vetoPoint_it.GetPx()
    vetoPoint_fPy = vetoPoint_it.GetPy()
    vetoPoint_fPz = vetoPoint_it.GetPz()
    vetoPoint_fPt = ROOT.TMath.Sqrt(vetoPoint_fPx**2 + vetoPoint_fPy**2)
    vetoPoint_fX = vetoPoint_it.GetX()
    vetoPoint_fY = vetoPoint_it.GetY()
    vetoPoint_fZ = vetoPoint_it.GetZ()
    h['vetoPoint_fPx'].Fill(vetoPoint_fPx)
    h['vetoPoint_fPy'].Fill(vetoPoint_fPy)
    h['vetoPoint_fPz'].Fill(vetoPoint_fPz)
    h['vetoPoint_fPt'].Fill(vetoPoint_fPt)
    h['vetoPoint_fX'].Fill(vetoPoint_fX)
    h['vetoPoint_fY'].Fill(vetoPoint_fY)
    h['vetoPoint_fZ'].Fill(vetoPoint_fZ)

  for strawtubesPoint_it in sTree.strawtubesPoint:
    strawtubesPoint_fPx = strawtubesPoint_it.GetPx()
    strawtubesPoint_fPy = strawtubesPoint_it.GetPy()
    strawtubesPoint_fPz = strawtubesPoint_it.GetPz()
    strawtubesPoint_fPt = ROOT.TMath.Sqrt(strawtubesPoint_fPx**2 + strawtubesPoint_fPy**2)
    strawtubesPoint_fX = strawtubesPoint_it.GetX()
    strawtubesPoint_fY = strawtubesPoint_it.GetY()
    strawtubesPoint_fZ = strawtubesPoint_it.GetZ()
    h['strawtubesPoint_fPx'].Fill(strawtubesPoint_fPx)
    h['strawtubesPoint_fPy'].Fill(strawtubesPoint_fPy)
    h['strawtubesPoint_fPz'].Fill(strawtubesPoint_fPz)
    h['strawtubesPoint_fPt'].Fill(strawtubesPoint_fPt)
    h['strawtubesPoint_fX'].Fill(strawtubesPoint_fX)
    h['strawtubesPoint_fY'].Fill(strawtubesPoint_fY)
    h['strawtubesPoint_fZ'].Fill(strawtubesPoint_fZ)

  for MCTrack_it in sTree.MCTrack:
    MCTrack_fPx = MCTrack_it.GetPx()
    MCTrack_fPy = MCTrack_it.GetPy()
    MCTrack_fPz = MCTrack_it.GetPz()
    MCTrack_fPt = ROOT.TMath.Sqrt(MCTrack_fPx**2 + MCTrack_fPy**2)
    MCTrack_fStartX = MCTrack_it.GetStartX()
    MCTrack_fStartY = MCTrack_it.GetStartY()
    MCTrack_fStartZ = MCTrack_it.GetStartZ()
    MCTrack_fM = MCTrack_it.GetMass()
    h['MCTrack_fPx'].Fill(MCTrack_fPx)
    h['MCTrack_fPy'].Fill(MCTrack_fPy)
    h['MCTrack_fPz'].Fill(MCTrack_fPz)
    h['MCTrack_fPt'].Fill(MCTrack_fPt)
    h['MCTrack_fStartX'].Fill(MCTrack_fStartX)
    h['MCTrack_fStartY'].Fill(MCTrack_fStartY)
    h['MCTrack_fStartZ'].Fill(MCTrack_fStartZ)
    h['MCTrack_fM'].Fill(MCTrack_fM)
    
sTree.GetEvent(0)
options.nEvents = min(sTree.GetEntries(),options.nEvents)

# import pi0Reco
#ut.bookHist(h,'pi0Mass','gamma gamma inv mass',100,0.,0.5)

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
