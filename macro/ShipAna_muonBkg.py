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

ut.bookHist(h,'Digi_SBTHits_YvsX', 'Y vs X SBT Hits', 400, -400., 400., 400, -400., 400.)
ut.bookHist(h,'Digi_SBTHits_Z', 'Z SBT Hits', 400, -3000., 3000)

def makePlots():
  ut.bookCanvas(h, key='Digi_SBTHits_YvsX_canvas', title='Digi SBT Y vs X hit map', nx=800, ny=600, cx=1, cy=1)
  h['Digi_SBTHits_YvsX'].Draw('colz')
  ut.bookCanvas(h, key='Digi_SBTHits_YvsX_canvas', title='Digi SBT Z hit map', nx=800, ny=600, cx=1, cy=1)
  h['Digi_SBTHits_Z'].Draw()
  
  print('finished making plots')

# start event loop
def myEventLoop(n):
  rc = sTree.GetEntry(n)
  sTree.GetEntry(n)
  for Digi_SBTHits_it in sTree.Digi_SBTHits:
    X_SBTHits = Digi_SBTHits_it.GetXYZ().X()
    Y_SBTHits = Digi_SBTHits_it.GetXYZ().Y()
    Z_SBTHits = Digi_SBTHits_it.GetXYZ().Z()

    h['Digi_SBTHits_YvsX'].Fill(X_SBTHits, Y_SBTHits)
    h['Digi_SBTHits_Z'].Fill(Z_SBTHits)
  
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
