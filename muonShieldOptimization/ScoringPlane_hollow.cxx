#include "ScoringPlane_hollow.h"
#include <math.h>
#include "vetoPoint.h"

#include "FairLogger.h"                 // for FairLogger, MESSAGE_ORIGIN
#include "FairVolume.h"
#include "FairGeoVolume.h"
#include "FairGeoNode.h"
#include "FairRootManager.h"
#include "FairGeoLoader.h"
#include "FairGeoInterface.h"
#include "FairGeoMedia.h"
#include "FairGeoBuilder.h"
#include "FairRun.h"
#include "FairRuntimeDb.h"
#include "ShipDetectorList.h"
#include "ShipStack.h"

#include "TClonesArray.h"
#include "TVirtualMC.h"
#include "TGeoManager.h"
#include "TGeoBBox.h"
#include "TGeoEltu.h"
#include "TGeoBoolNode.h"
#include "TGeoMaterial.h"
#include "TGeoCompositeShape.h"
#include "TGeoArb8.h"
#include "TParticle.h"
#include "TROOT.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TDatabasePDG.h"

#include <iostream>
using std::cout;
using std::endl;

ScoringPlane_hollow::ScoringPlane_hollow(): FairDetector("ScoringPlane_hollow", kTRUE, kStraw),
					    fTrackID(-1),
					    fVolumeID(-1),
					    fPos(),
					    fMom(),
					    fTime(-1.),
					    fLength(-1.),
					    fELoss(-1.),
					    fOnlyMuons(kFALSE),
					    //fSkipNeutrinos(kFALSE),
					    fxPos(0.0),
					    fyPos(0.0),
					    fzPos(3E8),
					    withNtuple(kFALSE),
					    fScoringPlanePointCollection(new TClonesArray("vetoPoint")),
					    fLastDetector(kFALSE),
					    fFastMuon(kFALSE),
					    fFollowMuon(kFALSE),
					    fSND(kFALSE),
					    fVetoName("veto"),
					    fLx(999.9), // cm
					    fLy(999.9), // cm
					    fLz(000.1),  // cm
					    fMediumName("vacuum"), // Initialize fMediumName to a default
					    fShapeType("Box"), // NEW: Default to "Box"
					    fD_Z_Arb8_hollow(0.0),
					    fArb8Corners_outer{{0.0}},
					    fArb8Corners_inner{{0.0}} {}

ScoringPlane_hollow::ScoringPlane_hollow(const char* name, Bool_t active, Bool_t islastdetector) : FairDetector(name, active, kStraw),
												   fTrackID(-1),
												   fVolumeID(-1),
												   fPos(),
												   fMom(),
												   fTime(-1.),
												   fLength(-1.),
												   fELoss(-1.),
												   fOnlyMuons(kFALSE),
												   //fSkipNeutrinos(kFALSE),
												   fxPos(0.0),
												   fyPos(0.0),
												   fzPos(3E8),
												   withNtuple(kFALSE),
												   fScoringPlanePointCollection(new TClonesArray("vetoPoint")),
												   fLastDetector(islastdetector),
												   fFastMuon(kFALSE),
												   fSND(kFALSE),
												   fFollowMuon(kFALSE),
												   fVetoName("veto"),
												   fLx(999.9), // cm
												   fLy(999.9), // cm
												   fLz(000.1),  // cm
												   fMediumName("vacuum"), // Initialize fMediumName to a default
												   fShapeType("Box"), // NEW: Default to "Box"
												   fD_Z_Arb8_hollow(0.0),
												   fArb8Corners_outer{{0.0}},
												   fArb8Corners_inner{{0.0}} {
}

ScoringPlane_hollow::ScoringPlane_hollow(const char* name, Bool_t active, Bool_t islastdetector,
					 Double_t Lx=999.9, Double_t Ly=999.9, Double_t Lz=0.1) : FairDetector(name, active, kStraw),
												  fTrackID(-1),
												  fVolumeID(-1),
												  fPos(),
												  fMom(),
												  fTime(-1.),
												  fLength(-1.),
												  fELoss(-1.),
												  fOnlyMuons(kFALSE),
												  //fSkipNeutrinos(kFALSE),
												  fxPos(0.0),
												  fyPos(0.0),
												  fzPos(3E8),
												  withNtuple(kFALSE),
												  fScoringPlanePointCollection(new TClonesArray("vetoPoint")),
												  fLastDetector(islastdetector),
												  fFastMuon(kFALSE),
												  fFollowMuon(kFALSE),
												  fSND(kFALSE),
												  fVetoName("veto"),
												  fMediumName("vacuum"), // Initialize fMediumName to a default
												  fShapeType("Box"), // NEW: Default to "Box"
												  fD_Z_Arb8_hollow(0.0),
												  fArb8Corners_outer{{0.0}},
												  fArb8Corners_inner{{0.0}} {
  fLx = Lx; // cm
  fLy = Ly; // cm
  fLz = Lz; // cm
  cout << this->GetName() << ", ConstructGeometry(): constructed ScoringPlane_hollow with X half width "<< fLx << "cm ";
  cout << " and Y half height " << fLy << "cm " << endl;
}

ScoringPlane_hollow::~ScoringPlane_hollow() {
  if (fScoringPlanePointCollection) {
    fScoringPlanePointCollection->Delete();
    delete fScoringPlanePointCollection;
  }
}

Bool_t  ScoringPlane_hollow::ProcessHits(FairVolume* vol) {
  /** This method is called from the MC stepping */
  //Set parameters at entrance of volume. Reset ELoss.
  if ( gMC->IsTrackEntering() ) {
    fELoss  = 0.;
    fTime   = gMC->TrackTime() * 1.0e09;
    fLength = gMC->TrackLength();
    gMC->TrackPosition(fPos);
    gMC->TrackMomentum(fMom);
  }
  // Sum energy loss for all steps in the active volume
  fELoss += gMC->Edep();
  
  // Create vetoPoint when exiting active volume
  if ( gMC->IsTrackExiting()    ||
       gMC->IsTrackStop()       || 
       gMC->IsTrackDisappeared()   ) {
    
    // if (fELoss == 0. ) { return kFALSE; } // if you do not want hits with zero eloss
    
    TParticle* p = gMC->GetStack()->GetCurrentTrack();
    Int_t pdgCode = p->GetPdgCode();
    if (!(fOnlyMuons && TMath::Abs(pdgCode)!=13)){ 
      fTrackID  = gMC->GetStack()->GetCurrentTrackNumber();
      Int_t detID;
      gMC->CurrentVolID(detID);
      TLorentzVector Pos;
      gMC->TrackPosition(Pos);
      TLorentzVector Mom;
      gMC->TrackMomentum(Mom);
      Double_t xmean = (fPos.X()+Pos.X())/2. ;
      Double_t ymean = (fPos.Y()+Pos.Y())/2. ;
      Double_t zmean = (fPos.Z()+Pos.Z())/2. ;
      
      AddHit(fTrackID, detID,
	     //TVector3(xmean, ymean, zmean), put entrance and exit instead
	     TVector3(fPos.X(), fPos.Y(), fPos.Z()),     // entrance position
	     TVector3(fMom.Px(), fMom.Py(), fMom.Pz()),  // entrance momentum
	     fTime, fLength, fELoss,pdgCode,
	     TVector3(Pos.X(),Pos.Y(),Pos.Z()),          // exit position
	     TVector3(Mom.Px(), Mom.Py(), Mom.Pz()) );   // exit momentum
      ShipStack* stack = (ShipStack*) gMC->GetStack();
      stack->AddPoint(kStraw);
    }
  }
  if (fLastDetector) gMC->StopTrack();
  return kTRUE;
}

void ScoringPlane_hollow::Initialize() {
  FairDetector::Initialize();
  TDatabasePDG* PDG = TDatabasePDG::Instance();
}

void ScoringPlane_hollow::EndOfEvent() {
  //std::cout << this->GetName() << this->GetName() << " EndOfEvent(): point collection has " << fScoringPlanePointCollection->GetEntries() << " entries" << std::endl;
  fScoringPlanePointCollection->Clear();
}

void ScoringPlane_hollow::PreTrack() {
  if (!fFastMuon){return;}
  if (TMath::Abs(gMC->TrackPid())!=13){
    gMC->StopTrack();
  }
}

void ScoringPlane_hollow::ConstructGeometry() {
  static FairGeoLoader *geoLoad=FairGeoLoader::Instance();
  static FairGeoInterface *geoFace=geoLoad->getGeoInterface();
  static FairGeoMedia *media=geoFace->getMedia();
  static FairGeoBuilder *geoBuild=geoLoad->getGeoBuilder();
  
  // Use the member variable fMediumName to get the tailored medium
  FairGeoMedium *shipMedium = media->getMedium(fMediumName);
  TGeoMedium* vac = gGeoManager->GetMedium(fMediumName);
  
  if (vac == NULL) {
    std::cout << this->GetName() << ", ConstructGeometry(): Creating medium '" << fMediumName << "'" << std::endl;
    geoBuild->createMedium(shipMedium);
  }
  vac = gGeoManager->GetMedium(fMediumName); // Re-get it to ensure it's loaded
  
  if (vac == NULL) {
    std::cerr << this->GetName() << ", ConstructGeometry(): ERROR - Medium '" << fMediumName << "' could not be created or retrieved!" << std::endl;
    // Handle error: perhaps fall back to a default medium or exit
    return;
  }
  
  TGeoVolume *top=gGeoManager->GetTopVolume();
  TGeoNavigator* nav = gGeoManager->GetCurrentNavigator();
  Double_t xLoc,yLoc,zLoc;

  xLoc = fxPos; // use external input
  yLoc = fyPos; // use external input
  zLoc = fzPos; // use external input
  
  // Declare sensPlane here so it's visible to the entire function scope
  TGeoVolume *sensPlane = nullptr;
  TString myname(this->GetName());
  TString shapename_prefix(myname); // Use a prefix for shapenames
  
  if (fShapeType == "Arb8_hollow") {
    TString outerName = shapename_prefix + "_outer";
    TString innerName = shapename_prefix + "_inner";
    
    // create shapes
    TGeoArb8* outerShape = new TGeoArb8(outerName, fD_Z_Arb8_hollow, fArb8Corners_outer.data());
    TGeoArb8* innerShape = new TGeoArb8(innerName, fD_Z_Arb8_hollow, fArb8Corners_inner.data());

    // register shapes
    gGeoManager->AddShape(outerShape);
    gGeoManager->AddShape(innerShape);

    // boolean subtraction
    TString boolExpr = outerName + "-" + innerName;

    TGeoCompositeShape* vesselShape = new TGeoCompositeShape(shapename_prefix + "_shape", boolExpr);

    sensPlane = new TGeoVolume(shapename_prefix, vesselShape, vac);

    cout << this->GetName()
         << ", ConstructGeometry(): Created hollow Arb8 decay vessel with dz=" << fD_Z_Arb8_hollow << " and custom corners." 
         << endl;
  }

  else {
    std::cerr << this->GetName() << ", ConstructGeometry(): ERROR - Unknown shape type '" << fShapeType << "'. Using default Box." << std::endl;
    sensPlane = gGeoManager->MakeBox(shapename_prefix, vac, fLx, fLy, fLz); // Fallback to Box
  }
  
  if (fMediumName.EqualTo("vacuum")) {
    sensPlane->SetLineColor(kBlue - 10);
  } else if (fMediumName.EqualTo("helium")) {
    sensPlane->SetLineColor(kGreen);
  } else {
    // Default color if no specific material match
    sensPlane->SetLineColor(kGray);
  }

  nav->GetCurrentNode()->GetVolume()->AddNode(sensPlane, 1, new TGeoTranslation(xLoc, yLoc, zLoc));
  AddSensitiveVolume(sensPlane);
  
  // only for fastMuon simulation, otherwise output becomes too big
  if (fFastMuon && fFollowMuon){
    const char* Vol  = "TGeoVolume";
    const char* Mag  = "Mag";
    //const char* Rock = "rock";
    const char* Cavern = "Cavern";
    //const char* Ooo = "Ooo"; //diluted magnets
    //const char* Shi  = "Shi"; // added by Massi, for shielding
    //const char* Coi  = "Coi"; // added by Massi, for coil
    //const char* Ram  = "Ram"; // added by Massi, for all Ram pieces, including Hadron Stopper
    //const char* Ain  = "AbsorberAdd";
    //const char* Aout = "AbsorberAddCore";
    TObjArray* volumelist = gGeoManager->GetListOfVolumes();
    int lastvolume = volumelist->GetLast();
    int volumeiterator=0;
    while ( volumeiterator<=lastvolume ) {
      const char* volumename = volumelist->At(volumeiterator)->GetName();
      const char* classname  = volumelist->At(volumeiterator)->ClassName();
      if (strstr(classname,Vol)){
         // Cavern is always stored in scoring plane
	if (strstr(volumename,Cavern)) {
           AddSensitiveVolume(gGeoManager->GetVolume(volumename));
           cout << this->GetName() << ", ConstructGeometry(): made sensitive for following muons: "<< volumename <<endl;
         }
         // Magnets are stored if there is no SND
         else if (strstr(volumename,Mag) && !fSND) {
           AddSensitiveVolume(gGeoManager->GetVolume(volumename));
           cout << this->GetName() << ", ConstructGeometry(): made sensitive for following muons: "<< volumename <<endl;
         }
      }
      volumeiterator++;
    }
  }
}

vetoPoint* ScoringPlane_hollow::AddHit(Int_t trackID, Int_t detID,
				TVector3 pos, TVector3 mom,
				Double_t time, Double_t length,
				Double_t eLoss, Int_t pdgCode,TVector3 Lpos, TVector3 Lmom) {
  TClonesArray& clref = *fScoringPlanePointCollection;
  Int_t size = clref.GetEntriesFast();
  return new(clref[size]) vetoPoint(trackID, detID, pos, mom,
				    time, length, eLoss, pdgCode,Lpos,Lmom);
}

void ScoringPlane_hollow::Register() {
  TString name  = fVetoName+"Point";
  TString title = fVetoName;
  FairRootManager::Instance()->Register(name, title, fScoringPlanePointCollection, kTRUE);
  // std::cout << this->GetName() << ",  Register() says: registered " << fVetoName <<" collection"<<std::endl;
}

TClonesArray* ScoringPlane_hollow::GetCollection(Int_t iColl) const {
  if (iColl == 0) { return fScoringPlanePointCollection; }
  else { return NULL; }
}

void ScoringPlane_hollow::Reset() {
  fScoringPlanePointCollection->Clear();
}

ClassImp(ScoringPlane_hollow)
