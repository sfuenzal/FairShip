#ifndef SCORINGPLANEHOLLOW_H
#define SCORINGPLANEHOLLOW_H

#include "FairDetector.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TGeoVolume.h"
#include "vetoPoint.h"
#include "TNtuple.h"
#include "TString.h"
#include <map>
#include <array>

class FairVolume;
class TClonesArray;

class ScoringPlane_hollow: public FairDetector {

public:

  /**      Name :  Detector Name
   *       Active: kTRUE for active detectors (ProcessHits() will be called)
   *               kFALSE for inactive detectors
   */
  ScoringPlane_hollow(const char* Name, Bool_t Active, Bool_t islastdetector);
  
  /**      Name :  Detector Name
   *       Active: kTRUE for active detectors (ProcessHits() will be called)
   *               kFALSE for inactive detectors
   *       geometry: full x, y and z extents in cm are Lx, Ly, Lz
   */
  ScoringPlane_hollow(const char* Name, Bool_t Active, Bool_t islastdetector, 
	       Double_t Lx, Double_t Ly, Double_t Lz);
  
  /**      default constructor    */
  ScoringPlane_hollow();
  
  /**       destructor     */
  virtual ~ScoringPlane_hollow();
  
  /**      Initialization of the detector is done here    */
  virtual void   Initialize();
  
  /**       this method is called for each step during simulation
   *       (see FairMCApplication::Stepping())
   */
  virtual Bool_t ProcessHits( FairVolume* v=0);

  /**       Registers the produced collections in FAIRRootManager.     */
  virtual void   Register();
  
  /** Gets the produced collections */
  virtual TClonesArray* GetCollection(Int_t iColl) const ;
  
  /**      has to be called after each event to reset the containers      */
  virtual void   Reset();
  
  /**      Create the detector geometry        */
  void ConstructGeometry();
  
  /** The following methods can be implemented if you need to make
   *  any optional action in your detector during the transport.
   */
  
  void SetMediumName(const char* name) {
    fMediumName = name;
  }
  
  // Setter for the shape type (e.g., "Box", "Trapezoid")
  void SetShapeType(const char* type) { fShapeType = type; }

  // This is intended for the veto counting plane
  void SetArb8Dimensions_hollow(Double_t dZ_arb8_hollow,
				const std::array<Double_t, 16>& corners_arb8_outer,
				const std::array<Double_t, 16>& corners_arb8_inner) {
    fD_Z_Arb8_hollow = dZ_arb8_hollow;
    fArb8Corners_outer = corners_arb8_outer;
    fArb8Corners_inner = corners_arb8_inner;
  }
  
  virtual void   CopyClones( TClonesArray* cl1,  TClonesArray* cl2 ,
			     Int_t offset) {;}
  virtual void   SetSpecialPhysicsCuts() {;}
  virtual void   EndOfEvent();
  virtual void   FinishPrimary() {;}
  virtual void   FinishRun() {;}
  virtual void   BeginPrimary() {;}
  virtual void   PostTrack() {;}
  virtual void   PreTrack();
  virtual void   BeginEvent() {;}
  
  vetoPoint* AddHit(Int_t trackID, Int_t detID,
		    TVector3 pos, TVector3 mom,
		    Double_t time, Double_t length,
		    Double_t eLoss,Int_t pdgcode,TVector3 Lpos, TVector3 Lmom);
  //inline void SetEnergyCut(Float_t emax) {EMax=emax;}// min energy to be copied to Geant4
  inline void SetOnlyMuons(){fOnlyMuons=kTRUE; std::cout<<"Massi ScoringPlane_hollow.SetOnlyMuons(): fOnlyMuons="<<fOnlyMuons<<std::endl;}
  inline void SetOpt4DP(){withNtuple=kTRUE;}
  //inline void SkipNeutrinos(){fSkipNeutrinos=kTRUE;}
  inline void SetXposition(Float_t x){fxPos=x;}
  inline void SetYposition(Float_t y){fyPos=y;}
  inline void SetZposition(Float_t z){fzPos=z;}
  inline void SetXYZposition(Float_t x,Float_t y,Float_t z){fxPos=x;fyPos=y;fzPos=z;}
  inline void SetIsLast(Bool_t islast){fLastDetector=islast;} // Added by Massi
  inline void SetVetoPointName(TString nam){fVetoName=nam;} // Added by Massi
  // kill all tracks except of muons:
  void SetFastMuon() {fFastMuon=true; std::cout<<"Massi ScoringPlane_hollow.SetFastMuon(): fFastMuon="<<fFastMuon<<std::endl;} 
  // make muon shield active to follow muons:
  void SetFollowMuon() {fFollowMuon=true;std::cout<<"Massi ScoringPlane_hollow.SetFollowMuon(): fFollowMuon="<<fFollowMuon<<std::endl;} 
  void SetSND() {fSND=true;std::cout<<"Guglielmo ScoringPlane_hollow.SetSND(): fSND="<<fSND<<std::endl;} 
private:
  /** Track information to be stored until the track leaves the
      active volume.
  */
  Bool_t     fFastMuon, fFollowMuon, fSND; 
  Int_t          fTrackID;           //!  track index
  Int_t          fVolumeID;          //!  volume id
  TLorentzVector fPos;               //!  position at entrance
  TLorentzVector fMom;               //!  momentum at entrance
  Double_t     fTime;              //!  time
  Double_t     fLength;            //!  length
  Double_t     fELoss;             //!  energy loss
  Double_t     fxPos;              //!  xPos, optional
  Double_t     fyPos;              //!  yPos, optional
  Double_t     fzPos;              //!  zPos, optional
  Bool_t withNtuple;               //! special option for Dark Photon physics studies
  TNtuple* fNtuple;               //!  
  //Float_t EMax;  //! max energy to transport
  Bool_t fOnlyMuons;      //! flag if only muons should be stored
  //Bool_t fSkipNeutrinos;  //! flag if neutrinos should be ignored
  TFile* fout; //!
  TClonesArray* fElectrons; //!
  Int_t index;

  // NEW: Member variable for the tailored medium name
  TString     fMediumName;
  
  /** container for data points */
  TClonesArray*  fScoringPlanePointCollection;
  ClassDef(ScoringPlane_hollow, 0)
  // massi, add this to control the stopMC:
  Bool_t fLastDetector;  //! if True then stop processing particles after this detector 
  Double_t     fLx;      //!  x full extent in cm 
  Double_t     fLy;      //!  y full extent in cm 
  Double_t     fLz;      //!  z full extent in cm 
  
  TString fVetoName;

  TString     fShapeType;
  Double_t fD_Z_Arb8_hollow;      // Half-length along Z for Arb8
  std::array<Double_t, 16> fArb8Corners_outer;  // Outer hollow decay vessel. 8 (x,y) corners for front and 8 (x,y) for back face
  std::array<Double_t, 16> fArb8Corners_inner;  // Inner hollow decay vessel. 8 (x,y) corners for front and 8 (x,y) for back face
};

#endif //SCORINGPLANEHOLLOW_H
