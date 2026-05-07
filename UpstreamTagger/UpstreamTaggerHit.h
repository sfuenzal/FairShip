#ifndef UpstreamTaggerHIT_H
#define UpstreamTaggerHIT_H 1
#include "FairVolume.h"
#include "ShipHit.h"
#include "UpstreamTaggerPoint.h"
#include "UpstreamTagger.h"
#include "TObject.h"
#include "TGeoShape.h"
#include "TGeoPhysicalNode.h"


class UpstreamTaggerHit : public ShipHit, UpstreamTagger
{
  public:

    /** Default constructor **/
    UpstreamTaggerHit();

    UpstreamTaggerHit(UpstreamTaggerPoint* p, UpstreamTagger* c, Double_t t0);

    /** Destructor **/
    virtual ~UpstreamTaggerHit();

    /** Accessors **/
    Double_t GetX();
    Double_t GetY();
    Double_t GetZ();
    TVector3 GetXYZ();


    TGeoNode* GetNode(Double_t &hit_final, Int_t &mod);
    std::vector<double> GetTime(Double_t x);
    std::vector<double> GetTime();
    std::vector<double> GetMeasurements();
    /** Modifier **/
    void SetPoint(Double_t p1, Double_t p2, Double_t p3){point_final[0]=p1;point_final[1]=p2;point_final[2]=p3;}
    /** Output to screen **/
    virtual void Print() const;

    void setInvalid() {flag = false;}
    void setIsValid() {flag = true;}

    //Rpc time is invalid if isValid returns False
    bool isValid() const {return flag;}
  private:
    UpstreamTaggerHit(const UpstreamTaggerHit& point);
    UpstreamTaggerHit operator=(const UpstreamTaggerHit& point);

    UpstreamTagger* c0;
    Double_t point_final[3];
    const Double_t * mom[3];

    Float_t flag;     ///< flag
    Float_t time;
    Double_t X, Y, Z;
  
    ClassDef(UpstreamTaggerHit,1);

};

#endif
