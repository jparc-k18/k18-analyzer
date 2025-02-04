#ifndef TPC_HITCLUSTER_HH
#define TPC_HITCLUSTER_HH

#include "TPCHit.hh"

#include "TVector3.h"
#include "TMatrixD.h"

#include <vector>

class TPCHitCluster
{

public:
  TPCHitCluster();
  TPCHitCluster(TPCHitCluster *cluster);
  virtual ~TPCHitCluster() {}

  void Clear(Option_t * = "");

  void SetCovMatrix(TMatrixD matrix);
  TMatrixD GetCovMatrix() const {return fCovMatrix;}
  
  Bool_t IsClusterised() const {return true;}

  Int_t GetNumHits() const {return fHits.size();}
  
  virtual void AddHit(TPCHit *hit);
  void ApplyModifiedHitInfo();
  
  void SetDFromCov(Double_t maxx, Double_t maxy, Double_t maxz, Bool_t setMin = false);

  void SetMeanTime(Double_t time) {fMeanTime = time;}
  Double_t GetMeanTime() const {return fMeanTime;}
  
  void SetLength(Double_t length) {fLength = length;}
  Double_t GetLength() const {return fLength;}

  void SetIsStable(Bool_t isStable) {fIsStable = isStable;}
  Int_t IsStable() const {
    if (fIsStable == true)
      return 1;
    else
      return 0;
  }
  
  void SetPOCA(TVector3 pos) {fPOCA = pos;}
  TVector3 GetPOCA() const {return fPOCA;}

  void SetType(Bool_t flag) {fType = flag;}
  
  Double_t GetType() const {
    if (fType == true)
      return 1.;
    else
      return -1.;
  }

  virtual void SetClusterID(Int_t clusterID);
  
  virtual void SetTrackID(Int_t trackID) {fTrackID = trackID ;}
  Int_t GetTrackID() const { return fTrackID;}
  
  void SetIsContinuousHits(Bool_t value) { fIsContinuousHits = value;}
  Bool_t IsContinuousHits() const {return fIsContinuousHits;}

  void ApplyCovLowLimit();

  void SetChi(Double_t chi) { fChi = chi;}
  Double_t GetChi() const {return fChi;}

  void SetLambda(Double_t dip) {fLambda = dip;}
  Double_t GetLambda() const {return fLambda;}

  void SetCharge(Double_t charge) {fCharge = charge;}
  void SetMeanCharge(Double_t charge) {fMeanCharge = charge;}
  Double_t GetCharge() const {return fCharge;}
  Double_t GetMeanCharge() const {return fMeanCharge;}
  Double_t GetChargeSqr() const {return fChargeSqr;}

  void SetDirection(TVector3 dir) {fDirection = dir;}
  TVector3 GetDirection() const {return fDirection;}
  
  void SetRow(Int_t row) {fRow = row;}
  Int_t GetRow() const {return fRow;}
  
  void SetLayer(Int_t layer) {fLayer = layer;}
  Int_t GetLayer() const {return fLayer;}
  
  std::vector<TPCHit*> *GetHits() {return &fHits;}
  TPCHit *GetHit(int i = 0) {return fHits.at(i);}
  
  Double_t GetX() const {return fX;}
  Double_t GetY() const {return fY;}
  Double_t GetZ() const {return fZ;}
  Double_t GetR() const {return sqrt(fX*fX + (fZ+143.)*(fZ+143.));}
  Double_t GetAlpha() const {return TMath::ATan2(fX, fZ+143.);}
  TVector3 GetPosition() const {return TVector3(fX,fY,fZ);}
  
  Double_t SetDx(Double_t dx)  {return fDx = dx;}
  Double_t SetDy(Double_t dy)  {return fDy = dy;}
  Double_t SetDz(Double_t dz)  {return fDz = dz;}

  Double_t GetDx() const {return fDx;}
  Double_t GetDy() const {return fDy;}
  Double_t GetDz() const {return fDz;}
  TVector3 GetCovariance() const { return TVector3(fDx, fDy, fDz);}

  inline Double_t operator[](int i) const {
    switch(i) {
    case 0:
      return fX;
    case 1:
      return fY;
    case 2:
      return fZ;
    default:
      Error("operator[](i)", "bad index (%d) returning 0",i);
    }
    return 0.;
  }

  inline Double_t GetVerticalRes(Double_t y, Double_t*par);
  inline Double_t GetTransverseRes(Double_t angle, Double_t y, Double_t*par);
  inline TVector3 GetResVector(TVector3 pos, Double_t *par_t, Double_t *par_y);
  
private:

  Int_t fClusterID;
  Int_t fTrackID;
  
  Bool_t fIsClusterised;
  Bool_t fIsStable;

  Double_t fX;
  Double_t fY;
  Double_t fZ;
  
  Double_t fDx;
  Double_t fDy;
  Double_t fDz;

  Double_t fCharge;
  Double_t fChargeSqr;
  Double_t fMeanCharge;
  
  Double_t fLength;

  Double_t fMeanTime;
  
  Int_t fRow;
  Int_t fLayer;

  TMatrixD fCovMatrix;
  Double_t xMean, yMean, zMean;
  Double_t xyMean, xzMean, yzMean;
  Double_t xxMean, yyMean, zzMean;
  
  //  std::vector<Int_t> fHitIDArray;
  std::vector<TPCHit*> fHits;

  TVector3 fPOCA;
  Double_t fPOCAx;
  Double_t fPOCAy;
  Double_t fPOCAz;

  Double_t fChi;
  Double_t fLambda;

  Double_t fS = 0;

  TVector3 fDirection;
  
  void CalculatePosition(TVector3 hitpos, Double_t charge);

  void CalculateCovMatrix(TVector3 hitpos, Double_t charge);

  Bool_t fIsContinuousHits;
  Bool_t fType;
  
};

#endif




