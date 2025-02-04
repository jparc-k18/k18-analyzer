#ifndef TPC_LINEAR_TRACK_HH
#define TPC_LINEAR_TRACK_HH

#include "TPCHit.hh"
#include "TPCHitCluster.hh"

#include <TObject.h>
#include <TVector3.h>

#include <vector>

class TPCLinearTrack
{

public:
  TPCLinearTrack();
  TPCLinearTrack(TPCLinearTrack *track);
  ~TPCLinearTrack(){}

  Bool_t Initialize();

  void SetIsFitted(Bool_t flag) {fIsFitted = flag;}
  Bool_t IsFitted() const {return fIsFitted;}

  void SetTrackID(Int_t id) {fTrackID = id;}
  Int_t GetTrackID() const {return fTrackID;}
  
  Int_t GetNumHits() const { return fHits.size();}
  TPCHit *GetHit(Int_t iHit) const {return fHits[iHit];}
  std::vector<TPCHit*> *GetHits() {return &fHits;}

  Int_t GetNumClusters() const {return fHitClusters.size();}
  Int_t GetNumStableClusters() const
  {
    Int_t num = 0.;
    for (auto cluster : fHitClusters) {
      if (cluster->IsStable())
	num++;
    }

    return num;
  }

  TPCHitCluster *GetCluster(Int_t iCl) const {return fHitClusters[iCl];}
  std::vector<TPCHitCluster*> *GetClusters() {return &fHitClusters;}

  void AddHit(TPCHit *hit);
  void RemoveHit(TPCHit *hit);
  void DeleteHits();


  void AddHitCluster(TPCHitCluster *cluster);

  void SetChargeSum(Double_t chargesum) {fChargeSum = chargesum;}
  Double_t GetChargeSum() const {return fChargeSum;}

  TVector3 GetMean() const {return TVector3(xMean, yMean, zMean);}

  Double_t GetMeanXX() const {return xxMean;}
  Double_t GetMeanYY() const {return yyMean;}
  Double_t GetMeanZZ() const {return zzMean;}

  Double_t GetMeanXY() const {return xyMean;}
  Double_t GetMeanYZ() const {return yzMean;}
  Double_t GetMeanZX() const {return zxMean;}


  Double_t xxCov() const {return fChargeSum * (xxMean - xMean * xMean);}
  Double_t yyCov() const {return fChargeSum * (yyMean - yMean * yMean);}
  Double_t zzCov() const {return fChargeSum * (zzMean - zMean * zMean);}

  Double_t xyCov() const { return fChargeSum * (xyMean - xMean * yMean);}
  Double_t yzCov() const { return fChargeSum * (yzMean - yMean * zMean);}
  Double_t zxCov() const { return fChargeSum * (zxMean - zMean * xMean);}

  Double_t GetXCov() const { return xxCov()/fChargeSum;}
  Double_t GetZCov() const { return zzCov()/fChargeSum;}

  void SetVertex(Int_t i, Double_t x, Double_t y, Double_t z) {
    fXVertex[i] = x;
    fYVertex[i] = y;
    fZVertex[i] = z;
  };
  
  void SetVertex(Int_t i, TVector3 vtx) {
    fXVertex[i] = vtx.X();
    fYVertex[i] = vtx.Y();
    fZVertex[i] = vtx.Z();
  };
  
  void SetXVertex(Int_t i, Double_t x) {fXVertex[i] = x;}
  void SetYVertex(Int_t i, Double_t y) {fYVertex[i] = y;}
  void SetZVertex(Int_t i, Double_t z) {fZVertex[i] = z;}

  TVector3 GetVertex(Int_t i) const {return TVector3(fXVertex[i], fYVertex[i], fZVertex[i]);}
  Double_t GetXVertex(Int_t i) const {return fXVertex[i];}
  Double_t GetYVertex(Int_t i) const {return fYVertex[i];}
  Double_t GetZVertex(Int_t i) const {return fZVertex[i];}
  
  void SetDirection(TVector3 vec) {fDirection = vec;}
  TVector3 GetDirection() const {return fDirection;}

  void SetNormal(TVector3 vec) {fNormal = vec;}
  TVector3 GetNormal() const { return fNormal;}

  void SetRMSLine(Double_t rms) { fRMSLine = rms;}
  void SetRMSPlane(Double_t rms) { fRMSPlane = rms;}

  Double_t GetRMSLine() const { return fRMSLine;}
  Double_t GetRMSPlane() const { return fRMSPlane;}

  
  Double_t GetSlope() const {return fSlope;}
  Double_t GetIntercept() const {return fIntercept;}

  Double_t GetChi2() const {return fChi2;}

  Double_t GetTrackLength() const {return fTrackLength;}
  Double_t GetTrackLengthY() const {return fTrackLengthY;}
  Double_t GetThetaY() const {return fThetaY;} //angle btw y and xz plane
  
  Double_t GetdEXZ() const {return fdEXZ;}

  TVector3 GetPointOnX(Double_t x);
  TVector3 GetPointOnY(Double_t y);
  TVector3 GetPointOnZ(Double_t z);

  TVector3 PerpLine(TVector3 pos) const;
  TVector3 PerpPlane(TVector3 pos) const;

  Double_t PerpDistLine(TVector3 pos);

  TVector3 GetClosestPointOnTrack(TVector3 pos) { return (pos + PerpLine(pos));}

  void FinalizeHits();
  
private:

  Int_t fTrackID;

  Bool_t fIsFitted;

  TVector3 fDirection;
  TVector3 fNormal;

  Double_t fChargeSum;
  Double_t xMean, yMean, zMean;
  Double_t xxMean, yyMean, zzMean;
  Double_t xyMean, yzMean, zxMean;
  
  Double_t fRMSLine, fRMSPlane;

  Double_t fSlope, fIntercept;

  Double_t fChi2;

  Double_t fTrackLength;
  Double_t fTrackLengthY;
  Double_t fThetaY;

  Double_t fdEXZ;

  Double_t fXVertex[2];
  Double_t fYVertex[2];
  Double_t fZVertex[2];
  
  std::vector<TPCHit*> fHits;
  std::vector<TPCHitCluster*> fHitClusters;
};


#endif
