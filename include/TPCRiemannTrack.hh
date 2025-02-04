#ifndef TPC_RIEMANN_TRACk_HH
#define TPC_RIEMANN_TRACk_HH

#include "TPCHit.hh"
#include "TPCHitCluster.hh"

#include <TObject.h>
#include <TVector3.h>

#include <vector>

class TPCRiemannTrack 
{
public:
  TPCRiemannTrack();
  TPCRiemannTrack(Int_t trackid);
  ~TPCRiemannTrack() {}

  Bool_t Initialize();

  void DeleteHits();


  void FinalizeClusters();

  Bool_t IsAccidental();
  
  Double_t GetRadius() const {return fRadius;}
  void SetRadius(Double_t r) {fRadius = r;}
  
  TVector3 GetCenter() const {return fCenter;}
  void SetCenter(Double_t x, Double_t z) {fCenter.SetXYZ(x,0,z);}
  
  Double_t GetAlpha() const {return fAlpha;}
  Double_t GetDip() const {
    if (fRadius <= 0)
      return -9999;
    
    return  TMath::ATan(fSlope/fRadius);
  }

  Bool_t IsFinished() {return fIsFinished;}
  Bool_t IsGood() {return fIsGood;}

  void SetTrackID(Int_t id) { fTrackID = id;}
  Int_t GetTrackID() const {return fTrackID;}
  
  void SetSlope(Double_t val) {fSlope = val;}
  Double_t GetSlope() const {return fSlope;}

  void SetOffset(Double_t val) {fOffset = val;}
  Double_t GetOffset() const {return fOffset;}

  Int_t GetNumHits() const {return fHits.size();}
  TPCHit *GetHit(Int_t iHit) const {return fHits[iHit];}
  TPCHit *GetLastHit() const {return fHits.back();}
  TPCHit *GetFirstHit() const {return fHits.front();}
  std::vector<TPCHit *> *GetHits() {return &fHits;}

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
  
  TPCHitCluster *GetCluster(Int_t iCluster) const {return fHitClusters[iCluster];}
  std::vector<TPCHitCluster *> *GetClusters() {return &fHitClusters;}

  TVector3 GetPosition(Double_t alpha) const {
    return TVector3(fRadius*TMath::Sin(alpha)+fCenter.X(), alpha*fSlope+fOffset, fRadius*TMath::Cos(alpha)+fCenter.Z());
  }  

  TVector3 GetDirection(Double_t alpha) const;

  void SetAlphaHead(Double_t val) {fAlphaHead = val;}
  Double_t GetAlphaHead() const {return fAlphaHead;}
  void SetAlphaTail(Double_t val) {fAlphaTail = val;}
  Double_t GetAlphaTail() const {return fAlphaTail;}

  void SetCharge(Bool_t flag) {fIsPositive = flag;}
  Int_t Charge() const {
    if (fIsPositive == true)
      return 1.;
    else
      return -1;
  }
  
  Int_t Helicity() const { return fSlope > 0 ? 1: -1;}
  
  Double_t GetMomentum(Bool_t isConst = false) const;
  Double_t GetMomResolution(Double_t ResT) const;

  void AddHit(TPCHit *hit);  
  void RemoveHit(TPCHit *hit);

  void AddHitCluster(TPCHitCluster *cluster);
  
  TVector3 ExtrapHead(Double_t length) const;
  TVector3 ExtrapTail(Double_t length) const;

  enum TPCFitStatus { kBad, kLine, kPlane, kHelix};

  void SetFitStatus(TPCFitStatus value) {fFitStatus = value;}
  TPCFitStatus GetFitStatus() const {return fFitStatus;}

  void SetIsBad() { fFitStatus = TPCRiemannTrack::kBad;}
  Bool_t IsBad() const { return fFitStatus == kBad ? true : false;}
  void SetIsLine() { fFitStatus = TPCRiemannTrack::kLine;}
  Bool_t IsLine() const { return fFitStatus == kLine ? true : false;}
  void SetIsPlane() { fFitStatus = TPCRiemannTrack::kPlane;}
  Bool_t IsPlane() const { return fFitStatus == kPlane ? true : false;}
  void SetIsHelix() { fFitStatus = TPCRiemannTrack::kHelix;}
  Bool_t IsHelix() const { return fFitStatus == kHelix ? true : false;}
  
  enum TPCRiemannTrackStatus { kInitialization, kExtension, kExtrapolation, kConfirmation};

  void SetTrackStatus(TPCRiemannTrackStatus value) {fTrackStatus = value;}
  TPCRiemannTrackStatus GetTrackStatus() const {return fTrackStatus;}

  void SetIsInitialized() {fTrackStatus = TPCRiemannTrack::kInitialization;}
  Bool_t IsInitialized() const {return fTrackStatus == kInitialization ? true : false;}
  void SetIsExtended() {fTrackStatus = TPCRiemannTrack::kExtension;}
  Bool_t IsExtended() const {return fTrackStatus == kExtension ? true : false;}
  void SetIsExtrapolated() {fTrackStatus = TPCRiemannTrack::kExtrapolation;}
  Bool_t IsExtrapolated() const {return fTrackStatus == kExtrapolation ? true : false;}
  void SetIsConfirmed() {fTrackStatus = TPCRiemannTrack::kConfirmation;}
  Bool_t IsConfirmed() const {return fTrackStatus == kConfirmation ? true : false;}

  enum TPCRiemannTrackProperty {kIsBeam, kIsKurama, kIsScat, kIsAccBeam};

  void SetTrackProperty(TPCRiemannTrackProperty value) {fTrackProperty = value;}
  TPCRiemannTrackProperty GetTrackProperty()  const {return fTrackProperty;}

  void SetIsBeam() {fTrackProperty = TPCRiemannTrackProperty::kIsBeam;}
  Bool_t IsBeam() const { return fTrackProperty == kIsBeam ? true : false;}
  void SetIsKurama() {fTrackProperty = TPCRiemannTrackProperty::kIsKurama;}
  Bool_t IsKurama() const { return fTrackProperty == kIsKurama ? true : false;}
  void SetIsScat() {fTrackProperty = TPCRiemannTrackProperty::kIsScat;}
  Bool_t IsScat() const { return fTrackProperty == kIsScat ? true : false;}
  void SetIsAccBeam() {fTrackProperty = TPCRiemannTrackProperty::kIsAccBeam;}
  Bool_t IsAccBeam() const { return fTrackProperty == kIsAccBeam ? true : false;}
  
  TVector3 PerpLine(TVector3 pos) const;
  TVector3 PerpPlane(TVector3 pos) const;

  void SetLineDirection(TVector3 dir){
    fRadius = dir.Z();
    fCenter.SetX(dir.X());
    fCenter.SetZ(dir.Y());
    
  } 
  TVector3 GetLineDirection() const {return TVector3(fCenter.X(),fCenter.Z(),fRadius);}

  void SetPlaneNormal(TVector3 normal) {
    fRadius = normal.Z();
    fCenter.SetX(normal.X());
    fCenter.SetZ(normal.Y());
  }
  
  TVector3 GetPlaneNormal() const {return TVector3(fCenter.X(),fCenter.Z(),fRadius);}

  TVector3 GetMean() const {return TVector3(xMean, yMean, zMean);}
  Double_t GetMeanXX() const { return xxMean;};
  Double_t GetMeanYY() const { return yyMean;};
  Double_t GetMeanZZ() const { return zzMean;};
  Double_t GetMeanXY() const { return xyMean;};
  Double_t GetMeanYZ() const { return yzMean;};
  Double_t GetMeanZX() const { return zxMean;};
  
  Double_t GetChargeSum()   const {return fChargeSum;}
  Double_t xxCov() const {return fChargeSum * (xxMean - xMean * xMean);}
  Double_t yyCov() const {return fChargeSum * (yyMean - yMean * yMean);}
  Double_t zzCov() const {return fChargeSum * (zzMean - zMean * zMean);}

  Double_t xyCov() const { return fChargeSum * (xyMean - xMean * yMean);}
  Double_t yzCov() const { return fChargeSum * (yzMean - yMean * zMean);}
  Double_t zxCov() const { return fChargeSum * (zxMean - zMean * xMean);}

  Double_t GetXCov() const { return xxCov()/fChargeSum;}
  Double_t GetYCov() const { return yyCov()/fChargeSum;}
  Double_t GetZCov() const { return zzCov()/fChargeSum;}
  
  void SetPDGEncoding(Int_t pid) { fPID = pid;}
  Int_t GetPDGEncoding() const {return fPID;}

  void SetRMSW(Double_t val) {fRMSW = val;}
  Double_t GetRMSW() const {return fRMSW;}
  void SetRMSH(Double_t val) {fRMSH = val;}
  Double_t GetRMSH() const {return fRMSH;}

  void SetPOCA(TVector3 val) {fPOCA = val;}
  TVector3 GetPOCA() const {return fPOCA;}

  void SetClosestDist(Double_t val) {fdistPOCA = val;}
  Double_t GetClosestDist() const {return fdistPOCA;}

  TVector3 PositionAtHead() const {
    return GetPosition(fAlphaHead);
  }
  
  TVector3 PositionAtTail() const {
    return GetPosition(fAlphaTail);
  }
  
  Double_t TrackLength() const {
    return std::abs(GetAlphaHead() - GetAlphaTail()) * fRadius / TMath::Cos(GetDip());
  }

  Double_t DistCircle(TVector3 pos) const {
    Double_t dx = pos.X() - fCenter.X();
    Double_t dz = pos.Z() - fCenter.Z();
    return sqrt(dx*dx + dz*dz) - fRadius;
  }
  
  TVector3 Map(TVector3 pos) const;
  Double_t ExtrapByMap(TVector3 pos, TVector3 &q, TVector3 &m) const;
  Double_t ExtrapToAlpha(Double_t alpha) const {
    return alpha * fRadius / TMath::Cos(GetDip());
  }
  Double_t ExtrapToAlpha(Double_t alpha, TVector3 &pt) const {
    pt.SetXYZ(fRadius*TMath::Sin(alpha) + fCenter.X(), alpha*fSlope + fOffset, fRadius*TMath::Cos(alpha) + fCenter.Z());
    Double_t length = alpha * fRadius / TMath::Cos(GetDip());
    return length;
  }

  Double_t AlphaAtPosition(TVector3 pos);
  
  Double_t ExtrapToPointAlpha(TVector3 pos, TVector3 &pt, Double_t &alpha) const;

  TVector3 InterpByLength(Double_t length) const;
  TVector3 InterpByRatio(Double_t ratio) const;

  Double_t Continuity();
  
  Bool_t ExtrapToZ(Double_t zl, TVector3 &pointOnHelix1, Double_t &alpha1,
		   TVector3 &pointOnHelix2, Double_t &alpha2) const;  
  Bool_t ExtrapToZ(Double_t z, Double_t alphaRef, TVector3 &pointOnHelix) const;
  Bool_t ExtrapToZ(Double_t z, TVector3 &pointOnHelix) const;
  Bool_t CheckExtrapToZ(Double_t z) const;


  Bool_t ExtrapToX(Double_t xl, TVector3 &pointOnHelix1, Double_t &alpha1,
		   TVector3 &pointOnHelix2, Double_t &alpha2) const;  
  Bool_t ExtrapToX(Double_t x, Double_t alphaRef, TVector3 &pointOnHelix) const;
  Bool_t ExtrapToX(Double_t x, TVector3 &pointOnHelix) const;
  Bool_t CheckExtrapToX(Double_t x) const;

  void DetermineCharge(TVector3 vertex);

  TVector3 Direction(Double_t alpha) const;
  
  TVector3 GetPrimaryVertex() const {
    return fVertex;
  }

  void SetDirectionAtTgt(TVector3 vec) {fdirTgt = vec;}
  TVector3 GetDirectionAtTgt() const {
    return fdirTgt;
  }

  Int_t GetWinding() const;

  Double_t GetTruncateddEdx(Double_t LL, Double_t HL);
  Double_t GetTruncatedMeandEdx(Double_t LL, Double_t HL);

  std::vector<Double_t> *GetdEdxArray() {return &fdEdxArray;}
  std::vector<Double_t> *GetMeandEdxArray() {return &fMeandEdxArray;}
  
private:
  
  Int_t fTrackID;
  
  TPCFitStatus fFitStatus;
  TPCRiemannTrackStatus fTrackStatus;
  TPCRiemannTrackProperty fTrackProperty;

  Int_t fPID;
  
  TVector3 fCenter; // center of helix in 3d (y = 0 )
  Double_t fRadius; // Radius of helix
  Double_t fAlpha;
  Double_t fAlphaHead, fAlphaTail;
  
  Double_t fSlope;
  Double_t fOffset;

  Double_t fRMS, fRMSPlane;

  Bool_t fIsFitted, fIsInitialized, fIsFinished, fIsGood;

  std::vector<TPCHit*> fHits;
  std::vector<TPCHitCluster*> fHitClusters;

  Double_t fChargeSum; // Weighting the average with hit error

  Double_t xMean, yMean, zMean;
  Double_t xxMean, yyMean, zzMean;
  Double_t xyMean, yzMean, zxMean;
  
  Bool_t fDoSort;

  Double_t fRMSW, fRMSH;

  Bool_t fIsPositive;
  TVector3 fVertex;

  Double_t fdistPOCA;
  TVector3 fPOCA;
  TVector3 fdirTgt;
  
  std::vector<Double_t> fdEdxArray;
  std::vector<Double_t> fMeandEdxArray;

};


class SortByX{
public:
  SortByX() {}
  Bool_t operator() (TPCHit *hit1, TPCHit *hit2) {
    return (hit1->GetPosition().X() > hit2->GetPosition().X());
  }
};

class SortByXInv{
public:
  SortByXInv() {}
  Bool_t operator() (TPCHit *hit1, TPCHit *hit2) {
    return (hit1->GetPosition().X() < hit2->GetPosition().X());
  }
};

class SortByY{
public:
  SortByY() {}
  Bool_t operator() (TPCHit *hit1, TPCHit *hit2) {
    return (hit1->GetPosition().Y() > hit2->GetPosition().Y());
  }
  Bool_t operator() (TPCHitCluster *hitcluster1, TPCHitCluster *hitcluster2) {
    return (hitcluster1->GetPosition().Y() > hitcluster2->GetPosition().Y());
  }
};

class SortByYInv{
public:
  SortByYInv() {}
  Bool_t operator() (TPCHit *hit1, TPCHit *hit2) {
    return (hit1->GetPosition().Y() < hit2->GetPosition().Y());
  }
  Bool_t operator() (TPCHitCluster *hitcluster1, TPCHitCluster *hitcluster2) {
    return (hitcluster1->GetPosition().Y() < hitcluster2->GetPosition().Y());
  }
};

class SortByYAbs{
public:
  SortByYAbs() {}
  Bool_t operator() (TPCHit *hit1, TPCHit *hit2) {
    return abs(hit1->GetPosition().Y()) > abs(hit2->GetPosition().Y());
  }
  Bool_t operator() (TPCHitCluster *hitcluster1, TPCHitCluster *hitcluster2) {
    return abs(hitcluster1->GetPosition().Y()) > abs(hitcluster2->GetPosition().Y());
  }
};

class SortByYAbsInv{
public:
  SortByYAbsInv() {}
  Bool_t operator() (TPCHit *hit1, TPCHit *hit2) {
    return abs(hit1->GetPosition().Y()) < abs(hit2->GetPosition().Y());
  }
  Bool_t operator() (TPCHitCluster *hitcluster1, TPCHitCluster *hitcluster2) {
    return abs(hitcluster1->GetPosition().Y()) < abs(hitcluster2->GetPosition().Y());
  }
};

class SortByZ{
public:
  SortByZ() {}
  Bool_t operator() (TPCHit *hit1, TPCHit *hit2) {
    return (hit1->GetPosition().Z() > hit2->GetPosition().Z());
  }
  Bool_t operator() (TPCHitCluster *hitcluster1, TPCHitCluster *hitcluster2) {
    return (hitcluster1->GetPosition().Z() > hitcluster2->GetPosition().Z());
  }
};

class SortByZInv{
public:
  SortByZInv() {}
  Bool_t operator() (TPCHit *hit1, TPCHit *hit2) {
    return (hit1->GetPosition().Z() < hit2->GetPosition().Z());
  }
  Bool_t operator() (TPCHitCluster *hitcluster1, TPCHitCluster *hitcluster2) {
    return (hitcluster1->GetPosition().Z() < hitcluster2->GetPosition().Z());
  }
};

class SortByR{
public:
  SortByR() {}
  Bool_t operator() (TPCHit *hit1, TPCHit *hit2) {
    return (hit1->GetR() > hit2->GetR());
  }
  Bool_t operator() (TPCHitCluster *hitcluster1, TPCHitCluster *hitcluster2) {
    return (hitcluster1->GetR() > hitcluster2->GetR());
  }
};

class SortByRInv{
public:
  SortByRInv() {}
  Bool_t operator() (TPCHit *hit1, TPCHit *hit2) {
    return (hit1->GetR() < hit2->GetR());
  }
  Bool_t operator() (TPCHitCluster *hitcluster1, TPCHitCluster *hitcluster2) {
    return (hitcluster1->GetR() < hitcluster2->GetR());
  }
};


class SortByDist{
private:
  TVector3 fVtx;
public:
  SortByDist(TVector3 vtx):fVtx(vtx) {}
  Bool_t operator() (TPCHit *hit1, TPCHit *hit2) {
    auto pos1 = hit1->GetPosition();
    auto pos2 = hit2->GetPosition();

    auto dist1 = (pos1-fVtx).Mag();
    auto dist2 = (pos2-fVtx).Mag();
    
    return dist1 > dist2;
  }
  Bool_t operator() (TPCHitCluster *hitcluster1, TPCHitCluster *hitcluster2) {
    auto pos1 = hitcluster1->GetPosition();
    auto pos2 = hitcluster2->GetPosition();

    auto dist1 = (pos1-fVtx).Mag();
    auto dist2 = (pos2-fVtx).Mag();
    
    return dist1 > dist2;
  }
};

class SortByDistInv{
  TVector3 fVtx;
public:
  SortByDistInv(TVector3 vtx):fVtx(vtx) {}
  Bool_t operator() (TPCHit *hit1, TPCHit *hit2) {
    auto pos1 = hit1->GetPosition();
    auto pos2 = hit2->GetPosition();

    auto dist1 = (pos1-fVtx).Mag();
    auto dist2 = (pos2-fVtx).Mag();
    
    return dist1 < dist2;
  }
  Bool_t operator() (TPCHitCluster *hitcluster1, TPCHitCluster *hitcluster2) {
    auto pos1 = hitcluster1->GetPosition();
    auto pos2 = hitcluster2->GetPosition();

    auto dist1 = (pos1-fVtx).Mag();
    auto dist2 = (pos2-fVtx).Mag();
    
    return dist1 < dist2;
  }
};


class SortByIncLength {
private:
  TPCRiemannTrack *trk;

public:
  SortByIncLength(TPCRiemannTrack *track):trk(track) {}
  bool operator() (TPCHit *hit1, TPCHit *hit2) {
    return trk->Map(hit1->GetPosition()).Z() < trk->Map(hit2->GetPosition()).Z();
    //    auto map1 = trk->Map(hit1->GetPosition()) - trk->Map(TVector3(0,0,-143));
    //    auto map2 = trk->Map(hit2->GetPosition()) - trk->Map(TVector3(0,0,-143));
    //    return map1.Z() < map2.Z();
  }
  bool operator() (TPCHitCluster *cluster1, TPCHitCluster *cluster2) {
    return trk->Map(cluster1->GetPosition()).Z() < trk->Map(cluster2->GetPosition()).Z();
    //    auto map1 = trk->Map(cluster1->GetPosition()) - trk->Map(TVector3(0,0,-143));
    //    auto map2 = trk->Map(cluster2->GetPosition()) - trk->Map(TVector3(0,0,-143));
    //    return map1.Z() < map2.Z();

  }
};

class SortByDecLength {
private:
  TPCRiemannTrack *trk;

public:
  SortByDecLength(TPCRiemannTrack *track):trk(track) {}
  bool operator() (TPCHit *hit1, TPCHit *hit2) {
    return trk->Map(hit1->GetPosition()).Z() > trk->Map(hit2->GetPosition()).Z();
    //    auto map1 = trk->Map(hit1->GetPosition()) - trk->Map(TVector3(0,0,-143));
    //    auto map2 = trk->Map(hit2->GetPosition()) - trk->Map(TVector3(0,0,-143));
    //    return map1.Z() > map2.Z();
  }

  bool operator() (TPCHitCluster *cluster1, TPCHitCluster *cluster2) {
    return trk->Map(cluster1->GetPosition()).Z() > trk->Map(cluster2->GetPosition()).Z();
    //    auto map1 = trk->Map(cluster1->GetPosition()) - trk->Map(TVector3(0,0,-143));
    //    auto map2 = trk->Map(cluster2->GetPosition()) - trk->Map(TVector3(0,0,-143));
    //    return map1.Z() > map2.Z();
  }
};

class SortByIncMeanTime {
public:
  SortByIncMeanTime(){}
  bool operator() (TPCHit *hit1, TPCHit *hit2) {
    return hit1->GetCTime() > hit2->GetCTime();
  }

  bool operator() (TPCHitCluster *cluster1, TPCHitCluster *cluster2) {
    return cluster1->GetMeanTime() > cluster2->GetMeanTime();
  }
};

class SortByDecMeanTime {
public:
  SortByDecMeanTime(){}
  bool operator() (TPCHit *hit1, TPCHit *hit2) {
    return hit1->GetCTime() < hit2->GetCTime();
  }

  bool operator() (TPCHitCluster *cluster1, TPCHitCluster *cluster2) {
    return cluster1->GetMeanTime() < cluster2->GetMeanTime();
  }
};

class SortBydEdx {
public:
  SortBydEdx(){}
  
  bool operator() (TPCHitCluster *cluster1, TPCHitCluster *cluster2) {
    return cluster1->GetCharge()/cluster1->GetLength() >  cluster2->GetCharge()/cluster2->GetLength();
  }
};

#endif

  
  
  
    
  
  
  
  
