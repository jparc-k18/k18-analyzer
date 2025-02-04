#include "TPCRiemannTrack.hh"
#include "TPCRiemannTrackSearch.hh"
#include "TPCRiemannFitter.hh"

#include "Kinematics.hh"
#include "DatabasePDG.hh"

#include <iostream>
#include <iterator>

#include <TMatrixD.h>
#include <TVectorD.h>
#include <TVector3.h>
#include <TVector2.h>
#include <TGraph.h>
#include <TF1.h>
#include <TMath.h>
#include <TPolyMarker3D.h>
#include <TCanvas.h>
#include <TApplication.h>
#include <TPolyLine3D.h>
#include <TSystem.h>

#include "ConfMan.hh"

namespace {
  static const TVector3 tgtPos(0,0,-143);
  //  const UserParamMan& gUser = UserParamMan::GetInstance();
  //  const DCGeomMan& gGeom = DCGeomMan::GetInstance();
}

TPCRiemannTrack::TPCRiemannTrack() {
  Initialize();
}

TPCRiemannTrack::TPCRiemannTrack(Int_t trackid)
{
  Initialize();
  fTrackID = trackid;
}

Bool_t TPCRiemannTrack::Initialize() {

  fTrackID = -9999;
  fPID = -9999;
  fCenter = TVector3(0.,0.,0.);
  fRadius = 0;
  fChargeSum = 0;
  fSlope = 0;
  fOffset = 0;
  fVertex = TVector3(-999,-999,-999);
  
  fIsFitted = false;
  fIsInitialized = false;
  fIsGood = false;
  fDoSort = true;
  fIsFinished = false;
  fIsPositive = true;
  
  fAlpha = -9999;
  fAlphaHead = -9999;
  fAlphaTail = -9999;

  xMean = 0; yMean = 0; zMean = 0;
  xxMean = 0; yyMean = 0; zzMean = 0;
  xyMean = 0; yzMean = 0; zxMean = 0;

  fRMSW = -9999;
  fRMSH = -9999;
  
  fFitStatus = kBad;
  
  return true;
}

void TPCRiemannTrack::AddHit(TPCHit *hit){

  Double_t x = hit->GetPosition().X();
  Double_t y = hit->GetPosition().Y();
  Double_t z = hit->GetPosition().Z();
  //  Double_t w = hit->GetDe();
  Double_t w = hit->GetCDe();
  
  Double_t W = fChargeSum + w;

  xMean = ( fChargeSum * xMean + w * x) / W;
  yMean = ( fChargeSum * yMean + w * y) / W;
  zMean = ( fChargeSum * zMean + w * z) / W;

  xxMean = ( fChargeSum * xxMean + w * x * x) / W;
  yyMean = ( fChargeSum * yyMean + w * y * y) / W;
  zzMean = ( fChargeSum * zzMean + w * z * z) / W;
  
  xyMean = ( fChargeSum * xyMean + w * x * y) / W;
  yzMean = ( fChargeSum * yzMean + w * y * z) / W;
  zxMean = ( fChargeSum * zxMean + w * z * x) / W;

  fChargeSum = W;

  fHits.push_back(hit);
}

void TPCRiemannTrack::RemoveHit(TPCHit *hit){

  Double_t x = hit->GetPosition().X();
  Double_t y = hit->GetPosition().Y();
  Double_t z = hit->GetPosition().Z();

  //  Double_t w = hit->GetDe();
  Double_t w = hit->GetCDe();

  Double_t W = fChargeSum - w;
  
  for (int i = 0 ; i < fHits.size() ; i++ ) {
    if (fHits[i] == hit){
      fHits.erase(fHits.begin() + i);
      break;
    }
  }

  xMean = ( fChargeSum * xMean - w * x) / W;
  yMean = ( fChargeSum * yMean - w * y) / W;
  zMean = ( fChargeSum * zMean - w * z) / W;
  
  xxMean = ( fChargeSum * xxMean - w * x * x) / W;
  yyMean = ( fChargeSum * yyMean - w * y * y) / W;
  zzMean = ( fChargeSum * zzMean - w * z * z) / W;
  
  xyMean = ( fChargeSum * xyMean - w * x * y) / W;
  yzMean = ( fChargeSum * yzMean - w * y * z) / W;
  zxMean = ( fChargeSum * zxMean - w * z * x) / W;

  fChargeSum = W;
}


void TPCRiemannTrack::FinalizeClusters() {

  /*
  Double_t alphaHead = GetAlphaHead();
  Double_t alphaTail = GetAlphaTail();
  
  TVector3 HeadPos = GetPosition(fAlphaHead);
  TVector3 TailPos = GetPosition(fAlphaTail);

  Double_t distHead = (HeadPos - tgtPos).Mag();
  Double_t distTail = (TailPos - tgtPos).Mag();
  
  if (distTail > distHead){
    Double_t val = (fOffset + fSlope*fAlphaHead) - (fOffset - fSlope*fAlphaTail);
    fAlphaHead = alphaTail;
    fAlphaTail = alphaHead;
    fSlope *= -1;
    fOffset += val;    
  }
  */

  TVector3 ptOnHelix, mappedPt;
  Double_t closestDist;
  Double_t poca;

  poca = ExtrapByMap(tgtPos, ptOnHelix, mappedPt);
  closestDist = (tgtPos - ptOnHelix).Mag();

  TVector3 dirTgt = Direction(AlphaAtPosition(ptOnHelix));
  //  dirTgt *= winding;
  SetPOCA(ptOnHelix);
  SetClosestDist(closestDist);
  SetDirectionAtTgt(dirTgt);
  

  auto maxdx = 0;
  auto maxdy = 0;
  auto maxdz = 0;
  
  for (auto cluster : fHitClusters) {

    TVector3 dir = Direction(AlphaAtPosition(cluster->GetPosition()));
    //    std::cout << dir.X() <<"\t" <<  dir.Y() <<"\t" <<  dir.Z() << std::endl;
    //    if (dirTgt.Z() < 0) dir *= -1;
    
    auto chi = TMath::ATan2((dir.X()), (dir.Z()));
    //    cluster->SetChi(TMath::PiOver2()-chi);
    cluster->SetChi(chi);
    cluster->SetLambda(GetDip());
    cluster->SetDirection(dir.Unit());

    //    std::cout << cluster->GetPosition().X() <<"\t" <<  cluster->GetPosition().Y() <<"\t" <<  cluster->GetPosition().Z() << std::endl;
    cluster->SetTrackID(fTrackID);

    //    SetClusterLength(cluster);
    //    if (cluster->IsStable()) {
    fdEdxArray.push_back(cluster->GetCharge()/cluster->GetLength());
    fMeandEdxArray.push_back(cluster->GetMeanCharge()/cluster->GetLength());

    auto cov = cluster->GetCovMatrix();
    auto charge = cluster->GetCharge();

    auto dx = sqrt(abs(cov(0,0)));
    auto dy = sqrt(abs(cov(1,1)));
    auto dz = sqrt(abs(cov(2,2)));


    //    if (dx > maxdx) maxdx = dx;
    //    if (dy > maxdx) maxdy = dy;
    //    if (dz > maxdx) maxdz = dz;
    auto nhit = cluster->GetNumHits();

   
    //    if (dx <= 0 || dy <= 0 || dz <= 0 || nhit <= 1)

    //    if (dx <= 0 || dy <= 0 || dz <= 0)
    
    //    if (nhit <= 1)
    cluster->SetDFromCov(maxdx, maxdy, maxdz);

  }

  Double_t alphaHead = GetAlphaHead();
  Double_t alphaTail = GetAlphaTail();
  
  TVector3 HeadPos = GetPosition(fAlphaHead);
  TVector3 TailPos = GetPosition(fAlphaTail);

  Double_t distHead = (HeadPos - tgtPos).Mag();
  Double_t distTail = (TailPos - tgtPos).Mag();
  

  //  if (distTail > distHead){
  //    std::sort(fHitClusters.begin(), fHitClusters.end(), SortByIncLength(this));
  //  }
  //  else {
  //    std::sort(fHitClusters.begin(), fHitClusters.end(), SortByDecLength(this));
  //  }

  
  std::sort(fHitClusters.begin(), fHitClusters.end(), SortByDistInv(tgtPos));  

    //  if (Helicity() == 1)
  //    //    std::sort(fHitClusters.begin(), fHitClusters.end(), SortByYInv());  
  //  else
  //    std::sort(fHitClusters.begin(), fHitClusters.end(), SortByIncLength(this));  
    //    std::sort(fHitClusters.begin(), fHitClusters.end(), SortByY());  

  //  Double_t dEdx = GetTruncateddEdx(0,0.7);
  //  Double_t dEdxcalc_pi = Kinematics::CalcTPCDedx(1,1000*pdg::PionMass(), GetMomentum());
  
}

Bool_t TPCRiemannTrack::IsAccidental() {

  Double_t momentum = GetMomentum();

  if (GetNumHits() > 20 && momentum > 1600 && TrackLength() > 400 )
    return true;

  return false;
}

void TPCRiemannTrack::AddHitCluster(TPCHitCluster *cluster) {
  fHitClusters.push_back(cluster);
}

TVector3 TPCRiemannTrack::PerpLine(TVector3 pos) const {

  TVector3 mean = GetMean();
  TVector3 dir = GetLineDirection();

  TVector3 subPos = pos - mean;
  TVector3 subPosUnit = subPos.Unit();
  Double_t cos = subPosUnit.Dot(dir);
  dir.SetMag(subPos.Mag()*cos);

  return dir - subPos;
}

TVector3 TPCRiemannTrack::PerpPlane(TVector3 pos) const {

  TVector3 normal = GetPlaneNormal();
  TVector3 mean = GetMean();

  Double_t perp = abs(normal * pos - normal * mean) / sqrt(normal * normal);

  return perp * normal;
}


TVector3 TPCRiemannTrack::Map(TVector3 pos) const {
  TVector3 q, m;
  ExtrapByMap(pos,q,m);

  return m;
}

Double_t TPCRiemannTrack::ExtrapByMap(TVector3 pos, TVector3 &q, TVector3 &m) const {

  Double_t head = ExtrapToAlpha(fAlphaHead);
  Double_t tail = ExtrapToAlpha(fAlphaTail);
  Double_t off = head;

  if (head > tail)
    off = tail;

  Double_t alpha;
  Double_t length = ExtrapToPointAlpha(pos,q,alpha);
  Double_t r = DistCircle(pos);
  Double_t y = pos.Y() - q.Y();

  m = TVector3(r, y/TMath::Cos(GetDip()), length + y*TMath::Sin(GetDip()) - off);

  return alpha * fRadius/TMath::Cos(GetDip());
}

Double_t TPCRiemannTrack::ExtrapToPointAlpha(TVector3 pos, TVector3 &pt, Double_t &alpha) const{

  Double_t alpha0 = TMath::ATan2(pos.X()-fCenter.X(), pos.Z()-fCenter.Z());
  TVector3 pt0(fRadius * TMath::Sin(alpha0) + fCenter.X(),
	       alpha0*fSlope + fOffset,
	       fRadius * TMath::Cos(alpha0) + fCenter.Z());
  Double_t y0 = std::abs(pt0.Y()-pos.Y());

  Double_t y1;
  Double_t alpha1 = alpha0;
  TVector3 pt1 = pt0;

  Double_t lengthY = std::abs(2*TMath::Pi() * fSlope);
  //  if (lengthY > 3 * fRMSH && lengthY > 5 && std::abs(GetDip()) < 1.5) {
  if (1) {

    Int_t count = 0;

    while(1) {
      alpha1 = alpha1 + 2*TMath::Pi();
      pt1.SetY(pt1.Y() + 2*TMath::Pi()*fSlope);
      y1 = std::abs(pt1.Y() - pos.Y());

      if (std::abs(y0) < std::abs(y1))
      //      if (y0-y1 < 1.0e+10)
	break;

      else {
	alpha0 = alpha1;
	pt0 = pt1;
	y0 = y1;
      }

      if (count++ > 50) break;
    }

    y1 = y0;
    alpha1 = alpha0;
    pt1 = pt0;

    count = 0;

    while(1) {
      alpha1 = alpha1 - 2*TMath::Pi();
      pt1.SetY(pt1.Y() - 2*TMath::Pi()*fSlope);
      y1 = std::abs(pt1.Y()-pos.Y());

      if (std::abs(y0) < std::abs(y1))
      //      if (y0-y1 < 1.0e+10)
	break;

      else {
	alpha0 = alpha1;
	pt0 = pt1;
	y0 = y1;
      }
      if (count++ > 50)
	break;
    }
  }

  pt = pt0;
  alpha = alpha0;
  Double_t length = alpha0 * fRadius / TMath::Cos(GetDip());

  return length;
}

Bool_t TPCRiemannTrack::ExtrapToZ(Double_t z, Double_t alphaRef, TVector3 &pointOnHelix) const {

  if (CheckExtrapToZ(z) == false)
    return false;

  Double_t xOff = sqrt(fRadius*fRadius - ( z - fCenter.Z()) * ( z - fCenter.Z()));

  Double_t x1 = fCenter.X() + xOff;
  Double_t x2 = fCenter.X() - xOff;

  Double_t alpha = TMath::ATan2(x1-fCenter.X(), z-fCenter.Z());
  Double_t alphaTemp = alpha;
  Double_t d1Cand = std::abs(alphaTemp-alphaRef);
  Double_t d1Temp = d1Cand;

  while(1) {
    alphaTemp = alpha + 2*TMath::Pi();
    d1Temp = std::abs(alphaTemp-alphaRef);
    if (d1Temp >= d1Cand)
      break;
    else {
      alpha = alphaTemp;
      d1Cand = d1Temp;
    }
  }

  while(1) {
    alphaTemp = alpha - 2*TMath::Pi();
    d1Temp = std::abs(alphaTemp-alphaRef);
    if (d1Temp >= d1Cand)
      break;
    else {
      alpha = alphaTemp;
      d1Cand = d1Temp;
    }
  }
  pointOnHelix = GetPosition(alpha);

  return true;  
}

Bool_t TPCRiemannTrack::ExtrapToZ(Double_t z, TVector3 &pointOnHelix) const {

  TVector3 position1, position2;
  Double_t alpha1, alpha2;

  if (ExtrapToZ(z, position1, alpha1, position2, alpha2) == false)
    return false;

  Double_t alphaMid = (fAlphaHead + fAlphaTail)/2.;

  if (std::abs(alpha1-alphaMid) < std::abs(alpha2-alphaMid))
    pointOnHelix = position1;
  else
    pointOnHelix = position2;

  return true;
}

Bool_t TPCRiemannTrack::ExtrapToZ(Double_t z, TVector3 &pointOnHelix1, Double_t &alpha1,
				  TVector3 &pointOnHelix2, Double_t &alpha2) const {
  if (CheckExtrapToZ(z) == false)
    return false;

  Double_t xOff = sqrt(fRadius*fRadius - ( z - fCenter.Z()) * ( z - fCenter.Z()));

  Double_t x1 = fCenter.X() + xOff;
  Double_t x2 = fCenter.X() - xOff;

  alpha1 = TMath::ATan2(x1-fCenter.X(), z-fCenter.Z());
  Double_t alpha1Temp = alpha1;
  Double_t d1Cand = std::abs(alpha1Temp-fAlphaHead);
  Double_t d1Temp = d1Cand;

  while(1) {
    alpha1Temp = alpha1 + 2*TMath::Pi();
    d1Temp = std::abs(alpha1Temp-fAlphaHead);
    if (d1Temp >= d1Cand)
      break;
    else {
      alpha1 = alpha1Temp;
      d1Cand = d1Temp;
    }
  }

  while(1) {
    alpha1Temp = alpha1 - 2*TMath::Pi();
    d1Temp = std::abs(alpha1Temp-fAlphaHead);
    if (d1Temp >= d1Cand)
      break;
    else {
      alpha1 = alpha1Temp;
      d1Cand = d1Temp;
    }
  }
  pointOnHelix1 = GetPosition(alpha1);

  alpha2 = TMath::ATan2(x2-fCenter.X(), z-fCenter.Z());
  Double_t alpha2Temp = alpha2;
  Double_t d2Cand = std::abs(alpha2Temp-fAlphaTail);
  Double_t d2Temp = d2Cand;

  while(1) {
    alpha2Temp = alpha2 + 2*TMath::Pi();
    d2Temp = std::abs(alpha2Temp-fAlphaTail);
    if (d2Temp >= d2Cand)
      break;
    else {
      alpha2 = alpha2Temp;
      d2Cand = d2Temp;
    }
  }

  while(1) {
    alpha2Temp = alpha2 - 2*TMath::Pi();
    d2Temp = std::abs(alpha2Temp-fAlphaTail);
    if (d2Temp >= d2Cand)
      break;
    else {
      alpha2 = alpha2Temp;
      d2Cand = d2Temp;
    }
  }
  pointOnHelix2 = GetPosition(alpha2);


  return true;  
}

Bool_t TPCRiemannTrack::CheckExtrapToZ(Double_t z) const {
  Double_t zRef = fCenter.Z() - z;
  Double_t mult = (zRef + fRadius) * (zRef - fRadius);

  if (mult > 0)
    return false;

  return true;
}

Bool_t TPCRiemannTrack::ExtrapToX(Double_t x, Double_t alphaRef, TVector3 &pointOnHelix) const {

  if (CheckExtrapToX(x) == false)
    return false;

  Double_t zOff = sqrt(fRadius*fRadius - ( x - fCenter.X()) * ( x - fCenter.X()));

  Double_t z1 = fCenter.Z() + zOff;
  Double_t z2 = fCenter.Z() - zOff;

  Double_t alpha = TMath::ATan2(x-fCenter.X(), z1-fCenter.Z());
  Double_t alphaTemp = alpha;
  Double_t d1Cand = std::abs(alphaTemp-alphaRef);
  Double_t d1Temp = d1Cand;

  while(1) {
    alphaTemp = alpha + 2*TMath::Pi();
    d1Temp = std::abs(alphaTemp-alphaRef);
    if (d1Temp >= d1Cand)
      break;
    else {
      alpha = alphaTemp;
      d1Cand = d1Temp;
    }
  }

  while(1) {
    alphaTemp = alpha - 2*TMath::Pi();
    d1Temp = std::abs(alphaTemp-alphaRef);
    if (d1Temp >= d1Cand)
      break;
    else {
      alpha = alphaTemp;
      d1Cand = d1Temp;
    }
  }
  pointOnHelix = GetPosition(alpha);

  return true;  
}

Bool_t TPCRiemannTrack::ExtrapToX(Double_t x, TVector3 &pointOnHelix) const {

  TVector3 position1, position2;
  Double_t alpha1, alpha2;

  if (ExtrapToX(x, position1, alpha1, position2, alpha2) == false)
    return false;

  Double_t alphaMid = (fAlphaHead + fAlphaTail)/2.;

  if (std::abs(alpha1-alphaMid) < std::abs(alpha2-alphaMid))
    pointOnHelix = position1;
  else
    pointOnHelix = position2;

  return true;
}

Bool_t TPCRiemannTrack::ExtrapToX(Double_t x, TVector3 &pointOnHelix1, Double_t &alpha1,
				  TVector3 &pointOnHelix2, Double_t &alpha2) const {
  if (CheckExtrapToX(x) == false)
    return false;

  Double_t zOff = sqrt(fRadius*fRadius - ( x - fCenter.X()) * ( x - fCenter.X()));

  Double_t z1 = fCenter.Z() + zOff;
  Double_t z2 = fCenter.Z() - zOff;

  alpha1 = TMath::ATan2(x-fCenter.X(), z1-fCenter.Z());
  Double_t alpha1Temp = alpha1;
  Double_t d1Cand = std::abs(alpha1Temp-fAlphaHead);
  Double_t d1Temp = d1Cand;

  while(1) {
    alpha1Temp = alpha1 + 2*TMath::Pi();
    d1Temp = std::abs(alpha1Temp-fAlphaHead);
    if (d1Temp >= d1Cand)
      break;
    else {
      alpha1 = alpha1Temp;
      d1Cand = d1Temp;
    }
  }

  while(1) {
    alpha1Temp = alpha1 - 2*TMath::Pi();
    d1Temp = std::abs(alpha1Temp-fAlphaHead);
    if (d1Temp >= d1Cand)
      break;
    else {
      alpha1 = alpha1Temp;
      d1Cand = d1Temp;
    }
  }
  pointOnHelix1 = GetPosition(alpha1);

  alpha2 = TMath::ATan2(x-fCenter.X(), z2-fCenter.Z());
  Double_t alpha2Temp = alpha2;
  Double_t d2Cand = std::abs(alpha2Temp-fAlphaTail);
  Double_t d2Temp = d2Cand;

  while(1) {
    alpha2Temp = alpha2 + 2*TMath::Pi();
    d2Temp = std::abs(alpha2Temp-fAlphaTail);
    if (d2Temp >= d2Cand)
      break;
    else {
      alpha2 = alpha2Temp;
      d2Cand = d2Temp;
    }
  }

  while(1) {
    alpha2Temp = alpha2 - 2*TMath::Pi();
    d2Temp = std::abs(alpha2Temp-fAlphaTail);
    if (d2Temp >= d2Cand)
      break;
    else {
      alpha2 = alpha2Temp;
      d2Cand = d2Temp;
    }
  }
  pointOnHelix2 = GetPosition(alpha2);


  return true;  
}

Bool_t TPCRiemannTrack::CheckExtrapToX(Double_t x) const {
  Double_t xRef = fCenter.X() - x;
  Double_t mult = (xRef + fRadius) * (xRef - fRadius);
  
  if (mult > 0)
    return false;
  
  return true;
}


Double_t TPCRiemannTrack::AlphaAtPosition(TVector3 pos) {

  Double_t alpha;
  TVector3 q(0,0,0);
  ExtrapToPointAlpha(pos, q, alpha);

  return alpha;
    
}

Double_t TPCRiemannTrack::GetMomentum(Bool_t isConst) const {

  if (fFitStatus != TPCRiemannTrack::kHelix)
    return -1;

  const Double_t &HSField0 = 0.986;
  //  const Double_t &HSField0 = 1.0;
  const Double_t &HSFieldCalc = ConfMan::Get<Double_t>("HSFLDCALC");
  const Double_t &HSFieldHall = ConfMan::Get<Double_t>("HSFLDHALL");
  const Double_t &HSFieldCalib = ConfMan::Get<Double_t>("HSFLDCALIB");
  
  Double_t MagneticField;
  if (!isConst){
    //    MagneticField = HSField0*(HSFieldHall/HSFieldCalc);}
    MagneticField = HSField0*(HSFieldHall*HSFieldCalib/HSFieldCalc);}
  else {
    MagneticField = 0.9;
  }
  
  Double_t cosDip = TMath::Cos(GetDip());
  Double_t coeff = 0.299792458;
  if (cosDip < 1e-02)
    return TMath::Abs(fRadius/1e-02 * coeff * MagneticField);

  return  TMath::Abs(fRadius/cosDip * coeff * MagneticField);
}


Double_t TPCRiemannTrack::GetMomResolution(Double_t ResT) const {

  if (fFitStatus != TPCRiemannTrack::kHelix)
    return -1;

  const Double_t &HSField0 = 0.986;
  //  const Double_t &HSField0 = 1.0;
  const Double_t &HSFieldCalc = ConfMan::Get<Double_t>("HSFLDCALC");
  const Double_t &HSFieldHall = ConfMan::Get<Double_t>("HSFLDHALL");
  const Double_t &HSFieldCalib = ConfMan::Get<Double_t>("HSFLDCALIB");

  ResT *= 1e-03;
  
  Double_t MagneticField;
  MagneticField = HSField0*(HSFieldHall*HSFieldCalib/HSFieldCalc);

  Double_t cosDip = TMath::Cos(GetDip());
  Double_t p = GetMomentum(false)/1e+3;
  Double_t pt;

  Double_t coeff = 0.299792458;

  if (cosDip < 1e-02)
    pt = p*1e-02;
  else
    pt = p*cosDip;

  Double_t L2 = TrackLength()*TrackLength()/1e+6;
  Double_t val = (ResT*pt*pt)/(coeff*MagneticField*L2) *  sqrt(720/(GetNumStableClusters()+4));
  
  return val;
    
}

TVector3 TPCRiemannTrack::ExtrapHead(Double_t length) const {

  Double_t alpha = fAlphaHead;
  Double_t dAlpha = std::abs(length*TMath::Cos(GetDip())/fRadius);

  if (fAlphaHead > fAlphaTail)
    alpha += dAlpha;
  else
    alpha -= dAlpha;

  return GetPosition(alpha);
}
			     
TVector3 TPCRiemannTrack::ExtrapTail(Double_t length) const {

  Double_t alpha = fAlphaTail;
  Double_t dAlpha = std::abs(length*TMath::Cos(GetDip())/fRadius);

  if (fAlphaHead < fAlphaTail)
    alpha += dAlpha;
  else
    alpha -= dAlpha;
  
  return GetPosition(alpha);
}
			       
void TPCRiemannTrack::DeleteHits() {
  for (Int_t iHit = 0 ; iHit < GetNumHits(); iHit++)
    delete fHits[iHit];
  fHits.clear();
}
	  
TVector3 TPCRiemannTrack::InterpByLength(Double_t length) const {

  TVector3 q;
  Double_t ratio = length/TrackLength();
  ExtrapToAlpha(ratio*fAlphaHead + (1-ratio)*fAlphaTail,q);

  return q;
}

TVector3 TPCRiemannTrack::InterpByRatio(Double_t ratio) const {

  TVector3 q;
  ExtrapToAlpha(ratio*fAlphaHead + (1-ratio)*fAlphaTail,q);

  return q;
}

Double_t TPCRiemannTrack::Continuity() {
  auto numHits = fHits.size();
  if (numHits < 2) return -1;

  std::sort(fHits.begin(), fHits.end(), SortByIncLength(this));

  Double_t total = 0 ;
  Double_t continuous = 0;
  TVector3 before = Map(fHits[0]->GetPosition());

  for (auto iHit = 1 ; iHit < numHits ; iHit++) {
    TVector3 current = Map(fHits[iHit]->GetPosition());
    auto length = std::abs(current.Z()-before.Z());
    std::cout << "Length : " << length << std::endl;
    total += length;
    if (length < 20)
      continuous += length;

    before = current;
  }

  return continuous/total;
}
  
void TPCRiemannTrack::DetermineCharge(TVector3 vtx) {
  
  fVertex = vtx;
  TVector3 mean = GetMean();
  
  TVector3 TrackHead = PositionAtHead();
  TVector3 TrackTail = PositionAtTail();
  
  Bool_t reverse = false;

  Double_t distHead = (TrackHead-fVertex).Mag();
  Double_t distTail = (TrackTail-fVertex).Mag();

  if (distTail > distHead) reverse = true;
  
  Double_t head, tail;

  //  if (fAlphaTail > fAlphaHead){
  head = ExtrapToAlpha(fAlphaHead);
  tail = ExtrapToAlpha(fAlphaTail);
  //  }
    //  else {
    //    head = ExtrapToAlpha(fAlphaTail);
    //    tail = ExtrapToAlpha(fAlphaHead);
    //  }    

  TVector3 q;
  Double_t alpha;
  Double_t Vertex = ExtrapToPointAlpha(fVertex,q,alpha);

    //  if (mean.Z() > -143.){
  if (abs(Vertex-tail) > abs(Vertex-head)){
    fIsPositive = true;
  }
  else {
    fIsPositive = false;
  }
  
}
/*
TVector3 TPCRiemannTrack::Direction(Double_t alpha) const
{
  
  TVector3 mean = GetMean();
  
  Double_t alphaTmp = alpha;

  if ((fAlphaTail) < (fAlphaHead))
    alphaTmp -= TMath::PiOver2();
  else 
    alphaTmp += TMath::PiOver2();

  TVector3 direction = GetPosition(alphaTmp) - GetCenter();
  direction.SetY(fSlope);

  if (direction.Z() < 0 && GetTrackProperty() == 1) {
    direction *= -1;
    return direction.Unit();
  }

  TVector3 TrackHead = PositionAtHead();
  TVector3 TrackTail = PositionAtTail();

  Double_t distHead = (TrackHead-tgtPos).Mag();
  Double_t distTail = (TrackTail-tgtPos).Mag();
  
  //  if (direction.Z() < 0) direction *= -1;
  
  if (distHead < distTail) {
    //    if (direction.Z() < 0) 
    //      direction *= -1;
  }

  Double_t head, tail;
  TVector3 q;
  head = ExtrapToAlpha(fAlphaHead);
  tail = ExtrapToAlpha(fAlphaTail);

  Double_t Vertex = ExtrapToPointAlpha(tgtPos,q,alphaTmp);

  //  if (abs(Vertex-tail) < abs(Vertex-head) && distHead > distTail)
  //    direction *= -1;
  
  if (abs(Vertex-tail) < abs(Vertex-head)){
    if (direction.Z() < 0)
      direction *= -1;
  }
  
  if (TrackHead.Z() < -143. || TrackTail.Z() < -143){
    if (GetTrackProperty() == 2 && mean.Z() < -135. ){
      //    if (GetTrackProperty() == 2){
      //    if (direction.Z() > 0) direction *= -1;
      //      direction *= -1;
    }
  }
  return direction.Unit();
  
}
*/

TVector3 TPCRiemannTrack::Direction(Double_t alpha) const
{
  
  TVector3 mean = GetMean();
  
  Double_t alphaTmp = alpha;

  if ((fAlphaTail) < (fAlphaHead))
    alphaTmp += TMath::PiOver2();
  else 
    alphaTmp -= TMath::PiOver2();

  TVector3 direction = GetPosition(alphaTmp) - GetCenter();
  auto directionZ = direction.Z();
  auto dirY = direction.Y();
  direction.SetY(0);
  //  direction.SetMag(fRadius*TMath::PiOver2());
  direction.SetMag(2.*fRadius*TMath::Pi());

  Double_t yLength = abs(2.*TMath::Pi()*fSlope);

  if (dirY > 0){
    //    direction *= -1;
    direction.SetY(yLength);
  }
  else {
    direction.SetY(-yLength);
  }

  if (direction.Z() < 0 && GetTrackProperty() == 1) {
    direction *= -1;
    return direction.Unit();
  }
  
  TVector3 TrackHead = PositionAtHead();
  TVector3 TrackTail = PositionAtTail();

  Double_t distHead = (TrackHead-tgtPos).Mag();
  Double_t distTail = (TrackTail-tgtPos).Mag();
  
  //  if (direction.Z() < 0) direction *= -1;
  
  if (distHead < distTail) {

    direction *= -1;
    
  }

  //  if ((fAlphaTail) < (fAlphaHead))
  //    direction *= -1;

  if (TrackHead.Z() < -143. && TrackTail.Z() < -143){
    if (GetTrackProperty() == 2 && mean.Z() < -135. ){
      //      if (direction.Z() > 0) direction *= -1;
    }
  }

  return direction.Unit();
  
}

  

Int_t TPCRiemannTrack::GetWinding() const {

  Double_t angleback = fHitClusters.back() -> GetAlpha();
  Double_t anglefront = fHitClusters.front() -> GetAlpha();
  Double_t angleav = 0.5*(angleback + anglefront);

  Double_t angle1;
  
  TVector3 hit(0,0,-143.);
  hit -= fCenter;
  //  angle1 = hit.Phi();
  angle1 = TMath::ATan2(hit.Z(), hit.X());

  if (angle1 > angleav)
    return 1;

  return -1;
}

Double_t TPCRiemannTrack::GetTruncateddEdx(Double_t LL, Double_t HL) {

  std::sort(fdEdxArray.begin(), fdEdxArray.end());

  auto numPoints = fdEdxArray.size();

  Int_t idxLow = Int_t(numPoints * LL);
  Int_t idxHigh = Int_t(numPoints * HL);

  numPoints = idxHigh - idxLow;
  if (numPoints < 3)
    return -1;

  Double_t dEdx = 0;
  for (Int_t idEdx = idxLow ; idEdx < idxHigh ; idEdx++) 
    dEdx += fdEdxArray[idEdx];
  dEdx = dEdx/numPoints;

  return dEdx;
}

Double_t TPCRiemannTrack::GetTruncatedMeandEdx(Double_t LL, Double_t HL) {

  std::sort(fMeandEdxArray.begin(), fMeandEdxArray.end());

  auto numPoints = fMeandEdxArray.size();

  Int_t idxLow = Int_t(numPoints * LL);
  Int_t idxHigh = Int_t(numPoints * HL);

  numPoints = idxHigh - idxLow;
  if (numPoints < 3)
    return -1;

  Double_t dEdx = 0;
  for (Int_t idEdx = idxLow ; idEdx < idxHigh ; idEdx++) 
    dEdx += fMeandEdxArray[idEdx];
  dEdx = dEdx/numPoints;

  return dEdx;
}
