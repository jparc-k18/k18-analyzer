#include "TPCLinearTrack.hh"
#include "TPCRiemannFitter.hh"

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

namespace {
  static const TVector3 tgtPos(0,0,-143);
}

TPCLinearTrack::TPCLinearTrack() {
  Initialize();
}

TPCLinearTrack::TPCLinearTrack(TPCLinearTrack *track) {

  Initialize();
  
  fTrackID = track->GetTrackID();

  fIsFitted = track->IsFitted();

  fRMSLine = track->GetRMSLine();
  fRMSPlane = track->GetRMSPlane();

  fDirection = track->GetDirection();
  fNormal = track->GetNormal();

  fSlope = track->GetSlope();
  fIntercept = track->GetIntercept();

  fChi2 = track->GetChi2();
  
  fTrackLength = track->GetTrackLength();
  fTrackLengthY = track->GetTrackLengthY();
  fThetaY = track->GetThetaY();

  fdEXZ = track->GetdEXZ();
  

  
  xMean = track->GetMean().X();
  yMean = track->GetMean().Y();
  zMean = track->GetMean().Z();
  
  xxMean = track->GetMeanXX();
  yyMean = track->GetMeanYY();
  zzMean = track->GetMeanZZ();
  
  xyMean = track->GetMeanXY();
  yzMean = track->GetMeanYZ();
  zxMean = track->GetMeanZX();

  fChargeSum = track->GetChargeSum();
  fHits = *(track->GetHits());
  
}

Bool_t TPCLinearTrack::Initialize() {

  fTrackID = -1;

  fIsFitted = false;

  xMean = 0; yMean = 0; zMean = 0;
  xxMean = 0; yyMean = 0; zzMean = 0;
  xyMean = 0; yzMean = 0; zxMean = 0;

  fChargeSum = 0;

  fRMSLine = -9999; fRMSPlane = -9999;

  fDirection = TVector3(-9999,-9999,-9999);
  fNormal = TVector3(-9999,-9999,-9999);



  
  return true;
}


void TPCLinearTrack::FinalizeHits() {

  auto fHitArray = GetHits();
  std::sort(fHitArray->begin(), fHitArray->end(), SortByZ());

  double sumX = 0, sumZ = 0, sumXX = 0, sumZZ = 0, sumXZ = 0;
  double slope, intercept ;
  //Linear Fitting
  // Calculate sums
  for (Int_t iHit = 0 ; iHit < GetNumHits(); iHit++){
    auto eachhit = fHitArray->at(iHit)->GetPosition();
    sumX += eachhit.x();
    sumZ += eachhit.z();
    sumXX += eachhit.x() * eachhit.x();
    sumZZ += eachhit.z() * eachhit.z();
    sumXZ += eachhit.x() * eachhit.z();
  }

  // Calculate means
  double meanX = sumX / GetNumHits();
  double meanZ = sumZ / GetNumHits();

  // Calculate variance
  double varX = 0, varZ = 0;
  double chi2 = 0;
  for (Int_t iHit = 0 ; iHit < GetNumHits(); iHit++){
    auto eachhit = fHitArray->at(iHit)->GetPosition();

    varX += (eachhit.x() - meanX) * (eachhit.x() - meanX);
    varZ += (eachhit.z() - meanZ) * (eachhit.z() - meanZ);


  }



  
  // Calculate the correlation coefficient
  double correlation = (GetNumHits() * sumXZ - sumX * sumZ) / sqrt((GetNumHits() * sumXX - sumX * sumX) * (GetNumHits() * sumZZ - sumZ * sumZ));

  // Calculate the Deming regression parameters
  slope = correlation * sqrt(varX) / sqrt(varZ);
  intercept = meanX - slope * meanZ;
  
  
  
  fSlope = slope;
  fIntercept = intercept;

  for (Int_t iHit = 0 ; iHit < GetNumHits(); iHit++){
    auto eachhit = fHitArray->at(iHit)->GetPosition();


    double predictedX = slope * eachhit.z() + intercept;
    chi2 += pow(eachhit.x() - predictedX, 2) / TMath::Abs(predictedX);
  }

  fChi2 = chi2;

  

  
  auto firstHit = fHitArray->at(0);
  auto firstHitPos = firstHit->GetPosition();
  double firstHitPosZ = firstHitPos.z();
  double firstHitPosX = slope * firstHitPosZ + intercept; //calculate from the result of linear fitting

  SetVertex(0, GetClosestPointOnTrack(firstHitPos));
  
  auto lastHit = fHitArray->back();
  auto lastHitPos = lastHit->GetPosition();
  double lastHitPosZ = lastHitPos.z();
  double lastHitPosX = slope * lastHitPosZ + intercept; //calculate from the result of linear fitting
  SetVertex(1, GetClosestPointOnTrack(lastHitPos));
  
  
  auto first_last = firstHitPos - lastHitPos;
  //fTrackLength = first_last.Mag();
  //Projection to xz plane  (still y is not clear)
  //fTrackLength = TMath::Sqrt(first_last.x()*first_last.x()+first_last.z()*first_last.z()); 
  fTrackLength = TMath::Sqrt( (lastHitPosX - firstHitPosX)*(lastHitPosX - firstHitPosX)+ (lastHitPosZ - firstHitPosZ)*(lastHitPosZ - firstHitPosZ) ); //Calculated

  fTrackLengthY = TMath::Abs(first_last.y());

  fThetaY = TMath::ATan(fTrackLengthY / fTrackLength)*180/TMath::Pi();

  fdEXZ = fChargeSum / fTrackLength; //sum de / xz tracklength

  

}

void TPCLinearTrack::AddHit(TPCHit *hit){

  Double_t x = hit->GetPosition().X();
  Double_t y = hit->GetPosition().Y();
  Double_t z = hit->GetPosition().Z();
  Double_t w = hit->GetDe();
  //Double_t w = hit->GetCDe();
  
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

void TPCLinearTrack::RemoveHit(TPCHit *hit){

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

void TPCLinearTrack::DeleteHits() {
  for (Int_t iHit = 0 ; iHit < GetNumHits(); iHit++)
    delete fHits[iHit];
  fHits.clear();
}

void TPCLinearTrack::AddHitCluster(TPCHitCluster *cluster) {
  fHitClusters.push_back(cluster);
}


TVector3 TPCLinearTrack::PerpLine(TVector3 pos) const {

  TVector3 mean = GetMean();
  TVector3 dir = GetDirection();
  
  TVector3 subPos = pos - mean;
  TVector3 subPosUnit = subPos.Unit();
  Double_t cos = subPosUnit.Dot(dir);
  dir.SetMag(subPos.Mag()*cos);

  return dir - subPos;
}

TVector3 TPCLinearTrack::PerpPlane(TVector3 pos) const {

  TVector3 normal = GetNormal();
  TVector3 mean = GetMean();

  Double_t perp = abs(normal * pos - normal * mean) / sqrt(normal * normal);

  return perp * normal;
}

Double_t TPCLinearTrack::PerpDistLine(TVector3 pos){

  TVector3 direction = GetDirection();
  TVector3 mean = GetMean();

  Double_t xC = (direction.Z() * (pos.Y() - mean.Y()) - direction.Y() * (pos.Z() - mean.Z()));
  Double_t yC = (direction.X() * (pos.Z() - mean.Z()) - direction.Z() * (pos.X() - mean.X()));
  Double_t zC = (direction.Y() * (pos.X() - mean.X()) - direction.X() * (pos.Y() - mean.Y()));

  return sqrt(xC*xC + yC*yC + zC*zC);
}


TVector3 TPCLinearTrack::GetPointOnX(Double_t x) {

  TVector3 direction = GetDirection();
  TVector3 mean = GetMean();

  Double_t d = (x-mean.X())/direction.X();

  Double_t y = d * direction.Y() + mean.Y();
  Double_t z = d * direction.Z() + mean.Z();

  return TVector3(x,y,z);
}
  
TVector3 TPCLinearTrack::GetPointOnY(Double_t y) {

  TVector3 direction = GetDirection();
  TVector3 mean = GetMean();

  Double_t d = (y-mean.Y())/direction.Y();

  Double_t x = d * direction.X() + mean.X();
  Double_t z = d * direction.Z() + mean.Z();

  return TVector3(x,y,z);
}
  
TVector3 TPCLinearTrack::GetPointOnZ(Double_t z) {

  TVector3 direction = GetDirection();
  TVector3 mean = GetMean();

  Double_t d = (z-mean.Z())/direction.Z();

  Double_t x = d * direction.X() + mean.X();
  Double_t y = d * direction.Y() + mean.Y();

  return TVector3(x,y,z);
}
  





