#include "TPCHitCluster.hh"
#include "TPCPadHelper.hh"
#include "TPCHit.hh"
#include "TMath.h"
#include "UserParamMan.hh"
#include "TPCParamMan.hh"
#include <iostream>

namespace {
  const UserParamMan& gUser = UserParamMan::GetInstance();
  const TPCParamMan& gTPC    = TPCParamMan::GetInstance();

  const auto& ResParamInnerLayerHSOn = gTPC.TPCResolutionParams(true, false); //B=1 T, Inner layers
  const auto& ResParamOuterLayerHSOn = gTPC.TPCResolutionParams(true, true); //B=1 T, Outer layers
  const auto& ResParamInnerLayerHSOff = gTPC.TPCResolutionParams(false, false); //B=0, Inner layers
  const auto& ResParamOuterLayerHSOff = gTPC.TPCResolutionParams(false, true); //B=0, Outer layers
}

TPCHitCluster::TPCHitCluster() {
  Clear();
}

TPCHitCluster::TPCHitCluster(TPCHitCluster *cluster) {
  fIsClusterised = cluster->IsClusterised();

  fX = cluster->GetX();
  fY = cluster->GetY();
  fZ = cluster->GetZ();

  fDx = cluster->GetDx();
  fDy = cluster->GetDy();
  fDz = cluster->GetDz();
 
  fCovMatrix.ResizeTo(3,3);
  fCovMatrix = cluster->GetCovMatrix();

  fCharge = cluster->GetCharge();
  fChargeSqr = cluster->GetChargeSqr();  
  fMeanCharge = cluster->GetMeanCharge();

  fHits = (*cluster->GetHits());
		  
  SetPOCA(cluster->GetPOCA());

  fIsContinuousHits = cluster -> IsContinuousHits();

  fMeanTime = cluster->GetMeanTime();

}

void TPCHitCluster::Clear(Option_t *)
{
  fX = 0;
  fY = 0;
  fZ = 0;

  fDx = 0;
  fDy = 0;
  fDz = 0;

  fCharge = 0;
  fMeanCharge = 0;
  fChargeSqr = 0;

  fCovMatrix.ResizeTo(3,3);
  for (Int_t i = 0 ; i < 9 ; i++)
    fCovMatrix(i/3,i%3) = 0;

  fPOCAx = -9999;
  fPOCAy = -9999;
  fPOCAz = -9999;

  fS = 0 ;

  fIsClusterised = kFALSE;
  fIsContinuousHits = kFALSE;

  xMean  =  0.;
  yMean  =  0.;
  zMean  =  0.;
  
  xxMean  =  0.;
  yyMean  =  0.;
  zzMean  =  0.;
  
  xyMean  =  0.;
  xzMean  =  0.;
  yzMean  =  0.;

  fHits.clear();

  fMeanTime = 0.;
}

void TPCHitCluster::SetCovMatrix(TMatrixD matrix) {
  fCovMatrix = matrix;
}

void TPCHitCluster::AddHit(TPCHit *hit) {

  //  auto charge = hit->GetDe();
  auto charge = hit->GetCDe();
  auto chargeSum = fCharge + charge;

  auto chargesqr = charge * charge;
  auto chargesqrSum = fChargeSqr + chargesqr;
  
  Double_t x  =  hit->GetPosition().X();
  Double_t y  =  hit->GetPosition().Y();
  Double_t z  =  hit->GetPosition().Z();
  
  Double_t t = hit->GetCTime();
  Double_t tSum = fMeanTime + t;

  fS = (fCharge * fS + charge * hit->GetS()) / chargeSum;

  Int_t numHits = GetNumHits();
    
  if (numHits == 0) {
    fX = x;
    fY = y;
    fZ = z;
    //    CalculateCovMatrix(hit->GetPosition(),charge);
    fCovMatrix(0,0) = hit->GetResolutionX() * hit->GetResolutionX();
    fCovMatrix(1,1) = hit->GetResolutionY() * hit->GetResolutionY();
    fCovMatrix(2,2) = hit->GetResolutionZ() * hit->GetResolutionZ();
  }  
  else {
    fMeanTime += t;
    
    CalculatePosition(hit->GetPosition(),charge);

    if (numHits == 1) {
      fCovMatrix(0,0) = 0;
      fCovMatrix(1,1) = 0;
      fCovMatrix(2,2) = 0;
      CalculateCovMatrix(hit->GetPosition(),charge);
    }
    else
      CalculateCovMatrix(hit->GetPosition(),charge);
  }
  
    if (fCovMatrix(0,0) == 0 || fCovMatrix(2,2) == 0){
      fCovMatrix(0,0) = 0.4;
      fCovMatrix(2,2) = 0.4;
    }

    fLayer = hit->GetLayer();
    
    TVector2 posPad(fX,fZ+143);
    Double_t meanRow = tpc::getMrow(fLayer,posPad.Phi());
    fRow = int(meanRow);
    
  fCharge = chargeSum;
  fChargeSqr = chargesqrSum;
  
  Int_t neff = fCharge * fCharge / fChargeSqr;
  
  fDx = sqrt(fCovMatrix(0,0));
  fDy = sqrt(fCovMatrix(1,1));
  fDz = sqrt(fCovMatrix(2,2));

  if (numHits > 1) {
    fDx = fDx/sqrt(neff);
    fDy = fDy/sqrt(neff);
    fDz = fDz/sqrt(neff);
  }
    
  fHits.push_back(hit);

  fMeanCharge = chargeSum/GetNumHits();
  fMeanTime = tSum/GetNumHits();
  
  hit->SetClusterID(fClusterID);

}

void TPCHitCluster::ApplyModifiedHitInfo() {

  fX = 0;
  fY = 0;
  fZ = 0;
  
  fDx = 0;
  fDy = 0;
  fDz = 0;

  fCharge = 0;
  fS = 0;

  xxMean =  0; yyMean = 0 ; zzMean = 0;
  xyMean =  0; yzMean = 0 ; xzMean = 0;
  xMean =  0;  yMean = 0 ; zMean = 0;
  
  fCovMatrix.ResizeTo(3,3);
  for (Int_t iElem = 0 ; iElem < 9; iElem++)
    fCovMatrix(iElem/3,iElem%3) = 0;

  Double_t default_cov[3] = {4*4, 1*1, 6*6};

  Int_t hit_idx = 0 ;
  
  for (auto hit : fHits) {
    //    auto charge = hit->GetDe();
    auto charge = hit->GetCDe();
    auto chargeSum = fCharge + charge;
  
    Double_t x  =  hit->GetPosition().X();
    Double_t y  =  hit->GetPosition().Y();
    Double_t z  =  hit->GetPosition().Z();
    
    xMean = (fCharge * xMean + charge * x)/chargeSum;
    yMean = (fCharge * yMean + charge * y)/chargeSum;
    zMean = (fCharge * zMean + charge * z)/chargeSum;
    
    xxMean = ( fCharge * xxMean + charge * x * x) / chargeSum;
    yyMean = ( fCharge * yyMean + charge * y * y) / chargeSum;
    zzMean = ( fCharge * zzMean + charge * z * z) / chargeSum;
    
    xyMean = ( fCharge * xyMean + charge * x * y) / chargeSum;
    xzMean = ( fCharge * xzMean + charge * x * z) / chargeSum;
    yzMean = ( fCharge * yzMean + charge * y * z) / chargeSum;
    
    Double_t xxCov = fCharge * (xxMean - xMean * xMean);
    Double_t yyCov = fCharge * (yyMean - yMean * yMean);
    Double_t zzCov = fCharge * (zzMean - zMean * zMean);
    
    Double_t xyCov = fCharge * (xyMean -  xMean * yMean);
    Double_t xzCov = fCharge * (xzMean -  xMean * zMean);
    Double_t yzCov = fCharge * (yzMean -  yMean * zMean);
  
    fS = (fCharge * fS + charge * hit->GetS()) / chargeSum;

    if (hit_idx == 0) {
      for (int i = 0 ; i < 3 ; i++)
	for (int j = 0 ; j < 3 ; j++)
	  fCovMatrix(i,j) = 0;

      for (int i = 0 ; i < 3 ; i++)
	fCovMatrix(i,i) = default_cov[i];

    }
    else {
      fCovMatrix(0,0) = xxCov;  fCovMatrix(0,1) = xyCov;  fCovMatrix(2,2) = xzCov;
      fCovMatrix(1,0) = xyCov;  fCovMatrix(1,1) = yyCov;  fCovMatrix(2,2) = yzCov;
      fCovMatrix(2,0) = xzCov;  fCovMatrix(2,1) = yzCov;  fCovMatrix(2,2) = zzCov;
    }
    
    fLayer = hit->GetLayer();      
    fRow = hit->GetRow();
    
    //    if (fLayer != -1) fCovMatrix(2,2) = default_cov[2];
    //    else if (fRow != -1) fCovMatrix(0,0) = default_cov[0];

    
    fCharge = chargeSum;
    
    fDx = sqrt(fCovMatrix(0,0));
    fDy = sqrt(fCovMatrix(1,1));
    fDz = sqrt(fCovMatrix(2,2));

    hit_idx++;
  }
}

void TPCHitCluster::CalculatePosition(TVector3 hitPos, Double_t charge) {
  TVector3 position(fX,fY,fZ);

  for (Int_t iPos = 0 ; iPos < 3 ; iPos++)
    position[iPos] += charge*(hitPos[iPos] - position[iPos])/(fCharge + charge);

  fX = position.X();
  fY = position.Y();
  fZ = position.Z();
}

void TPCHitCluster::CalculateCovMatrix(TVector3 hitPos, Double_t charge) {

  TVector3 position(fX,fY,fZ);

  for (Int_t iFirst = 0 ; iFirst < 3 ; iFirst++) {
    for (Int_t iSecond = 0; iSecond < iFirst+1 ; iSecond++) {
      fCovMatrix(iFirst, iSecond) = fCharge * fCovMatrix(iFirst, iSecond) / (fCharge + charge);
      fCovMatrix(iFirst, iSecond) += charge * (hitPos[iFirst] - position[iFirst]) * (hitPos[iSecond] - position[iSecond])/fCharge;
      fCovMatrix(iSecond, iFirst) = fCovMatrix(iFirst, iSecond);
    }
  }
  
  //  fDx = sqrt(fCovMatrix(0,0));
  //  fDy = sqrt(fCovMatrix(1,1));
  //  fDz = sqrt(fCovMatrix(2,2));

}

void TPCHitCluster::SetDFromCov(Double_t maxx, Double_t maxy, Double_t maxz, Bool_t setMin){
  /*
  std::cout<<"setd-1"<<std::endl;
  static const Double_t ResoXZ = gUser.GetParameter("ResoXZ");
  static const Double_t ResoY = gUser.GetParameter("ResoY");
  
  std::vector<Double_t> ResPar;
  std::cout<<"setd-2"<<std::endl;
  //  std::cout << fLayer << std::endl;
 if(fLayer < 10){
    ResPar = ResParamInnerLayerHSOn;
  }
  else{
    ResPar = ResParamOuterLayerHSOn;
  }

  std::cout<<"setd-3"<<std::endl;
  std::cout<<"data check : "<<ResPar[0]<<std::endl;
  double par_t[6]={
    ResPar[0],ResPar[1],ResPar[2],ResPar[3],ResPar[4],ResPar[5]};

  std::cout<<"setd-4"<<std::endl;

  double par_y[4] = {
    ResPar[6],ResPar[1],ResPar[7],ResPar[8]};

  std::cout<<"setd-5"<<std::endl;

  TVector3 resTPC = GetResVector(TVector3(fX,fY,fZ), par_t, par_y);
  
  std::cout<<"setd-6"<<std::endl;
  fDx = resTPC.X();
  fDy = resTPC.Y();
  fDz = resTPC.Z();
  std::cout<<"setd-7"<<std::endl;
  */
  fDx = 1.;
  fDy = 1.;
  fDz = 1.;


}


void TPCHitCluster::SetClusterID(Int_t clusterID) {
  //  TPCHit::SetClusterID(clusterID);
  for (auto hit: fHits)
    hit->SetClusterID(clusterID);
}


inline Double_t TPCHitCluster::GetVerticalRes(Double_t y, Double_t* par){

  Double_t p0 = par[0];
  Double_t p1 = par[1];
  Double_t p2 = par[2];
  Double_t p3 = par[3];
  
  Double_t dl = y + 300;

  Double_t val = sqrt(p0*p0 + p2*p2*dl/(p3*exp(-p1*dl)));

  
  return val;		      
}

inline Double_t TPCHitCluster::GetTransverseRes(Double_t angle, Double_t y, Double_t* par) {
  
  Double_t p0 = par[0];
  Double_t p1 = par[1];
  Double_t p2 = par[2];
  Double_t p3 = par[3];
  Double_t p4 = par[4];
  Double_t p5 = par[5];

  //  std::cout << p0 << "\t" << p1 << "\t" << p2 << "\t" << p3 << std::endl;

  Double_t dl = y + 300;
  
  Double_t val = sqrt(p0*p0 + p2*p2*dl/(p3*exp(-p1*dl)) + p4*p4*tan(angle)*tan(angle)/(12.*p5));
  //  std::cout << angle << "\t" << val << std::endl;
  return val;
}
      
    
inline TVector3 TPCHitCluster::GetResVector(TVector3 pos, Double_t *par_t, Double_t *par_y){

  Double_t theta = TMath::ATan2(pos.X(), pos.Z()+143.);
  //  Double_t theta = TMath::ATan2(pos.X(), pos.Z());
  Double_t chi = GetChi();
  //  std::cout << theta - chi << std::endl;
  Double_t angle = abs(theta-chi);
  
  Double_t resY = GetVerticalRes(pos.Y(), par_y);
  Double_t resT = GetTransverseRes(angle, pos.Y(), par_t);

  Double_t resX, resZ;
  if (TMath::Cos(angle) <1e-02){
    resX = 1e-02 * resT;
    resZ = sqrt(1.-1e-4)*resT;    
  }
  else {
    resX = abs(TMath::Cos(angle))*resT;
    resZ = abs(TMath::Sin(angle))*resT;
  }
  
  return TVector3(resX, resY, resZ);
}

