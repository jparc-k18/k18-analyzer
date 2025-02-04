#include "TPCRiemannFitter.hh"
#include "UserParamMan.hh"

#include <TVector2.h>
#include <iostream>
#include <TCanvas.h>
#include <TPolyMarker3D.h>
#include <TPolyLine3D.h>

namespace {
  const auto& gUser   = UserParamMan::GetInstance();
  const TVector3 tgtPos(0,0,-143);
  //  Double_t scale = 1.0;
}

Bool_t TPCRiemannFitter::FitLinear(TPCLinearTrack *track) {

  auto hitArray = track->GetHits();
  if (hitArray->size() < 4) { 
    //    std::cout <<"Not enough hits" << std::endl;
    return false;
  }

  fODRFitter->Reset();

  Double_t scale = 0.5;
  TVector3 mean = track->GetMean();
  Double_t xMean = mean.X();
  Double_t yMean = mean.Y();
  Double_t zMean = mean.Z();
  
  Double_t xCov = track->GetXCov();
  Double_t zCov = track->GetZCov();

  Double_t RSR = 2*sqrt(xCov+zCov);

  Double_t xMapMean = 0;
  Double_t yMapMean = 0;
  Double_t zMapMean = 0;

  Double_t x = 0 ; Double_t z = 0;

  for (auto hit: *hitArray) {

    x = hit->GetX() - xMean;
    z = hit->GetZ() - zMean;
    Double_t w = hit->GetCDe();
    //    Double_t w = hit->GetXDe();
    w = TMath::Power(w,scale);
    
    Double_t rEff = sqrt(x*x + z*z) / (2*RSR);
    Double_t denominator = 1 + rEff*rEff;

    Double_t xMap = x / denominator;
    Double_t zMap = z / denominator;
    Double_t yMap = 2 * RSR * rEff * rEff / denominator;
    //    Double_t yMap = (rEff*rEff/(2*RSR))/denominator;

    xMapMean += w * xMap;
    zMapMean += w * zMap;
    yMapMean += w * yMap;
  }

  Double_t weightSum = track->GetChargeSum();
  
  xMapMean = xMapMean / weightSum;
  yMapMean = yMapMean / weightSum;
  zMapMean = zMapMean / weightSum;

  fODRFitter -> SetCentroid(xMapMean, yMapMean, zMapMean);
  TVector3 mapMean(xMapMean, yMapMean, zMapMean);

  for (auto hit : *hitArray) {
    x = hit->GetX() - xMean;
    z = hit->GetZ() - zMean;
    Double_t w = hit->GetCDe();
    //    Double_t w = hit->GetDe();
    w = TMath::Power(w,scale);

    Double_t rEff = sqrt(x*x + z*z)/(2*RSR);
    Double_t denominator = 1 + rEff*rEff;

    Double_t xMap = x / denominator;
    Double_t zMap = z / denominator;
    Double_t yMap = 2 * RSR * rEff * rEff / denominator;
    //    Double_t yMap = (rEff*rEff/(2*RSR))/denominator;

    fODRFitter -> AddPoint(xMap, yMap, zMap, w);
  }
  
  fODRFitter->Solve();
    
  fODRFitter->ChooseEigenValue(0);
  TVector3 direction = fODRFitter->GetDirection();
  if (direction.Z() < 0) direction *= -1;
  track->SetDirection(direction);
  track->SetRMSLine(fODRFitter->GetRMSLine());

  fODRFitter->ChooseEigenValue(2);
  TVector3 normal = fODRFitter->GetNormal();
  track->SetNormal(normal);
  track->SetRMSPlane(fODRFitter->GetRMSPlane());

  track->SetIsFitted(true);

  return true;
}

Bool_t TPCRiemannFitter::FitPlane(TPCLinearTrack *track, TPCHit *hit,
				  Double_t &rmsL, Double_t &rmsP) {
  
  auto hitArray = track->GetHits();
  if (hitArray->size() < 4) { 
    //    std::cout <<"Not enough hits" << std::endl;
    return false;
  }

  fODRFitter->Reset();

  TPCLinearTrack *trackTmp = new TPCLinearTrack(track);
  trackTmp->AddHit(hit);

  FitPlane(trackTmp);

  fODRFitter->ChooseEigenValue(0);
  rmsL = fODRFitter->GetRMSLine();

  fODRFitter->ChooseEigenValue(2);
  rmsP = fODRFitter->GetRMSPlane();
  
  return true;
}

Bool_t TPCRiemannFitter::FitPlane(TPCLinearTrack *track) {

  if (track->GetNumHits() < 4)
    return false;
  
  fODRFitter->Reset();

  TVector3 mean = track->GetMean();
  Double_t xMean = mean.X(); Double_t yMean = mean.Y(); Double_t zMean = mean.Z();

  fODRFitter->SetCentroid(xMean,yMean,zMean);
  fODRFitter->SetMatrixA(track->xxCov(), track->xyCov(), track->zxCov(),
			 track->yyCov(), track->yzCov(),
			 track->zzCov());

  fODRFitter->SetWeightSum(track->GetChargeSum());
  fODRFitter->SetNumPoints(track->GetNumHits());

  fODRFitter->Solve();

  fODRFitter->ChooseEigenValue(0);
  TVector3 direction = fODRFitter->GetDirection();
  if (direction.Z() < 0) direction *= -1;
  track->SetDirection(direction);
  track->SetRMSLine(fODRFitter->GetRMSLine());
  
  fODRFitter->ChooseEigenValue(2);
  TVector3 normal = fODRFitter->GetNormal();
  track->SetNormal(normal);
  track->SetRMSPlane(fODRFitter->GetRMSPlane());

  track->SetIsFitted(true);
  return true;

}
Bool_t TPCRiemannFitter::FitPlane(TPCRiemannTrack *track) {
  

  if (track->GetNumHits() < 4)
    return false;
  
  fODRFitter->Reset();

  TVector3 mean = track->GetMean();
  Double_t xMean = mean.X(); Double_t yMean = mean.Y(); Double_t zMean = mean.Z();

  fODRFitter->SetCentroid(xMean,yMean,zMean);
  fODRFitter->SetMatrixA(track->xxCov(), track->xyCov(), track->zxCov(),
			 track->yyCov(), track->yzCov(),
			 track->zzCov());

  fODRFitter->SetWeightSum(track->GetChargeSum());
  fODRFitter->SetNumPoints(track->GetNumHits());

  if (fODRFitter->Solve() == false)
    return false;

  fODRFitter->ChooseEigenValue(2);
  TVector3 normal = fODRFitter->GetDirection();

  if (normal.Y() < 1e-10) {
    track->SetIsLine();
    fODRFitter->ChooseEigenValue(0);
    track->SetLineDirection(fODRFitter->GetDirection());
    track->SetRMSH(fODRFitter->GetRMSLine());

    return true;
  }

  track->SetIsPlane();
  track->SetPlaneNormal(normal);
  track->SetRMSH(fODRFitter->GetRMSPlane());


  return true;

  
}

Bool_t TPCRiemannFitter::Fit(TPCRiemannTrack *track){

  auto hitArray = track->GetHits();
  if (hitArray->size() < 4) { 
    //    std::cout <<"Not enough hits" << std::endl;
    return false;
  }

  Double_t scale = 1.0;
  const Double_t MaxTrackLength = gUser.GetParameter("MaxTrackLength");
  const Double_t MaxMeanCharge = gUser.GetParameter("MaxMeanCharge");

  const Bool_t DoScale = gUser.GetParameter("DoScale");

  if (DoScale) {  

    Double_t trackLength = track->TrackLength();
    Double_t meanCharge = track->GetChargeSum()/(track->GetNumHits());
    //    std::cout << meanCharge << std::endl;
    
    if (trackLength < MaxTrackLength){
      Double_t scaleTrackLength = (MaxTrackLength-trackLength)/MaxTrackLength;
      if (meanCharge < MaxMeanCharge)
	scaleTrackLength *= (1. - (MaxMeanCharge * MaxMeanCharge - meanCharge*meanCharge)/(MaxMeanCharge * MaxMeanCharge));
      scale += scaleTrackLength;
    }
  }
  
  fODRFitter->Reset();
  
  Double_t xMean = track->GetMean().X();
  Double_t zMean = track->GetMean().Z();

  Double_t xCov = track->GetXCov();
  Double_t zCov = track->GetZCov();

  Double_t RSR = 2*sqrt(xCov+zCov);

  Double_t xMapMean = 0;
  Double_t yMapMean = 0;
  Double_t zMapMean = 0;

  Double_t x = 0 ; Double_t z = 0;

  for (auto hit: *hitArray) {

    x = hit->GetX() - xMean;
    z = hit->GetZ() - zMean;
    
    Double_t w = hit->GetCDe();
    //    Double_t w = hit->GetDe();
    w = TMath::Power(w,scale); 
    
    Double_t rEff = sqrt(x*x + z*z) / (2*RSR);
    Double_t denominator = 1 + rEff*rEff;

    Double_t xMap = x / denominator;
    Double_t zMap = z / denominator;
    Double_t yMap = 2 * RSR * rEff * rEff / denominator;
    //    Double_t yMap = (rEff*rEff/(2*RSR))/denominator;

    xMapMean += w * xMap;
    zMapMean += w * zMap;
    yMapMean += w * yMap;
  }

  Double_t weightSum = track->GetChargeSum();
  
  xMapMean = xMapMean / weightSum;
  yMapMean = yMapMean / weightSum;
  zMapMean = zMapMean / weightSum;

  fODRFitter -> SetCentroid(xMapMean, yMapMean, zMapMean);
  TVector3 mapMean(xMapMean, yMapMean, zMapMean);

  for (auto hit : *hitArray) {
    x = hit->GetX() - xMean;
    z = hit->GetZ() - zMean;
    Double_t w = hit->GetCDe();
    w = TMath::Power(w,scale);

    Double_t rEff = sqrt(x*x + z*z)/(2*RSR);
    Double_t denominator = 1 + rEff*rEff;

    Double_t xMap = x / denominator;
    Double_t zMap = z / denominator;
    Double_t yMap = 2 * RSR * rEff * rEff / denominator;
    //    Double_t yMap = (rEff*rEff/(2*RSR))/denominator;

    fODRFitter -> AddPoint(xMap, yMap, zMap, w);
  }
  
  if (fODRFitter->Solve() == false)
    return false;
  
  fODRFitter->ChooseEigenValue(0);
  TVector3 uOnPlane = fODRFitter->GetDirection();
  fODRFitter->ChooseEigenValue(1);
  TVector3 vOnPlane = fODRFitter->GetDirection();
  fODRFitter->ChooseEigenValue(2);
  TVector3 nToPlane = fODRFitter->GetDirection();

  if (std::abs(nToPlane.Y()) < 1e-8) { 
    track->SetIsLine();
    return false;
  }
   
  TVector3 RSC(0,RSR,0);
  Double_t tRCC = nToPlane.Dot(mapMean - RSC)/nToPlane.Mag2();  
  TVector3 RCC = tRCC * nToPlane + RSC;
  Double_t CtoS = (RCC-RSC).Mag();
  Double_t RCR = sqrt(RSR*RSR-CtoS*CtoS);

  Double_t uConst = uOnPlane.Y()/sqrt(uOnPlane.Y()*uOnPlane.Y() + vOnPlane.Y()*vOnPlane.Y());
  Double_t vConst = sqrt(1-uConst*uConst);

  Double_t ref1 = uConst * uOnPlane.Y() + vConst * vOnPlane.Y();
  Double_t ref2 = uConst * uOnPlane.Y() - vConst * vOnPlane.Y();

  if (ref1 < 0 ) ref1 = -ref1;
  if (ref2 < 0 ) ref2 = -ref2;
  if (ref1 < ref2) vConst = -vConst;

  TVector3 toLouu = uConst *uOnPlane + vConst * vOnPlane;
  TVector3 louu = RCC + toLouu * RCR;
  TVector3 high = RCC - toLouu * RCR;

  TVector3 louuInvMap(louu.X()/(1-louu.Y()/(2*RSR)), 0, louu.Z()/(1-louu.Y()/(2*RSR)));
  TVector3 highInvMap(high.X()/(1-high.Y()/(2*RSR)), 0, high.Z()/(1-high.Y()/(2*RSR)));

  TVector3 FCC = 0.5*(louuInvMap + highInvMap);

  Double_t xCenter = FCC.X() + xMean;
  Double_t zCenter = FCC.Z() + zMean;
  Double_t radius = 0.5 * (louuInvMap - highInvMap).Mag();
  
  if (radius > 1e+8) {
    track->SetIsLine();
    return false;
  }

  track->SetCenter(xCenter, zCenter);
  track->SetRadius(radius);
  track->SetIsHelix();

  //  std::sort(hitArray->begin(), hitArray->end(), SortByX());
  std::sort(hitArray->begin(), hitArray->end(), SortByY());
  //  std::sort(hitArray->begin(), hitArray->end(), SortByZ());

  TVector3 position0 = hitArray->at(0)->GetPosition();

  x = position0.X()-xCenter;
  z = position0.Z()-zCenter;

  Double_t y = 0;

  Double_t alphaInit = TMath::ATan2(x,z);

  //  TVector2 xAxis(-x,z);
  //  TVector2 zAxis(z,x);

  TVector2 xAxis(-x,z);
  TVector2 zAxis(z,x);

  //  TVector2 zAxis(-z,x);
  //  TVector2 xAxis(x,z);

  zAxis = zAxis.Unit();
  xAxis = xAxis.Unit();

  Double_t expA= 0;
  Double_t expA2 = 0 ;
  Double_t expY = 0;
  Double_t expAY = 0;

  Double_t alphaStack = 0;
  Double_t alphaLast = 0;

  Double_t alphaMin = alphaInit;
  Double_t alphaMax = alphaInit;

  for (auto hit : *hitArray) {
    x = hit->GetX() - xCenter;
    y = hit->GetY();
    z = hit->GetZ() - zCenter;
    
    //    TVector2 v(x,z);
    TVector2 v(z,x);

    //    Double_t xRot = v*xAxis;
    //    Double_t zRot = v*zAxis;

    Double_t zRot = v*zAxis;
    Double_t xRot = v*xAxis;

    alphaLast = TMath::ATan2(xRot,zRot);
    //    alphaLast = TMath::ATan2(zRot,xRot);

    if (alphaLast > TMath::PiOver2() || alphaLast < -1 * TMath::PiOver2()) {
      Double_t t0 = alphaLast;

      alphaLast += alphaStack;
      alphaStack += t0;

      zAxis = v;
      xAxis = TVector2(-v.Y(),v.X());
      xAxis = xAxis.Unit();
      zAxis = zAxis.Unit();
    }

    else
      alphaLast += alphaStack;

    alphaLast = alphaLast + alphaInit;

    Double_t w = hit->GetCDe();
    //    Double_t w = hit->GetDe();

    expA += w * alphaLast;
    expA2 += w * alphaLast * alphaLast;
    expY += w * y;
    expAY += w * alphaLast * y;

    if (alphaLast < alphaMin)
      alphaMin = alphaLast;
    if (alphaLast > alphaMax)
      alphaMax = alphaLast;    
  }

  track->SetAlphaHead(alphaMin);
  track->SetAlphaTail(alphaMax);
  
  expA /= weightSum;
  expA2 /= weightSum;
  expY /= weightSum;
  expAY /= weightSum;

  Double_t slope = (expAY - expA * expY) / (expA2 - expA*expA);
  Double_t offset = (expA2*expY - expA*expAY) / (expA2 - expA*expA);

  if (std::isinf(slope)){
    track->SetIsLine();
    return false;
  }

  track->SetSlope(slope);
  track->SetOffset(offset);

  Double_t Sx = 0; Double_t Sy = 0;
  for (auto hit : *hitArray){
    TVector3 q = track->Map(hit->GetPosition());
    
    //    Double_t val = TMath::Sqrt(q.X() * q.X() + q.Z() * q.Z());

    //    Sx += hit->GetCDe() * val * val;
    Sx += hit->GetCDe() * q.X() * q.X();
    Sy += hit->GetCDe() * q.y() * q.y();
  }

  Double_t rmsR = sqrt(Sx / (track->GetChargeSum() * (1 - 3/hitArray->size())));
  Double_t rmsY = sqrt(Sy / (track->GetChargeSum() * (1 - 3/hitArray->size())));

  track->SetRMSW(rmsR);
  track->SetRMSH(rmsY);


  return true;

}

Bool_t TPCRiemannFitter::FitCluster(TPCLinearTrack *track) {

  if (track->GetNumClusters() < 4)
    return false;
  
  fODRFitter->Reset();

  auto clusterArray = track->GetClusters();

  Double_t xMean = 0;
  Double_t yMean = 0;
  Double_t zMean = 0;
  Double_t weightSum = 0;

  for (auto cluster : *clusterArray) {

    if (cluster->IsStable() == false)
      continue;

    Double_t w = cluster->GetCharge();

    xMean += w*cluster->GetX();
    yMean += w*cluster->GetY();
    zMean += w*cluster->GetZ();
    weightSum += w;
  }

  xMean = xMean/weightSum;
  yMean = yMean/weightSum;
  zMean = zMean/weightSum;

  Double_t xCov = track->GetXCov();
  Double_t zCov = track->GetZCov();

  Double_t RSR = 2*sqrt(xCov+zCov);

  Double_t xMapMean = 0;
  Double_t yMapMean = 0;
  Double_t zMapMean = 0;

  Double_t x = 0;
  Double_t z = 0;

  for (auto cluster : *clusterArray) {

    x = cluster->GetX() - xMean;
    z = cluster->GetZ() - zMean;
    Double_t w = cluster->GetCharge();

    Double_t rEff = sqrt(x*x + z*z)/(2*RSR);
    Double_t denominator = 1 + rEff*rEff;

    Double_t xMap = x/denominator;
    Double_t yMap = 2 * RSR * rEff * rEff / denominator;
    //    Double_t yMap = (rEff*rEff/(2*RSR))/denominator;
    Double_t zMap = z/denominator;

    xMapMean += w * xMap;
    yMapMean += w * yMap;
    zMapMean += w * zMap;
  }

  xMapMean = xMapMean/weightSum;
  yMapMean = yMapMean/weightSum;
  zMapMean = zMapMean/weightSum;

  fODRFitter -> SetCentroid(xMapMean, yMapMean, zMapMean);
  TVector3 mapMean(xMapMean, yMapMean, zMapMean);

  for (auto cluster : *clusterArray) {

    x = cluster->GetX() - xMean;
    z = cluster->GetZ() - zMean;
    Double_t w = cluster->GetCharge();

    Double_t rEff = sqrt(x*x + z*z)/(2*RSR);
    Double_t denominator = 1 + rEff*rEff;

    Double_t xMap = x/denominator;
    Double_t yMap = 2 * RSR * rEff * rEff / denominator;
    //    Double_t yMap = (rEff*rEff/(2*RSR))/denominator;
    Double_t zMap = z/denominator;

    fODRFitter->AddPoint(xMap,yMap,zMap,w);
  }

  if (fODRFitter->Solve() == false)
    return false;

  fODRFitter->ChooseEigenValue(0);
  TVector3 direction = fODRFitter->GetDirection();
  if (direction.Z() < 0) direction *= -1;
  track->SetDirection(direction);
  track->SetRMSLine(fODRFitter->GetRMSLine());
  
  fODRFitter->ChooseEigenValue(2);
  TVector3 normal = fODRFitter->GetNormal();
  track->SetNormal(normal);
  track->SetRMSPlane(fODRFitter->GetRMSPlane());

  track->SetIsFitted(true);

  return true;
}

Bool_t TPCRiemannFitter::FitCluster(TPCRiemannTrack *track) {


  Double_t scale = 0.5;

  //  Double_t scale = 1.0;
  const Double_t MaxTrackLength = gUser.GetParameter("MaxTrackLength");
  const Double_t MaxMeanCharge = gUser.GetParameter("MaxMeanCharge");

  const Bool_t DoScale = gUser.GetParameter("DoScale");

  if (DoScale) {  

    Double_t trackLength = track->TrackLength();
    Double_t meanCharge = track->GetChargeSum()/(track->GetNumHits());
    //    std::cout << meanCharge << std::endl;
    
    if (trackLength < MaxTrackLength){
      Double_t scaleTrackLength = (MaxTrackLength-trackLength)/MaxTrackLength;
      if (meanCharge < MaxMeanCharge)
	scaleTrackLength *= (1. - (MaxMeanCharge * MaxMeanCharge - meanCharge*meanCharge)/(MaxMeanCharge * MaxMeanCharge));
      scale += scaleTrackLength;
    }
  }
  
  if (track->GetNumClusters() < 4)
    return false;

  fODRFitter->Reset();

  auto clusterArray = track->GetClusters();

  Double_t xMean = 0;
  Double_t zMean = 0;
  Double_t weightSum = 0;

  for (auto cluster : *clusterArray) {

    if (cluster->IsStable() == false)
      continue;

    Double_t w = cluster->GetCharge();
    w = TMath::Power(w,scale);
    xMean += w*cluster->GetX();
    zMean += w*cluster->GetZ();
    weightSum += w;
  }
  
  xMean = xMean/weightSum;
  zMean = zMean/weightSum;

  Double_t xCov = track->GetXCov();
  Double_t zCov = track->GetZCov();

  Double_t RSR = 2*sqrt(xCov+zCov);

  Double_t xMapMean = 0;
  Double_t yMapMean = 0;
  Double_t zMapMean = 0;

  Double_t x = 0;
  Double_t z = 0;

  for (auto cluster : *clusterArray) {
    if (cluster->IsStable() == false)
      continue;

    x = cluster->GetX() - xMean;
    z = cluster->GetZ() - zMean;
    Double_t w = cluster->GetCharge();
    w = TMath::Power(w,scale);

    Double_t rEff = sqrt(x*x + z*z)/(2*RSR);
    Double_t denominator = 1 + rEff*rEff;

    Double_t xMap = x/denominator;
    Double_t yMap = 2 * RSR * rEff * rEff / denominator;
    //    Double_t yMap = (rEff*rEff/(2*RSR))/denominator;
    Double_t zMap = z/denominator;

    xMapMean += w * xMap;
    yMapMean += w * yMap;
    zMapMean += w * zMap;
  }

  xMapMean = xMapMean/weightSum;
  yMapMean = yMapMean/weightSum;
  zMapMean = zMapMean/weightSum;

  fODRFitter -> SetCentroid(xMapMean, yMapMean, zMapMean);
  TVector3 mapMean(xMapMean, yMapMean, zMapMean);

  for (auto cluster : *clusterArray) {

    if (cluster->IsStable() == false)
      continue;

    x = cluster->GetX() - xMean;
    z = cluster->GetZ() - zMean;
    Double_t w = cluster->GetCharge();
    w = TMath::Power(w,scale);

    Double_t rEff = sqrt(x*x + z*z)/(2*RSR);
    Double_t denominator = 1 + rEff*rEff;

    Double_t xMap = x/denominator;
    Double_t yMap = 2 * RSR * rEff * rEff / denominator;
    //    Double_t yMap = (rEff*rEff/(2*RSR))/denominator;
    Double_t zMap = z/denominator;

    fODRFitter->AddPoint(xMap,yMap,zMap,w);
  }

  if (fODRFitter->Solve() == false)
    return false;
  
  fODRFitter->ChooseEigenValue(0);
  TVector3 uOnPlane = fODRFitter->GetDirection();
  fODRFitter->ChooseEigenValue(1);
  TVector3 vOnPlane = fODRFitter->GetDirection();
  fODRFitter->ChooseEigenValue(2);
  TVector3 nToPlane = fODRFitter->GetDirection();

  if (std::abs(nToPlane.Y()) <1e-10) {
    track->SetIsLine();
    return false;
  }
   
  TVector3 RSC(0,RSR,0);
  Double_t tRCC = nToPlane.Dot(mapMean - RSC)/nToPlane.Mag2();  
  TVector3 RCC = tRCC*nToPlane + RSC;
  Double_t CtoS = (RCC-RSC).Mag();
  Double_t RCR = sqrt(RSR*RSR-CtoS*CtoS);

  Double_t uConst = uOnPlane.Y()/sqrt(uOnPlane.Y()*uOnPlane.Y() + vOnPlane.Y()*vOnPlane.Y());
  Double_t vConst = sqrt(1-uConst*uConst);

  Double_t ref1 = uConst*uOnPlane.Y() + vConst*vOnPlane.Y();
  Double_t ref2 = uConst*uOnPlane.Y() - vConst*vOnPlane.Y();

  if (ref1 < 0 ) ref1 = -ref1;
  if (ref2 < 0 ) ref2 = -ref2;
  if (ref1 < ref2) vConst = -vConst;

  TVector3 toLouu = uConst *uOnPlane + vConst * vOnPlane;
  TVector3 louu = RCC + toLouu * RCR;
  TVector3 high = RCC - toLouu * RCR;

  TVector3 louuInvMap(louu.X()/(1-louu.Y()/(2*RSR)), 0, louu.Z()/(1-louu.Y()/(2*RSR)));
  TVector3 highInvMap(high.X()/(1-high.Y()/(2*RSR)), 0, high.Z()/(1-high.Y()/(2*RSR)));

  TVector3 FCC = 0.5*(louuInvMap + highInvMap);

  Double_t xCenter = FCC.X() + xMean;
  Double_t zCenter = FCC.Z() + zMean;
  Double_t radius = 0.5 * (louuInvMap - highInvMap).Mag();
  
  if (radius > 1e+8) {
    track->SetIsLine();
    return false;
  }
  track->SetCenter(xCenter, zCenter);
  track->SetRadius(radius);
  track->SetIsHelix();

  TVector3 position0;
  Int_t firstIndex = 0;

  std::sort(clusterArray->begin(), clusterArray->end(), SortByY());

  for (auto cluster : *clusterArray) {
    if (cluster->IsStable()){
      position0 = cluster->GetPosition();
      break;
    }
   
    firstIndex++;
  }
  
  x = position0.X()-xCenter;
  z = position0.Z()-zCenter;

  Double_t y = 0;

  //  Double_t alphaInit = TMath::ATan2(x,z);
  //  TVector2 xAxis(z,x);
  //  TVector2 zAxis(-x,z);

  Double_t alphaInit = TMath::ATan2(x,z);
  TVector2 xAxis(-x,z);
  TVector2 zAxis(z,x);

  xAxis = xAxis.Unit();
  zAxis = zAxis.Unit();

  Double_t expA= 0;
  Double_t expA2 = 0 ;
  Double_t expY = 0;
  Double_t expAY = 0;

  Double_t alphaStack = 0;
  Double_t alphaLast = 0;

  Double_t alphaMin = alphaInit;
  Double_t alphaMax = alphaInit;

  auto numClusters = clusterArray->size();

  for (auto iCluster = firstIndex ; iCluster < numClusters ; iCluster++) {
    
    auto cluster = clusterArray->at(iCluster);
    
    if (cluster->IsStable() == false)
      continue;

    x = cluster->GetX() - xCenter;
    y = cluster->GetY();
    z = cluster->GetZ() - zCenter;
    
    TVector2 v(z,x);
    //    TVector2 v(x,z);

    Double_t xRot = v*xAxis;
    Double_t zRot = v*zAxis;

    alphaLast = TMath::ATan2(xRot,zRot);
    //    alphaLast = TMath::ATan2(zRot,xRot);

    if (alphaLast > TMath::PiOver2() || alphaLast < -1 * TMath::PiOver2()) {
    //    if (abs(alphaLast) < TMath::PiOver2()){
      Double_t t0 = alphaLast;

      alphaLast += alphaStack;
      alphaStack += t0;

      zAxis = v;
      xAxis = TVector2(-v.Y(),v.X());
      zAxis = zAxis.Unit();
      xAxis = xAxis.Unit();
    }

    else
      alphaLast += alphaStack;

    alphaLast = alphaLast + alphaInit;

    Double_t w = cluster->GetCharge();
    w = TMath::Power(w,scale);

    expA += w * alphaLast;
    expA2 += w * alphaLast * alphaLast;
    expY += w * y;
    expAY += w * alphaLast * y;

    if (alphaLast < alphaMin)
      alphaMin = alphaLast;
    if (alphaLast > alphaMax)
      alphaMax = alphaLast;    
  }

  track->SetAlphaHead(alphaMin);
  track->SetAlphaTail(alphaMax);
  
  expA /= weightSum;
  expA2 /= weightSum;
  expY /= weightSum;
  expAY /= weightSum;

  Double_t slope = (expAY - expA * expY) / (expA2 - expA*expA);
  Double_t offset = (expA2*expY - expA*expAY) / (expA2 - expA*expA);

  if (std::isinf(slope)){
    track->SetIsLine();
    return false;
  }

  track->SetSlope(slope);
  track->SetOffset(offset);

  Double_t Sx = 0; Double_t Sy = 0;
  for (auto cluster : *clusterArray){

    if (cluster->IsStable() == false)
      continue;
   
    TVector3 q = track->Map(cluster->GetPosition());

    Sx += cluster->GetCharge() * q.X() * q.X();
    Sy += cluster->GetCharge() * q.Y() * q.Y();
  }

  Double_t rmsR = sqrt(Sx / (track->GetChargeSum() * (1 - 3/clusterArray->size())));
  Double_t rmsY = sqrt(Sy / (track->GetChargeSum() * (1 - 3/clusterArray->size())));

  track->SetRMSW(rmsR);
  track->SetRMSH(rmsY);
  

  return true;
}



