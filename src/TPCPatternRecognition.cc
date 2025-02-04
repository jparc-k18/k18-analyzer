#include "TPCPatternRecognition.hh"
#include "TPCHit.hh"
#include "TPCHitCluster.hh"
#include "TPCPadHelper.hh"
#include "TPCRiemannTrack.hh"
//#include "TPCPIDAnalysis.hh"
#include "UserParamMan.hh"
#include "ConfMan.hh"
#include "Kinematics.hh"

#include <map>
#include <cmath>
#include <iostream>

#define DEBUG 0
#define scale 1.8

namespace {
  const auto& gUser   = UserParamMan::GetInstance();
  static const TVector3 tgtPos(0,0,-143);

}

TPCPatternRecognition::TPCPatternRecognition()
{}

TPCPatternRecognition::~TPCPatternRecognition()
{
  fTrackSearch = NULL;
}

Bool_t TPCPatternRecognition::Read(const std::vector<TPCHitContainer> &HitCont){

  fUsedHits = new std::vector<TPCHit*>;
  fUnusedHits = new std::vector<TPCHit*>;

  fCandHits = new std::vector<TPCHit*>;
  fBeamCandHits = new std::vector<TPCHit*>;
  fKuramaCandHits = new std::vector<TPCHit*>;
  
  fUsedHits->clear();
  fUnusedHits->clear();
  fCandHits->clear();  
  fBeamCandHits->clear();
  fKuramaCandHits->clear();  
  
  Initialize();

  for (Int_t layer = 0 ; layer < NumOfLayersTPC ; layer++) {
    for (Int_t ci = 0, n = HitCont[layer].size() ; ci < n ; ci++) {
      TPCHit* hit = HitCont[layer][ci];
      auto pos = hit->GetPosition();
      auto sigma = hit->GetSigma();

      if (fERadii != TVector3(-1,-1,-1)) {
	auto rVec = pos - fCutCenter;
	
	rVec =
	  TVector3(rVec.X()/fERadii.X(), rVec.Y()/fERadii.Y(), rVec.Z()/fERadii.Z());
	if (rVec.Mag() > 0.5) {
	
	  if (abs(sigma-3.9) < 1.5) {
	    if (!(tpc::Dead(layer, ci)) && hit->IsGood()){
	      if (hit->GetIsBeamCand())
		fBeamCandHits->push_back(hit);
	      else if (hit->GetIsKuramaCand())
		fKuramaCandHits->push_back(hit);
	      else 
		fCandHits->push_back(hit);
	    }
	  }
	}
      }
    }
  }

  return true;
}

Bool_t TPCPatternRecognition::ReadGeant4(const std::vector<TPCHitContainer> &HitCont){

  fUsedHits = new std::vector<TPCHit*>;
  fUnusedHits = new std::vector<TPCHit*>;

  fCandHits = new std::vector<TPCHit*>;
  fBeamCandHits = new std::vector<TPCHit*>;
  fKuramaCandHits = new std::vector<TPCHit*>;
  
  fUsedHits->clear();
  fUnusedHits->clear();
  fCandHits->clear();  
  fBeamCandHits->clear();
  fKuramaCandHits->clear();
  
  Initialize();

  for (Int_t layer = 0 ; layer < NumOfLayersTPC ; layer++) {
    for (Int_t ci = 0, n = HitCont[layer].size() ; ci < n ; ci++) {
      TPCHit* hit = HitCont[layer][ci];
      auto pos = hit->GetPosition();
      auto sigma = hit->GetSigma();

      if (fERadii != TVector3(-1,-1,-1)) {
	auto rVec = pos - fCutCenter;
	
	rVec =
	  TVector3(rVec.X()/fERadii.X(), rVec.Y()/fERadii.Y(), rVec.Z()/fERadii.Z());
	
	//if (rVec.Mag() > 1.0) {	  
	  if (!(tpc::Dead(layer, ci))) {
	    if(tpc::findPadID(pos.Z(),pos.X())!= -1000){
	      fCandHits->push_back(hit);
	    }
	  }
	  //}
      }
    }
  }

  return true;
}

Bool_t TPCPatternRecognition::ReadLinear(const std::vector<TPCHitContainer> &HitCont){

  fUsedHits = new std::vector<TPCHit*>;
  fUnusedHits = new std::vector<TPCHit*>;

  fLinearCandHits = new std::vector<TPCHit*>;
  fBeamCandHits = new std::vector<TPCHit*>;
  
  fUsedHits->clear();
  fUnusedHits->clear();

  fBeamCandHits->clear();
  fLinearCandHits->clear();  
  
  for (Int_t layer = 0 ; layer < NumOfLayersTPC ; layer++) {
    for (Int_t ci = 0, n = HitCont[layer].size() ; ci < n ; ci++) {
      TPCHit* hit = HitCont[layer][ci];
      if (!(tpc::Dead(layer, ci)) && hit->IsGood()){
	if (hit->GetIsBeamCand())
	  fBeamCandHits->push_back(hit);
	else 
	  fLinearCandHits->push_back(hit);
      }
    }
  }
  
  Initialize();

  return true;
}

Bool_t TPCPatternRecognition::Initialize() {

  fFitter = new TPCRiemannFitter();
  
  fHitArray = new std::vector<TPCHit*>;
  fRiemannTrackArray = new std::vector<TPCRiemannTrack*>;  
  fLinearTrackArray = new std::vector<TPCLinearTrack*>;  
  fHitClusterArray = new std::vector<TPCHitCluster*>;  

  fTrackSearch = new TPCRiemannTrackSearch();
  fLinearTrackSearch = new TPCLinearTrackSearch();
  
  fSortingMode = true;
  fSorting = TPCRiemannSort::kSortZ;
  
  fClusteringAngle = 60.;
  fClusteringMargin = 0.;
  
  fTrackSearch->SetSorting(fSorting);

  SetClusterCut(250,-250,250,-250);
  SetEllipsoidCut(TVector3(0,0,-143),TVector3(5,10,15),5);
    
  return true;				      
  
}
void TPCPatternRecognition::ExecuteLinear() {

  if (!fLinearTrackArray)
    std::cout << "Cannot Find RiemannTrackArray" << std::endl;
  fRiemannTrackArray -> clear();
  
  if (!fHitArray)
    std::cout << "Cannot find HitArray" << std::endl;
  fHitArray -> clear();

  std::vector<TPCLinearTrack *> LinearTemp;

  Int_t numUsedHits = 0;
  Int_t numHits;

  while (true) {
    fSorting = TPCRiemannSort::kSortDistance;
    fLinearTrackSearch->SetSorting(fSorting);
    BuildBeamTracks(fLinearTrackSearch, &LinearTemp, fSorting);      
    if (fBeamCandHits->size() == 0) break;
  }
  
  TPCLinearTrack *beamtrack;
  Int_t foundBeamTracks = LinearTemp.size();
  for (Int_t iTrack = 0 ; iTrack < foundBeamTracks ; iTrack++) {
      
    beamtrack = LinearTemp[iTrack];
    numHits = beamtrack->GetNumHits();
    if (numHits > 3 ) {
      numUsedHits += numHits;
      fLinearTrackArray -> push_back(beamtrack);
      //      beamtrack->SetIsBeam();
      for (Int_t iHit = 0 ; iHit < numHits; iHit++) {
	fHitArray->push_back(beamtrack->GetHit(iHit));
      }      
    }
    //    else {
    //      for (Int_t iHit = 0 ; iHit < numHits ; iHit++) 
    //	fCandHits->push_back(beamtrack->GetHit(iHit));      
    //    }
  }
  
  LinearTemp.clear();
  
  while (true) {

    fSorting = TPCRiemannSort::kSortDistance;
    fLinearTrackSearch->SetSorting(fSorting);
    BuildLinearTracks(fLinearTrackSearch, &LinearTemp, fSorting);
    if (fLinearCandHits->size() == 0) break;
  }

  TPCLinearTrack *lineartrack;
  Int_t foundLinearTracks = LinearTemp.size();
  for (Int_t iTrack = 0 ; iTrack < foundLinearTracks ; iTrack++) {

    lineartrack = LinearTemp[iTrack];
    numHits = lineartrack->GetNumHits();
    if (numHits > 6 && lineartrack->IsFitted()) {
      numUsedHits += numHits;
      fLinearTrackArray->push_back(lineartrack);
      for (Int_t iHit = 0 ; iHit < numHits ; iHit++) {
	fHitArray->push_back(lineartrack->GetHit(iHit));
      }
    }
  }
#if DEBUG
  std::cout << "Found Tracks : " << fLinearTrackArray->size() << std::endl;
#endif
  for (Int_t iTrack = 0 ; iTrack < fLinearTrackArray->size() ; iTrack++) {
    auto track = (TPCLinearTrack*)fLinearTrackArray->at(iTrack);
    track->FinalizeHits();
    HitClustering(track);    
  }  
}

void TPCPatternRecognition::Execute() {
  

  static const Bool_t DoPreTracking = gUser.GetParameter("DoPreTracking");

  static const Double_t minRMSW = gUser.GetParameter("MinRMSW");
  static const Double_t maxRMSW = gUser.GetParameter("MaxRMSW");
  static const Double_t minRMSH = gUser.GetParameter("MinRMSH");
  static const Double_t maxRMSH = gUser.GetParameter("MaxRMSH");

  static const Double_t ProxXCut = gUser.GetParameter("ProxXCut");
  static const Double_t ProxYCut = gUser.GetParameter("ProxYCut");
  static const Double_t DistCut = gUser.GetParameter("DistCut");



  if (!fRiemannTrackArray)
    std::cout << "Cannot Find RiemannTrackArray" << std::endl;
  fRiemannTrackArray -> clear();
  

  if (!fHitArray)
    std::cout << "Cannot find HitArray" << std::endl;
  fHitArray -> clear();


  std::vector<TPCRiemannTrack *> RiemannTemp;

  Int_t numBeamCandHits = fBeamCandHits->size();
  Int_t numKuramaCandHits = fKuramaCandHits->size();
  Int_t numCandHits = fCandHits->size();
  Int_t TotalHits = numBeamCandHits + numKuramaCandHits + numCandHits;
  Int_t numUsedHits = 0;
  Int_t numHits;

  
#if DEBUG
  std::cout << "" << std::endl;
  std::cout << "Starting pattern recognition..." << std::endl;
  std::cout << numBeamCandHits << "\t" << numKuramaCandHits << "\t" << numCandHits << std::endl;
  std::cout << " " << std::endl;
#endif

  while (true) {
    fSorting = TPCRiemannSort::kSortReverseZ;
    fTrackSearch->SetSorting(fSorting);
    //    fTrackSearch->SetRMSCut(1.5*minRMSW, 1.5*minRMSH, 1.5*maxRMSW, 1.5*maxRMSH);
    fTrackSearch->SetRMSCut(minRMSW, minRMSH, maxRMSW, maxRMSH);
    fTrackSearch->SetProxCut(ProxXCut,ProxYCut);
    fTrackSearch->SetDistCut(DistCut);
    //    fTrackSearch->SetProxCut(25,25);
    //    fTrackSearch->SetDistCut(75);
    BuildBeamTracks(fTrackSearch, &RiemannTemp, fSorting);      
    if (fBeamCandHits->size() == 0) break;
  }


  TPCRiemannTrack *beamtrack;
  Int_t foundbeamTracks = RiemannTemp.size();
  

  for (Int_t iTrack = 0 ; iTrack < foundbeamTracks ; iTrack++) {
      
    beamtrack = RiemannTemp[iTrack];
    numHits = beamtrack->GetNumHits();
    if (numHits > 4 ){
      numUsedHits += numHits;
      fRiemannTrackArray -> push_back(beamtrack);
      beamtrack->SetIsBeam();
      for (Int_t iHit = 0 ; iHit < numHits; iHit++) {
	fHitArray->push_back(beamtrack->GetHit(iHit));
      }      
    }
  }
  

  RiemannTemp.clear();


  TPCRiemannTrack *track;
  Int_t foundTracks;

  if (DoPreTracking) {
    while (true) {
      fSorting = TPCRiemannSort::kSortZ;
      fTrackSearch->SetSorting(fSorting);
      fTrackSearch->SetRMSCut(1.5*minRMSW, 1.5*minRMSH, 1.5*maxRMSW, 1.5*maxRMSH);
      //      fTrackSearch->SetRMSCut(minRMSW, minRMSH, maxRMSW, maxRMSH);
      fTrackSearch->SetProxCut(ProxXCut,ProxYCut);
      fTrackSearch->SetDistCut(DistCut);
      //      fTrackSearch->SetProxCut(20,15);
      //      fTrackSearch->SetDistCut(100);

      BuildTracks(fTrackSearch, &RiemannTemp, fSorting);      
      if (fCandHits->size() == 0) break;
    }

  
    foundTracks = RiemannTemp.size();

    for (Int_t iTrack = 0 ; iTrack < foundTracks ; iTrack++) {
      track = RiemannTemp[iTrack];
      numHits = track->GetNumHits();
      auto mean = track->GetMean();
      if (track->IsAccidental()){
	numUsedHits += numHits;
	fRiemannTrackArray -> push_back(track);
	track->SetIsAccBeam();
	for (Int_t iHit = 0 ; iHit < numHits; iHit++) {
	  fHitArray->push_back(track->GetHit(iHit));
	}            
      }
      //    else if (track->TrackLength() > 400 && abs(track->GetDip()) < 0.015 && track->GetRadius() > 1000){
      else if (track->TrackLength() > 400 && track->GetRadius() > 1000){
	numUsedHits += numHits;
	fRiemannTrackArray -> push_back(track);
	track->SetIsAccBeam();
	for (Int_t iHit = 0 ; iHit < numHits; iHit++) {
	  fHitArray->push_back(track->GetHit(iHit));
	}            
      }
      else {
	for (Int_t iHit = 0 ; iHit < numHits ; iHit++) 
	  fCandHits->push_back(track->GetHit(iHit));      
      }
    }
  
    RiemannTemp.clear();
  } //pretracking 
  //#endif

  while (true) {
    fSorting = TPCRiemannSort::kSortY;
    fTrackSearch->SetSorting(fSorting);
    fTrackSearch->SetRMSCut(minRMSW, minRMSH, maxRMSW, maxRMSH);
    fTrackSearch->SetProxCut(ProxXCut,ProxYCut);
    fTrackSearch->SetDistCut(DistCut);
      
    BuildTracks(fTrackSearch, &RiemannTemp, fSorting);      
    if (fCandHits->size() == 0) break;
  }
    

  for (Int_t iTrack = 0 ; iTrack < RiemannTemp.size() ; iTrack++) {
    track = RiemannTemp[iTrack];
    numHits = track->GetNumHits();
    auto mean = track->GetMean();
    if (numHits > 10 && track->TrackLength() > 200 && track->GetRadius() > 200) {
      if (abs(mean.Y()) < 180 ) {
	numUsedHits += numHits;
	fRiemannTrackArray -> push_back(track);
	track->SetIsScat();	
      }
      for (Int_t iHit = 0 ; iHit < numHits; iHit++) {
	fHitArray->push_back(track->GetHit(iHit));
      }            
    }
    else {
      for (Int_t iHit = 0 ; iHit < numHits ; iHit++) 
	fCandHits->push_back(track->GetHit(iHit));      
    }
  }

  RiemannTemp.clear();
  //  }
  while (true) {
    fSorting = TPCRiemannSort::kSortZ;
    fTrackSearch->SetSorting(fSorting);
    fTrackSearch->SetRMSCut(minRMSW, minRMSH, maxRMSW, maxRMSH);
    fTrackSearch->SetProxCut(ProxXCut,ProxYCut);
    fTrackSearch->SetDistCut(DistCut);
      
    BuildTracks(fTrackSearch, &RiemannTemp, fSorting);      
    if (fCandHits->size() == 0) break;
  }
    
  for (Int_t iTrack = 0 ; iTrack < RiemannTemp.size() ; iTrack++) {
    track = RiemannTemp[iTrack];
    numHits = track->GetNumHits();
    auto mean = track->GetMean();
    if (numHits > 10 && track->TrackLength() > 200 && track->GetRadius() > 100) {
      if (abs(mean.Y()) < 180 ) {
	numUsedHits += numHits;
	fRiemannTrackArray -> push_back(track);
	track->SetIsScat();	
      }
      for (Int_t iHit = 0 ; iHit < numHits; iHit++) {
	fHitArray->push_back(track->GetHit(iHit));
      }            
    }
    else {
      for (Int_t iHit = 0 ; iHit < numHits ; iHit++) 
	fCandHits->push_back(track->GetHit(iHit));      
    }
  }

  RiemannTemp.clear();

  while (true) {
    //    fSorting = TPCRiemannSort::kSortReverseAlpha;
    fSorting = TPCRiemannSort::kSortX;
    fTrackSearch->SetSorting(fSorting);
    fTrackSearch->SetRMSCut(minRMSW, minRMSH, maxRMSW, maxRMSH);
    fTrackSearch->SetProxCut(ProxXCut,ProxYCut);
    fTrackSearch->SetDistCut(DistCut);

    BuildTracks(fTrackSearch, &RiemannTemp, fSorting);      
    if (fCandHits->size() == 0) break;
  }
  
  for (Int_t iTrack = 0 ; iTrack < RiemannTemp.size() ; iTrack++) {
    track = RiemannTemp[iTrack];
    numHits = track->GetNumHits();
    auto mean = track->GetMean();
    if (numHits > 10 && track->TrackLength() > 200 && track->GetRadius() > 100) {
      if (abs(mean.Y()) < 180 ) {
	numUsedHits += numHits;
	fRiemannTrackArray -> push_back(track);
	track->SetIsScat();
      }
      for (Int_t iHit = 0 ; iHit < numHits; iHit++) {
	fHitArray->push_back(track->GetHit(iHit));
      }      
    }
    else {
      for (Int_t iHit = 0 ; iHit < numHits ; iHit++) 
	fCandHits->push_back(track->GetHit(iHit));      
    }
  }
  
  RiemannTemp.clear();


  while (true) {
    fSorting = TPCRiemannSort::kSortReverseAlpha;
    //    fSorting = TPCRiemannSort::kSortAlpha;
    fTrackSearch->SetSorting(fSorting);
    fTrackSearch->SetRMSCut(minRMSW, minRMSH, maxRMSW, maxRMSH);
    fTrackSearch->SetProxCut(ProxXCut,ProxYCut);
    fTrackSearch->SetDistCut(DistCut);

    BuildTracks(fTrackSearch, &RiemannTemp, fSorting);      
    if (fCandHits->size() == 0) break;
  }
  
  for (Int_t iTrack = 0 ; iTrack < RiemannTemp.size() ; iTrack++) {
    track = RiemannTemp[iTrack];
    numHits = track->GetNumHits();
    auto mean = track->GetMean();
    if (numHits > 10 && track->TrackLength() > 200 && track->GetRadius() > 50) {
      if (abs(mean.Y()) < 180 ) {
	numUsedHits += numHits;
	fRiemannTrackArray -> push_back(track);
	track->SetIsScat();
      }
      for (Int_t iHit = 0 ; iHit < numHits; iHit++) {
	fHitArray->push_back(track->GetHit(iHit));
      }      
    }
    else {
      for (Int_t iHit = 0 ; iHit < numHits ; iHit++) 
	fCandHits->push_back(track->GetHit(iHit));      
    }
  }
  
  RiemannTemp.clear();

  while (true) {
    fSorting = TPCRiemannSort::kSortPhi;
    fTrackSearch->SetSorting(fSorting);
    fTrackSearch->SetRMSCut(minRMSW, minRMSH, maxRMSW, maxRMSH);
    fTrackSearch->SetProxCut(ProxXCut,ProxYCut);
    fTrackSearch->SetDistCut(DistCut);

    BuildTracks(fTrackSearch, &RiemannTemp, fSorting);      
    if (fCandHits->size() == 0) break;
  }
  

  for (Int_t iTrack = 0 ; iTrack < RiemannTemp.size() ; iTrack++) {
    track = RiemannTemp[iTrack];
    numHits = track->GetNumHits();
    auto mean = track->GetMean();
    if (numHits > 10 && track->TrackLength() > 100 && track->GetRadius() > 50) {
      if (abs(mean.Y()) < 200 ) {
	numUsedHits += numHits;
	fRiemannTrackArray -> push_back(track);
	track->SetIsScat();	
      }
      for (Int_t iHit = 0 ; iHit < numHits; iHit++) {
	fHitArray->push_back(track->GetHit(iHit));
      }            
    }
    else {
      for (Int_t iHit = 0 ; iHit < numHits ; iHit++) 
	fCandHits->push_back(track->GetHit(iHit));      
    }
    
  }
  

  RiemannTemp.clear();

  while (true) {
    fSorting = TPCRiemannSort::kSortY;
    fTrackSearch->SetSorting(fSorting);
    fTrackSearch->SetRMSCut(scale*minRMSW, scale*minRMSH, scale*maxRMSW, scale*maxRMSH);
    fTrackSearch->SetProxCut(scale*ProxXCut,scale*ProxYCut);
    fTrackSearch->SetDistCut(scale*DistCut);

    BuildTracks(fTrackSearch, &RiemannTemp, fSorting);      
    if (fCandHits->size() == 0) break;
  }

  for (Int_t iTrack = 0 ; iTrack < RiemannTemp.size() ; iTrack++) {
    track = RiemannTemp[iTrack];
    numHits = track->GetNumHits();
    auto mean = track->GetMean();
    if (numHits > 6 && track->TrackLength() > 50 && track->GetRadius() > 50) {
      if (abs(mean.Y()) < 200 ) {
	numUsedHits += numHits;
	fRiemannTrackArray -> push_back(track);
	track->SetIsScat();	
      }
      for (Int_t iHit = 0 ; iHit < numHits; iHit++) {
	fHitArray->push_back(track->GetHit(iHit));
      }            
    }
    else {
      for (Int_t iHit = 0 ; iHit < numHits ; iHit++) 
	fCandHits->push_back(track->GetHit(iHit));      
    }
    
  }

  RiemannTemp.clear();

  while (true) {
    fSorting = TPCRiemannSort::kSortZ;
    fTrackSearch->SetSorting(fSorting);
    fTrackSearch->SetRMSCut(scale*minRMSW, scale*minRMSH, scale*maxRMSW, scale*maxRMSH);
    fTrackSearch->SetProxCut(scale*ProxXCut,scale*ProxYCut);
    fTrackSearch->SetDistCut(scale*DistCut);
    //    fTrackSearch->SetRMSCut(minRMSW, minRMSH, maxRMSW, maxRMSH);
    //    fTrackSearch->SetProxCut(ProxXCut,ProxYCut);
    //    fTrackSearch->SetDistCut(DistCut);

    BuildTracks(fTrackSearch, &RiemannTemp, fSorting);      
    if (fCandHits->size() == 0) break;
  }

  for (Int_t iTrack = 0 ; iTrack < RiemannTemp.size() ; iTrack++) {
    track = RiemannTemp[iTrack];
    numHits = track->GetNumHits();
    auto mean = track->GetMean();
    if (numHits > 6 && track->TrackLength() > 50 && track->GetRadius() > 50) {
      if (abs(mean.Y()) < 200 ) {
	numUsedHits += numHits;
	fRiemannTrackArray -> push_back(track);
	track->SetIsScat();	
      }
      for (Int_t iHit = 0 ; iHit < numHits; iHit++) {
	fHitArray->push_back(track->GetHit(iHit));
      }            
    }
    else {
      for (Int_t iHit = 0 ; iHit < numHits ; iHit++) 
	fCandHits->push_back(track->GetHit(iHit));      
    }
    
  }
  

  RiemannTemp.clear();

  
  while (true) {
    fSorting = TPCRiemannSort::kSortX;
    fTrackSearch->SetSorting(fSorting);
    fTrackSearch->SetRMSCut(scale*minRMSW, scale*minRMSH, scale*maxRMSW, scale*maxRMSH);
    fTrackSearch->SetProxCut(scale*ProxXCut,scale*ProxYCut);
    fTrackSearch->SetDistCut(scale*DistCut);
    //    fTrackSearch->SetRMSCut(minRMSW, minRMSH, maxRMSW, maxRMSH);
    //    fTrackSearch->SetProxCut(ProxXCut,ProxYCut);
    //    fTrackSearch->SetDistCut(DistCut);

    BuildTracks(fTrackSearch, &RiemannTemp, fSorting);      
    if (fCandHits->size() == 0) break;
  }

  for (Int_t iTrack = 0 ; iTrack < RiemannTemp.size() ; iTrack++) {
    track = RiemannTemp[iTrack];
    numHits = track->GetNumHits();
    auto mean = track->GetMean();
    if (numHits > 6 && track->TrackLength() > 50 && track->GetRadius() > 50) {
      if (abs(mean.Y()) < 200 ) {
	numUsedHits += numHits;
	fRiemannTrackArray -> push_back(track);
	track->SetIsScat();	
      }
      for (Int_t iHit = 0 ; iHit < numHits; iHit++) {
	fHitArray->push_back(track->GetHit(iHit));
      }            
    }
    else {
      for (Int_t iHit = 0 ; iHit < numHits ; iHit++) 
	fCandHits->push_back(track->GetHit(iHit));      
    }
    
  }
  

  RiemannTemp.clear();

  while (true) {
    fSorting = TPCRiemannSort::kSortAlpha;
    fTrackSearch->SetSorting(fSorting);
    fTrackSearch->SetRMSCut(scale*minRMSW, scale*minRMSH, scale*maxRMSW, scale*maxRMSH);
    fTrackSearch->SetProxCut(scale*ProxXCut,scale*ProxYCut);
    fTrackSearch->SetDistCut(scale*DistCut);
    //    fTrackSearch->SetRMSCut(minRMSW, minRMSH, maxRMSW, maxRMSH);
    //    fTrackSearch->SetProxCut(ProxXCut,ProxYCut);
    //    fTrackSearch->SetDistCut(DistCut);
    
    BuildTracks(fTrackSearch, &RiemannTemp, fSorting);      
    if (fCandHits->size() == 0) break;
  }

  for (Int_t iTrack = 0 ; iTrack < RiemannTemp.size() ; iTrack++) {
    track = RiemannTemp[iTrack];
    numHits = track->GetNumHits();
    auto mean = track->GetMean();
    if (numHits > 6 && track->TrackLength() > 50 && track->GetRadius() > 25) {
      if (abs(mean.Y()) < 180 ) {
	numUsedHits += numHits;
	fRiemannTrackArray -> push_back(track);
	track->SetIsScat();	
      }
      for (Int_t iHit = 0 ; iHit < numHits; iHit++) {
	fHitArray->push_back(track->GetHit(iHit));
      }            
    }
  }
  
  RiemannTemp.clear();

  
  


  foundTracks = fRiemannTrackArray->size();
  //  numUsedHits = fHitArray->size();
#if DEBUG
  std::cout << Form("Found %d Tracks, used %d(/%d) hits", foundTracks, numUsedHits, TotalHits) << std::endl;
  for (int i = 0; i < foundTracks; i++) {
    auto track = fRiemannTrackArray->at(i);
    std::cout << Form("%d Hits were used in Track %d", track->GetNumHits() , i) << std::endl;
  }
#endif

  TVector3 vertex = FindVertex(fRiemannTrackArray);

  for (Int_t iTrack = 0 ; iTrack < fRiemannTrackArray->size() ; iTrack++){
    
    auto track = (TPCRiemannTrack*)fRiemannTrackArray->at(iTrack);

    auto mean = track->GetMean();

    track->SetTrackID(iTrack);

    //    if (abs(vertex.Z()) < 250) track->DetermineCharge(vertex);
    //    else     track->DetermineCharge(tgtPos);

    track->DetermineCharge(tgtPos);

    if (!HitClustering(track))
      fRiemannTrackArray->erase(fRiemannTrackArray->begin()+iTrack);

    track->FinalizeClusters();

  }


  
  RiemannTemp.clear();

  delete fTrackSearch;
}

void TPCPatternRecognition::BuildTracks(TPCRiemannTrackSearch *TrackSearch,
					std::vector<TPCRiemannTrack*> *TrackletCont, Int_t sorting){

  fUnusedHits->clear();

  Int_t nHitIn = fCandHits->size();
  if (nHitIn == 0) return;
  
  TPCRiemannTrack *tracklet = TrackSearch->NewTrack(*fCandHits);  
  TrackSearch->InitTrack(*fCandHits, *fUnusedHits, tracklet);

#if DEBUG
  fUsedHits = tracklet->GetHits();
  std::cout << fUsedHits->size() << "\t" << fCandHits->size() << "\t" << fUnusedHits->size() << std::endl;
#endif
  
  //  if (tracklet->GetNumHits() > 6) {
  //  if (Initialized) {    
    Bool_t TrackExtended = TrackSearch->ExtendTrack(*fCandHits, *fUnusedHits, tracklet);    
    
#if DEBUG
    fUsedHits = tracklet->GetHits();
    std::cout << fUsedHits->size() << "\t" << fCandHits->size() << "\t" << fUnusedHits->size() << std::endl;
#endif
    //    if (tracklet->GetNumHits() > 16) {
    Bool_t TrackExtrap = TrackSearch->ExtrapTrack(*fCandHits, *fUnusedHits, tracklet);

#if DEBUG
    fUsedHits = tracklet->GetHits();
    std::cout << fUsedHits->size() << "\t" << fCandHits->size() << "\t" << fUnusedHits->size() << std::endl;
#endif
    //    }
    Bool_t TrackConfirm = TrackSearch->ConfirmTrack(*fCandHits, *fUnusedHits, tracklet);
    
    //  }

  TrackletCont->push_back(tracklet);
    
#if DEBUG
  std::cout << tracklet->TrackLength() << "\t" << tracklet->GetRadius() << std::endl;
#endif

}

void TPCPatternRecognition::BuildLinearTracks(TPCLinearTrackSearch *TrackSearch,
					      std::vector<TPCLinearTrack*> *TrackletCont, Int_t sorting) {

  Int_t nHitIn = fLinearCandHits->size();
  if (nHitIn == 0)
    return;

  TPCLinearTrack *tracklet = TrackSearch->NewTrack(*fLinearCandHits);
  TrackSearch->BuildTrack(*fLinearCandHits, *fUnusedHits, tracklet);

  TrackletCont->push_back(tracklet);

}

void TPCPatternRecognition::BuildBeamTracks(TPCLinearTrackSearch *TrackSearch,
					    std::vector<TPCLinearTrack*> *TrackletCont, Int_t sorting){
  Int_t nHitIn = fBeamCandHits->size();
  
  if (nHitIn == 0)
    return;
  
  TPCLinearTrack *tracklet = TrackSearch->NewTrack(*fBeamCandHits);
  TrackSearch->BuildTrack(*fBeamCandHits, *fUnusedHits, tracklet);

  TrackletCont->push_back(tracklet);
}

void TPCPatternRecognition::BuildBeamTracks(TPCRiemannTrackSearch *TrackSearch,
					    std::vector<TPCRiemannTrack*> *TrackletCont, Int_t sorting){
  Int_t nHitIn = fBeamCandHits->size();
  
  if (nHitIn == 0)
    return;
  
  TPCRiemannTrack *tracklet = TrackSearch->NewTrack(*fBeamCandHits);  
  
  TrackSearch->InitTrack(*fBeamCandHits, *fUnusedHits, tracklet);
  if (tracklet->GetNumHits() > 1) {     
    Bool_t TrackExtended = TrackSearch->ExtendTrack(*fBeamCandHits, *fUnusedHits, tracklet);
    Bool_t TrackExtrap = TrackSearch->ExtrapTrack(*fBeamCandHits, *fUnusedHits, tracklet);
    Bool_t TrackConfirm = TrackSearch->ConfirmTrack(*fBeamCandHits, *fUnusedHits, tracklet);
  }
  
  TrackletCont->push_back(tracklet);

}

void TPCPatternRecognition::BuildKuramaTracks(TPCRiemannTrackSearch *TrackSearch,
					    std::vector<TPCRiemannTrack*> *TrackletCont, Int_t sorting){
  Int_t nHitIn = fKuramaCandHits->size();
  
  if (nHitIn == 0)
    return;

  TPCRiemannTrack *tracklet = TrackSearch->NewTrack(*fKuramaCandHits);  
  
  TrackSearch->InitTrack(*fKuramaCandHits, *fUnusedHits, tracklet);
  if (tracklet->GetNumHits() > 2) { 
    
    Bool_t TrackExtended = TrackSearch->ExtendTrack(*fKuramaCandHits, *fUnusedHits, tracklet);   
    Bool_t TrackExtrap = TrackSearch->ExtrapTrack(*fKuramaCandHits, *fUnusedHits, tracklet);
    Bool_t TrackConfirm = TrackSearch->ConfirmTrack(*fKuramaCandHits, *fUnusedHits, tracklet);
  }
  
  TrackletCont->push_back(tracklet);
  
}


TVector3 TPCPatternRecognition::FindVertex(std::vector<TPCRiemannTrack *> *trackCont, Int_t nIter) {

  auto TestZ = [trackCont](TVector3 &v) {
    Double_t s = 0;
    Int_t numUsedTracks = 0;

    v.SetX(0);  v.SetY(0);

    auto numTracks = trackCont->size();

    for (auto iTrack = 0 ; iTrack < numTracks ; iTrack++) {
      TPCRiemannTrack *track = (TPCRiemannTrack*)trackCont->at(iTrack);
      if (track->GetTrackProperty() != 2) continue;
      
      TVector3 p,q;
      Bool_t extrapolated = track->ExtrapToZ(v.Z(),p);      
      if (extrapolated == false)
	continue;
      //      Double_t extrap = track->ExtrapByMap(v, p, q);
      
      v.SetX((numUsedTracks*v.X() + p.X())/(numUsedTracks+1));
      v.SetY((numUsedTracks*v.Y() + p.Y())/(numUsedTracks+1));
      
      if (numUsedTracks != 0)
	s = (double)numUsedTracks/(numUsedTracks+1)*s + (v-p).Mag()/numUsedTracks;

      if (s < 50)
	numUsedTracks++;
      
    }
    return s;
  };

  Double_t z0 = 0.;
  Double_t dz = 2.5;
  Double_t s0 = 1.0e+08;

  const Int_t numSamples = 200;
  Int_t halfOfSamples = (numSamples)/2;

  Double_t zArray[numSamples] = {0};
  for (Int_t iSample = 0 ; iSample < numSamples ; iSample++)
    zArray[iSample] = (iSample - halfOfSamples) * dz + z0;

  for (auto z : zArray) {
    TVector3 v(0,0,z);
    Double_t s = TestZ(v);

    if (s < s0) {
      s0 = s;
      z0 = z;
    }
  }

  while (nIter > 0) {
    dz = dz/halfOfSamples;
    for (Int_t iSample = 0 ; iSample < numSamples ; iSample++)
      zArray[iSample] = (iSample - halfOfSamples) * dz + z0;

    for (auto z : zArray) {
      TVector3 v(0,0,z);
      Double_t s = TestZ(v);

      if (s < s0) {
	s0 = s;
	z0 = z;
      }
    }

    nIter--;
  }
  
  TVector3 v(0,0,z0);
  TestZ(v);

  return v;
}


TPCHitCluster *TPCPatternRecognition::NewCluster(TPCHit *hit) {

  TPCHitCluster *cluster = new TPCHitCluster();
  cluster->AddHit(hit);
  cluster->SetClusterID(fHitClusterArray->size());

  return cluster;
}

Bool_t TPCPatternRecognition::HitClustering(TPCRiemannTrack *track){
    
  auto trackHits = track->GetHits();
  auto numHits = trackHits->size();

    
  std::sort(trackHits->begin(), trackHits->end(), SortByZ());
  //  if (track->Charge() > 0)
  //      std::sort(trackHits->begin(), trackHits->end(), SortByIncLength(track));
  //    else
  //      std::sort(trackHits->begin(), trackHits->end(), SortByDecLength(track));
    
  auto currentHit = trackHits->at(0);
  Bool_t clusterType = CheckClusterType(track, currentHit, nullptr);
  
  TPCHitCluster *lastCluster = nullptr;
  TPCHitCluster *frontCluster = nullptr;

  lastCluster = NewCluster(currentHit);

  Int_t rowMin = currentHit->GetRow();
  Int_t rowMax = currentHit->GetRow();

  Int_t layerMin = currentHit->GetLayer();
  Int_t layerMax = currentHit->GetLayer();

  //  lastCluster->SetRow(-1);
  lastCluster->SetLayer(currentHit->GetLayer());
  lastCluster->SetType(clusterType);
  
  std::vector<TPCHitCluster *> lastClusterBuilder;
  std::vector<TPCHitCluster *> frontClusterBuilder;

  //  SetClusterLength(lastCluster);
  lastClusterBuilder.push_back(lastCluster);

  for (auto iHit = 1 ; iHit < numHits; iHit++) {

    bool separateCluster = false;
    
    currentHit = trackHits->at(iHit);

    auto row = currentHit -> GetRow();
    auto layer = currentHit->GetLayer();

    auto phi = currentHit->GetPhi();
    //    auto alpha = currentHit->GetPhi();
    
    //    if (abs(phi) > TMath::PiOver2() && layer < 11) 
    if (abs(phi) > TMath::PiOver2()) 
      separateCluster = true;
    else
      separateCluster = false;

    if (!separateCluster) {
      
      bool createNewCluster = true;

      for (auto cluster : lastClusterBuilder) {
	if (layer == cluster->GetLayer()) {
	  if (CheckIsContinuousHits(cluster, currentHit)){
	    createNewCluster = false;
	    cluster->AddHit(currentHit);
	  }
	  else
	    createNewCluster = true;
	}
      }
      
      if (row < rowMin) rowMin = row;
      if (row > rowMax) rowMax = row;

      if (createNewCluster) {
	  
	clusterType = CheckClusterType(track, currentHit, (fClusteringMargin > 0 ? trackHits->at(iHit-1) : nullptr));	

	if (lastCluster != nullptr) {
	  track->AddHitCluster(lastCluster);
	  lastCluster->SetIsStable(true);
	}
	
	lastCluster = NewCluster(currentHit);
	//	lastCluster->SetRow(-1);
	lastCluster->SetLayer(layer);
	lastCluster->SetType(clusterType);
	//	SetClusterLength(lastCluster);
	lastClusterBuilder.push_back(lastCluster);  
      }
    }

    else {
      bool createNewCluster = true;

      for (auto cluster : frontClusterBuilder) {
	if (layer == cluster->GetLayer()) {
	  if (CheckIsContinuousHits(cluster, currentHit)){
	    createNewCluster = false;
	    cluster->AddHit(currentHit);
	  }
	  else
	    createNewCluster = true;
	}
      }
	  
      if (row < rowMin) rowMin = row;
      if (row > rowMax) rowMax = row;

      if (createNewCluster) {

	clusterType = CheckClusterType(track, currentHit, (fClusteringMargin > 0 ? trackHits->at(iHit-1) : nullptr));	

	if (frontCluster != nullptr) {
	  track->AddHitCluster(frontCluster);
	  frontCluster->SetIsStable(true);
	}
	frontCluster = NewCluster(currentHit);
	//	frontCluster->SetRow(-1);
	frontCluster->SetLayer(layer);
	frontCluster->SetType(clusterType);
	//	SetClusterLength(frontCluster);
	frontClusterBuilder.push_back(frontCluster);
      }
    }    
  }

  if (lastCluster != nullptr) {
    track->AddHitCluster(lastCluster);
    lastCluster->SetIsStable(true);
    lastCluster = nullptr;
  }

  if (frontCluster != nullptr) {
    track->AddHitCluster(frontCluster);
    frontCluster->SetIsStable(true);
    frontCluster = nullptr;
  }

  if (track->GetNumClusters() < 4)
    return false;

  auto clusterArray = track->GetClusters();
  for (auto cluster : *clusterArray) {

    auto pos = cluster->GetPosition();
    auto diff = (pos-tgtPos);
    if (diff.Mag() < 25. || cluster->GetLayer() < 2) cluster->SetIsStable(false);
  }


  if (fFitter->FitCluster(track) == false) {
    track->SetIsLine();
    
  }

  std::sort(clusterArray->begin(), clusterArray->end(), SortByDistInv(tgtPos));
  
  auto SetClusterLength = [track](TPCHitCluster *cluster) {

    auto row = cluster->GetRow();
    auto layer = cluster->GetLayer();
    auto pos = cluster->GetPosition();
    auto alpha = track->AlphaAtPosition(pos);
    auto dir = track->Direction(alpha);

    auto phi = (TMath::ATan2(pos.X() , pos.Z() + 143));
    auto theta = (TMath::ATan2((dir.X()), (dir.Z())));
    //    phi = TMath::PiOver2() - phi;
    
    Double_t angle = (theta-phi);

    //    if (abs(angle) > TMath::PiOver2())
    //      angle = angle + TMath::PiOver2();
    
    TVector3 prevPos;
    TVector3 postPos;
    
    Double_t length = 1;
    
    Double_t padlength = tpc::padParameter[layer][5];

    //    padlength = abs(padlength/TMath::Cos(angle));
    //    if (length > 2*tpc::padParameter[layer][5])
    //      padlength = tpc::padParameter[layer][5];
    //      length = tpc::padParameter[layer][5];
    
    if (abs(angle) < TMath::Pi()/4.) {
      Double_t prevZ = pos.Z() - padlength/2.;
      Double_t postZ = pos.Z() + padlength/2.;
      track->ExtrapToZ(prevZ, alpha, prevPos);
      track->ExtrapToZ(postZ, alpha, postPos);
      length = TMath::Abs(track->Map(prevPos).Z() - track->Map(postPos).Z());
    }
    else {
      Double_t prevX = pos.X() - padlength/2.;
      Double_t postX = pos.X() + padlength/2.;
      track->ExtrapToX(prevX, alpha, prevPos);
      track->ExtrapToX(postX, alpha, postPos);
      length = TMath::Abs(track->Map(prevPos).Z() - track->Map(postPos).Z());
    }

    if (length > sqrt(2)*tpc::padParameter[layer][5])
      length = sqrt(2)*tpc::padParameter[layer][5];
    
    //    std::cout << length << "\t" << angle << std::endl;
    cluster->SetLength(length);
					   
  };

  for (auto cluster : *track->GetClusters()) {

    SetClusterLength(cluster);
    //    std::cout << cluster->GetLength() << "\t" << cluster->GetType() << std::endl;
    
  }
  
  return true;
  
}

Bool_t TPCPatternRecognition::HitClustering(TPCLinearTrack *track){
    
  auto trackHits = track->GetHits();
  auto numHits = trackHits->size();

  std::sort(trackHits->begin(), trackHits->end(), SortByZ());
  
  auto currentHit = trackHits->at(0);
  TPCHitCluster *lastCluster = nullptr;
  TPCHitCluster *frontCluster = nullptr;

  lastCluster = NewCluster(currentHit);

  Int_t rowMin = currentHit->GetRow();
  Int_t rowMax = currentHit->GetRow();

  lastCluster->SetLayer(currentHit->GetLayer());
  
  std::vector<TPCHitCluster *> lastClusterBuilder;
  std::vector<TPCHitCluster *> frontClusterBuilder;

  lastClusterBuilder.push_back(lastCluster);

  for (auto iHit = 1 ; iHit < numHits; iHit++) {

    bool separateCluster = false;
    
    currentHit = trackHits->at(iHit);

    auto row = currentHit -> GetRow();
    auto layer = currentHit->GetLayer();

    auto phi = currentHit->GetPhi();
    
    if (abs(phi) > TMath::PiOver2()) 
      separateCluster = true;
    else
      separateCluster = false;

    if (!separateCluster) {
      
      bool createNewCluster = true;
      for (auto cluster : lastClusterBuilder) {
	if (layer == cluster->GetLayer()) {
	  if (CheckIsContinuousHits(cluster, currentHit)){
	    createNewCluster = false;
	    cluster->AddHit(currentHit);
	  }
	  else
	    createNewCluster = true;
	}
      }
      
      if (row < rowMin) rowMin = row;
      if (row > rowMax) rowMax = row;

      if (createNewCluster) {
	  
	if (lastCluster != nullptr) {
	  track->AddHitCluster(lastCluster);
	  lastCluster->SetIsStable(true);
	}
	
	lastCluster = NewCluster(currentHit);
	//	lastCluster->SetRow(-1);
	lastCluster->SetLayer(layer);
	lastClusterBuilder.push_back(lastCluster);  
      }
   } 

    else {
      bool createNewCluster = true;

      for (auto cluster : frontClusterBuilder) {
	if (layer == cluster->GetLayer()) {
	  if (CheckIsContinuousHits(cluster, currentHit)){
	    createNewCluster = false;
	    cluster->AddHit(currentHit);
	  }
	  else
	    createNewCluster = true;
	}
      }
	  
      if (row < rowMin) rowMin = row;
      if (row > rowMax) rowMax = row;

      if (createNewCluster) {

	if (frontCluster != nullptr) {
	  track->AddHitCluster(frontCluster);
	  frontCluster->SetIsStable(true);
	}
	frontCluster = NewCluster(currentHit);
	//	frontCluster->SetRow(-1);
	frontCluster->SetLayer(layer);
	frontClusterBuilder.push_back(frontCluster);
      }
    }    
  }
  
  if (lastCluster != nullptr) {
    track->AddHitCluster(lastCluster);
    lastCluster->SetIsStable(true);
    lastCluster = nullptr;
  }

  if (frontCluster != nullptr) {
    track->AddHitCluster(frontCluster);
    frontCluster->SetIsStable(true);
    frontCluster = nullptr;
  }

  Int_t numCluster = track->GetNumClusters();

  if (numCluster < 4)
    return false;

  if (fFitter->FitPlane(track) == false) 
    track->SetIsFitted(false);
  
  return true;
  
}

Bool_t TPCPatternRecognition::CheckIsContinuousHits(TPCHitCluster *cluster, TPCHit *hit) {

  auto hits = cluster -> GetHits();

  std::vector<Int_t> numbers;

  for (auto clhit : *hits)
    numbers.push_back(clhit->GetRow());

  numbers.push_back(hit->GetRow());

  std::sort(numbers.begin(), numbers.end());
  numbers.erase(unique(numbers.begin(), numbers.end()), numbers.end());

  if (numbers.size() < 2) {
    //    cluster->SetIsContinuousHits(true);
    return true;
  }

  for (auto i = 0 ; i < numbers.size() -1 ; i++)
    if (!(numbers[i]+1 == numbers[i+1]))
      return false;

  //  cluster->SetIsContinuousHits(true);

  return true;
}

Bool_t TPCPatternRecognition::CheckClusterType(TPCRiemannTrack *track, TPCHit *currentHit, TPCHit *prevHit){

  TVector3 q;
  Double_t alpha;
  track->ExtrapToPointAlpha(currentHit->GetPosition(), q, alpha);

  auto directionChangeAngle = fClusteringAngle * TMath::DegToRad();
  auto normAlpha = TMath::Abs(std::fmod(TMath::Abs(alpha), TMath::Pi()) - TMath::PiOver2());
  auto isPerp = (normAlpha > TMath::PiOver2() - directionChangeAngle);

  Double_t prevAlpha;
  if (prevHit != nullptr) {

    track->ExtrapToPointAlpha(prevHit->GetPosition(),q,prevAlpha);

    auto prevNormAlpha = TMath::Abs(std::fmod(TMath::Abs(alpha), TMath::Pi()) - TMath::PiOver2());
    auto isPerpPrev = (prevNormAlpha > TMath::PiOver2() - directionChangeAngle);

    if (isPerp == isPerpPrev)
      return isPerp;
    
    auto margin = fClusteringMargin * TMath::DegToRad();
    auto diffNormAlpha = TMath::Abs(normAlpha - TMath::PiOver2() + directionChangeAngle);
    
    if (diffNormAlpha < margin)
      return isPerpPrev;
  }
  
  return isPerp;
}
       
void TPCPatternRecognition::SetEllipsoidCut(TVector3 vector, TVector3 radii, Double_t margin){
  fCutCenter = vector;
  fERadii = radii;
  fCutMargin = margin;
}

void TPCPatternRecognition::SetClusterCut(Double_t left, Double_t right, Double_t top, Double_t bottom){
  fCCLeft = left;
  fCCRight = right;
  fCCTop = top;
  fCCBottom = bottom;
}


