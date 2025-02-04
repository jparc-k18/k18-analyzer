#include "TPCLinearTrackSearch.hh"
#include "TPCRiemannTrackSearch.hh"
#include "TPCHit.hh"
#include "TPCLinearTrack.hh"
#include "TPCRiemannSort.hh"
#include "TPCRiemannFitter.hh"
#include "TPCPadHelper.hh"
#include "DetectorID.hh"

#include <algorithm>
#include <iostream>
#include <iomanip>

#include <TMath.h>

#include "UserParamMan.hh"

namespace{
  const auto& gUser   = UserParamMan::GetInstance();
  const TVector3 tgtPos(0,0,-143);
}

TPCLinearTrackSearch::TPCLinearTrackSearch() {
  Initialize();
}

void TPCLinearTrackSearch::Initialize() {

  fFitter = new TPCRiemannFitter();

}

TPCLinearTrack *TPCLinearTrackSearch::NewTrack(std::vector<TPCHit*> &candHits) {

  TPCLinearTrack *track = new TPCLinearTrack();

  SetPosition(tgtPos);
  SortHits(candHits);
  track->AddHit(candHits.front());
  candHits.erase(candHits.begin() + 0);

  return track;
}

Bool_t TPCLinearTrackSearch::BuildTrack(std::vector<TPCHit*> &candHits, std::vector<TPCHit*> &unusedHits,
				       TPCLinearTrack *track) {
  
  std::vector<TPCHit*> candidates;
  std::vector<TPCHit*> badHits;

  GetAdjacentHits(track, candHits, candidates, unusedHits);

  Int_t numHits = candidates.size();
  Int_t count = 0;
  Double_t besttrackQuality = 0;
  
  while (numHits != 0) {

    fSorting = TPCRiemannSort::kSortZ;
    //fSorting = TPCRiemannSort::kSortY;
    SetPosition(track->GetMean());
    SortHits(candidates);

    for (Int_t iHit = 0 ; iHit < numHits ; iHit++) {

      TPCHit *hit = candidates.back();
      candidates.pop_back();

      Double_t trackQuality  = 0;

      Bool_t trackSurvive = CorrelateHT(track,hit, trackQuality);

      if (trackSurvive == true ) {
	trackQuality = besttrackQuality;
	track->AddHit(hit);
	//	std::cout << "Add Hits" << std::endl;
	//	if (track->GetNumHits() > 6)
	  //	  fFitter->FitLinear(track);

	fFitter->FitPlane(track);
	
      }
      
      else {
	badHits.push_back(hit);
      }
    }

    for (auto unusedHit : unusedHits)
      candHits.push_back(unusedHit);
    unusedHits.clear();

    GetAdjacentHits(track, candHits, candidates, unusedHits);

    if (count++ > 200) {
      for (auto cand : candidates)
	candHits.push_back(cand);
      candidates.clear();
      break;
    }

    numHits = candidates.size();      
  }

  for (int i = 0 ; i < unusedHits.size() ; i++)
    candHits.push_back(unusedHits.at(i));
  unusedHits.clear();

  for (auto badHit : badHits)
    candHits.push_back(badHit);
  badHits.clear();

  
  return true;
}

Bool_t TPCLinearTrackSearch::CorrelateHT(TPCLinearTrack *track, TPCHit *hit, Double_t &quality) {

  Bool_t surviveProx = false;
  Bool_t survivePerp = false;
  Bool_t surviveRMS = false;
  
  Bool_t CorrHTProx = CorrelateHTProx(track, hit, surviveProx, quality);
  Bool_t CorrHTPerp = CorrelateHTPerp(track, hit, survivePerp, quality);
  Bool_t CorrHTRMS = CorrelateHTRMS(track, hit, surviveRMS, quality);

  if ( survivePerp == false || surviveRMS == false) {
    //if (surviveProx == false || survivePerp == false || surviveRMS == false) {
    return false;
    //  if ( survivePerp == false || surviveRMS == false) {
    //    std::cout << "Failled in CorrHTRMS" << std::endl;
  }

  return true;
}

Bool_t TPCLinearTrackSearch::CorrelateHTProx(TPCLinearTrack *track, TPCHit *hit, Bool_t &survive, Double_t &quality){

  Int_t numHitsinTrack = track->GetNumHits();
  if (numHitsinTrack > 6){
    survive = true;
    return false;
  }
  
  survive = false;

  Double_t fProxCutX;
  Double_t fProxCutY;
  Double_t fProxCutZ;
  
  //Linear Track Adjacent cut condition
  if (hit->GetLayer() < 10){
    fProxCutX = 15.0;
    fProxCutZ = 12.0;
  }
  else {
    fProxCutX = 18.0;
    fProxCutZ = 15.0;
  }
  
  fProxCutY = 20;

  TVector3 hitPos = hit->GetPosition();

  for (int iHit = 0 ; iHit < numHitsinTrack ; iHit++) {

    TPCHit *trackhit = track->GetHit(iHit);
    TVector3 trackhitPos = trackhit->GetPosition();

    Double_t dx = abs(hitPos.X() - trackhitPos.X());
    Double_t dy = abs(hitPos.Y() - trackhitPos.Y());
    Double_t dz = abs(hitPos.Z() - trackhitPos.Z());

    Double_t xcorr = (fProxCutX - dx)/fProxCutX;
    Double_t ycorr = (fProxCutY - dy)/fProxCutY;
    Double_t zcorr = (fProxCutZ - dz)/fProxCutZ;

    //    std::cout << xcorr << "\t" << ycorr << "\t" << zcorr << std::endl;
    
    if ( xcorr >= 0 && ycorr >= 0 && zcorr >= 0){
    //if ( xcorr >= 0 && ycorr >= 0 ){

      survive = true;
      quality = TMath::Sqrt(xcorr*xcorr+ycorr*ycorr+zcorr*zcorr + quality)/(xcorr+ycorr+zcorr+quality);

      return true;
    }    
  }
  return true;
}

Bool_t TPCLinearTrackSearch::CorrelateHTPerp(TPCLinearTrack *track, TPCHit *hit, Bool_t &survive, Double_t &quality){

  Int_t numHitsinTrack = track->GetNumHits();

  if (numHitsinTrack < 4) {
    survive = true;
    return false;
  }
  
  survive = false;
  
  Double_t fPerpLineCut = 25;
  Double_t fPerpPlaneCut = 20;
  
  if (track->IsFitted() == false){
    fFitter->FitPlane(track);
  }

  //  TVector3 PerpLine = fFitter->PerpLine(track, hit);
  //  TVector3 PerpPlane = fFitter->PerpPlane(track, hit);
  TVector3 PerpLine = track->PerpLine(hit->GetPosition());
  TVector3 PerpPlane = track->PerpPlane(hit->GetPosition());
  PerpLine = PerpLine - PerpPlane;
  
  Double_t distLine = PerpLine.Mag();
  Double_t distPlane = PerpPlane.Mag();
  
  //  std::cout << distLine << "\t" << distPlane << std::endl;
  
  Double_t corrLine = (fPerpLineCut - distLine)/fPerpLineCut;
  Double_t corrPlane = (fPerpPlaneCut - distPlane)/fPerpPlaneCut;
  
  //  std::cout << corrLine << "\t" << corrPlane << std::endl;
  if (corrLine >= 0 && corrPlane >= 0) {

    survive = true;
    quality = TMath::Sqrt(quality*quality + corrLine*corrLine + corrPlane*corrPlane)/(corrLine + corrPlane +quality);
  }

  return true;
}

Bool_t TPCLinearTrackSearch::CorrelateHTRMS(TPCLinearTrack *track, TPCHit *hit, Bool_t &survive, Double_t &quality){

  Int_t numHitsinTrack = track->GetNumHits();

  if (numHitsinTrack < 4) {
    survive = true;
    return false;
  }

  survive = false;

  Double_t fRMSLineCut = 25;
  Double_t fRMSPlaneCut = 20;

  if (track->IsFitted() == false){
    fFitter->FitPlane(track);
  }

  Double_t RMSL, RMSP;
  fFitter->FitPlane(track,hit, RMSL, RMSP);
  /*
  RMSL = track->GetRMSLine();
  RMSP = track->GetRMSPlane();
  */
  
  //  std::cout << RMSL << "\t" << RMSP << std::endl;
  
  Double_t corrLine = (fRMSLineCut - RMSL)/fRMSLineCut;
  Double_t corrPlane = (fRMSPlaneCut - RMSP)/fRMSPlaneCut;

  //  std::cout << corrLine << "\t" << corrPlane << std::endl;
  
  if (corrLine >= 0 && corrPlane >= 0) {
    quality = TMath::Sqrt(quality*quality + corrLine*corrLine + corrPlane*corrPlane)/(corrLine + corrPlane + quality);
    survive = true;
  }

  return true;  
}

Bool_t TPCLinearTrackSearch::ExtrapTrack(std::vector<TPCHit*> &candHits, std::vector<TPCHit*> &unusedHits,TPCLinearTrack *track) {
  
  Int_t count = 0;
  Bool_t buildHead = true;
  Double_t extrapLength =5.;

  //track->SetIsExtrapolated();
  
    
  while (BuildByExtrap(candHits, unusedHits, track, buildHead, extrapLength)) {
    if (count++ > 200)
      break;
  }
  
  for (int i = 0 ; i < unusedHits.size() ; i++) 
    candHits.push_back(unusedHits.at(i));
  unusedHits.clear();

  count = 0;
  buildHead = !buildHead;
  extrapLength = 5.;

  while (BuildByExtrap(candHits, unusedHits, track, buildHead, extrapLength)) {
    if (count++ > 200)
      break;
  }
  
  for (int i = 0 ; i < unusedHits.size() ; i++) 
    candHits.push_back(unusedHits.at(i));
  unusedHits.clear();

  return true;
}


Bool_t TPCLinearTrackSearch::BuildByExtrap(std::vector<TPCHit*> &candHits, std::vector<TPCHit*> &unusedHits, TPCLinearTrack *track, Bool_t &buildHead, Double_t &extrapLength) {


  //Ensure there is no contents in unusedHits
  if (unusedHits.size() != 0) {
    for (auto hit : unusedHits)
      candHits.push_back(hit);
    unusedHits.clear();
  }


  //Get Hits from track
  std::vector<TPCHit*> *trackHits = track->GetHits();
  Int_t numTrackHits = track->GetNumHits();
  std::sort(trackHits->begin(),trackHits->end(),SortByZ()); 

  TVector3 head = trackHits->at(0)->GetPosition(); //Largest Z
  //std::cout<<"head : ("<<head.x()<<","<<head.y()<<","<<head.z()<<")"<<std::endl;

  

  TVector3 tail = trackHits->at(numTrackHits-1)->GetPosition(); //Smallest Z
  //std::cout<<"tail : ("<<tail.x()<<","<<tail.y()<<","<<tail.z()<<")"<<std::endl;

  TVector3 track_dir = head - tail;
  //std::cout<<"track direction : ("<<track_dir.x()<<","<<track_dir.y()<<","<<track_dir.z()<<")"<<std::endl;
  Double_t track_length = track_dir.Mag();
  //std::cout<<"track length : "<<track_length<<std::endl;

  if(track_length=0)return false;
  
  TVector3 track_ext = track_dir*(extrapLength/track_dir.Mag());
  //std::cout<<"ratio : "<<extrapLength/track_dir.Mag()<<std::endl;
  //std::cout<<"track extension : ("<<track_ext.x()<<","<<track_ext.y()<<","<<track_ext.z()<<")"<<std::endl;
  //std::cout<<"extention length : "<<track_ext.Mag()<<", Set value"<<extrapLength<<std::endl;
  TVector3 track_ext_rev = -1 * track_ext;

  TVector3 pos;

  if (buildHead)
    pos = head + track_ext;
  else
    pos = tail + track_ext_rev;

  
 
  
  //If the position of extrapolation is incorret, return false
  if (std::isnan(pos.X()) != 0 || std::isnan(pos.Z()) !=0){
    unusedHits = candHits;
    candHits.clear();
    return false;
  }
  
 
  //Double_t rms = track->GetRMSLine();
  //if (rms < 25) rms = 25.;
  
  Double_t rms = 30.;
  
  std::vector<TPCHit*> candidates;

  GetAdjacentHits(pos, rms, candHits, candidates, unusedHits);

  
  Int_t numHits = candidates.size();

  SetSorting(TPCRiemannSort::kSortZ);
  SortHits(candidates);
  
  Bool_t hitFound = false;
  
  for (Int_t iHit = 0 ; iHit < numHits; iHit++){
    
    TPCHit *hit = candidates[iHit];
    
    Double_t trackQuality  = 0;
    Bool_t trackSurvive = CorrelateHT(track,hit,trackQuality);

    //std::cout<<trackSurvive<<std::endl;
    if (trackSurvive) { // Update the values      
      track->AddHit(hit);
      hitFound = true;
      fFitter -> FitPlane(track);    
    }
    else
      unusedHits.push_back(hit);
    
    fFitter -> FitPlane(track);    
    
  }
  if (hitFound) {
    //std::cout<<"extrapolated hit found"<<std::endl;
    extrapLength = 5.;
  }
  else{
    extrapLength += 5.;
  }
  
  //extrapLength +=5.;
  //if (extrapLength > 1000) return false;
  if (pos.x()<-300||pos.x()>300) return false;
  if (pos.z()<-300||pos.z()>300) return false;
  /*
  if (hitFound) {
    //    std::cout << "Hit found" << std::endl;
    extrapLength = 5.;
    //if (helicity != track->Helicity());
    buildHead = !buildHead;
  }
  else {
    extrapLength += 5.;
    if (extrapLength > 200) return false;

  }
  */

  candidates.clear();
 
  for (int i = 0 ; i < unusedHits.size() ; i++) 
    candHits.push_back(unusedHits.at(i));
  unusedHits.clear();
  

  return true;
}




void TPCLinearTrackSearch::FinalizeTrack(TPCLinearTrack *track) {
  /*
  auto fHitClusterArray = track->GetClusters();
  std::sort(fHitClusterArray->begin(), fHitClusterArray->end(), SortByZ());
  
  auto firstCluster = fHitClusterArray->at(0);
  auto firstClusterPos = firstCluster->GetPosition();
  auto firstPerpLine = fFitter->PerpLine(track,firstClusterPos);
  track->SetVertex(0, firstClusterPos + firstPerpLine);
  
  auto lastCluster = fHitClusterArray->back();
  auto lastClusterPos = lastCluster->GetPosition();
  auto firstPerpLine = fFitter->PerpLine(track,lastClusterPos);
  track->SetVertex(1, lastClusterPos + lastPerpLine);
  */
}

void TPCLinearTrackSearch::GetAdjacentHits(TPCLinearTrack *track,std::vector<TPCHit*> &hitcont,
					   std::vector<TPCHit*> &candidates, std::vector<TPCHit*> &unusedHits){

  std::vector<TPCHit*> *trackHits = track->GetHits();

  for (int iHit = 0 ; iHit < hitcont.size() ; iHit++) {
    TPCHit *hit = hitcont.at(iHit);    
    auto hitpos = hit->GetPosition();
    
    Bool_t flag = false;
    
    for (auto trackhit : *trackHits) {
      auto trackhitpos = trackhit->GetPosition();
      
      auto resVec = trackhitpos-hitpos;
      if (resVec.Mag() < 50.)
	flag = true;
    }
    
    if (flag == true) 
      candidates.push_back(hit);    
    else 
      unusedHits.push_back(hit);    
  }
  
  hitcont.clear();
}

void TPCLinearTrackSearch::GetAdjacentHits(TVector3 pos, Double_t rms, std::vector<TPCHit*> &hitcont, std::vector<TPCHit*> &candidates, std::vector<TPCHit*> &unusedHits) {
  
  Int_t pad = tpc::findPadID(pos.Z(),pos.X());
  Int_t layer = tpc::getLayerID(pad);
  Int_t row = tpc::getRowID(pad);

  Int_t layerRange = (rms + 8)/12;
  Int_t rowRange = (rms + 6)/8;

  for (int iHit = 0 ; iHit < hitcont.size() ; iHit++) {

    Bool_t flag = false;
    
    Int_t hitlayer = hitcont.at(iHit)->GetLayer();
    Int_t hitrow = hitcont.at(iHit)->GetRow();
    
    Int_t numlayer, numrow;

    for (int ilayer  = -layerRange ; ilayer <= layerRange ; ilayer++) {
      Int_t maxrow = tpc::padParameter[layer][1];
      for (int irow = -rowRange ; irow <= rowRange ; irow++) {
	
	numlayer = layer + ilayer;
	numrow = row + irow;
	
	if (numlayer < 0 || numrow < 0 || numlayer > 31 || numrow >= maxrow)
	  continue;
		
	if (numlayer == hitlayer && numrow == hitrow) 
	  flag = true;
      }      
    }
    
    if (flag == true) 
      candidates.push_back(hitcont.at(iHit));
    else 
      unusedHits.push_back(hitcont.at(iHit));
   }
   //  std::cout << hitcont.size() << "\t" << candidates.size() << "\t" << unusedHits.size() << std::endl;

  hitcont.clear();
}

void TPCLinearTrackSearch::SortHits(std::vector<TPCHit*> &hits) {
  
  SortHitClass sortHit;
  sortHit.SetSorting(fSorting);
  sortHit.SetPosition(fPos);
  std::sort(hits.begin(), hits.end(), sortHit);
}
