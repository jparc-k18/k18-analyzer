#include "TPCRiemannTrackSearch.hh"
#include "TPCHit.hh"
#include "TPCRiemannTrack.hh"
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

TPCRiemannTrackSearch::TPCRiemannTrackSearch(){
  fMinRMSH = -9999;
  fMinRMSW = -9999;
  
  fMaxRMSH = -9999;
  fMaxRMSW = -9999;

  fProxXCut = -9999;
  fProxYCut = -9999;
  
  Initialize();
}

TPCRiemannTrackSearch::~TPCRiemannTrackSearch(){
}

void TPCRiemannTrackSearch::Initialize() {

  fFitter = new TPCRiemannFitter();
    
}

TPCRiemannTrack *TPCRiemannTrackSearch::NewTrack(std::vector<TPCHit*> &candHits) {

  TPCRiemannTrack *track = new TPCRiemannTrack();
  SetPosition(tgtPos);
  SortHits(candHits);
  track->AddHit(candHits.front());
  candHits.erase(candHits.begin() + 0);
  
  return track;
}

Bool_t TPCRiemannTrackSearch::InitTrack(std::vector<TPCHit*> &candHits, std::vector<TPCHit*> &unusedHits,TPCRiemannTrack *track)
{

  track->SetIsInitialized();

  std::vector<TPCHit*> candidates;
  std::vector<TPCHit*> badHits;
  
  GetAdjacentHits(track, candHits, candidates, unusedHits);
  
  Int_t numHits = candidates.size();
  Int_t count = 0;
  Bool_t terminate = false;
  while (numHits != 0) {
    
    fSorting = TPCRiemannSort::kSortDistance;
    SetPosition(track->GetMean());
    SortHits(candidates);
    
    for (Int_t iHit = 0 ; iHit < numHits; iHit++){
      //      TPCHit *hit = candHits.at(iHit);
      TPCHit *hit = candidates.back();
      candidates.pop_back();
      
      Bool_t trackSurvive = PreCorrelate(track,hit);

      if (trackSurvive) { // Update the values
	
	track->AddHit(hit);      	
	//	std::cout << "Added" << std::endl;
	if (track->GetNumHits() > 6) {	  
	  
	  if (track->GetNumHits() > 12) {
	    for (auto cand : candidates)
	      candHits.push_back(cand);
	    candidates.clear();
	    terminate = true;
	    break;
	  }
	  fFitter -> Fit(track);	   
	}

	fFitter->FitPlane(track);	

      }
      else {
	badHits.push_back(hit);
      }
    }
    
    for (auto unusedhit : unusedHits)
      candHits.push_back(unusedhit);
    unusedHits.clear();

    if (terminate == false) {
      GetAdjacentHits(track, candHits, candidates, unusedHits);       
      
      if (count++ > 200) {
	for (auto cand : candidates)
	  candHits.push_back(cand);
	candidates.clear();
	break;
      }
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

Bool_t TPCRiemannTrackSearch::ExtendTrack(std::vector<TPCHit*> &candHits, std::vector<TPCHit*> &unusedHits,
					  TPCRiemannTrack *track) {
  
  track->SetIsExtended();  

  std::vector<TPCHit*> candidates;
  std::vector<TPCHit*> badHits;

  GetAdjacentHits(track, candHits, candidates, unusedHits);

  Int_t numHits = candidates.size();
  Int_t count = 0;
  Bool_t terminate = false;
  
  while (numHits != 0 ) {

    fSorting = TPCRiemannSort::kSortReverseCharge;
    //    fSorting = TPCRiemannSort::kSortCharge;
    SortHits(candidates);

    fFitter -> Fit(track);
    
    for (Int_t iHit = 0 ; iHit < numHits; iHit++){
      
      TPCHit *hit = candidates.back();
      candidates.pop_back();
      
      Bool_t trackSurvive = Correlate(track,hit);
      
      if (trackSurvive){
	track->AddHit(hit);      
	fFitter -> Fit(track);
	/*
	if (track->GetNumHits() > 60) {
	  for (auto cand : candidates)
	    candHits.push_back(cand);
	  candidates.clear();
	  terminate = true;
	  break;
	}
	*/
      }
      else 
	badHits.push_back(hit);

    }//Loop over hits

    for (auto unusedhit : unusedHits)
      candHits.push_back(unusedhit);
    unusedHits.clear();

    //    if (terminate == false){
      GetAdjacentHits(track, candHits, candidates, unusedHits);       
      //    numHits = candidates.size();
      
      if (count++ > 200) {
	for (auto cand : candidates)
	  candHits.push_back(cand);
	candidates.clear();
	break;
      }
      
      //    }
    
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

Bool_t TPCRiemannTrackSearch::ExtrapTrack(std::vector<TPCHit*> &candHits, std::vector<TPCHit*> &unusedHits,
					  TPCRiemannTrack *track) {
  
  Int_t count = 0;
  Bool_t buildHead = true;
  Double_t extrapLength = 1.5; 

  track->SetIsExtrapolated();
    
  while (BuildByExtrap(candHits, unusedHits, track, buildHead, extrapLength)) {
    if (count++ > 200)
      break;
  }
  
  for (int i = 0 ; i < unusedHits.size() ; i++) 
    candHits.push_back(unusedHits.at(i));
  unusedHits.clear();

  count = 0;
  buildHead = !buildHead;
  extrapLength = 1.5;

  while (BuildByExtrap(candHits, unusedHits, track, buildHead, extrapLength)) {
    if (count++ > 200)
      break;
  }
  
  for (int i = 0 ; i < unusedHits.size() ; i++) 
    candHits.push_back(unusedHits.at(i));
  unusedHits.clear();

  return true;
}

Bool_t TPCRiemannTrackSearch::BuildByExtrap(std::vector<TPCHit*> &candHits, std::vector<TPCHit*> &unusedHits, TPCRiemannTrack *track, Bool_t &buildHead, Double_t &extrapLength) {
    
  const Double_t scale = gUser.GetParameter("Scale");

  if (unusedHits.size() != 0) {
    for (auto hit : unusedHits)
      candHits.push_back(hit);
    unusedHits.clear();
  }
  auto helicity = track->Helicity();

  TVector3 pos;
  if (buildHead)
    pos = track->ExtrapHead(extrapLength);
  else
    pos = track->ExtrapTail(extrapLength);

  //  pos.Print();
  
  if (std::isnan(pos.X()) != 0 || std::isnan(pos.Z()) !=0){
    unusedHits = candHits;
    candHits.clear();
    return false;
  }

  Double_t rms = scale * track->GetRMSW();
  //  std::cout <<  rms << std::endl;
  if (rms < 10) rms = 10.;
  
  std::vector<TPCHit*> candidates;
  GetAdjacentHits(pos, rms, candHits, candidates, unusedHits);
  
  Int_t numHits = candidates.size();
  //  std::cout << "NumHits : " << numHits << std::endl;
  //  SetSorting(TPCRiemannSort::kSortCharge);
  //  SortHits(candidates);
  
  Bool_t hitFound = false;
  
  for (Int_t iHit = 0 ; iHit < numHits; iHit++){
    
    TPCHit *hit = candidates[iHit];
    
    Bool_t trackSurvive = Correlate(track,hit);
      
    if (trackSurvive) { // Update the values      
      track->AddHit(hit);
      //      std::cout << "Hit found" << std::endl;
      hitFound = true;
      fFitter -> Fit(track);    
    }
    else
      unusedHits.push_back(hit);
    
    fFitter -> Fit(track);    
    
  }
  if (hitFound) {
    extrapLength = 1.5; 
    if (helicity != track->Helicity());
    buildHead = !buildHead;
  }
  else {
    extrapLength += 1.5;
    //    if (extrapLength > 2.0 * track->TrackLength()) {
    if (extrapLength >  0.5 * track->TrackLength()) {
      return false;
    }
  }

  candidates.clear();
 
  for (int i = 0 ; i < unusedHits.size() ; i++) 
    candHits.push_back(unusedHits.at(i));
  unusedHits.clear();
  

  return true;
}

Bool_t TPCRiemannTrackSearch::ConfirmTrack(std::vector<TPCHit*> &candHits, std::vector<TPCHit*> &unusedHits,
					   TPCRiemannTrack *track) {

  auto tailToHead = false;
  if (track->PositionAtTail().Z() > track->PositionAtHead().Z())
    tailToHead = true;

  track->SetIsConfirmed();
  
  //  if (candHits.size() > 0)
  ConfirmTrackHits(candHits, unusedHits, track, tailToHead);
  
  for (int i = 0 ; i < unusedHits.size() ; i++) 
    candHits.push_back(unusedHits.at(i));
  unusedHits.clear();
  
  tailToHead = !tailToHead;
  
  //  if (candHits.size() > 0)
  ConfirmTrackHits(candHits, unusedHits, track, tailToHead);

  for (int i = 0 ; i < unusedHits.size() ; i++) 
    candHits.push_back(unusedHits.at(i));
  unusedHits.clear();
  
  return true;
}


Bool_t TPCRiemannTrackSearch::ConfirmTrackHits(std::vector<TPCHit*> &candHits, std::vector<TPCHit*> &unusedHits,
					       TPCRiemannTrack *track, Bool_t &tailToHead){
  
  Int_t numHits = track->GetNumHits();

  std::vector<TPCHit*> *trackHits = track->GetHits();
  
  if (tailToHead)
    std::sort(trackHits->begin(), trackHits->end(), SortByIncLength(track));
  else
    std::sort(trackHits->begin(), trackHits->end(), SortByDecLength(track));
  
  
  TVector3 q,m;
  auto pre = track->ExtrapByMap(trackHits->at(numHits-1)->GetPosition(),q,m);

  auto extrapLength = 1.5;
  
  for (int iHit = 1; iHit < numHits ; iHit++) {

    TPCHit *trackHit = trackHits->at(numHits-iHit-1);
    auto cur = track->ExtrapByMap(trackHit->GetPosition(),q,m);
    Bool_t trackSurvive = Correlate(track,trackHit);

    if (!trackSurvive) {
      track->RemoveHit(trackHit);
      unusedHits.push_back(trackHit);
      auto helicity = track->Helicity();
      fFitter->Fit(track);
      if (helicity != track->Helicity())
	tailToHead = !tailToHead;
    }

    auto dLength = abs(cur-pre);
    extrapLength = 1.5; //TODO
    while (dLength > 0 && BuildByInterp(candHits, unusedHits, track, tailToHead, extrapLength)) { 
      dLength -= 1.5;
    }  
  }
  
  //  extrapLength = 1.5;  
  //  while (BuildByExtrap(candHits, unusedHits, track,tailToHead, extrapLength)) {}    

  return true;

}

Bool_t TPCRiemannTrackSearch::BuildByInterp(std::vector<TPCHit*> &candHits, std::vector<TPCHit*> &unusedHits,
					    TPCRiemannTrack *track, Bool_t &buildHead, Double_t &extrapLength) {

  const Double_t scale = gUser.GetParameter("Scale");

  if (unusedHits.size() != 0) {
    for (auto hit : unusedHits)
      candHits.push_back(hit);
    unusedHits.clear();
  }
  
  auto helicity = track->Helicity();
  TVector3 pos;


  if (buildHead)
    pos = track->InterpByLength(extrapLength);
  else
    pos = track->InterpByLength(track->TrackLength() - extrapLength);
  
  if (std::isnan(pos.X()) != 0 || std::isnan(pos.Z()) !=0){
    candHits.clear();
    return false;
  }
  
  Double_t rms = scale * track->GetRMSW();
  if (rms < 10) rms = 10.;
  
  std::vector<TPCHit*> candidates;
  GetAdjacentHits(pos, rms, candHits, candidates, unusedHits);

  Int_t numHits = candidates.size();    
  SortHits(candidates);
  
  Bool_t hitFound = false;
  
  for (Int_t iHit = 0 ; iHit < numHits; iHit++){

    TPCHit *hit = candidates[iHit];
    
    Bool_t trackSurvive = Correlate(track,hit);
      
    if (trackSurvive) { // Update the values

      track->AddHit(hit);
      hitFound = true;
      fFitter -> Fit(track);    

    }
    else
      unusedHits.push_back(hit);
    
  }
  if (hitFound) {
    //    std::cout << "Hit found" << std::endl;
    extrapLength = 0.5;
    if (helicity != track->Helicity());
    buildHead = !buildHead;
  }
  else {
    extrapLength += 0.5;
    if (extrapLength > track->TrackLength()) {
      return false;
    }
  }

  candidates.clear();

  return true;
}

void TPCRiemannTrackSearch::GetAdjacentHits(TPCRiemannTrack *track,std::vector<TPCHit*> &hitcont,
					    std::vector<TPCHit*> &candidates, std::vector<TPCHit*> &unusedHits){

  const Double_t scale = gUser.GetParameter("Scale");

  auto trackhits = track->GetHits();

  if (track->IsExtended()) {
  
      for (int iHit = 0 ; iHit < hitcont.size() ; iHit++) {
	TPCHit *hit = hitcont.at(iHit);    
	auto hitpos = hit->GetPosition();

	Bool_t flag = false;

	for (auto trackhit : *trackhits) {
	  auto trackhitpos = trackhit->GetPosition();
	
	  auto resVec = trackhitpos-hitpos;

	  if (resVec.Mag() < scale*fDistCut){
	    //	    if (fabs(resVec.Y()) < 1.2*fProxYCut) {
	      flag = true;
	      //	    }
	  }
	}
    
	if (flag == true) 
	  candidates.push_back(hit);    
	else 
	  unusedHits.push_back(hit);    
      }
    }

    else {
      for (int iHit = 0 ; iHit < hitcont.size() ; iHit++) {
	TPCHit *hit = hitcont.at(iHit);    
	auto hitpos = hit->GetPosition();

	Bool_t flag = false;

	for (auto trackhit : *trackhits) {
	  auto trackhitpos = trackhit->GetPosition();
	
	  auto resVec = trackhitpos-hitpos;
	  if (resVec.Mag() < fDistCut)
	    flag = true;
	}
	
	if (flag == true) 
	  candidates.push_back(hit);    
	else 
	  unusedHits.push_back(hit);    
      }
    }

  //  std::cout << candidates.size() << std::endl;
  hitcont.clear();
}

void TPCRiemannTrackSearch::GetAdjacentHits(TVector3 pos, Double_t rms,
					    std::vector<TPCHit*> &hitcont, std::vector<TPCHit*> &candidates,
					    std::vector<TPCHit*> &unusedHits) {
  
  const Double_t scale = gUser.GetParameter("Scale");
  //  std::cout << hitcont.size() << "\t" << candidates.size() << "\t" << unusedHits.size() << std::endl;

  Int_t pad = tpc::findPadID(pos.Z(),pos.X());
  Int_t layer = tpc::getLayerID(pad);
  Int_t row = tpc::getRowID(pad);

  Int_t layerRange = 1;
  //  Int_t rowRange = rms/scale;

  Int_t rowRange = (rms)/2;
  //    Int_t rowRange = (rms + 4)/2;
  //  Int_t rowRange = (rms + 3)/6;
  //  Int_t rowRange = (rms + 2)/4;

  //  std::  cout << layerRange << "\t" << rowRange << std::endl;
  
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
		

	if (numlayer == hitlayer && numrow == hitrow) {
	  //	  std::cout << ilayer << "\t" << irow << "\t" << numlayer << "\t" << numrow << std::endl;
	  //	  std::cout << "Hit " << "\t" << hitlayer << "\t" << hitrow << std::endl;
	  flag = true;
	}
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

Bool_t TPCRiemannTrackSearch::TrackQualityCheck(TPCRiemannTrack *track){

  Double_t continuity = track->Continuity();
  std::cout << "Continuity : " << continuity << "\t" << track->TrackLength() << std::endl;
  if (continuity < 0.6) {
    if (track->TrackLength() * continuity < 500)
      return false;
  }

  if (track->GetRadius() < 30)
    return false;

  return true;
}


Bool_t TPCRiemannTrackSearch::PreCorrelate(TPCRiemannTrack *track, TPCHit* hit) {

  const Double_t scale = gUser.GetParameter("Scale");

    Double_t quality = 0;
  
  Int_t layer = hit->GetLayer();
  Int_t row = hit->GetRow();

  auto hitpos = hit->GetPosition();

  auto trackHits = track->GetHits();  

  Bool_t isNear = false;
  for (auto trackhit : *trackHits) {
    
    if (hit == trackhit){
      std::cout << "Same Hits" << std::endl;
      return false;
    }

    if (layer == trackhit->GetLayer() && row == trackhit->GetRow())
      return false;
    
    auto pos = trackhit->GetPosition();    
    auto dx = abs(hitpos.X()-pos.X());
    auto dy = abs(hitpos.Y()-pos.Y());

    if (dx < fProxXCut && dy < fProxYCut)
      isNear = true;   
  }
  
  if (isNear == false)
    return false;
  
  if (track->IsBad()) {
    quality = 1;
  }
  
  else if (track->IsLine()) {
    auto perp = track->PerpLine(hit->GetPosition());
    
    Double_t rmsCut = track->GetRMSH();
    if (rmsCut < fMinRMSH) rmsCut = fMinRMSH;
    if (rmsCut > fMaxRMSH) rmsCut = fMaxRMSH;
    rmsCut = scale *rmsCut;
    
    if (perp.Y() > rmsCut){
      quality = 0;
    }
    else {
      perp.SetY(0);
      // if (perp.Mag() < 15) quality = 1;
      if (perp.Mag() < fProxYCut) quality = 1;
    }
  }
  
  else if (track->IsPlane()) {
    Double_t dist = (track->PerpPlane(hit->GetPosition())).Mag();

    Double_t rmsCut = track->GetRMSH();

    if (rmsCut < fMinRMSH) rmsCut = fMinRMSH;
    if (rmsCut > fMaxRMSH) rmsCut = fMaxRMSH;
    
    rmsCut = scale*rmsCut;

    if (dist < rmsCut ) quality = 1;
    //    std::cout << dist << "\t" << rmsCut << std::endl;
  }
  
  if (quality > 0)
    return true;
  else
    return false;

}

Bool_t TPCRiemannTrackSearch::Correlate(TPCRiemannTrack* track, TPCHit* hit){

  Double_t scale = gUser.GetParameter("Scale");

  const Double_t MaxTrackLength = gUser.GetParameter("MaxTrackLength");
  Double_t trackLength = track->TrackLength();
  if (trackLength < MaxTrackLength)
    scale = scale + (MaxTrackLength - trackLength)/MaxTrackLength;
  
  Double_t rmsWCut = track->GetRMSW();
  if (rmsWCut < fMinRMSW) rmsWCut = fMinRMSW;
  if (rmsWCut > fMaxRMSW) rmsWCut = fMaxRMSW;
  rmsWCut = scale * rmsWCut;

  Double_t rmsHCut = track->GetRMSH();
  if (rmsHCut < fMinRMSH) rmsHCut = fMinRMSH;
  if (rmsHCut > fMaxRMSH) rmsHCut = fMaxRMSH;
  rmsHCut = scale * rmsHCut;

  auto qHead = track->Map(track->PositionAtHead());
  auto qTail = track->Map(track->PositionAtTail());
  auto q = track->Map(hit->GetPosition());  

  Double_t quality = 0;

  /*

    auto LengthAlphaCut = [track](Double_t dLength) {
    Double_t dip = track->GetDip();
    Double_t radius = track->GetRadius();
    Double_t val = dLength * TMath::Cos(dip)/radius;
    if (dLength > 0 ) {
      if (dLength > 0.5 * track->TrackLength()){
	if (abs(val) > TMath::Pi()/2.)
	  return true;
      }
    }    
    return false;
  };

  if (qHead.Z() > qTail.Z()) {
    if (LengthAlphaCut(q.Z() - qHead.Z())) {
      quality = 0.;
      return false;
    }
    if (LengthAlphaCut(qTail.Z() - q.Z())) {
      quality = 0.;
      return false;
    }
  }
  else {
    if (LengthAlphaCut(q.Z() - qTail.Z())) {
      quality = 0.;
      return false;
    }
    if (LengthAlphaCut(qHead.Z() - q.Z())) {
      quality = 0.;
      return false;
    }
  }
  */
  Double_t dr = std::abs(q.X());
  if (dr < rmsWCut && std::abs(q.Y()) < rmsHCut)
    quality = sqrt((dr-rmsWCut)*(dr-rmsWCut))/rmsWCut;

  if (quality > 0)
    return true;
  
  else
    return false;
  
}

void TPCRiemannTrackSearch::SortHits(std::vector<TPCHit*> &hits) {

  SortHitClass sortHit;
  sortHit.SetSorting(fSorting);
  sortHit.SetPosition(fPos);
  sortHit.SetTrack(fTrk);
  std::sort(hits.begin(), hits.end(), sortHit);
}

Bool_t SortHitClass::operator() (TPCHit *hit1, TPCHit *hit2){
  Double_t a1, a2;
  TVector3 d1;
  TVector3 d2;

  switch (fSorting) {

  case 0:
    a1 = hit1->GetPosition().X();
    a2 = hit2->GetPosition().X();
    return a1 > a2;
    break;

  case 1:
    a1 = hit1->GetPosition().Y();
    a2 = hit2->GetPosition().Y();
    return (a1) > (a2);
    break;

  case -1:
    a1 = hit1->GetPosition().Y();
    a2 = hit2->GetPosition().Y();
    return (a1) < (a2);
    break;

  case 2:
    a1 = hit1->GetPosition().Z();
    a2 = hit2->GetPosition().Z();
    return a1 > a2;
    break;

  case -2:
    a1 = hit1->GetPosition().Z();
    a2 = hit2->GetPosition().Z();
    return a1 < a2;
    break;
   
  case 4:
    d1 = hit1->GetPosition();
    d1 -= fPosition;
    
    d2 = hit2->GetPosition();
    d2 -= fPosition;
    
    a1 = d1.Mag(); a2 = d2.Mag();

    if (a1 == a2 ) return d1.Z() > d2.Z();
    return a1 > a2;
    break;

  case -4:
    d1 = hit1->GetPosition();
    d1 -= fPosition;
    
    d2 = hit2->GetPosition();
    d2 -= fPosition;
    
    a1 = d1.Mag(); a2 = d2.Mag();

    if (a1 == a2 ) return d1.Z() < d2.Z();
    return a1 < a2;
    break;

  case 5:

    a1 = hit1->GetAlpha();
    a2 = hit2->GetAlpha(); 
    
    if (a1 < -1*TMath::PiOver2()) a1 += TMath::TwoPi();
    if (a2 < -1*TMath::PiOver2()) a2 += TMath::TwoPi();
    return a1 > a2;
    break;

  case -5:
	      
    a1 = hit1->GetAlpha();
    a2 = hit2->GetAlpha(); 
    
    if (a1 < -1*TMath::PiOver2()) a1 += TMath::TwoPi();
    if (a2 < -1*TMath::PiOver2()) a2 += TMath::TwoPi();
    return a1 < a2;
    break;

  case 6:

    d1 = fTrk->Map(hit1->GetPosition());
    a1 = d1.Z();

    d2 = fTrk->Map(hit2->GetPosition());
    a2 = d2.Z();

    return a1 > a2;
    break;
    
  case -6:

    d1 = fTrk->Map(hit1->GetPosition());
    a1 = d1.Z();

    d2 = fTrk->Map(hit2->GetPosition());
    a2 = d2.Z();

    return a1 < a2;
    break;

  case -3:

    a1 = hit1->GetR();
    a2 = hit2->GetR(); 
    return a1 < a2;
    break;

  case 7:

    d1 = hit1->GetPosition();
    //    a1 = hit1->GetDe();
    a1 = hit1->GetCDe();

    d2 = hit2->GetPosition();
    //    a2 = hit2->GetDe();
    a2 = hit2->GetCDe();

    if ( a1 == a2)
      return d1.Y() > d2.Y();

    return a1 > a2;

    break;
    
  case -7:

    d1 = hit1->GetPosition();
    a1 = hit1->GetCDe();

    d2 = hit2->GetPosition();
    a2 = hit2->GetCDe();

    if ( a1 == a2)
      return d1.Y() < d2.Y();

    return a1 < a2;

    break;

  case 8:

    a1 = hit1->GetTheta();
    a2 = hit2->GetTheta();

    return a1 > a2;

    break;

  case -8:
    a1 = hit1->GetTheta();
    a2 = hit2->GetTheta();

    return a1 < a2;

    break;
    
  case 3:
  default:

    a1 = hit1->GetR();
    a2 = hit2->GetR(); 
    return a1 > a2;

  }
}



