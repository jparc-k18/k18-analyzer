#ifndef TPC_RIEMANN_TRACK_SEARCH_HH
#define TPC_RIEMANN_TRACK_SEARCH_HH

#include "TPCRiemannFitter.hh"

#include <vector>

#include <TString.h>
#include <TMath.h>

class TPCHit;
class TPCRiemannTrack;

typedef std::vector<TPCHit*> TPCHitContainer;

class TPCRiemannTrackSearch {

public:

  TPCRiemannTrackSearch();
  ~TPCRiemannTrackSearch();

  void SetMinHitsForFit(Int_t numHits) {fMinHitsForFit = numHits;}
  
  void SetSorting(Int_t sorting) {fSorting = sorting;} //true : use internal sorting whein adding hits to track candidates
  void SetPosition(TVector3 pos) {fPos = pos;} 
  void SetTrack(TPCRiemannTrack *track) {fTrk = track;} 

  Int_t GetSorting() {return fSorting;}

  Int_t BuildTracks(std::vector<TPCHit*> &hits, std::vector<TPCRiemannTrack *> &TrackletCand);

  TPCRiemannTrack *NewTrack(std::vector<TPCHit*> &hits);
  
  Bool_t InitTrack(std::vector<TPCHit*> &hits, std::vector<TPCHit*> &unusedhits, TPCRiemannTrack *track);
  Bool_t ExtendTrack(std::vector<TPCHit*> &hits, std::vector<TPCHit*> &unusedhits, TPCRiemannTrack *track);
  Bool_t ExtrapTrack(std::vector<TPCHit*> &hits, std::vector<TPCHit*> &unusedhits, TPCRiemannTrack *track);
  Bool_t ConfirmTrack(std::vector<TPCHit*> &candHits, std::vector<TPCHit*> &unusedhits, TPCRiemannTrack *track);

  Bool_t ConfirmTrackHits(std::vector<TPCHit*> &candHits, std::vector<TPCHit*> &unusedhits, TPCRiemannTrack *track, Bool_t &tailToHead);
  
  void SortHits(std::vector<TPCHit *> &hits);
  void SortClusters(std::vector<TPCHitCluster *> &clusters);

  void Initialize();

  Bool_t BuildByExtrap(std::vector<TPCHit*> &hits, std::vector<TPCHit*> &unusedhits,
		       TPCRiemannTrack *track, Bool_t &buildHead, Double_t &extrapLength);
  Bool_t BuildByInterp(std::vector<TPCHit*> &hits, std::vector<TPCHit*> &unusedhits,
		       TPCRiemannTrack *track, Bool_t &buildHead, Double_t &extrapLength);
  
  
  void GetAdjacentHits(TPCRiemannTrack *track,
		       std::vector<TPCHit*> &hitcont, std::vector<TPCHit*> &candidates, std::vector<TPCHit*> &unusedHits);
  void GetAdjacentHits(TVector3 pos, Double_t rms,
		       std::vector<TPCHit*> &hitcont, std::vector<TPCHit*> &candidates, std::vector<TPCHit*> &unusedHits);
  Bool_t TrackQualityCheck(TPCRiemannTrack *track);

  Bool_t PreCorrelate(TPCRiemannTrack *track, TPCHit* hit);
  Bool_t Correlate(TPCRiemannTrack *track, TPCHit* hit);

  void SetRMSCut(Double_t MinRMSW, Double_t MinRMSH, Double_t MaxRMSW, Double_t MaxRMSH){
    fMinRMSW = MinRMSW;
    fMinRMSH = MinRMSH;
    
    fMaxRMSW = MaxRMSW;
    fMaxRMSH = MaxRMSH;    
  }

  void SetProxCut(Double_t dx, Double_t dy){
    fProxXCut = dx;
    fProxYCut = dy;
  }
  
  void SetDistCut(Double_t dist){
    fDistCut = dist;
  }
  
private:

  TPCRiemannFitter *fFitter;

  Int_t fSorting;
  Int_t fMinHitsForFit;
  TVector3 fPos;

  TPCRiemannTrack *fTrk;

  Double_t fMinRMSH;
  Double_t fMinRMSW;
  
  Double_t fMaxRMSH;
  Double_t fMaxRMSW;

  Double_t fProxXCut;
  Double_t fProxYCut;

  Double_t fDistCut;
};

class SortHitClass {
public:
  Bool_t operator() (TPCHit *hit1, TPCHit *hit2);

  void SetSorting(Int_t sorting) {fSorting = sorting;}
  void SetPosition(TVector3 p) {fPosition = p;}
  void SetTrack(TPCRiemannTrack *track) {fTrk = track;}

private:
  Int_t fSorting;
  TVector3 fPosition;
  TPCRiemannTrack *fTrk;
};


#endif
		     
  

