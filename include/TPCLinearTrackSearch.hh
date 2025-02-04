#ifndef TPC_LINEAR_TRACK_SEARCH_HH
#define TPC_LINEAR_TRACK_SEARCH_HH

#include "TPCHit.hh"
#include "TPCRiemannFitter.hh"
#include "TPCLinearTrack.hh"

#include <vector>

#include <TVector3.h>
#include <TMath.h>

class TPCLinearTrackSearch
{

public:


  TPCLinearTrackSearch() ;
  ~TPCLinearTrackSearch() {};

  Int_t BuildTracks(std::vector<TPCHit*> &hits, std::vector<TPCLinearTrack *> &TrackletCand);

  TPCLinearTrack *NewTrack(std::vector<TPCHit*> &hits);

  Bool_t BuildTrack(std::vector<TPCHit*> &hits, std::vector<TPCHit*> &unusedhits, TPCLinearTrack *track);

  Bool_t ExtrapTrack(std::vector<TPCHit*> &hits, std::vector<TPCHit*> &unusedhits, TPCLinearTrack *track);

  void SetSorting(Int_t sorting) {fSorting = sorting;} //true : use internal sorting whein adding hits to track candidates
  void SetPosition(TVector3 pos) {fPos = pos;} 

  void SortHits(std::vector<TPCHit *> &hits);
  void SortClusters(std::vector<TPCHitCluster *> &clusters);

  void Initialize();



  void FinalizeTrack(TPCLinearTrack *track);

  Bool_t BuildByExtrap(std::vector<TPCHit*> &hits, std::vector<TPCHit*> &unusedhits,TPCLinearTrack *track, Bool_t &buildHead, Double_t &extrapLength);

  Bool_t CorrelateHT(TPCLinearTrack *track, TPCHit *hit, Double_t &quality);

  Bool_t CorrelateHTPerp(TPCLinearTrack *track, TPCHit *hit, Bool_t &survive, Double_t &quality);
  Bool_t CorrelateHTProx(TPCLinearTrack *track, TPCHit *hit, Bool_t &survive, Double_t &quality);
  Bool_t CorrelateHTRMS(TPCLinearTrack *track, TPCHit *hit, Bool_t &survive, Double_t &quality);

  void GetAdjacentHits(TPCLinearTrack *track,std::vector<TPCHit*> &hitcont,
		       std::vector<TPCHit*> &candidates, std::vector<TPCHit*> &unusedHits);
void GetAdjacentHits(TVector3 pos, Double_t rms,
		       std::vector<TPCHit*> &hitcont, std::vector<TPCHit*> &candidates, std::vector<TPCHit*> &unusedHits);
private:

  TPCRiemannFitter *fFitter;

  Int_t fSorting;
  Bool_t fSortingMode;
  TVector3 fPos;
};

#endif
