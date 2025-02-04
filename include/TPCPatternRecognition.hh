// -*- C++ -*-

#ifndef TPC_PATTERN_RECOGNITION_HH
#define TPC_PATTERN_RECOGNITION_HH

#include "TPCHit.hh"
#include "TPCHitCluster.hh"

#include "TPCRiemannSort.hh"
#include "DCAnalyzer.hh"
#include "TPCRiemannTrackSearch.hh"
#include "TPCRiemannTrack.hh"

#include "TPCLinearTrackSearch.hh"
#include "TPCLinearTrack.hh"

#include "TPCRiemannFitter.hh"
//#include "TPCTrack.hh"


#include <TClonesArray.h>
#include <TString.h>

class TPCHitCluster;
class TPCRiemannTrack;
class TPCLinearTrack;
class TPCHit;

typedef std::vector<TPCHit*> TPCHitContainer;
typedef std::vector<TPCRiemannTrack*> TPCRiemannTrackContainer;

class TPCPatternRecognition
{
public:
  static const TString &ClassName();
  TPCPatternRecognition();
  virtual ~TPCPatternRecognition();

  Bool_t Initialize();
  Bool_t Read(const std::vector<TPCHitContainer> &HitCont);
  Bool_t ReadGeant4(const std::vector<TPCHitContainer> &HitCont);
  Bool_t ReadLinear(const std::vector<TPCHitContainer> &HitCont);

  virtual void Execute();
  virtual void ExecuteLinear();

  void BuildTracks(TPCRiemannTrackSearch *TrackSearch,
		   std::vector<TPCRiemannTrack* > *TrackletCont, Int_t sorting);
  void BuildBeamTracks(TPCLinearTrackSearch *TrackSearch,
		   std::vector<TPCLinearTrack* > *TrackletCont, Int_t sorting);
  void BuildBeamTracks(TPCRiemannTrackSearch *TrackSearch,
		   std::vector<TPCRiemannTrack* > *TrackletCont, Int_t sorting);
  void BuildKuramaTracks(TPCRiemannTrackSearch *TrackSearch,
		   std::vector<TPCRiemannTrack* > *TrackletCont, Int_t sorting);
  void BuildLinearTracks(TPCLinearTrackSearch *TrackSearch,
		   std::vector<TPCLinearTrack* > *TrackletCont, Int_t sorting);
  
  TVector3 FindVertex(std::vector<TPCRiemannTrack* > *trackCont, Int_t nIter = 1);

  TPCHitCluster *NewCluster(TPCHit *hit);
  Bool_t HitClustering(TPCRiemannTrack *track);
  Bool_t HitClustering(TPCLinearTrack *track);

  Bool_t CheckIsContinuousHits(TPCHitCluster *cluster, TPCHit *hit);
  //  Bool_t CheckIsContinuousHits(TPCHitCluster *cluster);
  Bool_t CheckClusterType(TPCRiemannTrack *track, TPCHit *currentHit, TPCHit *prevHit);
  
  TPCRiemannTrackContainer GetFinalTracks() const {return *fRiemannTrackArray;}
  TPCLinearTrackContainer GetFinalLinearTracks() const {return *fLinearTrackArray;}

  void SetEllipsoidCut(TVector3 vector, TVector3 radii, Double_t margin);
  void SetClusterCut(Double_t left, Double_t right, Double_t top, Double_t bottom);

  inline Double_t TransverseDistance(Double_t x_center, Double_t z_center, Double_t x, Double_t z);
  
private:
  
  std::vector <TPCRiemannTrack*> *fRiemannTrackArray;
  std::vector <TPCLinearTrack*> *fLinearTrackArray;
  
  std::vector <TPCHitCluster*> *fHitClusterArray;
  std::vector <TPCHit*> *fHitArray;  

  std::vector<TPCHit*> *fUsedHits;
  std::vector<TPCHit*> *fUnusedHits;
  std::vector<TPCHit*> *fCandHits;
  std::vector<TPCHit*> *fBeamCandHits;
  std::vector<TPCHit*> *fKuramaCandHits;
  std::vector<TPCHit*> *fLinearCandHits;
  
  TPCRiemannTrackSearch *fTrackSearch;
  TPCLinearTrackSearch *fLinearTrackSearch;
  TPCRiemannFitter *fFitter;
  
  Bool_t fSortingMode;
  Int_t fSorting;

  Double_t fCCLeft;
  Double_t fCCRight;
  Double_t fCCTop;
  Double_t fCCBottom;

  TVector3 fCutCenter;
  Double_t fCRadius;
  Double_t fZLength;
  Double_t fSRadius;
  TVector3 fERadii;
  Double_t fCutMargin;

  Double_t fClusteringAngle;
  Double_t fClusteringMargin;

};

inline const TString& TPCPatternRecognition::ClassName()
{
  static TString s_name("TPCPatternRecognition");
  return s_name;    
}

#endif

