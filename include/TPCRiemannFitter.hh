#ifndef TPC_RIEMANN_FITTER
#define TPC_RIEMANN_FITTER

#include <vector>
#include "Rtypes.h"

#include "ODRFitter.hh"
#include "TPCHit.hh"
#include "TPCRiemannTrack.hh"
#include "TPCLinearTrack.hh"

class TPCRiemannFitter
{
public:
  TPCRiemannFitter() : fODRFitter(new ODRFitter()) {}
  ~TPCRiemannFitter() {};

  Bool_t FitPlane(TPCRiemannTrack *track);
  Bool_t FitPlane(TPCLinearTrack *track);
  Bool_t FitPlane(TPCLinearTrack *track, TPCHit *hit,
		  Double_t &rmsL, Double_t &rmsP);

  Bool_t FitCluster(TPCRiemannTrack *track);
  Bool_t FitCluster(TPCLinearTrack *track);

  Bool_t FitLinear(TPCLinearTrack *track);

  Bool_t Fit(TPCRiemannTrack *track);

private:
  ODRFitter *fODRFitter;

};

#endif
