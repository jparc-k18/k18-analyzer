/*
  DCPairHitCluster.cc

  2012/1/24
*/

#include "DCPairHitCluster.hh"

#ifdef MemoryLeak
debug::Counter DCPairHitCluster::sm_counter("DCPariHitCluster");
#endif

DCPairHitCluster::DCPairHitCluster( DCLTrackHit *hitA, DCLTrackHit *hitB )
  : hitA_(hitA), hitB_(hitB), nhits_(0)
{
  if(hitA_) ++nhits_;
  if(hitB_) ++nhits_;

#ifdef MemoryLeak
  ++sm_counter;
#endif
}

DCPairHitCluster::~DCPairHitCluster()
{
#ifdef MemoryLeak
  --sm_counter;
#endif
}
