/*
  DCPairHitCluster.hh

  2012/1/24
*/

#ifndef DCPairHitCluster_h
#define DCPairHitCluster_h 1

#ifdef MemoryLeak
#include "DebugCounter.hh"
#endif

class DCLTrackHit;

class DCPairHitCluster
{
public:
  DCPairHitCluster( DCLTrackHit *hitA, DCLTrackHit *hitB=0 );
  ~DCPairHitCluster();

private:
  DCLTrackHit *hitA_, *hitB_;
  int nhits_;

#ifdef MemoryLeak
  static debug::Counter sm_counter;
#endif

public:
  int NumberOfHits( void ) const { return nhits_; }
  DCLTrackHit *GetHit( int i ) const
  {
    if(i==0)      return hitA_;
    else if(i==1) return hitB_;
    else          return 0;
  }

};

#endif
