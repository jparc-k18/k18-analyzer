/*
  BH2Cluster.cc

  2012/1/24
*/

#include "BH2Cluster.hh"
#include "BH2Hit.hh"
#include "HodoAnalyzer.hh"

#include <cstring>
#include <string>
#include <iostream>
#include <stdexcept>
#include <sstream>


BH2Cluster::BH2Cluster( BH2Hit *hitA, BH2Hit *hitB, BH2Hit *hitC )
  : hitA_(hitA), hitB_(hitB), hitC_(hitC), csize_(0), gfastatus_(true)
{
  if(hitA) ++csize_;
  if(hitB) ++csize_;
  if(hitC) ++csize_;

  calculate();
}

BH2Hit * BH2Cluster::GetHit( int i ) const
{
  if(i==0) return hitA_;
  else if(i==1) return hitB_;
  else if(i==2) return hitC_;
  else return 0;
}

void BH2Cluster::calculate( void )
{
  double ms=0., mt=0., de=0., t0=0., dt=0;
  if( hitA_ ){
    ms += hitA_->SegmentId();
    mt += hitA_->CTime();
    de += hitA_->DeltaE();
    t0 += hitA_->CTime0();
    dt += hitA_->GetT();
  }
  if( hitB_ ){
    ms += hitB_->SegmentId();
    mt += hitB_->CTime();
    de += hitB_->DeltaE();
    t0 += hitB_->CTime0();
    dt += hitB_->GetT();
  }
  if( hitC_ ){
    ms += hitC_->SegmentId();
    mt += hitC_->CTime();
    de += hitC_->DeltaE();
    t0 += hitC_->CTime0();
    dt += hitC_->GetT();
  }
  ms /= double(csize_);
  mt /= double(csize_);
  t0 /= double(csize_); 
  dt /= double(csize_);
  MeanTime_=mt; dE_=de; MeanSeg_=ms; time0_=t0; Tdiff_=dt;
}

bool BH2Cluster::ReCalc( bool applyRecursively )
{
  if( applyRecursively ){
    if(hitA_) hitA_->ReCalc(applyRecursively);
    if(hitB_) hitB_->ReCalc(applyRecursively);
    if(hitC_) hitC_->ReCalc(applyRecursively);
  }
  calculate();
  return true;
}
