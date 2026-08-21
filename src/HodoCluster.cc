/*
  HodoCluster.cc

  2012/1/24
*/

#include "HodoCluster.hh"
#include "Hodo2Hit.hh"
#include "HodoAnalyzer.hh"

#include <cstring>
#include <string>
#include <iostream>
#include <stdexcept>
#include <sstream>

HodoCluster::HodoCluster( Hodo2Hit *hitA, Hodo2Hit *hitB, Hodo2Hit *hitC )
  : hitA_(hitA), hitB_(hitB), hitC_(hitC), csize_(0), gfastatus_(true)
{
  if(hitA) ++csize_;
  if(hitB) ++csize_;
  if(hitC) ++csize_;

  calculate();
}

Hodo2Hit * HodoCluster::GetHit( int i ) const
{
  if(i==0) return hitA_;
  else if(i==1) return hitB_;
  else if(i==2) return hitC_;
  else return 0;
}

void HodoCluster::calculate( void )
{
  double ms=0., mt=0., de=0., dt=0.;
  if( hitA_ ){
    ms += hitA_->SegmentId();
    mt += hitA_->CMeanTime();
    de += hitA_->DeltaE();
    dt += (hitA_->GetTDown()-hitA_->GetTUp());
  }
  if( hitB_ ){
    ms += hitB_->SegmentId();
    mt += hitB_->CMeanTime();
    de += hitB_->DeltaE();
    dt += (hitB_->GetTDown()-hitB_->GetTUp());
  }
  if( hitC_ ){
    ms += hitC_->SegmentId();
    mt += hitC_->CMeanTime();
    de += hitC_->DeltaE();
    dt += (hitC_->GetTDown()-hitC_->GetTUp());
  }
  ms /= double(csize_);
  mt /= double(csize_);
  dt /= double(csize_);

  MeanTime_=mt; dE_=de; MeanSeg_=ms, Tdiff_=dt;

} 

bool HodoCluster::ReCalc( bool applyRecursively )
{
  if( applyRecursively ){
    if(hitA_) hitA_->ReCalc(applyRecursively);
    if(hitB_) hitB_->ReCalc(applyRecursively);
    if(hitC_) hitC_->ReCalc(applyRecursively);
  }
  calculate();

  return true;
}
