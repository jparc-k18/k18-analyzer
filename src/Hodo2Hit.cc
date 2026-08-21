/*
  Hodo2Hit.cc

  2012/1/24
*/

#include "Hodo2Hit.hh"
#include "ConfMan.hh"
#include "HodoParamMan.hh"
#include "HodoPHCMan.hh"
#include "RawData.hh"

#include <iostream>
#include <iomanip>
#include <string>
#include <cmath>
#include <stdexcept>
#include <sstream>

Hodo2Hit::Hodo2Hit( HodoRawHit *rhit, int index)
  : raw_(rhit), Status_(false), index_(index)
{
}

Hodo2Hit::~Hodo2Hit()
{
}

double Hodo2Hit::DeltaE( void ) const
{
  return sqrt(fabs(a1_*a2_));
}

bool Hodo2Hit::calculate( void )
{
  static const std::string funcname = "[Hodo2Hit::calculate]";

  Status_=false;
  if( raw_->GetNumOfTdcHits()!=2 ) return Status_;

  int tdc1=raw_->GetTdc1(index_), tdc2=raw_->GetTdc2(index_);
  if( tdc1<0 || tdc2<0 ) return Status_;

  ConfMan *confMan=ConfMan::GetConfManager();
  if(!confMan){
    std::cerr << funcname << ": cannot get confManager" << std::endl;
 
    return Status_;
  }
  HodoParamMan *hodoMan = confMan->GetHodoParamManager();
  HodoPHCMan   *phcMan  = confMan->GetHodoPHCManager();
  if(!hodoMan){
    std::cerr << funcname << ": cannot get HodoParamManager" << std::endl; 
    return Status_;
  }

  int cid=raw_->DetectorId(), plid=raw_->PlaneId(),
    seg=raw_->SegmentId();
  int adc1=raw_->GetAdc1(), adc2=raw_->GetAdc2();

  if( !hodoMan->GetTime(cid,plid,seg,0,tdc1,t1_) ||
      !hodoMan->GetTime(cid,plid,seg,1,tdc2,t2_) ) return Status_;
  
  if( adc1>=0 ){
    if( !hodoMan->GetDe(cid,plid,seg,0,adc1,a1_) ) return Status_;
  }
  else
    a1_=0.;

  if( adc2>=0 ){
    if( !hodoMan->GetDe(cid,plid,seg,1,adc2,a2_) ) return Status_;
  }
  else
    a2_=0.;
  
  ct1_=t1_; ct2_=t2_;
  
  if( phcMan ){
    phcMan->doCorrection(cid,plid,seg,0,t1_,a1_,ct1_ );
    phcMan->doCorrection(cid,plid,seg,1,t2_,a2_,ct2_ );
  }

  return Status_=true;
} 
