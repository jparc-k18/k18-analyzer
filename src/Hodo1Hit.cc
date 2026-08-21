/*
  Hodo1Hit.cc

  2012/1/24
*/

#include "Hodo1Hit.hh"

#include "ConfMan.hh"
#include "HodoParamMan.hh"
#include "HodoPHCMan.hh"
#include "RawData.hh"

#include <iostream>
#include <iomanip>
#include <string>
#include <cstring>
#include <stdexcept>
#include <sstream>


Hodo1Hit::Hodo1Hit( HodoRawHit *rhit )
  : raw_(rhit), Status_(false)
{}

Hodo1Hit::~Hodo1Hit()
{}

bool Hodo1Hit::calculate( void )
{
  static const std::string funcname = "[Hodo1Hit::calculate]";

  Status_=false;
  if( raw_->GetNumOfTdcHits()!=1 ) return Status_;

  // detector information
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

  int cid=raw_->DetectorId(), plid=raw_->PlaneId(), seg=raw_->SegmentId();

  // hit information
  multi_hit = raw_->SizeTdc1();
  int UorD=0;
  int adc = raw_->GetAdc1(0);
  if( 0 > raw_->GetTdc1(0) ){
    UorD=1;
    adc = raw_->GetAdc2(0);
  }

  double dE = 0;
  if( adc>=0 ){
    //std::cout<<adc<<std::endl;
    //std::cout<<cid<<"  "<<plid<<"  "<<seg<<"  "<<UorD<<"  "<<adc<<"  "<<dE<<"  "<<std::endl;
    if( !hodoMan->GetDe(cid, plid, seg, UorD, adc, dE) ){
      std::cerr << funcname << " : something wrong at GetDe" << std::endl;
      return Status_;      
    }
  }
  else{
    dE=0.;
  }
  a_.push_back(dE);

  for(int m = 0; m<multi_hit; ++m){
    int    tdc  = 0;
    double time = 0;
    if(0 == UorD){
      tdc = raw_->GetTdc1(m);
    }else{
      tdc = raw_->GetTdc2(m);
    }
    
    if(tdc<0){continue;}

    if( !hodoMan->GetTime(cid, plid, seg, UorD, tdc, time) ){
      std::cerr << funcname << " : something wrong at GetTime" << std::endl;
      return Status_;
    }
    t_.push_back(time);
    
    double ctime = time;
    if( phcMan ){
      phcMan->doCorrection(cid, plid, seg, UorD, time, dE, ctime );
    }
    ct_.push_back(ctime);
    Status_ = true;
  }

  return Status_;
}
