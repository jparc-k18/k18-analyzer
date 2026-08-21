/*
  BH2Hit.cc

  2012/1/24
*/

#include "BH2Hit.hh"
#include "ConfMan.hh"
#include "HodoParamMan.hh"
#include "HodoPHCMan.hh"
#include "RawData.hh"

#include <iostream>
#include <iomanip>
#include <string>
#include <cstring>
#include <cmath>
#include <stdexcept>
#include <sstream>

bool BH2Hit::calculate( void )
{
  static const std::string funcname = "[BH2Hit::calculate]";
  
  Status_=false;
  if( !this->Hodo1Hit::calculate()) return Status_;

  ConfMan *confMan=ConfMan::GetConfManager();
  if(!confMan){
    std::cerr << funcname << ": cannot get confManager" << std::endl; 
    return Status_;
  }
  HodoParamMan *hodoMan = confMan->GetHodoParamManager();
  if(!hodoMan){
    std::cerr << funcname << ": cannot get HodoParamManager" << std::endl; 
    return Status_;
  }

  int cid=raw_->DetectorId(), plid=raw_->PlaneId(),
    seg=raw_->SegmentId();

  hodoMan->GetTime(cid,plid,seg,2,0,Tofs_);

  return Status_=true;

}
