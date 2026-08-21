/*
  DCHit.cc

  2012/1/24
*/

#include <cmath>
#include <iostream>
#include <string>
#include <stdexcept>
#include <sstream>
#include <iterator>
#include <algorithm>

#include "DCHit.hh"
#include "ConfMan.hh"
#include "DCGeomMan.hh"
#include "DCTdcCalibMan.hh"
#include "DCDriftParamMan.hh"
#include "DCLTrackHit.hh"
#include "DCParameters.hh"

const double Deg2Rad = acos(-1)/180.;
const double Rad2Deg = 180./acos(-1);

#ifdef MemoryLeak
debug::Counter DCHit::sm_counter("DCHit");
#endif

DCHit::DCHit()
  : layer_(-1), wire_(-1)
{
#ifdef MemoryLeak
  ++sm_counter;
#endif
}

DCHit::DCHit( int layer )
  : layer_( layer ), wire_(-1), mwpcflag_(false)
{
#ifdef MemoryLeak
  ++sm_counter;
#endif
}

DCHit::DCHit( int layer, double wire )
  : layer_(layer), wire_(wire), mwpcflag_(false)
{
#ifdef MemoryLeak
  ++sm_counter;
#endif
}

DCHit::~DCHit()
{
  clearRegisteredHits();
#ifdef MemoryLeak
  --sm_counter;
#endif
}

void DCHit::SetTdcVal( int tdc )
{
  tdc_.push_back(tdc); 
  belongTrack_.push_back(false);
  dlRange_.push_back(false);
}

void DCHit::SetAdcVal(int adc )
{
  adc_.push_back(adc);
}

void DCHit::clearRegisteredHits( void )
{
  int n=Cont_.size();
  for(int i=0; i<n; ++i)
    delete Cont_[i];
}

bool DCHit::CalcDCObservables( void )
{
  static const std::string funcname="[DCHit::CalcObservables]";
  
  ConfMan *confMan=ConfMan::GetConfManager();
  if(!confMan) return false;
  DCGeomMan *geomMan=confMan->GetDCGeomManager();
  if(!geomMan) return false;
  DCTdcCalibMan *calibMan=confMan->GetDCTdcCalibManager();
  if(!calibMan) return false;
  DCDriftParamMan *driftMan=confMan->GetDCDriftParamManager();
  if(!driftMan) return false;

  wpos_=geomMan->calcWirePosition(layer_,wire_);
  angle_=geomMan->GetTiltAngle(layer_);
  
  bool Status = true;
  int nhitdc  = tdc_.size();
  for (int i=0; i<nhitdc; i++) {
    double ctime;
    if(!calibMan->GetTime( layer_, wire_, tdc_[i], ctime ))
      return false;
    
    double dtime, dlength;
    bool status=driftMan->calcDrift( layer_, wire_, ctime, dtime, dlength );
    
    if (status == false) Status = status;
    
    dt_.push_back(dtime);
    dl_.push_back(dlength);
    
    if(layer_>=100){
      if( dl_[i]>MinDLBc[layer_-100] && dl_[i]<MaxDLBc[layer_-100] )
	dlRange_[i]=true;
    }
    else{
      if( dl_[i]>MinDLSdc[layer_] && dl_[i]<MaxDLSdc[layer_] )
	dlRange_[i]=true;
    }
  }

  return Status;
}

bool DCHit::CalcMWPCObservables()
{
  static const std::string funcname="[DCHit::CalcObservables]";
  
  ConfMan *confMan=ConfMan::GetConfManager();
  if(!confMan) return false;
  DCGeomMan *geomMan=confMan->GetDCGeomManager();
  if(!geomMan) return false;
  DCTdcCalibMan *calibMan=confMan->GetDCTdcCalibManager();
  if(!calibMan) return false;

  angle_=geomMan->GetTiltAngle(layer_);
  wpos_=geomMan->calcWirePosition(layer_,wire_);

  bool Status = true;
  int nhitdc  = tdc_.size();

  for (int i=0; i<nhitdc; i++) {
    if (dt_.size()<nhitdc){
      double ctime;
      if(!calibMan->GetTime( layer_, wire_, tdc_[i], ctime ))
	return false;
  
      dt_.push_back( ctime );
      dl_.push_back( 0.);

      double ctt; // tdc value of trailing edge --> time
      if(!calibMan->GetTime( layer_, wire_, trailing_[i], ctt ))
	return false;
      trailingTime_.push_back(ctt);
    }
    
    if( dt_[i]>MinDLBc[layer_-100] && dt_[i]<MaxDLBc[layer_-100] ){
      dlRange_[i]=true;
      Status=true; 
    }
    else Status=false; 
  }

  return Status;
}

bool DCHit::CalcFiberObservables()
{
  static const std::string funcname="[DCHit::CalcFiberObservables]";
  
  ConfMan *confMan=ConfMan::GetConfManager();
  if(!confMan) return false;
  DCGeomMan *geomMan=confMan->GetDCGeomManager();
  if(!geomMan) return false;

  angle_=geomMan->GetTiltAngle(layer_);
  
  bool Status = true;
  int nhitdc  = tdc_.size();
  
  for (int i=0; i<nhitdc; i++) {
    if (dt_.size()<nhitdc){
      dt_.push_back( tdc_[i] );
      dl_.push_back( 0.);
    }
    
    dlRange_[i]=true;
  }
  
  return Status;
}

bool DCHit::CalcSsdObservables( void )
{
  static const std::string funcname = "[DCHit::CalcSsdObservables]";
  ConfMan *confMan = ConfMan::GetConfManager();
  if(!confMan) return false;
  DCGeomMan *geomMan = confMan->GetDCGeomManager();
  if(!geomMan) return false;
  wpos_ = geomMan->calcWirePosition(layer_,wire_);
  bool Status = false;
  std::vector<int> adc_peek;
  if(adc_.size()>0)
    adc_peek = adc_;
  else
    return false;
  sort(adc_peek.begin(),adc_peek.end()-1,std::greater<int>());
  if(adc_peek[0]==adc_[0]) return false;
  if(adc_.size()>5) Status = true;
  return Status;
}
// bool DCHit::CalcObservablesSimulation( double dlength)
// {
//   static const std::string funcname="[DCHit::CalcObservablesSimulation]";

//   ConfMan *confMan=ConfMan::GetConfManager();
//   if(!confMan) return false;
//   DCGeomMan *geomMan=confMan->GetDCGeomManager();
//   if(!geomMan) return false;

//   wpos_=geomMan->calcWirePosition(layer_,wire_);
//   angle_=geomMan->GetTiltAngle(layer_);

//   dl_ = dlength;
//   bool status=true;

//   if(layer_>=100){
//     if( dl_>MinDLBc[layer_-100] && dl_<MaxDLBc[layer_-100] )
//       dlRange_=true;
//   }
//   else {
//     if( dl_>MinDLSdc[layer_] && dl_<MaxDLSdc[layer_] )
//       dlRange_=true;
//   }
//   return status;
// }

//______________________________________________________________________________
void
DCHit::print(const std::string& arg) const
{
  std::cout << "DCHit::print " << arg << std::endl;
  std::cout << "* layer = " << layer_ << "\n"
	    << "* wire  = " << wire_ << "\n"
	    << "* tdc .size() = " << tdc_.size() << "\n   ";
  std::copy(tdc_.begin(), tdc_.end(),
	    std::ostream_iterator<int>(std::cout, " "));
  std::cout << "\n" 
	    << "* trailing .size() = " << trailing_.size() << "\n val = ";
  std::copy(trailing_.begin(), trailing_.end(),
	    std::ostream_iterator<int>(std::cout, " "));
  std::cout << "\n"
	    << "* drift time .size() = " << dt_.size() << "\n val = ";
  std::copy(dt_.begin(), dt_.end(),
	    std::ostream_iterator<double>(std::cout, " "));
  std::cout << "\n"
	    << "* drift length .size() = " << dl_.size() << "\n val = ";
  std::copy(dl_.begin(), dl_.end(),
	    std::ostream_iterator<double>(std::cout, " "));
  std::cout << "\n"
	    << "* trailing time .size() = " << trailingTime_.size() << "\n val = ";
  std::copy(trailingTime_.begin(), trailingTime_.end(),
	    std::ostream_iterator<double>(std::cout, " "));
  std::cout << "\n"
	    << "* wpos = " << wpos_ << "\n";
  std::cout << "* angle = " << angle_ << "\n";
  std::cout << "* belongTrack .size() = " << belongTrack_.size() << "\n val = ";
  std::copy(belongTrack_.begin(), belongTrack_.end(),
	    std::ostream_iterator<bool>(std::cout, " "));
  std::cout << "\n* dlRange .size() = " << dlRange_.size() << "\n val = ";
  std::copy(dlRange_.begin(), dlRange_.end(),
	    std::ostream_iterator<bool>(std::cout, " "));
  std::cout << "\n "
	    << "* clsize = " << clsize_ << "\n"
	    << "* mwpcflag = " << mwpcflag_ << "\n"
	    << "* mean wire = " << mwire_ << "\n"
	    << "* mean wire pos = " << mwpos_
	    << "\n-----------------" << std::endl;
  return;
}

