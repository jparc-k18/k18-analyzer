/**
 *  file: Hodo1Hit.cc
 *  date: 2017.04.10
 *
 */

#include "HodoWave1Hit.hh"

#include <cstring>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>

#include <std_ostream.hh>

#include "DebugCounter.hh"
#include "HodoParamMan.hh"
#include "HodoPHCMan.hh"
#include "RawData.hh"

namespace
{
  const std::string& class_name("HodoWave1Hit");
  const HodoParamMan& gHodo = HodoParamMan::GetInstance();
  const HodoPHCMan&   gPHC  = HodoPHCMan::GetInstance();

}


//______________________________________________________________________________
HodoWave1Hit::HodoWave1Hit( HodoRawHit *rhit, int ref_sample )
  : HodoHit(), m_raw(rhit), m_is_calculated(false), m_ref_sample(ref_sample)
{
  debug::ObjectCounter::increase(class_name);
}

//______________________________________________________________________________
HodoWave1Hit::~HodoWave1Hit( void )
{
  debug::ObjectCounter::decrease(class_name);
}

//______________________________________________________________________________
bool
//HodoWave1Hit::Calculate( void )
HodoWave1Hit::Calculate( bool tdc_flag )
{
  static const std::string func_name("["+class_name+"::"+__func__+"()]");

  if( m_is_calculated ){
    hddaq::cout << func_name << " already calculated" << std::endl;
    return false;
  }

  //if( tdc_flag && m_raw->GetNumOfTdcHits()!=1 )
  //return false;

  if( !gHodo.IsReady() ){
    hddaq::cout << func_name << " HodoParamMan must be initialized" << std::endl;
    return false;
  }

  if( !gPHC.IsReady() ){
    hddaq::cout << func_name << " HodoPHCMan must be initialized" << std::endl;
    return false;
  }

  int cid  = m_raw->DetectorId();
  int plid = m_raw->PlaneId();
  int seg  = m_raw->SegmentId();

  // hit information
  m_multi_hit_l = m_raw->SizeTdc1()  > m_raw->SizeTdc2()?  m_raw->SizeTdc1()  : m_raw->SizeTdc2();
  m_multi_hit_t = m_raw->SizeTdcT1() > m_raw->SizeTdcT2()? m_raw->SizeTdcT1() : m_raw->SizeTdcT2();
  int UorD=0;
  int nh_fadc = m_raw->SizeAdc1();

  double charge = 0;
  double charge_pedestal = 0;
  for (int i=0; i<nh_fadc; i++) {
    int adc = m_raw->GetAdc1(i);
    double pedestal = gHodo.GetP0( cid, i, seg, UorD+2);
    double V_per_ch = gHodo.GetP1( cid, i, seg, UorD+2);
    if (pedestal < 0)
      pedestal = gHodo.GetP0( cid, 0, seg, UorD+2);
    if (V_per_ch < 0)
      V_per_ch = gHodo.GetP1( cid, 0, seg, UorD+2);

    double V = (adc-pedestal)*V_per_ch;

    m_a_wf.push_back(V);

    int tdc = (i-m_ref_sample);
    double T_per_ch = gHodo.GetGain( cid, plid, seg, UorD+2);
    double T_offset = gHodo.GetOffset( cid, plid, seg, UorD+2);
    double time = tdc*T_per_ch + T_offset;
    /*
    if( !gHodo.GetTime( cid, plid, seg, UorD+2, tdc, time ) ){
      hddaq::cerr << "#E " << func_name
		  << " something is wrong at GetTime("
		  << cid  << ", " << plid << ", " << seg  << ", "
		  << UorD << ", " << tdc  << ", " << time << ")" << std::endl;
      return false;
    }
    */
    m_t_wf.push_back(time);

    if (time>-100 && time<200) {
      double T_per_ch = gHodo.GetGain( cid, 0, seg, UorD+2);
      double q = V/50.*T_per_ch*1000.; //pC
      charge += q;
    }
    
    if (time>-400 && time<-100) {
      double T_per_ch = gHodo.GetGain( cid, 0 , seg, UorD+2);
      double q = V/50.*T_per_ch*1000.;
      charge_pedestal += q;
    }
  }

  //if (charge<0)
  //charge *= -1.0;

  m_q.push_back(charge);
  m_q_pedestal.push_back(charge_pedestal);
  


  double dE = 0.;
  if( !gHodo.GetDe2( cid, plid, seg, UorD, charge, dE ) ){
    /*
    hddaq::cerr << "#E " << func_name
		<< " something is wrong at GetDe("
		<< cid  << ", " << plid << ", " << seg << ", "
		<< UorD << ", " << charge  << ", " << dE  << ")" << std::endl;
    return false;
    */
  }


  m_a.push_back(dE);

  int mhit = m_multi_hit_l;
  for( int m=0; m<mhit; ++m ){
    int    tdc  = -999;
    double time = -999.;
    if(0 == UorD){
      tdc = m_raw->GetTdc1(m);
    }else{
      tdc = m_raw->GetTdc2(m);
    }

    if( tdc<0 ) continue;

    if( !gHodo.GetTime( cid, plid, seg, UorD, tdc, time ) ){
      hddaq::cerr << "#E " << func_name
		  << " something is wrong at GetTime("
		  << cid  << ", " << plid << ", " << seg  << ", "
		  << UorD << ", " << tdc  << ", " << time << ")" << std::endl;
      return false;
    }
    m_t.push_back(time);

    double ctime = time;
    //gPHC.DoCorrection(cid, plid, seg, UorD, time, dE, ctime );
    m_ct.push_back(ctime);

    m_flag_join.push_back(false);
  }


  m_is_calculated = true;
  return true;
}

// ____________________________________________________________
bool
HodoWave1Hit::JoinedAllMhit()
{
  bool ret = true;
  for(int i = 0; i<m_flag_join.size(); ++i){
    ret = ret & m_flag_join[i];
  }// for(i)
  return ret;
}
