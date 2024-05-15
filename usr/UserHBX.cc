/**
 *  file: UserHBX.cc
 *  Ref. file: UserSkeleton.cc
 *  date: 2023.06.13
 *
 */

#include <iostream>
#include <sstream>
#include <cmath>
#include <iomanip>
#include <string>

#include "ConfMan.hh"
#include "UserParamMan.hh"
#include "HodoParamMan.hh"
#include "DCGeomMan.hh"
#include "DetectorID.hh"
#include "RMAnalyzer.hh"
#include "MathTools.hh"
#include "RawData.hh"
#include "HodoRawHit.hh"
#include "RootHelper.hh"
#include "UnpackerManager.hh"
#include "VEvent.hh"
#include "Unpacker.hh"

///////////include each Ge gain (ADC ch/keV)/////////////

namespace
{
  using namespace root;
  const auto qnan = TMath::QuietNaN();
  const std::string& classname("EventHBX");
  RMAnalyzer& gRM = RMAnalyzer::GetInstance();
  const UserParamMan& gUser = UserParamMan::GetInstance();
}

namespace
{
  using namespace hddaq::unpacker;
  const UnpackerManager& gUnpacker = GUnpacker::get_instance();
}

//______________________________________________________________________________
struct Event
{
  int runnum;
  int evnum;
  int spill;
  int rm1evnum;
  int rm2evnum;
  int rm1spill;
  int rm2spill;

  //Event Sync Clock
  Int_t hrtdc[NumOfChHRTDC][MaxDepth];
  Int_t syncclock[MaxDepth];

  //Trig flag
  int trignhits;
  int trigpat[NumOfSegTrig];
  int trigflag[NumOfSegTrig];

  //Ge Reset
  int gereset[NumOfSegGe];

  //Ge ADC
  //double geadc[NumOfSegGe];
  double geadc[NumOfSegGe*2];
  //Ge TFA
  int genhits;
  int gehitpat[MaxHits];
  double getfa[NumOfSegGe][MaxDepth];

  //Ge CRM
  int gecrmnhits;
  int gecrmhitpat[MaxHits];
  double gecrm[NumOfSegGe][MaxDepth];

  //Ge BGO
  int gebgonhits;
  int gebgohitpat[MaxHits];
  double gebgot[NumOfSegBGO][MaxDepth];

  //HUL-Scaler and HUL-RM for HBX and E70
  int scaler[NumOfPlaneScaler][NumOfSegScaler];
  int evnum_hulrm[NumOfPlaneHulRm];
  int spill_hulrm[NumOfPlaneHulRm];

  void clear();
};

//______________________________________________________________________________
void Event::clear()
{
  evnum      = 0;
  spill      = 0;
  trignhits  = 0;
  genhits    = 0;
  gecrmnhits = 0;
  gebgonhits = 0;
  rm1evnum   = 0;
  rm2evnum   = 0;
  rm1spill   = 0;
  rm2spill   = 0;
  for(int i=0; i<NumOfChHRTDC; ++i){
    for(int j=0; j<MaxDepth; ++j){
      hrtdc[i][j]      = -999;
    }
  }
  for(int i=0; i<MaxDepth; ++i){
    syncclock[i] = -999;
  }
  for(int i=0; i<NumOfPlaneScaler; ++i){
    for(int j=0; j<NumOfSegScaler; ++j){
      scaler[i][j] = -1;      
    }
  }
  for(int i=0; i<NumOfPlaneHulRm; ++i){
    evnum_hulrm[i] = -1;
    spill_hulrm[i] = -1;
  }
  for(int i=0; i<NumOfSegTrig; ++i){
    trigpat[i]  = -1;
    trigflag[i] = -1;
  }
  for(int i=0; i<NumOfSegGe; ++i){
    gereset[i] = qnan;
    //geadc[i]   = qnan;
  }
  for(int i=0; i<NumOfSegGe*2; i++){
    geadc[i] = qnan;
  }
  for(int m=0; m<MaxDepth; ++m){
    for(int i=0; i<NumOfSegGe; ++i){
      getfa[i][m] = qnan;
      gecrm[i][m] = qnan;
    }
    for(int i=0; i<NumOfSegBGO; ++i){
      gebgot[i][m] = qnan ;
    }
  }
  for(int i=0; i<MaxHits; ++i){
    gebgohitpat[i]  = qnan;
    gehitpat[i]     = qnan;
    gecrmhitpat[i]  = qnan;
  }
}

//______________________________________________________________________________
struct Dst
{
  int evnum;
  int spill;
  int rm1spill;
  int rm2spill;
  int rm1evnum;
  int rm2evnum;
  
  //  Double_t hrtdc[NumOfChHRTDC][MaxDepth];
  Int_t hrtdc[NumOfChHRTDC][MaxDepth];
  Int_t syncclock[MaxDepth];

  //Trig flag
  int trignhits;
  int trigpat[NumOfSegTrig];
  int trigflag[NumOfSegTrig];

  //Ge Reset
  int geReset[NumOfSegGe];
  //Ge ADC
  //double geAdc[NumOfSegGe];
  double geAdc[NumOfSegGe*2];
  //Ge TFA
  int nhGeTfa;
  double geTfa[NumOfSegGe][MaxDepth];

  //Ge CRM
  int nhGeCrm;
  double geCrm[NumOfSegGe][MaxDepth];

  //Ge BGO
  int nhBgo;
  double bgoTdc[NumOfSegBGO][MaxDepth];

  //HUL-Scaler and HUL-RM for HBX and E70
  int scaler[NumOfPlaneScaler][NumOfSegScaler];
  int evnum_hulrm[NumOfPlaneHulRm];
  int spill_hulrm[NumOfPlaneHulRm];

  void clear();
};

//______________________________________________________________________________
void
Dst::clear()
{
  evnum     = 0;
  spill     = 0;
  trignhits = 0;
  nhGeTfa   = 0;
  nhGeCrm   = 0;
  nhBgo     = 0;
  rm1spill  = 0;
  rm2spill  = 0;
  rm1evnum  = 0;
  rm2evnum  = 0;


  for(int i=0; i<NumOfPlaneScaler; ++i){
    for(int j=0; j<NumOfSegScaler; ++j){
      scaler[i][j] = -1;
    }
  }
  for(int i=0; i<NumOfPlaneHulRm; ++i){
    evnum_hulrm[i] = -1;
    spill_hulrm[i] = -1;
  }

  for(int i=0; i<NumOfChHRTDC; ++i){
    for(int j=0; j<MaxDepth; ++j){
      hrtdc[i][j] = -999;
    }
  }
  for(int i=0; i<MaxDepth; ++i){
    syncclock[i] = -999;
  }
  for(int i=0; i<NumOfSegTrig; ++i){
    trigpat[i]=-1;
    trigflag[i]=-1;
  }
  for(int i=0; i<NumOfSegGe; ++i){
    geReset[i]=qnan;
    //geAdc[i]=qnan;
    for(int m=0; m<MaxDepth; ++m){
      geTfa[i][m]=qnan;
      geCrm[i][m]=qnan;
    }
  }
  for(int i=0; i<NumOfSegGe; i++){
    geAdc[i] = qnan;
  }
  for(int i=0; i<NumOfSegBGO; ++i){
    for(int m=0; m<MaxDepth; ++m){
      bgoTdc[i][m]=qnan;
    }
  }
}
//______________________________________________________________________________
namespace root
{
Event  event;
Dst    dst;
TH1   *h[MaxHist];
TTree *tree;
TTree *hbx;

enum eDetHid{
  GeHid    = 270000,
  BGOHid = 1140000,
};
}

//______________________________________________________________________________
bool
ProcessingBegin()
{
  event.clear();
  dst.clear();
  return true;
}

//______________________________________________________________________________
bool
ProcessingNormal()
{
  // static const double MinTFACUT = gUser.GetParameter("TCut", 0);
  // static const double MaxTFACUT = gUser.GetParameter("TCut", 1);
  static const double MinTFACUT = 3800;
  static const double MaxTFACUT = 4300;//decided by run4028 tmpdata
  // static const double MinCRMCUT = gUser.GetParameter("CCut", 0);
  // static const double MaxCRMCUT = gUser.GetParameter("CCut", 1);
  static const double MinCRMCUT = 3800;
  static const double MaxCRMCUT = 4300;//decided by run4028 2024/4/22 tmp data
  //  static const int MinBGOCUT = gUser.GetParameter("BCut", 0);
  //  static const int MaxBGOCUT = gUser.GetParameter("BCut", 1);

  // rawData = new RawData;
  RawData rawData;


  //  gRM.Decode();

  // event.runnum = gRM.RunNumber();
  // event.evnum  = gRM.EventNumber();
  // event.spill  = gRM.SpillNumber();
  // dst.evnum    = gRM.EventNumber();
  // dst.spill    = gRM.SpillNumber();
  // std::cout << "evnum: " << event.evnum << std::endl;
  // std::cout << "spill num: " << event.spill << std::endl;
    
  rawData.DecodeHits("HBXTFlag");
  std::bitset<NumOfSegTrig> trigger_flag;
  //   for(const auto& hit: rawData.GetHodoRawHitContainer("TFlag")){
  for(const auto& hit: rawData.GetHodoRawHitContainer("HBXTFlag")){
    Int_t seg = hit->SegmentId();
    Int_t tdc = hit->GetTdc();
    //tdc_clock = hit->GetTdc();
    //tdc_dead = hit->GetTdc;
    if( tdc>0 ){
      event.trigpat[trigger_flag.count()] = seg;
      event.trigflag[seg] = tdc;
      dst.trigpat[trigger_flag.count()] = seg;
      dst.trigflag[seg] = tdc;
      trigger_flag.set(seg);
      HF1(10, seg);
      HF1(10+seg, tdc);
    }
  }
  //  if(trigger_flag[trigger::kSpillEnd]) return true;
  //int seg_bgo[NumOfSegGe];
    //int seg_bgo0[4]
    int seg_bgo[NumOfSegBGO+16];
    
    int seg_bgo1[4] = {0, 1, 10, 11};//slot1-1
    int seg_bgo2[4] = {1, 2, 3, 4};
    int seg_bgo3[4] = {4, 5, 6, 7};
    int seg_bgo4[4] = {7, 8, 9, 10};
    int seg_bgo5[4] = {16, 17, 18, 19};//slot2-1
    int seg_bgo6[4] = {19, 20, 21, 22};
    int seg_bgo7[4] = {12, 13, 22, 23};
    int seg_bgo8[4] = {13, 14, 15, 16};
    int seg_bgo9[4] = {24, 25, 34, 35};//slot3-1
    int seg_bgo10[4] = {25, 26, 27, 28};
    int seg_bgo11[4] = {28, 29, 30, 31};
    int seg_bgo12[4] = {31, 32, 33, 34};
    int seg_bgo13[4] = {36, 37, 46, 47};//slot4-1
    int seg_bgo14[4] = {37, 38, 39, 40};
    int seg_bgo15[4] = {40, 41, 42, 43};
    int seg_bgo16[4] = {43, 44, 45, 46};
    for(int i=0; i<64; i++){
      if(0<=i&&i<=3)seg_bgo[i] = seg_bgo1[i];
      if(4<=i&&i<=7)seg_bgo[i] = seg_bgo2[i-4];
      if(8<=i&&i<=11)seg_bgo[i] = seg_bgo3[i-8];
      if(12<=i&&i<=15)seg_bgo[i] = seg_bgo4[i-12];
      if(16<=i&&i<=19)seg_bgo[i] = seg_bgo5[i-16];
      if(20<=i&&i<=23)seg_bgo[i] = seg_bgo6[i-20];
      if(24<=i&&i<=27)seg_bgo[i] = seg_bgo7[i-24];
      if(28<=i&&i<=31)seg_bgo[i] = seg_bgo8[i-28];
      if(32<=i&&i<=35)seg_bgo[i] = seg_bgo9[i-32];
      if(36<=i&&i<=39)seg_bgo[i] = seg_bgo10[i-36];
      if(40<=i&&i<=43)seg_bgo[i] = seg_bgo11[i-40];
      if(44<=i&&i<=47)seg_bgo[i] = seg_bgo12[i-44];
      if(48<=i&&i<=51)seg_bgo[i] = seg_bgo13[i-48];
      if(52<=i&&i<=55)seg_bgo[i] = seg_bgo14[i-52];
      if(56<=i&&i<=59)seg_bgo[i] = seg_bgo15[i-56];
      if(60<=i&&i<=63)seg_bgo[i] = seg_bgo16[i-60];
    }
    
    int tdc_bgo[NumOfSegBGO];
    int flag_bgo[NumOfSegGe];
    for(int seg_b = 0; seg_b<NumOfSegBGO; seg_b++){
      int nhit_bgo =  gUnpacker.get_entries( DetIdBGO, 0, seg_b, 0, 0 );
      if(nhit_bgo>0){
	tdc_bgo[seg_b] = gUnpacker.get( DetIdBGO, 0, seg_b, 0, 0, 0 );
	HF1 (GeHid + 74, double(seg_b) );
	//tdc_bgo[seg] = 100;
	}
    }
    for(int k=0; k<NumOfSegGe; k++){
      flag_bgo[k]=0;
      if(tdc_bgo[seg_bgo[4*k+3]]>500||tdc_bgo[seg_bgo[4*k+2]]>500||tdc_bgo[seg_bgo[4*k+1]]>500||tdc_bgo[seg_bgo[4*k]]>500){
	flag_bgo[k] = 1;
      }
      int nhit_adc[NumOfSegGe];
      int nhit_tdc[NumOfSegGe];
      nhit_adc[k] = gUnpacker.get_entries( DetIdGe, 0, k, 0, 0 ); //adc
      nhit_tdc[k] = gUnpacker.get_entries( DetIdGe, 0, k, 0, 1 ); //adc
      //std::cout<<"kama chan"<<std::endl;
      if(nhit_adc[k]>0&&flag_bgo[k]==0&&nhit_tdc[k]>0){
	int adc = gUnpacker.get( DetIdGe, 0, k, 0, 0, 0);
	int tdc = gUnpacker.get( DetIdGe, 0, k, 0, 1, 0); //tdc(first hit)
	double energy = adc*0.21;
	if(MinTFACUT<tdc&&tdc<MaxTFACUT){
	HF1 ( GeHid + 100*(k+1) +16, double(adc) ); 
	HF1 ( GeHid + 100*(k+1) +50, double(energy) );
	if(dst.trigflag[2]>0&&dst.trigflag[6]>0)HF1 ( GeHid + 100*(k+1) +17, double(adc) ); //lsoxge adc(spill off)
	if(dst.trigflag[1]>0&&dst.trigflag[5]>0){
	  HF1 ( GeHid + 100*(k+1) +18, double(adc) ); //l1 and apill onadc
	  HF1 ( GeHid + 100*(k+1) +51, double(energy) ); //l1 and apill onadc
	}
	//if(dst.trigflag[1]>0)HF1 ( GeHid + 100*(k+1) +19, double(adc) ); //l1
	}//tfa cut range
      }//nhitadc&nhittdc&flagbgo
    }//numofsegge k
    
  //---Ge ADC & TDC-----------------------------------------------------
    for( int seg=0; seg<NumOfSegGe*2; ++seg ){
    int nhit_a = gUnpacker.get_entries( DetIdGe, 0, seg, 0, 0 ); //adc
    int nhit_t = gUnpacker.get_entries( DetIdGe, 0, seg, 0, 1 ); //tfa_leading
    int nhit_c = gUnpacker.get_entries( DetIdGe, 0, seg, 0, 4 ); //crm_leading
    int nhit_r = gUnpacker.get_entries( DetIdGe, 0, seg, 0, 7 ); //reset_time
    int nhit_bgo[NumOfSegBGO];
    nhit_bgo[seg]= gUnpacker.get_entries( DetIdGe, 0, seg, 0, 7 ); //bgo 
    if( nhit_a>0 ){
      int adc = gUnpacker.get( DetIdGe, 0, seg, 0, 0, 0);
      double energy = adc*0.20;
      event.geadc[seg] = adc;
      dst.geAdc[seg] = adc;
      ///adc wo tdc cut///////
      if(adc>0){
      HF1( GeHid+100*(seg+1)+10, double(adc) );
      HF1( GeHid+100*(seg+1)+52, double(energy) );
      //////adc w973U crm cut///
      if( nhit_t>0){
        int tfa  = gUnpacker.get( DetIdGe, 0, seg, 0, 1, 0);
        if( tfa>0 ){
          HF1( GeHid+100*(seg+1)+21, double(adc) );
          HF2( GeHid+100*(seg+1)+24, double(tfa), double(adc) );
	  if(100<adc&&adc<300)HF1( GeHid+100*(seg+1)+61, double(tfa) ); //20keV-60keV
	  if(300<adc&&adc<500)HF1( GeHid+100*(seg+1)+62, double(tfa) ); //60keV-100keV
	  if(500<adc&&adc<700)HF1( GeHid+100*(seg+1)+63, double(tfa) ); //100keV-140keV
	  if(700<adc&&adc<900)HF1( GeHid+100*(seg+1)+64, double(tfa) ); //140keV-180keV
	  if(900<adc&&adc<1100)HF1( GeHid+100*(seg+1)+65, double(tfa) ); //180keV-220keV
	  if(1100<adc&&adc<1300)HF1( GeHid+100*(seg+1)+66, double(tfa) ); //220keV-260keV
	  if(1300<adc&&adc<1500)HF1( GeHid+100*(seg+1)+67, double(tfa) ); //260keV-300keV
	  ////ge hitpat//////
	  HF1( GeHid+73, double(seg) );
        }
        else{
          HF1( GeHid+100*(seg+1)+22, double(adc) );
        }
        //adc w/ tfa
        if( MinTFACUT<tfa && tfa<MaxTFACUT ) HF1( GeHid+100*(seg+1)+23, double(adc) );
	//only L1 and spill on trig adc before tdccut///
	if(dst.trigflag[1]>0&&dst.trigflag[5]>0){
	  HF1 ( GeHid + 100*(seg+1) +19, double(adc) );
	  //HF1 ( GeHid + 100*(seg+1) +55, double(energy) );
	}
      //------adc(w crm)wTrigger flag-----//
	if(MinTFACUT<tfa && tfa<MaxTFACUT){ 
	  if(dst.trigflag[2]>0){//lsoxge
	    if(dst.trigflag[5]>0){//spill on 
	      HF1( GeHid+100*(seg+1)+11, double(adc) );
	      HF1( GeHid+100*(seg+1)+12, double(adc) );
	      HF1( GeHid+100*(seg+1)+53, double(energy) );
	    }
	    if(dst.trigflag[6]>0){//spill off
	      HF1( GeHid+100*(seg+1)+13, double(adc) );
	      HF1( GeHid+100*(seg+1)+14, double(adc) );
	      HF1( GeHid+100*(seg+1)+54, double(energy) );
	    }
	  }//LSOxGe trig
	  //////L1 and spill on adc///////
	  if(dst.trigflag[1]>0&&dst.trigflag[5]>0){
	    HF1 ( GeHid + 100*(seg+1) +15, double(adc) );
	    HF1 ( GeHid + 100*(seg+1) +55, double(energy) );
	  }
	}//tfa cut
      }//for nhit_t
      HF2( GeHid+100 +0, seg+0.5, double(adc) );
      }//for adc
      /*
      //adc w/ crm (no need at e96?)
      if( nhit_t>0){
	int crm  = gUnpacker.get( DetIdGe, 0, seg, 0, 1, 0);
	if( crm>0 ){
	  HF1( GeHid+100*(seg+1)+31, double(adc) );
	  HF2( GeHid+100*(seg+1)+34, double(crm), double(adc) );
	}
	else{
	  HF1( GeHid+100*(seg+1)+32, double(adc) );
	}
	//adc w/ crm cut
	if( MinCRMCUT<crm && crm<MaxCRMCUT ) HF1( GeHid+100*(seg+1)+33, double(adc) );
      }
      */
      //adc w/ reset
      if( nhit_r>0){
	int reset_time = gUnpacker.get( DetIdGe, 0, seg, 0, 7);
	HF2( GeHid+100*(seg+1)+44, double(reset_time), double(adc) );
      }//nhit_r
    }
    //-------tdc hist------------------------------------------------
    //tfa
    if( nhit_t>0 ){
      HF1( GeHid+1, seg+0.5); //0 origin
      for(int i = 0; i<nhit_t; ++i){
	int tfa  = gUnpacker.get( DetIdGe, 0, seg, 0, 1, i );
	if( tfa>0 ){
	  HF1( GeHid+100*(seg+1)+20, double(tfa) ); //0 origin
	  //------tfa w/adc range select------//
	  HF2( GeHid+100+2, seg+0.5, double(tfa) );
	}
	event.getfa[seg][i] = tfa;
	dst.geTfa[seg][i] = tfa;
      }//for nhit
      /////973U first hit tdc //////
      int tfa_first = gUnpacker.get( DetIdGe, 0, seg, 0, 1, 0 );
      if(tfa_first>0){
	HF1( GeHid+100*(seg+1)+25, double(tfa_first) ); //0 origin
      }
    }

    //crm
    if( nhit_c>0 ){
      for(int i = 0; i<nhit_c; ++i){
	int crm  = gUnpacker.get( DetIdGe, 0, seg, 0, 4, i );
	if(crm > 0){
	  HF1( GeHid+100*(seg+1)+30, double(crm) ); //0 origin
	  HF2( GeHid+100+3, seg+0.5, double(crm) );
	  event.gecrm[seg][i] = crm;
	  dst.geCrm[seg][i] = crm;
	}
      }
      ////////////671 first hit tdc/////////////
      int crm_first = gUnpacker.get( DetIdGe, 0, seg, 0, 4, 0 );
      if(crm_first>0){
        HF1( GeHid+100*(seg+1)+35, double(crm_first) ); //0 origin
      }
    }//for nhit_c

    //reset
    if( nhit_r>0 ){
      int reset_time = gUnpacker.get( DetIdGe, 0, seg, 0, 7);
      HF1( GeHid+100*(seg+1)+40, double(reset_time) );
      HF2( GeHid+100+4, seg+0.5, double(reset_time) );
      event.gereset[seg] = reset_time;
      dst.geReset[seg] = reset_time;
    }
  }

  //---BGO--------------------------------------------------------------
  for(int seg = 0; seg<NumOfSegBGO; ++seg){
    int nhit = gUnpacker.get_entries( DetIdBGO, 0, seg, 0, 0 );
    if( nhit>0 ){
      HF1( BGOHid+1, seg+0.5); //0 origin
      for(int i = 0; i<nhit; ++i){
	int tdc = gUnpacker.get( DetIdBGO, 0, seg, 0, 0, i )  ;
	if(tdc > 0){
	  HF1( BGOHid+100*(seg+1)+0, double(tdc) );
	  HF2( BGOHid+0, seg+0.5, double(tdc) );
	}
	event.gebgot[seg][i] = tdc;
	dst.bgoTdc[seg][i] = tdc;
      }//for i
      /////BGO first hit/////
      int tdc = gUnpacker.get( DetIdBGO, 0, seg, 0, 0, 0 );
      if(tdc > 0){
	HF1( BGOHid+100*(seg+1)+10, double(tdc) );
	///////////////////////////////////
      }
    }
  }
  // need to edit for multi-hit depth? 
  // and debug because scaler data is not at all decoded 
  // int scaler_0;
  // int scaler_15;
  // int daq_live;
  // for( int seg=0; seg<NumOfSegScaler; ++seg ){
  //   int nhit = 0;
  //   //    nhit = gUnpacker.get_entries( DetIdScaler, 2, 0, seg, 0 ); //original
  //   nhit = gUnpacker.get_entries( DetIdScaler, 1, 0, seg, 0 ); // confirm if this line is right 
  //   // especially in the order of the number 
  //   // (this order seems to be different from other modules)
  //   if( nhit>0 ){
  //     // std::cout << " / / / / / /  / / / / / / / / / / / / / / " << std::endl;
  //     // std::cout << "nhit: " << std::endl;
  //     // std::cout << " / / / / / /  / / / / / / / / / / / / / / " << std::endl;
  //     //      int data = gUnpacker.get( DetIdScaler, 2, 0, seg, 0 ); //original
  //     HF1 (GeHid +75, double(seg) ); //all scaler hitpat
  //     if(0<=seg&&seg<=15)HF1 (GeHid +76, double(seg) ); //trigflag scaler hitpat
  //     if(0<=seg&&seg<=15)HF1 (GeHid +77, double(seg) ); //973u crm scaler hitpat
  //     if(0<=seg&&seg<=15)HF1 (GeHid +78, double(seg) ); //671 crm scaler hitpat
  //     if(0<=seg&&seg<=15)HF1 (GeHid +79, double(seg) ); //reset scaler hitpat
  //     int data = gUnpacker.get( DetIdScaler, 1, 0, seg, 0 );
  //     event.scaler[seg] = data;
  //     scaler_0 = event.scaler[0];
  //     scaler_15 = event.scaler[15];
  //     daq_live = scaler_0 - scaler_15;
  //   }
  // }
  // if(dst.trigflag[3]>0)HF1 (GeHid+71, double(daq_live));//spill on end daqlivetime                           
  // if(dst.trigflag[4]>0)HF1 (GeHid+72, double(daq_live));//spill off end daqlivetime  

  // // ---HUL-RM--------------------------------------------------------
  std::cout << "hul-rm" << std::endl;
  for(int plane=0; plane<NumOfPlaneHulRm; ++plane){
    static const int device_id = gUnpacker.get_device_id("HUL-RM");
    int dtype = 0; // evnum
    int nhit = gUnpacker.get_entries( device_id , plane, 0, 0, dtype );
    if( nhit<=0 ) continue;
    int data = 0;
    data = gUnpacker.get( device_id, plane, 0, 0, 0 );
    event.evnum_hulrm[plane] = data;
    dst.evnum_hulrm[plane] = data;
    data = gUnpacker.get( device_id, plane, 0, 0, 1 );
    event.spill_hulrm[plane] = data;
    dst.spill_hulrm[plane] = data;
  }
  // // ---HUL-Scaler----------------------------------------------------
  std::cout << "hul-scaler" << std::endl;
  for(int plane=0; plane<NumOfPlaneScaler; ++plane){
    for(int seg=0; seg<NumOfSegScaler; ++seg ){
      int nhit = 0;
      nhit = gUnpacker.get_entries( DetIdScaler, plane, 0, seg, 0 ); 
      if( nhit>0 ){
  	int data = gUnpacker.get( DetIdScaler, plane, 0, seg, 0);
  	event.scaler[plane][seg] = data;
      }
    }
  }

 
  if(trigger_flag[trigger::kSpillOnEnd]) return true;
  
  return true;
}

//______________________________________________________________________________
bool
ProcessingEnd()
{
  tree->Fill();
  hbx->Fill();
  return true;
}
//______________________________________________________________________________
void
InitializeEvent( void )
{
  event.evnum = 0;
  event.spill = 0;
  event.rm1evnum = 0;
  event.rm2evnum = 0;
  event.rm1spill = 0;
  event.rm2spill = 0;
  event.gebgonhits  = 0;
  event.trignhits = 0;

  for( int it=0; it<NumOfChHRTDC; ++it){
    for (int m=0; m<MaxDepth; ++m){
      event.hrtdc[it][m] = -9999;
    }
  }
  for (int m=0; m<MaxDepth; ++m){
    event.syncclock[m] = -9999;
  }
  for( int it=0; it<NumOfSegTrig; it++){
    event.trigpat[it] = -1;
    event.trigflag[it] = -1;
    dst.trigpat[it]  = -1;
    dst.trigflag[it] = -1;
  }

  //Ge Reset
  for( int it=0; it<NumOfSegGe; it++){
    event.gereset[it] = 0;
  }

  //Ge ADC
  for( int it=0; it<NumOfSegGe*2; it++){
    event.geadc[it] = -9999.;
  }

  //Ge TDC
  for( int it=0; it<MaxHits; ++it ){
    event.gehitpat[it]  = -1;
  }
  for( int it=0; it<NumOfSegGe; it++){
    for(int m = 0; m<MaxDepth; ++m){
      event.getfa[it][m] = -9999.;
      event.gecrm[it][m] = -9999.;
    }
  }

  //Ge BGO
  for( int it=0; it<MaxHits; ++it ){
    event.gebgohitpat[it]  = -1;
  }
  for( int it=0; it<NumOfSegBGO; it++){
    for(int m = 0; m<MaxDepth; ++m){
      event.gebgot[it][m] = -9999.;
    }
  }

  // HBX and E70 scaler 
  for(int i=0; i<NumOfPlaneScaler; ++i){
    for(int j=0; j<NumOfSegScaler; ++j){
      event.scaler[i][j] = -1;      
    }
  }

  ////Dst////////////////////////
  dst.evnum = 0;
  dst.spill = 0;
  dst.rm1evnum = 0;
  dst.rm2evnum = 0;
  dst.rm1spill = 0;
  dst.rm2spill = 0;
  dst.trignhits = 0;

  dst.nhGeTfa = 0;
  dst.nhGeCrm = 0;
  dst.nhBgo = 0;

  for( int it=0; it<NumOfChHRTDC; ++it){
    for (int m=0; m<MaxDepth; ++m){
      dst.hrtdc[it][m] = -9999;
    }
  }
  for (int m=0; m<MaxDepth; ++m){
    dst.syncclock[m] = -9999;
  }

  for( int it=0; it<NumOfSegGe; ++it ){
    dst.geReset[it] = -9999;
    //dst.geAdc[it] = -9999.;
    for( int m=0; m<MaxDepth; ++m ){
      dst.geTfa[it][m] = -9999.;
      dst.geCrm[it][m] = -9999.;
      }
  }
  for(int it=0; it<NumOfSegGe*2; ++it){
    dst.geAdc[it] = -9999;
  }
  for( int it=0; it<NumOfSegBGO; ++it ){
    for( int m=0; m<MaxDepth; ++m ){
      dst.bgoTdc[it][m] = -9999.;
    }
  }

  // HUL-Scaler and HUL-RM for HBX and E70
  for(int i=0; i<NumOfPlaneScaler; ++i){
    for(int j=0; j<NumOfSegScaler; ++j){
      dst.scaler[i][j] = -1;      
    }
  }

}


//______________________________________________________________________________
namespace
{

  const int    NbinAdc = 8192;
  const double MinAdc  =   0.;
  const double MaxAdc  = 8192.;
  /*
  const int    NbinGeTdc = 16384;
  const double MinGeTdc  =   0.;
  const double MaxGeTdc  = 16384.;
  */
  const int    NbinGeTdc = 8192;
  const double MinGeTdc  =   0.;
  const double MaxGeTdc  = 8192.;
  
  const int    NbinBGOTdc = 8192;
  const double MinBGOTdc  =   0.;
  const double MaxBGOTdc  = 8192.;

}

//______________________________________________________________________________
bool
ConfMan::InitializeHistograms( void )
{
  HBTree("tree","tree of hbx");
  tree->Branch("runnum", &event.runnum, "runnum/I");
  tree->Branch("evnum",  &event.evnum,  "evnum/I");
  tree->Branch("rm1evnum",  &event.rm1evnum,  "rm1evnum/I");
  tree->Branch("rm2evnum",  &event.rm2evnum,  "rm2evnum/I");
  tree->Branch("rm1spill",  &event.rm1spill,  "rm1spill/I");
  tree->Branch("rm2spill",  &event.rm1spill,  "rm2spill/I");
  //  tree->Branch("hrtdc",  &event.hrtdc,  "hrtdc/I");
  tree->Branch("spill",  &event.spill,  "spill/I");
  tree->Branch("hrtdc",  event.hrtdc,  Form("hrtdc[%d][%d]/I", NumOfChHRTDC, MaxDepth));
  tree->Branch("syncclock",  event.syncclock,  Form("syncclock[%d]/I", MaxDepth));
  //Trig
  tree->Branch("trignhits", &event.trignhits, "trignhits/I");
  HB1(  1, "Status", 20, 0., 20. );
  HB1( 10, "Trigger HitPat", NumOfSegTrig, 0., Double_t(NumOfSegTrig) );
  for(Int_t i=0; i<NumOfSegTrig; ++i){
    HB1( 10+i+1, Form("Trigger Trig %d", i+1), 0x1000, 0, 0x1000 );
  }

  //---Ge ADC & TDC-----------------------------------------------------
  for( int i=1; i<=NumOfSegGe*2; ++i ){
    TString title10  = Form("Ge_671_%d Adc", i);
    TString title11  = Form("Ge_973U_%d Adc", i);
    TString title20  = Form("Ge-%d 973U CRM", i);
    TString title21  = Form("Ge-%d Adc(w 973CRM)", i);
    TString title22  = Form("Ge-%d Adc(w/o 973CRM)", i);
    TString title23  = Form("Ge-%d Adc(973CRM Cut)", i);
    TString title24  = Form("Adc%%973CRM-%d", i);
    TString title30  = Form("Ge-%d 671Crm", i);
    TString title31  = Form("Ge-%d Adc(w 671Crm)", i);
    TString title32  = Form("Ge-%d Adc(w/o 671Crm)", i);
    TString title33  = Form("Ge-%d Adc(671Crm Cut)", i);
    TString title34  = Form("Adc%%671Crm-%d", i);
    TString title40  = Form("Ge-%d Reset", i);
    TString title44  = Form("Adc%%Reset-%d", i);
    TString title45  = Form("Ge-%d 671crm first-%d", i);
    TString title46  = Form("Ge-%d 973crm first-%d", i);
    //TString title46  = Form("Ge-%d 973crm first-%d", i);
    /////////for fit/////////
    TString title12  = Form("Ge_671_%d Adc (200keV spill on) wtdc", i);
    TString title13  = Form("Ge_671_%d Adc (300keV spill on) wtdc", i);
    TString title14  = Form("Ge_671_%d Adc (200keV spill off) wtdc", i);
    TString title15  = Form("Ge_671_%d Adc (300keV spill off) wtdc", i);
    TString title16  = Form("Ge_671_%d Adc (spillon&L1) wtdc", i);
    TString title17  = Form("Ge_671_%d Adc all (wo/BGO) wtdc", i);
    TString title18  = Form("Ge_671_%d Adc lsoxge spilloff (wo/BGO) wtdc", i);
    TString title19  = Form("Ge_671_%d Adc L1 spill on (wo/BGO) wtdc", i);
    TString title19_2  = Form("Ge_671_%d Adc (L1&spillon) (no bgo select)", i);
    TString title51  = Form("Ge_671_%d energy L1&spillon(wo/BGO)", i);
    TString title55  = Form("Ge_671_%d energy L1&spillon(all)", i);
    TString title61  = Form("Ge_973CRM_%d_20-60keV", i);
    TString title62  = Form("Ge_973CRM_%d_60-100keV", i);
    TString title63  = Form("Ge_973CRM_%d_100-140keV", i);
    TString title64  = Form("Ge_973CRM_%d_140-180keV", i);
    TString title65  = Form("Ge_973CRM_%d_180-220keV", i);
    TString title66  = Form("Ge_973CRM_%d_220-260keV", i);
    TString title67  = Form("Ge_973CRM_%d_260-300keV", i);
    ////////////////////////////////////////////////////////////////
    if(1<=i&&i<=16){
      HB1( GeHid +100*i +10, title10, NbinAdc, MinAdc, MaxAdc );
      ////for fit ////////
      HB1( GeHid +100*i +11, title12, NbinAdc, MinAdc, MaxAdc );
      HB1( GeHid +100*i +12, title13, NbinAdc, MinAdc, MaxAdc );
      HB1( GeHid +100*i +13, title14, NbinAdc, MinAdc, MaxAdc );
      HB1( GeHid +100*i +14, title15, NbinAdc, MinAdc, MaxAdc );
      HB1( GeHid +100*i +15, title16, NbinAdc, MinAdc, MaxAdc );
      HB1( GeHid +100*i +16, title17, NbinAdc, MinAdc, MaxAdc );
      HB1( GeHid +100*i +17, title18, NbinAdc, MinAdc, MaxAdc );
      HB1( GeHid +100*i +18, title19, NbinAdc, MinAdc, MaxAdc );
      HB1( GeHid +100*i +19, title19_2, NbinAdc, MinAdc, MaxAdc );
      HB1( GeHid +100*i +50, title51, 1000, 0, 1000 );
      HB1( GeHid +100*i +55, title55, 1000, 0, 1000 );
      HB1( GeHid +100*i +61, title61, 3000, 0, 5000 );
      HB1( GeHid +100*i +62, title62, 3000, 0, 5000 );
      HB1( GeHid +100*i +63, title63, 3000, 0, 5000 );
      HB1( GeHid +100*i +64, title64, 3000, 0, 5000 );
      HB1( GeHid +100*i +65, title65, 3000, 0, 5000 );
      HB1( GeHid +100*i +66, title66, 3000, 0, 5000 );
      HB1( GeHid +100*i +67, title67, 3000, 0, 5000 );

      /////////////////////
    }
    if(17<=i&&i<=32){
      HB1( GeHid +100*i +10, title11, NbinAdc, MinAdc, MaxAdc );
    }

    HB1( GeHid +100*i +20, title20, NbinGeTdc, MinGeTdc, MaxGeTdc );
    HB1( GeHid +100*i +21, title21, NbinAdc, MinAdc, MaxAdc );
    HB1( GeHid +100*i +22, title22, NbinAdc, MinAdc, MaxAdc );
    HB1( GeHid +100*i +23, title23, NbinAdc, MinAdc, MaxAdc );
    HB2( GeHid +100*i +24, title24,
    	 NbinGeTdc/16, MinGeTdc, MaxGeTdc, NbinAdc/8, MinAdc, MaxAdc);
    HB1( GeHid +100*i +25, title45, NbinGeTdc, MinGeTdc, MaxGeTdc );
    
    HB1( GeHid +100*i +30, title30, NbinGeTdc, MinGeTdc, MaxGeTdc );
    HB1( GeHid +100*i +31, title31, NbinAdc, MinAdc, MaxAdc );
    HB1( GeHid +100*i +32, title32, NbinAdc, MinAdc, MaxAdc );
    HB1( GeHid +100*i +33, title33, NbinAdc, MinAdc, MaxAdc );
    HB2( GeHid +100*i +34, title34,
    	 NbinGeTdc/16, MinGeTdc, MaxGeTdc, NbinAdc/8, MinAdc, MaxAdc);
    HB1( GeHid +100*i +35, title46, NbinGeTdc, MinGeTdc, MaxGeTdc );


    HB1( GeHid +100*i +40, title40, NbinGeTdc, MinGeTdc, MaxGeTdc );
    HB2( GeHid +100*i +44, title44,
    	 NbinGeTdc/8, MinGeTdc, MaxGeTdc, NbinAdc/8, MinAdc, MaxAdc);

  }

  //Hits
  HB1( GeHid +0, "#Hits Ge",        NumOfSegGe+1, 0., double(NumOfSegGe+1) );
  HB1( GeHid +1, "Hitpat Ge",       NumOfSegGe,   0., double(NumOfSegGe)   );
  //daq livetime
  HB1( GeHid +71, "DAQlivetime(spillon)", 2000,21000000,24000000 );
  HB1( GeHid +72, "DAQlivetime(spilloff)", 2000, 110000000, 13000000);
  //hitpat
  HB1( GeHid +73, "Ge Hitpat",  NumOfSegGe+1, 0., double(NumOfSegGe+1) );
  HB1( GeHid +74, "BGO Hitpat",  NumOfSegBGO+1, 0., double(NumOfSegBGO+1) );
  //scaler hitpat
  HB1( GeHid +75, "ALL Scaler Hitpat",  NumOfSegGe+1, 0., double(NumOfSegGe+1) );
  HB1( GeHid +76, "TrigFlag scaler Hitpat",  NumOfSegGe+1, 0., double(NumOfSegGe+1) );
  HB1( GeHid +77, "973U CRM Hitpat",  NumOfSegGe+1, 0., double(NumOfSegGe+1) );
  HB1( GeHid +78, "671 CRM Hitpat",  NumOfSegGe+1, 0., double(NumOfSegGe+1) );
  HB1( GeHid +79, "Ge Reset Scaler Hitpat",  NumOfSegGe+1, 0., double(NumOfSegGe+1) );

  /*
  HB2( GeHid +100 +1, "Adc%Ch",
       NumOfSegGe, 0, NumOfSegGe, NbinAdc, MinAdc, MaxAdc);
  HB2( GeHid +100 +2, "Tfa%Ch",
       NumOfSegGe, 0, NumOfSegGe, NbinGeTdc, MinGeTdc, MaxGeTdc);
  HB2( GeHid +100 +3, "Crm%Ch",
       NumOfSegGe, 0, NumOfSegGe, NbinGeTdc, MinGeTdc, MaxGeTdc);
  HB2( GeHid +100 +4, "Reset%Ch",
       NumOfSegGe, 0, NumOfSegGe, NbinGeTdc, MinGeTdc, MaxGeTdc);
  */
  /*
  HB2( GeHid +9999, "Ge_segvsBGO_seg",
       NumOfSegGe, 0, NumOfSegGe, NumOfSegBGO, 0, NumOfSegBGO);
  */
  //---BGO--------------------------------------------------------------
  //HB1( BGOHid +0, "#Hits BGO",        NumOfSegBGO+1, 0., double(NumOfSegBGO+1) );
  //HB1( BGOHid +1, "Hitpat BGO",       NumOfSegBGO,   0., double(NumOfSegBGO)   );

  for( int i=1; i<=NumOfSegBGO; ++i ){
    TString title_1 = Form("BGO-%d Tdc", i);
    TString title_2 = Form("BGO-%d Tdc first", i);
    HB1( BGOHid +100*i +0, title_1, NbinBGOTdc, MinBGOTdc, MaxBGOTdc );
    HB1( BGOHid +100*i +10, title_2, NbinBGOTdc, MinBGOTdc, MaxBGOTdc );
  }


  /////tree/////////////////////////////////
  //Flag
  tree->Branch("trigpat", event.trigpat, Form("trigpat[%d]/I", NumOfSegTrig));
  tree->Branch("trigflag", event.trigflag, Form("trigflag[%d]/I", NumOfSegTrig));

  //Ge ADC
  tree->Branch("geadc",   event.geadc,  Form("geadc[%d]/D", NumOfSegGe*2));

  //Ge TDC
  tree->Branch("genhits",   &event.genhits,   "genhits/I");
  tree->Branch("gehitpat",   event.gehitpat,  Form("gehitpat[%d]/I", NumOfSegGe));
  tree->Branch("ge973crm",        event.getfa,       Form("ge973crm[%d][%d]/D", NumOfSegGe, MaxDepth));
  tree->Branch("ge671crm",        event.gecrm,       Form("ge671crm[%d][%d]/D", NumOfSegGe, MaxDepth));
  tree->Branch("gereset", event.gereset, Form("gereset[%d]/I", NumOfSegGe));

  //BGO
  tree->Branch("gebgonhits",   &event.gebgonhits,   "gebgonhits/I");
  tree->Branch("gebgohitpat",   event.gebgohitpat,  Form("gebgohitpat[%d]/I", NumOfSegBGO));
  tree->Branch("gebgot",        event.gebgot,       Form("gebgot[%d][%d]/D", NumOfSegBGO, MaxDepth));

  //Ge Scaler
  tree->Branch("scaler", event.scaler, Form("scaler[%d]/I", NumOfSegScaler));

  /////Dst/////////////////////////////////
  hbx = new TTree( "hbx", "Data Summary Table of hbx" );
  hbx->Branch("evnum", &dst.evnum, "evnum/I");
  hbx->Branch("spill", &dst.spill, "spill/I");
  hbx->Branch("rm1evnum", &dst.rm1evnum, "rm1evnum/I");
  hbx->Branch("rm2evnum", &dst.rm2evnum, "rm2evnum/I");
  hbx->Branch("rm1spill", &dst.rm1spill, "rm1spill/I");
  hbx->Branch("rm2spill", &dst.rm2spill, "rm2spill/I");
  hbx->Branch("hrtdc",  dst.hrtdc,  Form("hrtdc[%d][%d]/I", NumOfChHRTDC, MaxDepth));
  hbx->Branch("syncclock",  dst.syncclock,  Form("syncclock[%d]/I", MaxDepth)); 

  hbx->Branch("trignhits", &dst.trignhits, "trignhits/I");
  hbx->Branch("trigpat",    dst.trigpat,   "trigpat[trignhits]/I");
  hbx->Branch("trigflag",   dst.trigflag,  Form("trigflag[%d]/I", NumOfSegTrig));

  hbx->Branch("geReset", dst.geReset,  Form("geReset[%d]/I", NumOfSegGe));
  hbx->Branch("geAdc",   dst.geAdc,  Form("geAdc[%d]/D", NumOfSegGe*2));

  hbx->Branch("nhGeTfa", &dst.nhGeTfa, "nhGeTfa/I");
  hbx->Branch("nhGeCrm", &dst.nhGeCrm, "nhGeCrm/I");
  hbx->Branch("nhBgo", &dst.nhBgo, "nhBgo/I");

  hbx->Branch("geTfa",  dst.geTfa,   Form("geTfa[%d][%d]/D", NumOfSegGe, MaxDepth));
  hbx->Branch("geCrm",  dst.geCrm,   Form("geCrm[%d][%d]/D", NumOfSegGe, MaxDepth));
  hbx->Branch("bgoTdc", dst.bgoTdc,  Form("bgoTdc[%d][%d]/D", NumOfSegBGO, MaxDepth));

  HPrint();
  return true;
}

//______________________________________________________________________________
bool
ConfMan::InitializeParameterFiles( void )
{
  return
    ( InitializeParameter<DCGeomMan>("DCGEO")        &&
      InitializeParameter<HodoParamMan>("HDPRM")     );
      // InitializeParameter<HodoParamMan>("HDPRM")     &&
      // InitializeParameter<UserParamMan>("USER")      );
}

//______________________________________________________________________________
bool
ConfMan::FinalizeProcess( void )
{
  return true;
}
