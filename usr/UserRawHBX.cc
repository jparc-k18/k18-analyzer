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
#include "DetectorID.hh"
#include "RMAnalyzer.hh"
#include "MathTools.hh"
#include "RawData.hh"
#include "HodoRawHit.hh"
#include "RootHelper.hh"
#include "UnpackerManager.hh"
#include "VEvent.hh"
#include "Unpacker.hh"
#include "UserParamMan.hh"

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
  runnum     = 0;
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
  int runnum;
  int evnum;
  int spill;
  int rm1spill;
  int rm2spill;
  int rm1evnum;
  int rm2evnum;

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
  runnum    = 0;
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
  TTree *maindaq;

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
  static const double MinTFACUT = 0;
  static const double MaxTFACUT = 100;
  static const double MinCRMCUT = 0;
  static const double MaxCRMCUT = 100;

  RawData rawData;

  //  gRM.Decode();
  // //  event.runnum = gRM.RunNumber();
  // event.evnum  = gRM.EventNumber();
  // event.spill  = gRM.SpillNumber();
  // //  dst.runnum   = gRM.RunNumber();
  // dst.evnum    = gRM.EventNumber();
  // dst.spill    = gRM.SpillNumber();

  rawData.DecodeHits("HBXTFlag");
  std::bitset<NumOfSegTrig> trigger_flag;
  for(const auto& hit: rawData.GetHodoRawHitContainer("HBXTFlag")){
    Int_t seg = hit->SegmentId();
    Int_t tdc = hit->GetTdc();
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
  

  //  ---Ge ADC & TDC-----------------------------------------------------
  // for( int seg=0; seg<NumOfSegGe*2; ++seg ){
  //   int nhit_a = gUnpacker.get_entries( DetIdGe, 0, seg, 0, 0 ); //adc
  //   int nhit_t = gUnpacker.get_entries( DetIdGe, 0, seg, 0, 1 ); //tfa_leading
  //   int nhit_c = gUnpacker.get_entries( DetIdGe, 0, seg, 0, 4 ); //crm_leading
  //   int nhit_r = gUnpacker.get_entries( DetIdGe, 0, seg, 0, 7 ); //reset_time
  //   if( nhit_a>0 ){
  //     int adc = gUnpacker.get( DetIdGe, 0, seg, 0, 0, 0);
  //     event.geadc[seg] = adc;
  //     dst.geAdc[seg] = adc;
  //     HF1( GeHid+100*(seg+1)+10, double(adc) );
  //     HF2( GeHid+100 +0, seg+0.5, double(adc) );
  //     //adc w/ tfa
  //     if( nhit_t>0){
  // 	int tfa  = gUnpacker.get( DetIdGe, 0, seg, 0, 1, 0);
  // 	if( tfa>0 ){
  // 	  HF1( GeHid+100*(seg+1)+21, double(adc) );
  // 	  HF2( GeHid+100*(seg+1)+24, double(tfa), double(adc) );
  // 	}
  // 	else{
  // 	  HF1( GeHid+100*(seg+1)+22, double(adc) );
  // 	}
  // 	//adc w/ tfa cut
  // 	if( MinTFACUT<tfa && tfa<MaxTFACUT ) HF1( GeHid+100*(seg+1)+23, double(adc) );
  //     }
  //     //adc w/ crm
  //     if( nhit_t>0){
  // 	int crm  = gUnpacker.get( DetIdGe, 0, seg, 0, 1, 0);
  // 	if( crm>0 ){
  // 	  HF1( GeHid+100*(seg+1)+31, double(adc) );
  // 	  HF2( GeHid+100*(seg+1)+34, double(crm), double(adc) );
  // 	}
  // 	else{
  // 	  HF1( GeHid+100*(seg+1)+32, double(adc) );
  // 	}
  // 	//adc w/ crm cut
  // 	if( MinCRMCUT<crm && crm<MaxCRMCUT ) HF1( GeHid+100*(seg+1)+33, double(adc) );
  //     }
  //     //adc w/ reset
  //     if( nhit_r>0){
  // 	int reset_time = gUnpacker.get( DetIdGe, 0, seg, 0, 7);
  // 	HF2( GeHid+100*(seg+1)+44, double(reset_time), double(adc) );
  //     }
  //   }
  //   //tfa
  //   if( nhit_t>0 ){
  //     HF1( GeHid+1, seg+0.5); //0 origin
  //     for(int i = 0; i<nhit_t; ++i){
  // 	int tfa  = gUnpacker.get( DetIdGe, 0, seg, 0, 1, i );
  // 	if( tfa>0 ){
  // 	  HF1( GeHid+100*(seg+1)+20, double(tfa) ); //0 origin
  // 	  HF2( GeHid+100+2, seg+0.5, double(tfa) );
  // 	}
  // 	event.getfa[seg][i] = tfa;
  // 	dst.geTfa[seg][i] = tfa;
  //     }
  //   }
  //   //crm
  //   if( nhit_c>0 ){
  //     for(int i = 0; i<nhit_c; ++i){
  // 	int crm  = gUnpacker.get( DetIdGe, 0, seg, 0, 4, i );
  // 	if(crm > 0){
  // 	  HF1( GeHid+100*(seg+1)+30, double(crm) ); //0 origin
  // 	  HF2( GeHid+100+3, seg+0.5, double(crm) );
  // 	  event.gecrm[seg][i] = crm;
  // 	  dst.geCrm[seg][i] = crm;
  // 	}
  //     }
  //   }
  //   //reset
  //   if( nhit_r>0 ){
  //     int reset_time = gUnpacker.get( DetIdGe, 0, seg, 0, 7);
  //     HF1( GeHid+100*(seg+1)+40, double(reset_time) );
  //     HF2( GeHid+100+4, seg+0.5, double(reset_time) );
  //     event.gereset[seg] = reset_time;
  //     dst.geReset[seg] = reset_time;
  //   }
  // }
  // //---BGO--------------------------------------------------------------
  // for(int seg = 0; seg<NumOfSegBGO; ++seg){
  //   int nhit = gUnpacker.get_entries( DetIdBGO, 0, seg, 0, 0 );
  //   if( nhit>0 ){
  //     HF1( BGOHid+1, seg+0.5); //0 origin
  //     for(int i = 0; i<nhit; ++i){
  // 	int tdc = gUnpacker.get( DetIdBGO, 0, seg, 0, 0, i )  ;
  // 	if(tdc > 0){
  // 	  HF1( BGOHid+100*(seg+1)+0, double(tdc) );
  // 	  HF2( BGOHid+0, seg+0.5, double(tdc) );
  // 	}
  // 	event.gebgot[seg][i] = tdc;
  // 	dst.bgoTdc[seg][i] = tdc;
  //     }
  //   }
  // }
  // ---HUL-RM--------------------------------------------------------
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
  maindaq->Fill();
  return true;
}
//______________________________________________________________________________
void
InitializeEvent( void )
{
  event.runnum = 0;
  event.evnum = 0;
  event.spill = 0;
  event.rm1evnum = 0;
  event.rm2evnum = 0;
  event.rm1spill = 0;
  event.rm2spill = 0;
  event.gebgonhits  = 0;
  event.trignhits = 0;

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
  for(int i=0; i<NumOfPlaneHulRm; ++i){
    event.evnum_hulrm[i] = -1;
    event.spill_hulrm[i] = -1;
  }


  ////Dst////////////////////////
  dst.runnum = 0;
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
  for(int i=0; i<NumOfPlaneHulRm; ++i){
    dst.evnum_hulrm[i] = -1;
    dst.spill_hulrm[i] = -1;
  }
}

//______________________________________________________________________________
namespace
{

  const int    NbinAdc = 8192;
  const double MinAdc  =   0.;
  const double MaxAdc  = 8192.;

  const int    NbinGeTdc = 16384;
  const double MinGeTdc  =   0.;
  const double MaxGeTdc  = 16384.;

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
  tree->Branch("spill",  &event.spill,  "spill/I");
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

    if(1<=i&&i<=16){
    HB1( GeHid +100*i +10, title10, NbinAdc, MinAdc, MaxAdc );
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

    HB1( GeHid +100*i +30, title30, NbinGeTdc, MinGeTdc, MaxGeTdc );
    HB1( GeHid +100*i +31, title31, NbinAdc, MinAdc, MaxAdc );
    HB1( GeHid +100*i +32, title32, NbinAdc, MinAdc, MaxAdc );
    HB1( GeHid +100*i +33, title33, NbinAdc, MinAdc, MaxAdc );
    HB2( GeHid +100*i +34, title34,
    	 NbinGeTdc/16, MinGeTdc, MaxGeTdc, NbinAdc/8, MinAdc, MaxAdc);

    HB1( GeHid +100*i +40, title40, NbinGeTdc, MinGeTdc, MaxGeTdc );
    HB2( GeHid +100*i +44, title44,
    	 NbinGeTdc/16, MinGeTdc, MaxGeTdc, NbinAdc/8, MinAdc, MaxAdc);

  }

  //Hits
  HB1( GeHid +0, "#Hits Ge",        NumOfSegGe+1, 0., double(NumOfSegGe+1) );
  HB1( GeHid +1, "Hitpat Ge",       NumOfSegGe,   0., double(NumOfSegGe)   );
  HB2( GeHid +100 +1, "Adc%Ch",
       NumOfSegGe, 0, NumOfSegGe, NbinAdc, MinAdc, MaxAdc);
  HB2( GeHid +100 +2, "Tfa%Ch",
       NumOfSegGe, 0, NumOfSegGe, NbinGeTdc, MinGeTdc, MaxGeTdc);
  HB2( GeHid +100 +3, "Crm%Ch",
       NumOfSegGe, 0, NumOfSegGe, NbinGeTdc, MinGeTdc, MaxGeTdc);
  HB2( GeHid +100 +4, "Reset%Ch",
       NumOfSegGe, 0, NumOfSegGe, NbinGeTdc, MinGeTdc, MaxGeTdc);


  //---BGO--------------------------------------------------------------
  HB1( BGOHid +0, "#Hits BGO",        NumOfSegBGO+1, 0., double(NumOfSegBGO+1) );
  HB1( BGOHid +1, "Hitpat BGO",       NumOfSegBGO,   0., double(NumOfSegBGO)   );

  for( int i=1; i<=NumOfSegBGO; ++i ){
    TString title = Form("BGO-%d Tdc", i);
    HB1( BGOHid +100*i +0, title, NbinBGOTdc, MinBGOTdc, MaxBGOTdc );
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
  hbx->Branch("runnum", &dst.runnum, "runnum/I");
  hbx->Branch("evnum", &dst.evnum, "evnum/I");
  hbx->Branch("spill", &dst.spill, "spill/I");
  hbx->Branch("rm1evnum", &dst.rm1evnum, "rm1evnum/I");
  hbx->Branch("rm2evnum", &dst.rm2evnum, "rm2evnum/I");
  hbx->Branch("rm1spill", &dst.rm1spill, "rm1spill/I");
  hbx->Branch("rm2spill", &dst.rm2spill, "rm2spill/I");
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

  maindaq = new TTree( "maindaq", "Data Summary Table of maindata" );
  maindaq->Branch("runnum", &dst.runnum, "runnum/I");
  maindaq->Branch("evnum",  &dst.evnum,  "evnum/I");
  maindaq->Branch("spill",  &dst.spill,  "spill/I");
  maindaq->Branch("scaler",  dst.scaler,  Form("scaler[%d][%d]", NumOfPlaneScaler, NumOfSegScaler));
  maindaq->Branch("evnum_hulrm", dst.evnum_hulrm, Form("evnum_hulrm[%d]", NumOfPlaneHulRm));
  maindaq->Branch("spill_hulrm", dst.evnum_hulrm, Form("evnum_hulrm[%d]", NumOfPlaneHulRm));

  HPrint();
  return true;
}

//______________________________________________________________________________
bool
ConfMan::InitializeParameterFiles( void )
{
  return
    ( InitializeParameter<UserParamMan>("USER") );
}

//______________________________________________________________________________
bool
ConfMan::FinalizeProcess( void )
{
  return true;
}
