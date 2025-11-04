/**
 *  file: UserHBX.cc
 *  Ref. file: UserSkeleton.cc
 *  date: 2025.10.14
 *  author imamoto
 */

#include <iostream>
#include <sstream>
#include <cmath>
#include <iomanip>
#include <string>
#include <vector>

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
#include "Utility.hh"
#include "UserParamMan.hh"
#include "HodoParamMan.hh"
#include "HodoParamMan.hh"
#include "DCGeomMan.hh"

namespace
{
  using namespace root;
  using namespace trigger;
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
  int unixtime;

  //Event Sync Clock
  Int_t hrtdc[NumOfChHRTDC][MaxDepth];
  Int_t syncclock[MaxDepth];

  //HBXX Trig flag
  int HBXtrignhits;
  int HBXtrigpat[NumOfSegTrig];
  int HBXtrigflag[NumOfSegTrig];

  //Trig flag
  int trignhits;
  int trigpat[NumOfSegTrig];
  int trigflag[NumOfSegTrig];

  //Ge Reset
  double gereset[NumOfSegGe];

  //Ge ADC
  double geadc[NumOfSegGe];

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

  //Ge Scaler
  int scaler1[NumOfSegScaler];
  int scaler2[NumOfSegScaler];

  void clear();
};

//______________________________________________________________________________
void Event::clear()
{
  evnum      = 0;
  spill      = 0;
  HBXtrignhits  = 0;
  trignhits  = 0;
  genhits    = 0;
  gecrmnhits = 0;
  gebgonhits = 0;
  rm1evnum   = 0;
  rm2evnum   = 0;
  rm1spill   = 0;
  rm2spill   = 0;
  unixtime   = 0;
  for(int i=0; i<NumOfChHRTDC; ++i){
    for(int j=0; j<MaxDepth; ++j){
      hrtdc[i][j]      = -999;
    }
  }
  for(int i=0; i<MaxDepth; ++i){
    syncclock[i] = -999;
  }
  for(int i=0; i<NumOfSegScaler; ++i){
    scaler1[i] =-1;
    scaler2[i] =-1;
  }
  for(int i=0; i<NumOfSegTrig; ++i){
    HBXtrigpat[i]  = -1;
    HBXtrigflag[i] = -1;
    trigpat[i]  = -1;
    trigflag[i] = -1;
  }
  for(int i=0; i<NumOfSegGe; ++i){
    gereset[i] = qnan;
    geadc[i]   = qnan;
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

  //HBX Trig flag
  int HBXtrignhits;
  int HBXtrigpat[NumOfSegTrig];
  int HBXtrigflag[NumOfSegTrig];

  //Trig flag
  int trignhits;
  int trigpat[NumOfSegTrig];
  int trigflag[NumOfSegTrig];

  //Ge Reset
  double geReset[NumOfSegGe];
  //Ge ADC
  double geAdc[NumOfSegGe];

  //Ge TFA
  int nhGeTfa;
  double geTfa[NumOfSegGe][MaxDepth];

  //Ge CRM
  int nhGeCrm;
  double geCrm[NumOfSegGe][MaxDepth];

  //Ge BGO
  int nhBgo;
  double bgoTdc[NumOfSegBGO][MaxDepth];

  void clear();
};

//______________________________________________________________________________
void
Dst::clear()
{
  evnum     = 0;
  spill     = 0;
  HBXtrignhits = 0;
  trignhits = 0;
  nhGeTfa   = 0;
  nhGeCrm   = 0;
  nhBgo     = 0;
  rm1spill  = 0;
  rm2spill  = 0;
  rm1evnum  = 0;
  rm2evnum  = 0;

  for(int i=0; i<NumOfChHRTDC; ++i){
    for(int j=0; j<MaxDepth; ++j){
      hrtdc[i][j] = -999;
    }
  }
  for(int i=0; i<MaxDepth; ++i){
    syncclock[i] = -999;
  }
  for(int i=0; i<NumOfSegTrig; ++i){
    HBXtrigpat[i]=-1;
    HBXtrigflag[i]=-1;
    trigpat[i]=-1;
    trigflag[i]=-1;
  }
  for(int i=0; i<NumOfSegGe; ++i){
    geReset[i]=qnan;
    geAdc[i]=qnan;
    for(int m=0; m<MaxDepth; ++m){
      geTfa[i][m]=qnan;
      geCrm[i][m]=qnan;
    }
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
    GeHidraw    = 100000,
    GeHidcalib    = 200000,
    TwodimHid = 300000,
    TDCHid = 400000,
    BGOHid = 500000,
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
  static const double MinTFACUT = gUser.GetParameter("TCut", 0);
  static const double MaxTFACUT = gUser.GetParameter("TCut", 1);

  static const double MinCRMCUT = gUser.GetParameter("CCut", 0);
  static const double MaxCRMCUT = gUser.GetParameter("CCut", 1);

  static const int MinBGOCUT = gUser.GetParameter("BCut", 0);
  static const int MaxBGOCUT = gUser.GetParameter("BCut", 1);

  static const int MinRESETCUT = gUser.GetParameter("RCut", 0);
  static const int MaxRESETCUT = gUser.GetParameter("RCut", 1);

  static double Calib_slo[NumOfSegGe] = {0};
  static double Calib_int[NumOfSegGe] = {0};


  for(int i=0; i<NumOfSegGe; i++){
    Calib_slo[i] = gUser.GetParameter(Form("CALIB_%d",i), 0);
    Calib_int[i] = gUser.GetParameter(Form("CALIB_%d",i), 1);
  }

  //flag condition

  struct flagsel {
    std::function<bool(std::vector<bool> &)> condition;
    std::function<void(Int_t,int,double,double)> action;
  };

  std::vector<flagsel> flagsels = {
    {[](std::vector<bool> &v){return  true;},
    [](Int_t type, int seg, double adc, double calib){
      HF1(GeHidraw + type*100 + seg, adc);
      HF1(GeHidcalib + type*100 + seg, calib);}},

    {[](std::vector<bool> &v){return  v[kTrigCPS];},
    [](Int_t type, int seg, double adc, double calib){
      HF1(GeHidraw + FlagIdKBeamTOF + type*100 + seg, adc);
      HF1(GeHidcalib + FlagIdKBeamTOF + type*100 + seg, calib);}},

    {[](std::vector<bool> &v){return  v[kTrigEPS];},
     [](Int_t type, int seg, double adc, double calib){
       HF1(GeHidraw + FlagIdLSOxGe + type*100 + seg, adc);
       HF1(GeHidcalib + FlagIdLSOxGe + type*100 + seg, calib);}},

    {[](std::vector<bool> &v){return  v[kTrigFPS];},
     [](Int_t type, int seg, double adc, double calib){
       HF1(GeHidraw + FlagIdGeself + type*100 + seg, adc);
       HF1(GeHidcalib + FlagIdGeself + type*100 + seg, calib);}},

    {[](std::vector<bool> &v){return  v[kTrigEPS] && v[kL1SpillOn];},
     [](Int_t type, int seg, double adc, double calib){
       HF1(GeHidraw + FlagIdLSOxGexSpillOn + type*100 + seg, adc);
       HF1(GeHidcalib + FlagIdLSOxGexSpillOn + type*100 + seg, calib);}},

    {[](std::vector<bool> &v){return  v[kTrigEPS] && v[kL1SpillOff];},
     [](Int_t type, int seg, double adc, double calib){
       HF1(GeHidraw + FlagIdLSOxGexSpillOff + type*100 + seg, adc);
       HF1(GeHidcalib + FlagIdLSOxGexSpillOff + type*100 + seg, calib);}}
  };


  RawData rawData;

  //  gRM.Decode();

  // event.runnum = gRM.RunNumber();
  // event.evnum  = gRM.EventNumber();
  // event.spill  = gRM.SpillNumber();
  // dst.evnum    = gRM.EventNumber();
  // dst.spill    = gRM.SpillNumber();
  // std::cout << "evnum: " << event.evnum << std::endl;
  // std::cout << "spill num: " << event.spill << std::endl;

  rawData.DecodeHits("TFlag");
  rawData.DecodeHits("HBXTFlag");
  std::bitset<NumOfSegTrig> trigger_flag;
  std::bitset<NumOfSegTrig> HBXtrigger_flag;

  std::vector<bool> tflag(NTriggerFlag,false);

  for(const auto& hit: rawData.GetHodoRawHitContainer("TFlag")){
    Int_t seg = hit->SegmentId();
    Int_t tdc = hit->GetTdc();
    if( tdc>0 ){
      event.trigpat[trigger_flag.count()] = seg;
      event.trigflag[seg] = tdc;
      dst.trigpat[trigger_flag.count()] = seg;
      dst.trigflag[seg] = tdc;
      trigger_flag.set(seg);
      HF1(10, seg-1);
      HF1(10+seg, tdc);
      tflag[seg] = true;
    }
  }

  for(const auto& hit: rawData.GetHodoRawHitContainer("HBXTFlag")){
    Int_t seg = hit->SegmentId();
    Int_t tdc = hit->GetTdc();
    if( tdc>0 ){
      event.HBXtrigpat[HBXtrigger_flag.count()] = seg;
      event.HBXtrigflag[seg] = tdc;
      dst.HBXtrigpat[HBXtrigger_flag.count()] = seg;
      dst.HBXtrigflag[seg] = tdc;
      HBXtrigger_flag.set(seg);
      HF1(100, seg-1);
      HF1(100+seg, tdc);
    }
  }

    int nhit_trig[NumOfSegTrig] = {0};
    int trigtdc[NumOfSegTrig] = {0};

    for(int seg=0; seg<NumOfSegTrig; seg++){
      nhit_trig[seg] = gUnpacker.get_entries( DetIdHBXTrig, 0, seg, 0, 1 );

      if(nhit_trig[seg]>0){
	trigtdc[seg] = gUnpacker.get( DetIdHBXTrig, 0, seg, 0, 1, 0);
      }
    }
    // --unixtime--
    int time = 0;
    time = utility::UnixTime();
    event.unixtime = time;
    //

  //---Ge ADC & TDC-----------------------------------------------------
  for( int seg=0; seg<NumOfSegGe; ++seg ){
    int nhit_a = gUnpacker.get_entries( DetIdGe, 0, seg, 0, 0 ); //adc
    int nhit_t = gUnpacker.get_entries( DetIdGe, 0, seg, 0, 1 ); //tfa_leading
    int nhit_r = gUnpacker.get_entries( DetIdGe, 0, seg, 0, 7 ); //reset_time
    int nhit_b[BtoG[seg].size()] = {};

    //adc
    if( nhit_a>0 ){
      double adc = gUnpacker.get( DetIdGe, 0, seg, 0, 0, 0);
      double calib = adc*Calib_slo[seg] + Calib_int[seg]; // a*x + b

      event.geadc[seg] = adc;
      dst.geAdc[seg] = adc;

       for(auto& f : flagsels){
	if (f.condition(tflag)) f.action(DataIdAdc,seg,adc,calib);
      }

      bool tfaflag = false;
      bool tfacutflag = false;
      //adc w/ tfa
      if( nhit_t>0){
	for(int i=0; i<nhit_t; i++){
	  int tfa = gUnpacker.get(DetIdGe, 0, seg, 0, 1, i);
	  HF2(TwodimHid+100*DataIdTfa + seg, adc, tfa);
	  if(tfa>0 )tfaflag = true;
	  if(MinTFACUT<tfa && tfa<MaxTFACUT )tfacutflag = true;
	}

	if(tfaflag){
	  for(auto& f : flagsels){
	    if (f.condition(tflag)) f.action(DataIdAdcwtfa,seg,adc,calib);
	  }
	}
	if(tfacutflag){
	  for(auto& f : flagsels){
	    if (f.condition(tflag)) f.action(DataIdAdcwtfacut,seg,adc,calib);
	  }
	}
      }

      //adc w/ bgo cut

      bool bgohitflag = false;
      bool bgocutflag = false;
      int crystal = 0;
      for(auto& v: BtoG[seg]){
	nhit_b[crystal] = gUnpacker.get_entries( DetIdBGO, 0, BtoG[seg][crystal], 0, 0 );
	if(nhit_b[crystal]>0){
      	  for(int j=0; j<nhit_b[crystal]; j++){
      	    int bgotdc = gUnpacker.get( DetIdBGO, 0, BtoG[seg][crystal], 0, 0, j );
      	    if(bgotdc>MinBGOCUT && bgotdc<MaxBGOCUT){
      	      bgohitflag = bgohitflag + 1;
      	    }else{
      	      bgohitflag = bgohitflag;
      	    }
      	  }
      	}
      }

      if(bgohitflag==0) bgocutflag = true;

      if(bgocutflag){
	for(auto& f : flagsels){
	  if (f.condition(tflag)) f.action(DataIdAdcwbgocut,seg,adc,calib);
	}
      }

      bool resetcutflag = false;
      //adc w/ reset
      if( nhit_r>0){
	for(int i=0; i<nhit_r; i++){
	  double reset_time = gUnpacker.get( DetIdGe, 0, seg, 0, 7);
	  HF2(TwodimHid+100*DataIdReset + seg, reset_time, adc);
	  if(MinRESETCUT<reset_time && reset_time<MaxRESETCUT )resetcutflag = true;

	  if(resetcutflag){
	    for(auto& f : flagsels){
	      if (f.condition(tflag)) f.action(DataIdAdcwresetcut,seg,adc,calib);
	    }
	  }

	  if(resetcutflag && bgocutflag && tfacutflag){
	    for(auto& f : flagsels){
	      if (f.condition(tflag)) f.action(DataIdAdcwbgotfaresetcut,seg,adc,calib);
	    }
	  }
	}
      }
    }


      //tfa
      if( nhit_t>0 ){
	for(int i = 0; i<nhit_t; ++i){
	  double tfa  = gUnpacker.get( DetIdGe, 0, seg, 0, 1, i );
	  if( tfa>0 ){
	    HF1( TDCHid+100*DataIdTfa+seg, tfa ); //0 origin
	  }
	  event.getfa[seg][i] = tfa;
	  dst.geTfa[seg][i] = tfa;
	}
      }

      //reset
      if( nhit_r>0 ){
	double reset_time = gUnpacker.get( DetIdGe, 0, seg, 0, 7);
	HF1( TDCHid+100*DataIdReset+seg, reset_time );
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
	double tdc = gUnpacker.get( DetIdBGO, 0, seg, 0, 0, i )  ;
	if(tdc > 0){
	  HF1( BGOHid+100+seg, tdc );
	  HF2( BGOHid+0, seg+0.5, tdc );
	}
	event.gebgot[seg][i] = tdc;
	dst.bgoTdc[seg][i] = tdc;
      }
    }
  }


  for( int seg=0; seg<NumOfSegScaler; ++seg ){
    int nhit = gUnpacker.get_entries( DetIdScaler, 1, 0, seg, 0 );
    if( nhit>0 ){
      int data = gUnpacker.get( DetIdScaler, 1, 0, seg, 0 );
      event.scaler1[seg] = data;
    }
  }

  for( int seg=0; seg<NumOfSegScaler; ++seg ){
    int nhit = gUnpacker.get_entries( DetIdScaler, 0, 0, seg, 0 );
    if( nhit>0 ){
      int data = gUnpacker.get( DetIdScaler, 0, 0, seg, 0 );
      event.scaler2[seg] = data;
    }
  }


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
  event.unixtime = 0;

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
  for( int it=0; it<NumOfSegGe; it++){
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

  //Scaler
  for( int it=0; it<NumOfSegScaler; it++){
    event.scaler1[it] = 0;
    event.scaler2[it] = 0;
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
    dst.geAdc[it] = -9999.;
    for( int m=0; m<MaxDepth; ++m ){
      dst.geTfa[it][m] = -9999.;
      dst.geCrm[it][m] = -9999.;
      }
  }

  for( int it=0; it<NumOfSegBGO; ++it ){
    for( int m=0; m<MaxDepth; ++m ){
      dst.bgoTdc[it][m] = -9999.;
    }
  }
}

//______________________________________________________________________________
namespace
{

  const int    NbinAdc = 8192;
  const double MinAdc  =   0.;
  const double MaxAdc  = 8192.;

  const int    NbinAdcCalib = 1600;
  const double MinAdcCalib  =   0.;
  const double MaxAdcCalib  = 800.;

  const int    NbinGeTdc = 8192;
  const double MinGeTdc  =   0.;
  const double MaxGeTdc  = 8192.;

  const int    NbinBGOTdc = 8192;
  const double MinBGOTdc  =   0.;
  const double MaxBGOTdc  = 8192.;

  const int    NbinGeReset = 4000;
  const double MinGeReset  =   0.;
  const double MaxGeReset  = 20000;


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
  tree->Branch("unixtime", &event.unixtime, "unixtime/I");
  //  tree->Branch("hrtdc",  &event.hrtdc,  "hrtdc/I");
  tree->Branch("spill",  &event.spill,  "spill/I");
  tree->Branch("hrtdc",  event.hrtdc,  Form("hrtdc[%d][%d]/I", NumOfChHRTDC, MaxDepth));
  tree->Branch("syncclock",  event.syncclock,  Form("syncclock[%d]/I", MaxDepth));
  //Trig
  tree->Branch("trignhits", &event.trignhits, "trignhits/I");

  HB1(  1, "Status", 20, 0., 20. );
  HB1(10, "Trigger HitPat", NumOfSegTrig, 0., Double_t(NumOfSegTrig));
  for(Int_t i=0; i<NumOfSegTrig; ++i){
    HB1(10+i+1, Form("Trigger Flag %d", i+1), 0x1000, 0, 0x1000);
  }
  HB1( 100, "HBXX Trigger HitPat", NumOfSegTrig, 0., Double_t(NumOfSegTrig) );
  for(Int_t i=0; i<NumOfSegTrig; ++i){
    HB1( 100+i+1, Form("HBXX Trigger Trig %d", i+1), 0x1000, 0, 0x1000 );
  }


  //---Ge ADC & TDC-----------------------------------------------------
  for( int i=1; i<=NumOfSegGe; ++i ){

    TString title0_0  = Form("Ge-%d Adc", i);
    TString title1_0  = Form("Ge-%d Adc(w Tfa)", i);
    TString title2_0  = Form("Ge-%d Adc(Tfa Cut)", i);
    TString title3_0  = Form("Ge-%d Adc(Bgo Cut)", i);
    TString title4_0  = Form("Ge-%d Adc(Reset Cut)", i);
    TString title5_0  = Form("Ge-%d Adc(Tfa BGO Reset Cut)", i);

    TString title0_1  = Form("Ge-%d Adc KBeamxTOF ", i);
    TString title1_1  = Form("Ge-%d Adc(w Tfa KBeamxTOF)", i);
    TString title2_1  = Form("Ge-%d Adc(Tfa Cut KBeamxTOF)", i);
    TString title3_1  = Form("Ge-%d Adc(Bgo Cut KBeamxTOF)", i);
    TString title4_1  = Form("Ge-%d Adc(Reset Cut KBeamxTOF)", i);
    TString title5_1  = Form("Ge-%d Adc(Tfa BGO Reset CutKBeamxTOF)", i);

    TString title0_2  = Form("Ge-%d Adc LSOxGe", i);
    TString title1_2  = Form("Ge-%d Adc(w Tfa LSOxGe)", i);
    TString title2_2  = Form("Ge-%d Adc(Tfa Cut LSOxGe)", i);
    TString title3_2  = Form("Ge-%d Adc(Bgo Cut LSOxGe)", i);
    TString title4_2  = Form("Ge-%d Adc(Reset Cut LSOxGe)", i);
    TString title5_2  = Form("Ge-%d Adc(Tfa BGO Reset Cut LSOxGe)", i);

    TString title0_3  = Form("Ge-%d Adc Ge self", i);
    TString title1_3  = Form("Ge-%d Adc(w Tfa Ge self)", i);
    TString title2_3  = Form("Ge-%d Adc(Tfa Cut Ge self)", i);
    TString title3_3  = Form("Ge-%d Adc(Bgo Cut Ge self)", i);
    TString title4_3  = Form("Ge-%d Adc(Reset Cut Ge self)", i);
    TString title5_3  = Form("Ge-%d Adc(Tfa BGO Reset Cut Ge self)", i);

    TString title0_4  = Form("Ge-%d Adc LSOxGexSpillOn", i);
    TString title1_4  = Form("Ge-%d Adc(w Tfa LSOxGexSpillOn)", i);
    TString title2_4  = Form("Ge-%d Adc(Tfa Cut LSOxGexSpillOn)", i);
    TString title3_4  = Form("Ge-%d Adc(Bgo Cut LSOxGexSpillOn)", i);
    TString title4_4  = Form("Ge-%d Adc(Reset Cut LSOxGexSpillOn)", i);
    TString title5_4  = Form("Ge-%d Adc(Tfa BGO Reset Cut LSOxGexSpillOn)", i);

    TString title0_5  = Form("Ge-%d Adc LSOxGexSpillOff", i);
    TString title1_5  = Form("Ge-%d Adc(w Tfa LSOxGexSpillOff)", i);
    TString title2_5  = Form("Ge-%d Adc(Tfa Cut LSOxGexSpillOff)", i);
    TString title3_5  = Form("Ge-%d Adc(Bgo Cut LSOxGexSpillOff)", i);
    TString title4_5  = Form("Ge-%d Adc(Reset Cut LSOxGexSpillOff)", i);
    TString title5_5  = Form("Ge-%d Adc(Tfa BGO Reset Cut LSOxGexSpillOff)", i);

    TString title0_6  = Form("Ge-%d Calib", i);
    TString title1_6  = Form("Ge-%d Calib(w Tfa)", i);
    TString title2_6  = Form("Ge-%d Calib(Tfa Cut)", i);
    TString title3_6  = Form("Ge-%d Calib(Bgo Cut)", i);
    TString title4_6  = Form("Ge-%d Calib(Reset Cut)", i);
    TString title5_6  = Form("Ge-%d Calib(Tfa BGO Reset Cut)", i);

    TString title0_7  = Form("Ge-%d Calib KBeamxTOF ", i);
    TString title1_7  = Form("Ge-%d Calib(w Tfa KBeamxTOF)", i);
    TString title2_7  = Form("Ge-%d Calib(Tfa Cut KBeamxTOF)", i);
    TString title3_7  = Form("Ge-%d Calib(Bgo Cut KBeamxTOF)", i);
    TString title4_7  = Form("Ge-%d Calib(Reset Cut KBeamxTOF)", i);
    TString title5_7  = Form("Ge-%d Calib(Tfa BGO Reset CutKBeamxTOF)", i);

    TString title0_8  = Form("Ge-%d Calib LSOxGe", i);
    TString title1_8  = Form("Ge-%d Calib(w Tfa LSOxGe)", i);
    TString title2_8  = Form("Ge-%d Calib(Tfa Cut LSOxGe)", i);
    TString title3_8  = Form("Ge-%d Calib(Bgo Cut LSOxGe)", i);
    TString title4_8  = Form("Ge-%d Calib(Reset Cut LSOxGe)", i);
    TString title5_8  = Form("Ge-%d Calib(Tfa BGO Reset Cut LSOxGe)", i);

    TString title0_9  = Form("Ge-%d Calib Ge self", i);
    TString title1_9  = Form("Ge-%d Calib(w Tfa Ge self)", i);
    TString title2_9  = Form("Ge-%d Calib(Tfa Cut Ge self)", i);
    TString title3_9  = Form("Ge-%d Calib(Bgo Cut Ge self)", i);
    TString title4_9  = Form("Ge-%d Calib(Reset Cut Ge self)", i);
    TString title5_9  = Form("Ge-%d Calib(Tfa BGO Reset Cut Ge self)", i);

    TString title0_10  = Form("Ge-%d Calib LSOxGexSpillOn", i);
    TString title1_10  = Form("Ge-%d Calib(w Tfa LSOxGexSpillOn)", i);
    TString title2_10  = Form("Ge-%d Calib(Tfa Cut LSOxGexSpillOn)", i);
    TString title3_10  = Form("Ge-%d Calib(Bgo Cut LSOxGexSpillOn)", i);
    TString title4_10  = Form("Ge-%d Calib(Reset Cut LSOxGexSpillOn)", i);
    TString title5_10  = Form("Ge-%d Calib(Tfa BGO Reset Cut LSOxGexSpillOn)", i);

    TString title0_11  = Form("Ge-%d Calib LSOxGexSpillOff", i);
    TString title1_11  = Form("Ge-%d Calib(w Tfa LSOxGexSpillOff)", i);
    TString title2_11  = Form("Ge-%d Calib(Tfa Cut LSOxGexSpillOff)", i);
    TString title3_11  = Form("Ge-%d Calib(Bgo Cut LSOxGexSpillOff)", i);
    TString title4_11  = Form("Ge-%d Calib(Reset Cut LSOxGexSpillOff)", i);
    TString title5_11  = Form("Ge-%d Calib(Tfa BGO Reset Cut LSOxGexSpillOff)", i);

    TString title6  = Form("Ge-%d Tfa vs adc", i);
    TString title7  = Form("Ge-%d Reset vs adc", i);

    TString title8  = Form("Ge-%d Tfa", i);
    TString title9  = Form("Ge-%d Crm", i);
    TString title10  = Form("Ge-%d Reset", i);


    int  m = i - 1;

    //no selecting trigflag//

    HB1( GeHidraw + DataIdAdc*100 + m, title0_0, NbinAdc, MinAdc, MaxAdc );
    HB1( GeHidraw + DataIdAdcwtfa*100 + m, title1_0, NbinAdc, MinAdc, MaxAdc );
    HB1( GeHidraw + DataIdAdcwtfacut*100 + m, title2_0, NbinAdc, MinAdc, MaxAdc );
    HB1( GeHidraw + DataIdAdcwbgocut*100 + m, title3_0, NbinAdc, MinAdc, MaxAdc );
    HB1( GeHidraw + DataIdAdcwresetcut*100 + m, title4_0, NbinAdc, MinAdc, MaxAdc );
    HB1( GeHidraw + DataIdAdcwbgotfaresetcut*100 + m, title5_0, NbinAdc, MinAdc, MaxAdc );
    HB1( GeHidcalib + DataIdAdc*100 + m, title0_6, NbinAdc, MinAdc, MaxAdc );
    HB1( GeHidcalib + DataIdAdcwtfa*100 + m, title1_6, NbinAdc, MinAdc, MaxAdc );
    HB1( GeHidcalib + DataIdAdcwtfacut*100 + m, title2_6, NbinAdc, MinAdc, MaxAdc );
    HB1( GeHidcalib + DataIdAdcwbgocut*100 + m, title3_6, NbinAdc, MinAdc, MaxAdc );
    HB1( GeHidcalib + DataIdAdcwresetcut*100 + m, title4_6, NbinAdc, MinAdc, MaxAdc );
    HB1( GeHidcalib + DataIdAdcwbgotfaresetcut*100 + m, title5_6, NbinAdc, MinAdc, MaxAdc );

    HB1( GeHidraw + FlagIdKBeamTOF + m, title0_1, NbinAdc, MinAdc, MaxAdc );
    HB1( GeHidraw + FlagIdKBeamTOF + DataIdAdcwtfa*100 + m, title1_0, NbinAdc, MinAdc, MaxAdc );
    HB1( GeHidraw + FlagIdKBeamTOF + DataIdAdcwtfacut*100 + m, title2_0, NbinAdc, MinAdc, MaxAdc );
    HB1( GeHidraw + FlagIdKBeamTOF + DataIdAdcwbgocut*100 + m, title3_0, NbinAdc, MinAdc, MaxAdc );
    HB1( GeHidraw + FlagIdKBeamTOF + DataIdAdcwresetcut*100 + m, title4_0, NbinAdc, MinAdc, MaxAdc );
    HB1( GeHidraw + FlagIdKBeamTOF + DataIdAdcwbgotfaresetcut*100 + m, title5_0, NbinAdc, MinAdc, MaxAdc );
    HB1( GeHidcalib + FlagIdKBeamTOF + DataIdAdc*100 + m, title0_7, NbinAdc, MinAdc, MaxAdc );
    HB1( GeHidcalib + FlagIdKBeamTOF + DataIdAdcwtfa*100 + m, title1_7, NbinAdc, MinAdc, MaxAdc );
    HB1( GeHidcalib + FlagIdKBeamTOF + DataIdAdcwtfacut*100 + m, title2_7, NbinAdc, MinAdc, MaxAdc );
    HB1( GeHidcalib + FlagIdKBeamTOF + DataIdAdcwbgocut*100 + m, title3_7, NbinAdc, MinAdc, MaxAdc );
    HB1( GeHidcalib + FlagIdKBeamTOF + DataIdAdcwresetcut*100 + m, title4_7, NbinAdc, MinAdc, MaxAdc );
    HB1( GeHidcalib + FlagIdKBeamTOF + DataIdAdcwbgotfaresetcut*100 + m, title5_7, NbinAdc, MinAdc, MaxAdc );

    HB1( GeHidraw + FlagIdLSOxGe + DataIdAdc*100 + m, title0_2, NbinAdc, MinAdc, MaxAdc );
    HB1( GeHidraw + FlagIdLSOxGe + DataIdAdcwtfa*100 + m, title1_2, NbinAdc, MinAdc, MaxAdc );
    HB1( GeHidraw + FlagIdLSOxGe + DataIdAdcwtfacut*100 + m, title2_2, NbinAdc, MinAdc, MaxAdc );
    HB1( GeHidraw + FlagIdLSOxGe + DataIdAdcwbgocut*100 + m, title3_2, NbinAdc, MinAdc, MaxAdc );
    HB1( GeHidraw + FlagIdLSOxGe + DataIdAdcwresetcut*100 + m, title4_2, NbinAdc, MinAdc, MaxAdc );
    HB1( GeHidraw + FlagIdLSOxGe + DataIdAdcwbgotfaresetcut*100 + m, title5_2, NbinAdc, MinAdc, MaxAdc );
    HB1( GeHidcalib + FlagIdLSOxGe + DataIdAdc*100 + m, title0_8, NbinAdc, MinAdc, MaxAdc );
    HB1( GeHidcalib + FlagIdLSOxGe + DataIdAdcwtfa*100 + m, title1_8, NbinAdc, MinAdc, MaxAdc );
    HB1( GeHidcalib + FlagIdLSOxGe + DataIdAdcwtfacut*100 + m, title2_8, NbinAdc, MinAdc, MaxAdc );
    HB1( GeHidcalib + FlagIdLSOxGe + DataIdAdcwbgocut*100 + m, title3_8, NbinAdc, MinAdc, MaxAdc );
    HB1( GeHidcalib + FlagIdLSOxGe + DataIdAdcwresetcut*100 + m, title4_8, NbinAdc, MinAdc, MaxAdc );
    HB1( GeHidcalib + FlagIdLSOxGe + DataIdAdcwbgotfaresetcut*100 + m, title5_8, NbinAdc, MinAdc, MaxAdc );

    HB1( GeHidraw + FlagIdGeself + DataIdAdc*100 + m, title0_3, NbinAdc, MinAdc, MaxAdc );
    HB1( GeHidraw + FlagIdGeself + DataIdAdcwtfa*100 + m, title1_3, NbinAdc, MinAdc, MaxAdc );
    HB1( GeHidraw + FlagIdGeself + DataIdAdcwtfacut*100 + m, title2_3, NbinAdc, MinAdc, MaxAdc );
    HB1( GeHidraw + FlagIdGeself + DataIdAdcwbgocut*100 + m, title3_3, NbinAdc, MinAdc, MaxAdc );
    HB1( GeHidraw + FlagIdGeself + DataIdAdcwresetcut*100 + m, title4_3, NbinAdc, MinAdc, MaxAdc );
    HB1( GeHidraw + FlagIdGeself + DataIdAdcwbgotfaresetcut*100 + m, title5_3, NbinAdc, MinAdc, MaxAdc );
    HB1( GeHidcalib + FlagIdGeself + DataIdAdc*100 + m, title0_9, NbinAdc, MinAdc, MaxAdc );
    HB1( GeHidcalib + FlagIdGeself + DataIdAdcwtfa*100 + m, title1_9, NbinAdc, MinAdc, MaxAdc );
    HB1( GeHidcalib + FlagIdGeself + DataIdAdcwtfacut*100 + m, title2_9, NbinAdc, MinAdc, MaxAdc );
    HB1( GeHidcalib + FlagIdGeself + DataIdAdcwbgocut*100 + m, title3_9, NbinAdc, MinAdc, MaxAdc );
    HB1( GeHidcalib + FlagIdGeself + DataIdAdcwresetcut*100 + m, title4_9, NbinAdc, MinAdc, MaxAdc );
    HB1( GeHidcalib + FlagIdGeself + DataIdAdcwbgotfaresetcut*100 + m, title5_9, NbinAdc, MinAdc, MaxAdc );

    HB1( GeHidraw + FlagIdLSOxGexSpillOn + DataIdAdc*100 + m, title0_4, NbinAdc, MinAdc, MaxAdc );
    HB1( GeHidraw + FlagIdLSOxGexSpillOn + DataIdAdcwtfa*100 + m, title1_4, NbinAdc, MinAdc, MaxAdc );
    HB1( GeHidraw + FlagIdLSOxGexSpillOn + DataIdAdcwtfacut*100 + m, title2_4, NbinAdc, MinAdc, MaxAdc );
    HB1( GeHidraw + FlagIdLSOxGexSpillOn + DataIdAdcwbgocut*100 + m, title3_4, NbinAdc, MinAdc, MaxAdc );
    HB1( GeHidraw + FlagIdLSOxGexSpillOn + DataIdAdcwresetcut*100 + m, title4_4, NbinAdc, MinAdc, MaxAdc );
    HB1( GeHidraw + FlagIdLSOxGexSpillOn + DataIdAdcwbgotfaresetcut*100 + m, title5_4, NbinAdc, MinAdc, MaxAdc );
    HB1( GeHidcalib + FlagIdLSOxGexSpillOn + DataIdAdc*100 + m, title0_10, NbinAdc, MinAdc, MaxAdc );
    HB1( GeHidcalib + FlagIdLSOxGexSpillOn + DataIdAdcwtfa*100 + m, title1_10, NbinAdc, MinAdc, MaxAdc );
    HB1( GeHidcalib + FlagIdLSOxGexSpillOn + DataIdAdcwtfacut*100 + m, title2_10, NbinAdc, MinAdc, MaxAdc );
    HB1( GeHidcalib + FlagIdLSOxGexSpillOn + DataIdAdcwbgocut*100 + m, title3_10, NbinAdc, MinAdc, MaxAdc );
    HB1( GeHidcalib + FlagIdLSOxGexSpillOn + DataIdAdcwresetcut*100 + m, title4_10, NbinAdc, MinAdc, MaxAdc );
    HB1( GeHidcalib + FlagIdLSOxGexSpillOn + DataIdAdcwbgotfaresetcut*100 + m, title5_10, NbinAdc, MinAdc, MaxAdc );

    HB1( GeHidraw + FlagIdLSOxGexSpillOff + DataIdAdc*100 + m, title0_5, NbinAdc, MinAdc, MaxAdc );
    HB1( GeHidraw + FlagIdLSOxGexSpillOff + DataIdAdcwtfa*100 + m, title1_5, NbinAdc, MinAdc, MaxAdc );
    HB1( GeHidraw + FlagIdLSOxGexSpillOff + DataIdAdcwtfacut*100 + m, title2_5, NbinAdc, MinAdc, MaxAdc );
    HB1( GeHidraw + FlagIdLSOxGexSpillOff + DataIdAdcwbgocut*100 + m, title3_5, NbinAdc, MinAdc, MaxAdc );
    HB1( GeHidraw + FlagIdLSOxGexSpillOff + DataIdAdcwresetcut*100 + m, title4_5, NbinAdc, MinAdc, MaxAdc );
    HB1( GeHidraw + FlagIdLSOxGexSpillOff + DataIdAdcwbgotfaresetcut*100 + m, title5_5, NbinAdc, MinAdc, MaxAdc );
    HB1( GeHidcalib + FlagIdLSOxGexSpillOff + DataIdAdc*100 + m, title0_11, NbinAdc, MinAdc, MaxAdc );
    HB1( GeHidcalib + FlagIdLSOxGexSpillOff + DataIdAdcwtfa*100 + m, title1_11, NbinAdc, MinAdc, MaxAdc );
    HB1( GeHidcalib + FlagIdLSOxGexSpillOff + DataIdAdcwtfacut*100 + m, title2_11, NbinAdc, MinAdc, MaxAdc );
    HB1( GeHidcalib + FlagIdLSOxGexSpillOff + DataIdAdcwbgocut*100 + m, title3_11, NbinAdc, MinAdc, MaxAdc );
    HB1( GeHidcalib + FlagIdLSOxGexSpillOff + DataIdAdcwresetcut*100 + m, title4_11, NbinAdc, MinAdc, MaxAdc );
    HB1( GeHidcalib + FlagIdLSOxGexSpillOff + DataIdAdcwbgotfaresetcut*100 + m, title5_11, NbinAdc, MinAdc, MaxAdc );

    HB2(TwodimHid + 100*DataIdTfa + m, title6, NbinAdc, MinAdc, MaxAdc, NbinGeTdc, MinGeTdc, MaxGeTdc);
    HB2(TwodimHid + 100*DataIdReset + m, title7, NbinGeReset, MinGeReset, MaxGeReset, NbinAdc, MinAdc, MaxAdc);

    HB1(TDCHid + 100*DataIdTfa + m,title8 , NbinGeTdc, MinGeTdc, MaxGeTdc);
    HB1(TDCHid + 100*DataIdReset + m, title10, NbinGeReset, MinGeReset, MaxGeReset);

  }

    TString title2_30 = "Ge vs BGO hitpat Spill On LSOxGe";
    TString title3_30 = "Ge vs BGO hitpat Spill Off LSOxGe";
    TString title4_30 = "Ge vs BGO hitpat Spill On L1";

  //Hits
    // HB1( GeHid +30*Offset, "#Hits Ge",        NumOfSegGe+1, 0., double(NumOfSegGe+1) );
    // HB1( GeHid +30*Offset+1, "Hitpat Ge",       NumOfSegGe,   0., double(NumOfSegGe)   );
    // HB2( GeHid +30*Offset+2, title2_30,
    // 	 NumOfSegBGO, 0, NumOfSegBGO, NumOfSegGe, 0, NumOfSegGe);
    // HB2( GeHid +30*Offset+3, title3_30,
    // 	 NumOfSegBGO, 0, NumOfSegBGO, NumOfSegGe, 0, NumOfSegGe);
    // HB2( GeHid +30*Offset+4, title4_30,
    // 	 NumOfSegBGO, 0, NumOfSegBGO, NumOfSegGe, 0, NumOfSegGe);


    // HB2( GeHid +30*Offset+3, "Tfa%Ch",
    // 	 NumOfSegGe, 0, NumOfSegGe, NbinGeTdc, MinGeTdc, MaxGeTdc);
    // HB2( GeHid +30*Offset+4, "Crm%Ch",
    // 	 NumOfSegGe, 0, NumOfSegGe, NbinGeTdc, MinGeTdc, MaxGeTdc);
    // HB2( GeHid +30*Offset+5, "Reset%Ch",
    // 	 NumOfSegGe, 0, NumOfSegGe, NbinGeTdc, MinGeTdc, MaxGeTdc);


  //---BGO--------------------------------------------------------------
  HB1( BGOHid +0, "#Hits BGO",        NumOfSegBGO+1, 0., double(NumOfSegBGO+1) );
  HB1( BGOHid +1, "Hitpat BGO",       NumOfSegBGO,   0., double(NumOfSegBGO)   );

  for( int i=1; i<=NumOfSegBGO; ++i ){
    TString title = Form("BGO-%d Tdc", i);
    int m=i-1;
    HB1( BGOHid+100+m, title, NbinBGOTdc, MinBGOTdc, MaxBGOTdc );
  }



/////tree/////////////////////////////////
  //Flag
  tree->Branch("HBXtrigpat", event.HBXtrigpat, Form("HBXtrigpat[%d]/I", NumOfSegTrig));
  tree->Branch("HBXtrigflag", event.HBXtrigflag, Form("HBXtrigflag[%d]/I", NumOfSegTrig));
  tree->Branch("trigpat", event.trigpat, Form("trigpat[%d]/I", NumOfSegTrig));
  tree->Branch("trigflag", event.trigflag, Form("trigflag[%d]/I", NumOfSegTrig));

  //Ge ADC
  tree->Branch("geadc",   event.geadc,  Form("geadc[%d]/D", NumOfSegGe));

  //Ge TDC
  tree->Branch("genhits",   &event.genhits,   "genhits/I");
  tree->Branch("gehitpat",   event.gehitpat,  Form("gehitpat[%d]/I", NumOfSegGe));
  tree->Branch("getfa",        event.getfa,       Form("getfa[%d][%d]/D", NumOfSegGe, MaxDepth));
  tree->Branch("gecrm",        event.gecrm,       Form("gecrm[%d][%d]/D", NumOfSegGe, MaxDepth));
  tree->Branch("gereset", event.gereset, Form("gereset[%d]/I", NumOfSegGe));

  //BGO
  tree->Branch("gebgonhits",   &event.gebgonhits,   "gebgonhits/I");
  tree->Branch("gebgohitpat",   event.gebgohitpat,  Form("gebgohitpat[%d]/I", NumOfSegBGO));
  tree->Branch("gebgot",        event.gebgot,       Form("gebgot[%d][%d]/D", NumOfSegBGO, MaxDepth));

  //Ge Scaler
  tree->Branch("scaler1", event.scaler1, Form("scaler1[%d]/I", NumOfSegScaler));
  tree->Branch("scaler2", event.scaler2, Form("scaler2[%d]/I", NumOfSegScaler));

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
  hbx->Branch("geAdc",   dst.geAdc,  Form("geAdc[%d]/D", NumOfSegGe));
  hbx->Branch("geAdc",   dst.geAdc,  Form("geAdc[%d]/D", NumOfSegGe));

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
Bool_t
ConfMan::InitializeParameterFiles( void )
{
  return
    (InitializeParameter<DCGeomMan>("DCGEO")        &&
     InitializeParameter<HodoParamMan>("HDPRM")     &&
     InitializeParameter<UserParamMan>("USER"));
}

//______________________________________________________________________________
bool
ConfMan::FinalizeProcess( void )
{
  return true;
}
