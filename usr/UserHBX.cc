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
  rm1evnum   = 0;
  rm2evnum   = 0;
  rm1spill   = 0;
  rm2spill   = 0;
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

  for( int it=0; it<NumOfChHRTDC; ++it){
    for (int m=0; m<MaxDepth; ++m){
      dst.hrtdc[it][m] = -9999;
    }
  }
  for (int m=0; m<MaxDepth; ++m){
    dst.syncclock[m] = -9999;
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


  /////tree/////////////////////////////////
  //Flag
  tree->Branch("trigpat", event.trigpat, Form("trigpat[%d]/I", NumOfSegTrig));
  tree->Branch("trigflag", event.trigflag, Form("trigflag[%d]/I", NumOfSegTrig));

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
