/**
 *  file: UseHBC.cc
 *  Ref. file: UserSkeleton.cc
 *  date: 2020.12.09
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
#include "KuramaLib.hh"
#include "MathTools.hh"
#include "RawData.hh"
#include "RootHelper.hh"
#include "UnpackerManager.hh"
#include "VEvent.hh"

namespace
{
  using namespace root;
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
VEvent::VEvent( void )
{
}

//______________________________________________________________________________
VEvent::~VEvent( void )
{
}

//______________________________________________________________________________
class EventHBX : public VEvent
{
private:
  RawData *rawData;

public:
        EventHBX( void );
       ~EventHBX( void );
  bool  ProcessingBegin( void );
  bool  ProcessingEnd( void );
  bool  ProcessingNormal( void );
  bool  InitializeHistograms( void );
  void  InitializeEvent( void );
};

//______________________________________________________________________________
EventHBX::EventHBX( void )
  : VEvent(),
    rawData(0)
{
}

//______________________________________________________________________________
EventHBX::~EventHBX( void )
{
  if (rawData) delete rawData;
}

//______________________________________________________________________________
struct Event
{
  int runnum;
  int evnum;
  int spill;

  //Trig flag
  int trignhits;  
  int trigpat[NumOfSegTrig];
  int trigflag[NumOfSegTrig];
  int hbxtrignhits;  
  int hbxtrigpat[NumOfSegHbxTrig];
  int hbxtrigflag[NumOfSegHbxTrig];

  //Ge Reset
  int gereset[NumOfSegGe];

  //Ge ADC
  double geadc[NumOfSegGe];

  double lsogeadc[NumOfSegGe];
  double coingeadc[NumOfSegGe];
  double spillongeadc[NumOfSegGe];
  double spilloffgeadc[NumOfSegGe];

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
  int scaler[NumOfSegScaler];
};
//______________________________________________________________________________
struct Dst
{
  int evnum;
  int spill;
};

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
EventHBX::ProcessingBegin( void )
{
  InitializeEvent();
  return true;
}

//______________________________________________________________________________
bool
EventHBX::ProcessingNormal( void )
{
  static const double MinTFACUT = gUser.GetParameter("TCut", 0);
  static const double MaxTFACUT = gUser.GetParameter("TCut", 1);
  static const double MinCRMCUT = gUser.GetParameter("CCut", 0);
  static const double MaxCRMCUT = gUser.GetParameter("CCut", 1);
  //  static const int MinBGOCUT = gUser.GetParameter("BCut", 0);
  //  static const int MaxBGOCUT = gUser.GetParameter("BCut", 1);

  gRM.Decode();   

  event.runnum = gRM.RunNumber();
  event.evnum  = gRM.EventNumber();
  event.spill  = gRM.SpillNumber();
  dst.evnum    = gRM.EventNumber();
  dst.spill    = gRM.SpillNumber();
  
  // Trigger Flag
  {
    int trignhits = 0;
    for(int seg=0; seg<NumOfSegTrig; seg++){
      int nhits_trig = gUnpacker.get_entries( DetIdTrig, 0, seg, 0, 1 );
      for(int hit=0; hit<nhits_trig; hit++){
	int trig = gUnpacker.get( DetIdTrig, 0, seg, 0, 1, hit );
	if(trig>0){
	  event.trigpat[seg] = seg; 
	  event.trigflag[seg] = trig; 
	  HF1( 10, seg);
	  HF1( 10+seg+1, trig);
	  trignhits++;
	}
      }
    } 
    int hbxtrignhits = 0;
    for(int seg=0; seg<NumOfSegHbxTrig; seg++){
      int nhits_hbxtrig = gUnpacker.get_entries( DetIdHbxTrig, 0, seg, 0, 1 );
      for(int hit=0; hit<nhits_hbxtrig; hit++){
	int hbxtrig = gUnpacker.get( DetIdHbxTrig, 0, seg, 0, 1, hit );
	if(hbxtrig>0){
	  event.hbxtrigpat[seg] = seg; 
	  event.hbxtrigflag[seg] = hbxtrig; 
	  HF1( 100, seg);
	  HF1( 100+seg+1, hbxtrig);
	  hbxtrignhits++;
	}
      }
    } 
    event.trignhits = trignhits;
    event.hbxtrignhits = hbxtrignhits;
  }
   
  //  if( trigflag[SpillEndFlag] ) return true;   
   
  //---Ge ADC & TDC-----------------------------------------------------
  for( int seg=0; seg<NumOfSegGe; ++seg ){
    int nhit_a = gUnpacker.get_entries( DetIdGe, 0, seg, 0, 0 ); //adc
    int nhit_t = gUnpacker.get_entries( DetIdGe, 0, seg, 0, 1 ); //tfa
    int nhit_c = gUnpacker.get_entries( DetIdGe, 0, seg, 0, 4 ); //crm
    int nhit_r = gUnpacker.get_entries( DetIdGe, 0, seg, 0, 7 ); //reset

    if( nhit_a>0 ){
      int adc = gUnpacker.get( DetIdGe, 0, seg, 0, 0 );
      event.geadc[seg] = adc;
      HF1( GeHid+100*(seg+1)+10, double(adc) );
      HF2( GeHid+100 +0, seg+0.5, double(adc) );
      if( event.hbxtrigflag[LSOGeFlag]>0  ) event.lsogeadc[seg]      = adc;	 
      if( event.hbxtrigflag[GeCoinFlag]>0 ) event.coingeadc[seg]     = adc;	 
      if( event.hbxtrigflag[SpillOnFlag]>0) event.spillongeadc[seg]  = adc;
      if( event.hbxtrigflag[SpillOffFlag]>0  ) event.spilloffgeadc[seg] = adc;

      //adc w/ tfa
      if( nhit_t>0){
	int tfa  = gUnpacker.get( DetIdGe, 0, seg, 0, 1, 0);
	if( tfa>0 ){
	  HF1( GeHid+100*(seg+1)+21, double(adc) );
	  HF2( GeHid+100*(seg+1)+24, double(tfa), double(adc) );
	}
	else{
	  HF1( GeHid+100*(seg+1)+22, double(adc) );
	}
	//adc w/ tfa cut
	if( MinTFACUT<tfa && tfa<MaxTFACUT ) HF1( GeHid+100*(seg+1)+23, double(adc) );
      }

      //adc w/ crm
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

      //adc w/ reset
      if( nhit_r>0){
	int reset_time = gUnpacker.get( DetIdGe, 0, seg, 0, 7);
	HF2( GeHid+100*(seg+1)+44, double(reset_time), double(adc) );
      }
    }
   
    //tfa
    if( nhit_t>0 ){
      HF1( GeHid+1, seg+0.5); //0 origin
      for(int i = 0; i<nhit_t; ++i){
	int tfa  = gUnpacker.get( DetIdGe, 0, seg, 0, 1, i );
	if( tfa>0 ){
	  HF1( GeHid+100*(seg+1)+20, double(tfa) ); //0 origin
	  HF2( GeHid+100+2, seg+0.5, double(tfa) );
	}
	event.getfa[seg][i] = tfa;
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
	}
      }
    }

    //reset
    if( nhit_r>0 ){
      int reset_time = gUnpacker.get( DetIdGe, 0, seg, 0, 7);
      HF1( GeHid+100*(seg+1)+40, double(reset_time) );
      HF2( GeHid+100+4, seg+0.5, double(reset_time) );
      event.gereset[seg] = reset_time;
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
	event.gebgot[seg][nhit] = tdc;
      }
    }
  }
  
  return true;
}

//______________________________________________________________________________
bool
EventHBX::ProcessingEnd( void )
{
  tree->Fill();
  return true;

}

//______________________________________________________________________________
void
EventHBX::InitializeEvent( void )
{
  event.evnum = 0;
  event.spill = 0;
  event.gebgonhits  = 0;
  event.trignhits = 0;
  event.hbxtrignhits = 0;

  for( int it=0; it<NumOfSegTrig; it++){
    event.trigpat[it] = -1;
    event.trigflag[it] = -1;
  }
  for( int it=0; it<NumOfSegHbxTrig; it++){
    event.hbxtrigpat[it] = -1;
    event.hbxtrigflag[it] = -1;
  }
 
  //Ge Reset
  for( int it=0; it<NumOfSegGe; it++){
    event.gereset[it] = 0;
  }

  //Ge ADC
  for( int it=0; it<NumOfSegGe; it++){
    event.geadc[it] = -9999.;

    event.lsogeadc[it] = -9999.;	 
    event.coingeadc[it] = -9999.;	
    event.spillongeadc[it] = -9999.;
    event.spilloffgeadc[it] = -9999.;
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
    event.scaler[it] = 0;
  }
}

//______________________________________________________________________________
VEvent*
ConfMan::EventAllocator( void )
{
  return new EventHBX;
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
  tree->Branch("spill",  &event.spill,  "spill/I");
  //Trig
  tree->Branch("trignhits", &event.trignhits, "trignhits/I");
  tree->Branch("hbxtrignhits", &event.hbxtrignhits, "hbxtrignhits/I");

  HB1(  1, "Status", 20, 0., 20. );
  HB1( 10, "Trigger HitPat", NumOfSegTrig, 0., Double_t(NumOfSegTrig) );
  for(Int_t i=0; i<NumOfSegTrig; ++i){
    HB1( 10+i+1, Form("Trigger Trig %d", i+1), 0x1000, 0, 0x1000 );
  }
  HB1( 100, "HBX Trigger HitPat", NumOfSegHbxTrig, 0., Double_t(NumOfSegHbxTrig) );
  for(Int_t i=0; i<NumOfSegHbxTrig; ++i){
    HB1( 100+i+1, Form("HBX Trigger Trig %d", i+1), 0x1000, 0, 0x1000 );
  }

  
  //---Ge ADC & TDC-----------------------------------------------------
  for( int i=1; i<=NumOfSegGe; ++i ){
    TString title10  = Form("Ge-%d Adc", i);

    TString title20  = Form("Ge-%d Tfa", i);	     
    TString title21  = Form("Ge-%d Adc(w Tfa)", i);  
    TString title22  = Form("Ge-%d Adc(w/o Tfa)", i);
    TString title23  = Form("Ge-%d Adc(Tfa Cut)", i);
    TString title24  = Form("Adc%%Tfa-%d", i); 

    TString title30  = Form("Ge-%d Crm", i);	     
    TString title31  = Form("Ge-%d Adc(w Crm)", i);  
    TString title32  = Form("Ge-%d Adc(w/o Crm)", i);
    TString title33  = Form("Ge-%d Adc(Crm Cut)", i);
    TString title34  = Form("Adc%%Crm-%d", i); 

    TString title40  = Form("Ge-%d Reset", i);	     
    TString title44  = Form("Adc%%Reset-%d", i); 

    //TString title41 = Form("Ge-%d Adc [LSO*Ge]", i);
    //TString title42 = Form("Ge-%d Adc [GeCoin]", i);
    //TString title43 = Form("Ge-%d Adc [Spill ON]", i);
    //TString title44 = Form("Ge-%d Adc [Spill OFF]", i);

    HB1( GeHid +100*i +10, title10, NbinAdc, MinAdc, MaxAdc );

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

  //HBXFlag
  tree->Branch("hbxtrigpat", event.hbxtrigpat, Form("hbxtrigpat[%d]/I", NumOfSegHbxTrig)); 
  tree->Branch("hbxtrigflag", event.hbxtrigflag, Form("hbxtrigflag[%d]/I", NumOfSegHbxTrig)); 

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
  tree->Branch("scaler", event.scaler, Form("scaler[%d]/I", NumOfSegScaler)); 

  //Ge adc  w/ Flag
  tree->Branch("lsogeadc"     , event.lsogeadc,      Form("lsogeadc[%d]/I",      NumOfSegGe)); 
  tree->Branch("coingeadc"    , event.coingeadc,     Form("coingeadc[%d]/I",     NumOfSegGe)); 
  tree->Branch("spillongeadc" , event.spillongeadc,  Form("spillongeadc[%d]/I",  NumOfSegGe)); 
  tree->Branch("spilloffgeadc", event.spilloffgeadc, Form("spilloffgeadc[%d]/I", NumOfSegGe)); 

  /////Dst/////////////////////////////////
  hbx = new TTree( "hbx", "Data Summary Table of hbx" );


  HPrint();
  return true;
}

//______________________________________________________________________________
bool
ConfMan::InitializeParameterFiles( void )
{
  return
    ( InitializeParameter<DCGeomMan>("DCGEO")        &&
      InitializeParameter<HodoParamMan>("HDPRM")     &&
      InitializeParameter<UserParamMan>("USER")      );
}

//______________________________________________________________________________
bool
ConfMan::FinalizeProcess( void )
{
  return true;
}
