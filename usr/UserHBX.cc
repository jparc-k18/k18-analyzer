/**
 *  file: UserSkeleton.cc
 *  date: 2017.04.10
 *
 */

#include <iostream>
#include <sstream>
#include <cmath>

#include "ConfMan.hh"
#include "DetectorID.hh"
#include "RMAnalyzer.hh"
#include "RootHelper.hh"
#include "HodoRawHit.hh"
#include "KuramaLib.hh"
#include "RawData.hh"
#include "UnpackerManager.hh"
#include "VEvent.hh"

namespace
{
  using namespace root;
  const std::string& classname("EventSkeleton");
  RMAnalyzer& gRM = RMAnalyzer::GetInstance();
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
class EventSkeleton : public VEvent
{
private:
  RawData *rawData;

public:
        EventSkeleton( void );
       ~EventSkeleton( void );
  bool  ProcessingBegin( void );
  bool  ProcessingEnd( void );
  bool  ProcessingNormal( void );
  bool  InitializeHistograms( void );
  void  InitializeEvent( void );
};

//______________________________________________________________________________
EventSkeleton::EventSkeleton( void )
  : VEvent(),
    rawData(0)
{
}

//______________________________________________________________________________
EventSkeleton::~EventSkeleton( void )
{
  if (rawData) delete rawData;
}

//______________________________________________________________________________
struct Event
{
  int runnum;
  int evnum;
  int spill;

  int trignhits;
  //int trigpat[NumOfSegTrig];
  //int trigflag[NumOfSegTrig];

  //Ge Reset
  int gereset[NumOfSegGe];

  //Ge ADC
  double geadc[NumOfSegGe];

  //Ge TFA
  int genhits;
  int gehitpat[MaxHits];
  double getfa[NumOfSegGe][MaxDepth];
  double gecrm[NumOfSegGe][MaxDepth];
  
  //Ge CRM
  int gecrmnhits;
  int gecrmhitpat[MaxHits];

  
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
EventSkeleton::ProcessingBegin( void )
{
  InitializeEvent();
  return true;
}

//______________________________________________________________________________
bool
EventSkeleton::ProcessingNormal( void )
{

   gRM.Decode();   

   event.runnum = gRM.RunNumber();
   event.evnum  = gRM.EventNumber();
   event.spill  = gRM.SpillNumber();
   dst.evnum    = gRM.EventNumber();
   dst.spill    = gRM.SpillNumber();
  
   //Trig
   
   
   //---Ge ADC & TDC-----------------------------------------------------
   for( int seg=0; seg<NumOfSegGe; ++seg ){
     int nhit_a = gUnpacker.get_entries( DetIdGe, 0, seg, 0, 0 ); //adc
     int nhit_t = gUnpacker.get_entries( DetIdGe, 0, seg, 0, 1 ); //tfa
     int nhit_c = gUnpacker.get_entries( DetIdGe, 0, seg, 0, 4 ); //crm
     int nhit_r = gUnpacker.get_entries( DetIdGe, 0, seg, 0, 7 ); //reset

     int adc = 0;
     //int tfa = 0;
     //int crm = 0;

     if( nhit_a>0 ){
       adc = gUnpacker.get( DetIdGe, 0, seg, 0, 0 );
       event.geadc[seg] = adc;
       HF1( GeHid+100*(seg+1)+1, double(adc) );
       HF2( GeHid+21, seg+0.5, double(adc) );

       //adc w/ tfa
       if( nhit_t>0){
	 int tfa  = gUnpacker.get( DetIdGe, 0, seg, 0, 1, 0);
	 if( tfa>0 ){
	   HF1( GeHid+100*(seg+1)+3, double(adc) );
	   HF2( GeHid+100*(seg+1)+31, double(tfa), double(adc) );
	 }
	 else{
	   HF1( GeHid+100*(seg+1)+4, double(adc) );
	 }
       }

       //adc w/ crm
       if( nhit_t>0){
	 int crm  = gUnpacker.get( DetIdGe, 0, seg, 0, 1, 0);
	 if( crm>0 ){
	   HF1( GeHid+100*(seg+1)+6, double(adc) );
	   HF2( GeHid+100*(seg+1)+32, double(crm), double(adc) );
	 }
	 else{
	   HF1( GeHid+100*(seg+1)+7, double(adc) );
	 }
       }
 
     }
   
     //tfa
     if( nhit_t>0 ){
       HF1( GeHid+1, seg+0.5); //0 origin
       for(int i = 0; i<nhit_t; ++i){
	 int tfa  = gUnpacker.get( DetIdGe, 0, seg, 0, 1, i );
	 if( tfa>0 ){
	   HF1( GeHid+100*(seg+1)+2, double(tfa) ); //0 origin
	   HF2( GeHid+22, seg+0.5, double(tfa) );
	 }
	 event.getfa[seg][i] = tfa;
       }
       
       //if(event.gebgot[seg-1][0] > 0){
       //event.gebgohitpat[gebgonhits++] = seg;
       //}
       
     }

     //crm
     if( nhit_c>0 ){
       for(int i = 0; i<nhit_c; ++i){
	 //HF1( GeHid+1, seg+0.5); //0 origin
	 int crm  = gUnpacker.get( DetIdGe, 0, seg, 0, 4, i )  ;
	 if(crm > 0){
	   HF1( GeHid+100*(seg+1)+5, double(crm) ); //0 origin
	   HF2( GeHid+23, seg+0.5, double(crm) );
	 event.gecrm[seg][i] = crm;

       }
       /*
       if(event.gebgot[seg-1][0] > 0){
	 event.gebgohitpat[gebgonhits++] = seg;
	 }
       */
       }
     }

     //reset
     if( nhit_r>0 ){
       int reset_time = gUnpacker.get( DetIdGe, 0, seg, 0, 7 );
       HF1( GeHid+100*(seg+1)+8, double(reset_time) );
       HF2( GeHid+24, seg+0.5, double(reset_time) );
       HF2( GeHid+100*(seg+1)+33, double(reset_time), double(adc) );
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
	   HF1( BGOHid+100*(seg+1)+3, double(tdc) );
	   HF2( BGOHid+21, seg+0.5, double(tdc) );
	 }
	 event.gebgot[seg][nhit] = tdc;
	 
	 //int trailing = gUnpacker.get( DetIdBGO, 0, seg, 0, 1, i )  ;
	 //std::cout << leading << std::endl; 
       }
       /*
	 if(event.gebgot[seg-1][0] > 0){
	 event.gebgohitpat[gebgonhits++] = seg;
	 }
       */
     }
   }

   
#if 0
   // Ge Scaler
   
   for( int seg=0; seg<NumOfSegScaler; ++seg ){
     int nhit = gUnpacker.get_entries( DetIdScaler, 3, 0, seg, 0 );
     if( nhit>0 ){
       int data = gUnpacker.get( DetIdScaler, 3, 0, seg, 0 );
       event.scaler[seg] = data;
     }
   }
   // for(l)
#endif




#if 0
   //Ge ADC
   for( int seg=0; seg<NumOfSegGe; ++seg ){
     int nhit = gUnpacker.get_entries( DetIdGe, 0, seg, 0, 0 );
     if( nhit>0 ){
       int adc = gUnpacker.get( DetIdGe, 0, seg, 0, 0 );
       HF1( GeHid+100*(seg+1)+1, double(adc) );
       event.geadc[seg] = adc;
     }
   } 
#endif   
 
#if 0
   //Ge TFA
   for(int seg = 0; seg<NumOfSegGe; ++seg){
     int nhit = gUnpacker.get_entries( DetIdGe, 0, seg, 0, 1 );
     if( nhit>0 ){
       //std::cout << "nhit" << nhit << std::endl;
       for(int i = 0; i<nhit; ++i){
	 HF1( GeHid+1, seg+0.5); //0 origin
	 int leading  = gUnpacker.get( DetIdGe, 0, seg, 0, 1, i )  ;
	 if(leading > 0) HF1( GeHid+100*(seg+1)+3, double(leading) ); //0 origin
	 event.getfa[seg][i] = leading;
	 //int trailing = gUnpacker.get( DetIdGe, 0, seg, 0, 1, i )  ;
       }
       /*
       if(event.gebgot[seg-1][0] > 0){
	 event.gebgohitpat[gebgonhits++] = seg;
       }
       */
     }
   }
#endif

#if 0
   
   //Ge CRM
   for(int seg = 0; seg<NumOfSegGe; ++seg){
     int nhit = gUnpacker.get_entries( DetIdGe, 0, seg, 0, 4 );
     if( nhit>0 ){
       for(int i = 0; i<nhit; ++i){

	 //HF1( GeHid+1, seg+0.5); //0 origin
	 int tdc  = gUnpacker.get( DetIdGe, 0, seg, 0, 4, i )  ;

	 if(tdc > 0) HF1( GeHid+100*(seg+1)+4, double(tdc) ); //0 origin
	 event.gecrm[seg][i] = tdc;

       }
       /*
       if(event.gebgot[seg-1][0] > 0){
	 event.gebgohitpat[gebgonhits++] = seg;
	 }
       */
     }
   }


#endif

#if 0
   //Ge Reset
   for( int seg=0; seg<NumOfSegGe; ++seg ){
     int nhit = gUnpacker.get_entries( DetIdGe, 0, seg, 0, 7 );

     if( nhit>0 ){
       int reset_time = gUnpacker.get( DetIdGe, 0, seg, 0, 7 );
       HF2( GeHid+24, seg+0.5, double(reset_time) );
       event.gereset[seg] = reset_time;
     }

   } 
#endif

  

   return true;
}

//______________________________________________________________________________
bool
EventSkeleton::ProcessingEnd( void )
{
  tree->Fill();
  return true;

}

//______________________________________________________________________________
void
EventSkeleton::InitializeEvent( void )
{
  event.evnum = 0;
  event.spill = 0;
  event.gebgonhits  = 0;

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
    event.scaler[it] = 0;
  }
}

//______________________________________________________________________________
VEvent*
ConfMan::EventAllocator( void )
{
  return new EventSkeleton;
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
  
  //---Ge ADC & TDC-----------------------------------------------------
  for( int i=1; i<=NumOfSegGe; ++i ){
    TString title1  = Form("Ge-%d Adc", i);
    TString title2  = Form("Ge-%d Tfa", i);
    TString title3  = Form("Ge-%d Adc(w Tfa)", i);
    TString title4  = Form("Ge-%d Adc(w/o Tfa)", i);
    TString title5  = Form("Ge-%d Crm", i);
    TString title6  = Form("Ge-%d Adc(w Crm)", i);
    TString title7  = Form("Ge-%d Adc(w/o Crm)", i);
    TString title8  = Form("Ge-%d Reset", i);

    TString title31 = Form("Adc%%Tfa-%d", i);
    TString title32 = Form("Adc%%Crm-%d", i);
    TString title33 = Form("Adc%%Reset-%d", i);

    HB1( GeHid +100*i +1, title1, NbinAdc, MinAdc, MaxAdc );
    HB1( GeHid +100*i +2, title2, NbinGeTdc, MinGeTdc, MaxGeTdc );
    HB1( GeHid +100*i +3, title3, NbinAdc, MinAdc, MaxAdc );
    HB1( GeHid +100*i +4, title4, NbinAdc, MinAdc, MaxAdc );
    HB1( GeHid +100*i +5, title5, NbinGeTdc, MinGeTdc, MaxGeTdc );
    HB1( GeHid +100*i +6, title6, NbinAdc, MinAdc, MaxAdc );
    HB1( GeHid +100*i +7, title7, NbinAdc, MinAdc, MaxAdc );
    HB1( GeHid +100*i +8, title8, NbinGeTdc, MinGeTdc, MaxGeTdc );

    HB2( GeHid +100*i +31, title31, 
    	 NbinGeTdc/16, MinGeTdc, MaxGeTdc, NbinAdc/8, MinAdc, MaxAdc);
    HB2( GeHid +100*i +32, title32,
    	 NbinGeTdc/16, MinGeTdc, MaxGeTdc, NbinAdc/8, MinAdc, MaxAdc);
    HB2( GeHid +100*i +33, title33, 
	 NbinGeTdc/16, MinGeTdc, MaxGeTdc, NbinAdc/8, MinAdc, MaxAdc);
  }

  //2D Hist
  HB2( GeHid +21, "Adc%Ch",
       NumOfSegGe, 0, NumOfSegGe, NbinAdc, MinAdc, MaxAdc);
  HB2( GeHid +22, "Tdc%Ch",
       NumOfSegGe, 0, NumOfSegGe, NbinGeTdc, MinGeTdc, MaxGeTdc);
  HB2( GeHid +23, "Tdc%Ch",
       NumOfSegGe, 0, NumOfSegGe, NbinGeTdc, MinGeTdc, MaxGeTdc);
  HB2( GeHid +24, "Reset%Ch",
       NumOfSegGe, 0, NumOfSegGe, NbinGeTdc, MinGeTdc, MaxGeTdc);
  HB2( BGOHid +21, "BGOTdc%Ch",
       NumOfSegBGO, 0, NumOfSegBGO, NbinBGOTdc, MinBGOTdc, MaxBGOTdc);


  //Hits
  HB1( GeHid +0, "#Hits Ge",        NumOfSegGe+1, 0., double(NumOfSegGe+1) );
  HB1( GeHid +1, "Hitpat Ge",       NumOfSegGe,   0., double(NumOfSegGe)   );


  //---BGO--------------------------------------------------------------
  HB1( BGOHid +0, "#Hits BGO",        NumOfSegBGO+1, 0., double(NumOfSegBGO+1) );
  HB1( BGOHid +1, "Hitpat BGO",       NumOfSegBGO,   0., double(NumOfSegBGO)   );

  for( int i=1; i<=NumOfSegBGO; ++i ){
    TString title13 = Form("BGO-%d Tdc", i);
    HB1( BGOHid +100*i +3, title13, NbinBGOTdc, MinBGOTdc, MaxBGOTdc );
  }


  /////tree/////////////////////////////////
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
    ( InitializeParameter<DCGeomMan>("DCGEO")    &&
      InitializeParameter<HodoParamMan>("HDPRM") );
}

//______________________________________________________________________________
bool
ConfMan::FinalizeProcess( void )
{
  return true;
}
