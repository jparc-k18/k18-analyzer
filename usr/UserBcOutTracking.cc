/*
  UserBcOutTracking.cc

  2012/1/24
*/

#include "VEvent.hh"

#include <iostream>
#include <sstream>
#include <cmath>
#include <string>
#include <iomanip>
#include <set>

#include "UnpackerManager.hh"
#include "ConfMan.hh"
#include "HistHelper.hh"
#include "DetectorID.hh"
#include "RawData.hh"
#include "DCRawHit.hh"
#include "SksLib.hh"

#include "FileManager.hh"
#include "ParameterManager.hh"

const int BEAM_TRIG  = 1;
const int PIK_TRIG   = 2;
const int PIPI_TRIG  = 11;
const int PIP_TRIG   = 12;
const int BH2_TRIG   = 4;

const double Deg2Rad = acos(-1.)/180.;
const double Rad2Deg = 180./acos(-1.);

const double MaxMultiHitBcOut = 50.0;

const double TdcLow  = 100.;
const double TdcHigh = 1500.;

#define HodoCut 1

 VEvent::VEvent()
 {
 }

 VEvent::~VEvent()
 {
 }

class EventBcOutTracking 
  : public VEvent
{

private:
  RawData *rawData;
  DCAnalyzer *DCAna;
  HodoAnalyzer *hodoAna;

public:
  EventBcOutTracking();
  ~EventBcOutTracking();

  bool ProcessingBegin();
  bool ProcessingEnd();
  bool ProcessingNormal();
  bool InitializeHistograms();
  void InitializeEvent();
};

EventBcOutTracking::EventBcOutTracking()
  : VEvent(),
    rawData(0), 
    DCAna(new DCAnalyzer()),
    hodoAna(new HodoAnalyzer())
{
}

EventBcOutTracking::~EventBcOutTracking()
{
  if (DCAna)   delete DCAna;
  if (hodoAna) delete hodoAna;
  if (rawData) delete rawData;
}

bool EventBcOutTracking::ProcessingBegin()
{
  //  ConfMan::InitializeHistograms(); 
 return true;
}

#ifndef MaxHits 
#define MaxHits 20
#endif

#ifndef MaxDepth 
#define MaxDepth 100
#endif

//For Tree
struct Event{  
  int runnum;
  int evnum;
  int trigpat[NumOfMisc];
  int trigflag[NumOfMisc];

  //BH2
  int bh2ut[NumOfSegBH2];
  int bh2dt[NumOfSegBH2];
  int rawnhBh2;
  int rawBh2Seg[NumOfSegBH2];

  //Beam
  double btof[MaxHits];
  int multip_bh1;
  int multip_bh2;

  int nhBh2;
  double T0Seg;
  double Bh2Seg[MaxHits];
  double tBh2[MaxHits];
  double deBh2[MaxHits];

  int nhBh1;
  double Bh1Seg[MaxHits];
  double tBh1[MaxHits];
  double deBh1[MaxHits];

  int ntrack;
  double chisqr[MaxHits];  
  double x_Bh2[MaxHits];
  double x0[MaxHits];  
  double y0[MaxHits];  
  double u0[MaxHits];  
  double v0[MaxHits]; 
  double pos[12][50]; 

  // BC34
  int BCOUT_TDC[NumOfLayersBcOut][MaxDepth];
  int NHits_wF[NumOfSegBH2][12];
};
static Event event;

// BH2 Cut Parameters
const double MinDeltaEBH2 = 0.5; //for Pion
const double MaxDeltaEBH2 = 3.0; //for Pion

// BH1 Cut Parameters
const double MinDeltaEBH1 = 0.5; //for Pion
const double MaxDeltaEBH1 = 3.0; //for Pion

// BH1-BH2
const double MinBeamToF  = -1.0;
const double MaxBeamToF  =  1.0;

bool EventBcOutTracking::ProcessingNormal()
{
  rawData = new RawData;
  rawData->DecodeHits();

  //Tree
  TTree *tree = dynamic_cast<TTree *>(gFile->Get("tree"));
  InitializeEvent();
  event.runnum =
    hddaq::unpacker::GUnpacker::get_instance().get_run_number();
  event.evnum =
    hddaq::unpacker::GUnpacker::get_instance().get_event_number();

  //Misc
  /*
  bool fl_trigger = false;
  {
    const HodoRHitContainer &cont=rawData->GetMiscRawHC();
    ParamMan& Trig = FileManager::GetInstance().GetParamMan("TriggerFlag");
    Trig.stChangeParamSet("TriggerFlag");
    static int TargetTrig = (int)Trig.GetParameter(0, "TARGET");
    int nh=cont.size();
    for( int i=0; i<nh; ++i ){
      HodoRawHit *hit=cont[i];
      int seg=hit->SegmentId()+1;
      int T=hit->GetTdc1();

      if( T>0 ){
	event.trigpat[i] = seg;
	if(TargetTrig == seg){ fl_trigger = true;}
      } 
    }
  }
  */
  HF1( 1, 0. );

  //------------------------Cut
  //if(!fl_trigger){return true;}
  HF1( 1, 1. );

#if HodoCut
 //////////////BH2 time 0
  hodoAna->DecodeBH2Hits(rawData);
  int ncBh2=hodoAna->GetNClustersBH2();
  event.nhBh2=ncBh2;
  if( ncBh2==0 ){ tree->Fill(); return true; }

  HF1( 1, 2. );

  //////////////BH2 Analysis
  BH2Cluster *clBH2Time0=hodoAna->GetClusterBH2(0);
  double time0=clBH2Time0->CTime0();

  std::vector<BH2Cluster*> valid_bh2;
  {
    int ncOk=0;
    double mint=clBH2Time0->CMeanTime();
    for( int i=0; i<ncBh2; ++i ){
      BH2Cluster *cl=hodoAna->GetClusterBH2(i);
      double t=cl->CMeanTime();
      double dEbh2 = cl->DeltaE();
      //------------------------Cut
      if( dEbh2<MinDeltaEBH2 || MaxDeltaEBH2<dEbh2 ) continue;
      event.Bh2Seg[i] = cl->MeanSeg()+1;
      event.tBh2[i]   = t;
      event.deBh2[i]  = cl->DeltaE();
      ++ncOk;
      if( fabs(t)<fabs(mint) ){
	clBH2Time0 = cl;
	mint=t; time0=clBH2Time0->CTime0();
      }

      if(-1.<t && t < 1.){
	valid_bh2.push_back(cl);
      }
    }

    double trig = clBH2Time0->CMeanTime();
    if(-0.5 > trig || trig > 0.5 || ncOk<1){
      //return true;
    }
    
    if(false
       || valid_bh2.size() > 1
       || clBH2Time0->ClusterSize() > 1
       ){
      event.multip_bh2 = 1;
    }
    event.T0Seg = clBH2Time0->MeanSeg() +1;
  }

 //////////////BH2 Analysis
  {
    const HodoRHitContainer &cont=rawData->GetBH2RawHC();
    int nh=cont.size();
    for( int i=0; i<nh; ++i ){
      HodoRawHit *hit=cont[i];
      int seg=hit->SegmentId();
      event.bh2ut[seg] = hit->GetTdcUp();
      event.bh2dt[seg] = hit->GetTdcDown();
    }

    event.rawnhBh2 = hodoAna->GetNHitsBH2();
    for(int i = 0; i<event.rawnhBh2; ++i){
      const BH2Hit* const h = hodoAna->GetHitBH2(i);
      if(!h){
	continue;
      }
      event.rawBh2Seg[i] = h->SegmentId();
    }
  }

  HF1( 1, 3. );

 //////////////BH1 Analysis
  hodoAna->DecodeBH1Hits(rawData); 
  int ncBh1=hodoAna->GetNClustersBH1();
  event.nhBh1=ncBh1;
  if(0 == ncBh1){ tree->Fill(); return true; };
  HF1( 1, 4. );

  std::vector<HodoCluster*> valid_bh1;
  HodoCluster* clBH1Time0 = hodoAna->GetClusterBH1(0);  
  {
    int ncOk=0;
    double min_tof = clBH1Time0->CMeanTime() - time0;
    for( int i=0; i<ncBh1; ++i ){
      HodoCluster *cl=hodoAna->GetClusterBH1(i);
      double btof= (cl->CMeanTime())-time0;
      double dEbh1 = cl->DeltaE();
      //------------------------Cut
      if( dEbh1<MinDeltaEBH1 || MaxDeltaEBH1<dEbh1 ) continue;
      event.Bh1Seg[i]=cl->MeanSeg()+1;
      event.tBh1[i]=cl->CMeanTime();
      event.deBh1[i]=cl->DeltaE();
      //------------------------Cut
      if( MinBeamToF<btof && btof<MaxBeamToF ){
	valid_bh1.push_back(cl);
	++ncOk;
	if( fabs(btof)<fabs(min_tof) ){
	  clBH1Time0 = cl;
	  min_tof = btof;
	}
      }
      else{
	cl->GoodForAnalysis( false );
      }
    }
    //------------------------Cut
    //    if( ncOk<1 ) return true;

    if(false
       || valid_bh1.size() > 1
       || clBH1Time0->ClusterSize() > 2
       ){
      event.multip_bh1 = 1;
    }
  }
  HF1( 1, 5. );
#endif

  HF1( 1, 10. );

  //////////////BC3&4 number of hit layer
  DCAna->DecodeRawHits( rawData );  

  //BC3&BC4

  double multi_BcOut=0.;
  {
    for( int layer=1; layer<=NumOfLayersBcOut; ++layer ){
      const DCHitContainer &contOut =DCAna->GetBcOutHC(layer);
      int nhOut=contOut.size();
      multi_BcOut += double(nhOut);

      HF1( 100*layer, nhOut );

      int plane_eff = (layer-1)*3;
      bool fl_valid_sig = false; 
      int Current = 0;
      for( int i=0; i<nhOut; ++i ){
	DCHit *hit=contOut[i];
	double wire=hit->GetWire();
	
	HF1( 100*layer+1, wire-0.5 );

	int nhtdc = hit->GetTdcSize();
	int tdc1st = -1;
	for( int k=0; k<nhtdc; k++ ){
	  int tdc = hit->GetTdcVal(k);
	  if( (TdcLow< tdc && tdc < TdcHigh) && tdc > tdc1st ){
	    tdc1st=tdc;
	    fl_valid_sig = true;
	  }
	}

	if(nhOut<8){
	  if( TdcLow< tdc1st && tdc1st < TdcHigh ){
	    event.pos[layer-1][Current++] = hit->GetWirePosition();
	  }
	}

	if(i < MaxDepth){
	  event.BCOUT_TDC[layer-1][i] = tdc1st;
	}else{
	  std::cout << "over flow" << std::endl;
	}
	HF1( 100*layer+2, tdc1st );
	HF1( 10000*layer+int(wire), tdc1st );

	int nhdt = hit->GetDriftTimeSize();
	for( int k=0; k<nhdt; k++ ){
	  double dt = hit->GetDriftTime(k);
	  HF1( 100*layer+3, dt );
	  HF1( 10000*layer+1000+int(wire), dt );
	}
	int nhdl = hit->GetDriftTimeSize();
	for( int k=0; k<nhdl; k++ ){
	  double dl = hit->GetDriftLength(k);
	  HF1( 100*layer+4, dl );
	}
      }
      
      if(fl_valid_sig){
	++plane_eff;
	//	std::cout << "!";
      }else{
	//	std::cout << "." ;
      }
      HF1(38, plane_eff);
    }    
  }
  //  std::cout << std::endl;
  if( multi_BcOut/double(NumOfLayersBcOut) > MaxMultiHitBcOut )
    { tree->Fill(); return true; }


  
  HF1( 1, 11. );
  
#if 1
  // Bc Out
  //  std::cout << "==========TrackSearch BcOut============" << std::endl;
  DCAna->TrackSearchBcOut( );
  //  DCAna->ChiSqrCutBcOut(10);
  //DCAna->TrackSearchBcOut();

  //DCAna->ApplyBh2SegmentCut(event.T0Seg);
  int nt=DCAna->GetNtracksBcOut();
  HF1( 10, double(nt) );
  // the tree arrays are dimensioned to MaxHits
  if( nt>MaxHits ) nt=MaxHits;
  event.ntrack=nt;
  int valid_nt = 0;
  for( int it=0; it<nt; ++it ){
    DCLocalTrack *tp=DCAna->GetTrackBcOut(it);

    int nh=tp->GetNHit();
    double chisqr=tp->GetChiSquare();
    double x0=tp->GetX0(), y0=tp->GetY0();
    double u0=tp->GetU0(), v0=tp->GetV0();

    double cost = 1./sqrt(1.+u0*u0+v0*v0);
    double theta = acos(cost)*Rad2Deg;

    if(chisqr < 10){++valid_nt;}

    event.chisqr[it]=chisqr;
    event.x0[it]=x0;
    event.y0[it]=y0;
    event.u0[it]=u0;
    event.v0[it]=v0;

    //    double x_bh2 = tp->GetX(603.63);
    //event.x_Bh2[it] = x_bh2;

    HF1( 11, double(nh) );
    HF1( 12, chisqr ); 
    HF1( 14, x0 ); HF1( 15, y0 );
    HF1( 16, u0 ); HF1( 17, v0 );
    HF2( 18, x0, u0 ); HF2( 19, y0, v0 );
    HF2( 20, x0, y0 );

    double xtgt=tp->GetX(1800.), ytgt=tp->GetY(1800.);
    double utgt=u0, vtgt=v0;
    HF1( 21, xtgt ); HF1( 22, ytgt );
    HF1( 23, utgt ); HF1( 24, vtgt );
    HF2( 25, xtgt, utgt ); HF2( 26, ytgt, vtgt );
    HF2( 27, xtgt, ytgt );

    double xbac=tp->GetX(603.), ybac=tp->GetY(885.);
    double ubac=u0, vbac=v0;
    HF1( 31, xbac ); HF1( 32, ybac );
    HF1( 33, ubac ); HF1( 34, vbac );
    HF2( 35, xbac, ubac ); HF2( 36, ybac, vbac );
    HF2( 37, xbac, ybac );

    double Xangle = -1000*std::atan(u0);

    HF2(51, -tp->GetX(245.), Xangle);
    HF2(52, -tp->GetX(600.), Xangle);
    HF2(53, -tp->GetX(1200.), Xangle);
    HF2(54, -tp->GetX(1600.), Xangle);

    HF2(61, tp->GetX(600.), tp->GetX(245.));
    HF2(62, tp->GetX(1200.), tp->GetX(245.));
    HF2(63, tp->GetX(1600.), tp->GetX(245.));

    for( int ih=0; ih<nh; ++ih ){
      DCLTrackHit *hit=tp->GetHit(ih);
      int layerId=hit->GetLayer()-112; 

      HF1( 13, layerId );
      double wire=hit->GetWire();
      double dt=hit->GetDriftTime(), dl=hit->GetDriftLength();
      HF1( 100*layerId+11, wire-0.5 );
      HF1( 100*layerId+12, dt );
      HF1( 100*layerId+13, dl );
      HF1( 10000*layerId+ 5000 +int(wire), dt);
      double xcal=hit->GetXcal(), ycal=hit->GetYcal();
      double pos=hit->GetLocalHitPos(), res=hit->GetResidual();
      HF1( 100*layerId+14, pos );
      HF1( 100*layerId+15, res );
      HF2( 100*layerId+16, pos, res );
      HF2( 100*layerId+17, xcal, ycal);
      //      HF1( 100000*layerId+50000+wire, res);      
      double wp=hit->GetWirePosition();
      double sign=1.;
      if( pos-wp<0. ) sign=-1;
      HF2( 100*layerId+18, sign*dl, res );
      double xlcal=hit->GetLocalCalPos();
      HF2( 100*layerId+19, dt, xlcal-wp);

      if (theta>=0 && theta<15) 
	HF1( 100*layerId+71, res );
      else if (theta>=15 && theta<30) 
	HF1( 100*layerId+72, res );
      else if (theta>=30 && theta<45) 
	HF1( 100*layerId+73, res );
      else if (theta>=45) 
	HF1( 100*layerId+74, res );
      
      if (fabs(dl-fabs(xlcal-wp))<2.0) {      
	HFProf( 100*layerId+20, dt, fabs(xlcal-wp));
	HF2( 100*layerId+22, dt, fabs(xlcal-wp));
	HFProf( 100000*layerId+3000+int(wire), xlcal-wp,dt);
	HF2( 100000*layerId+4000+int(wire), xlcal-wp,dt);
      }
    }
  }
  //  if(event.multip_bh2 == 0){
  //  HF1(40, valid_nt);
  //}

#endif

  HF1( 1, 12. );

  tree->Fill();

  return true;
}

void EventBcOutTracking::InitializeEvent( void )
{
  event.runnum = -1;
  event.evnum = -1;
  // size branches must never stay negative at Fill
  event.ntrack   = 0;
  event.nhBh2  = 0;
  event.rawnhBh2  = 0;
  event.nhBh1  = 0;
  event.T0Seg = -100;
  event.multip_bh1 = 0;
  event.multip_bh2 = 0;

  for( int it=0; it<MaxHits; it++){
    event.Bh2Seg[it] = -1;
    event.tBh2[it] = -9999.0;
    event.deBh2[it] = -9999.0;

    event.Bh1Seg[it] = -1;
    event.tBh1[it] = -9999.0;
    event.deBh1[it] = -9999.0;
    event.btof[it]  = -1;
  }

  for( int it=0; it<MaxHits; it++){
    event.chisqr[it] = -1.0;
    event.x0[it] = -9999.0;
    event.y0[it] = -9999.0;
    event.u0[it] = -9999.0;
    event.v0[it] = -9999.0;
    event.x_Bh2[it] = -9999.0;
  }

  for( int it=0; it<NumOfMisc; it++){ 
    event.trigpat[it] = -1;
    event.trigflag[it] = -1;
  }

  for(int i = 0; i<NumOfSegBH2; ++i){
    event.bh2ut[i] = -1;
    event.bh2dt[i] = -1;
    event.rawBh2Seg[i] = -100;
    for(int layer = 0; layer<12; ++layer){
      event.NHits_wF[i][layer] = -10;
    }
  }

  for(int layer = 0; layer<NumOfLayersBcOut; ++layer){
    for(int m = 0; m<MaxDepth; ++m){
      event.BCOUT_TDC[layer][m] = -1;
    }
  }

  for(int layer = 0; layer<12; ++layer){
    for(int i = 0; i<50; ++i){
      event.pos[layer][i] = -9999.;
    }
  }
}

bool EventBcOutTracking::ProcessingEnd()
{
  return true;
}

VEvent *ConfMan::EventAllocator()
{
  return new EventBcOutTracking;
}


const int NbinBcOutTdc   = 2000;
const double MinBcOutTdc =    0.;
const double MaxBcOutTdc = 2000.;

// 1 ns/bin with integer edges: the drift time is quantized at the
// TDC LSB (1 ns), any other bin width makes a moire pattern
const int NbinBcOutDT   =   70;
const double MinBcOutDT =  -10.;
const double MaxBcOutDT =   60.;

// drift length is clamped to [0, 1.5] mm (DriftLength6)
const int NbinBcOutDL   = 100;
const double MinBcOutDL =  -0.25;
const double MaxBcOutDL =   1.75;

bool ConfMan:: InitializeHistograms()
{  
  HB1( 1, "Status", 20, 0., 20. );

  //***********************Chamber
  // BC3
  for( int i=1; i<=NumOfLayersBc; ++i ){
    std::ostringstream title1, title2, title3, title4, title5;
    title1 << "#Hits BC3#" << std::setw(2) << i;
    title2 << "Hitpat BC3#" << std::setw(2) << i;
    title3 << "Tdc BC3#" << std::setw(2) << i;
    title4 << "Drift Time BC3#" << std::setw(2) << i;
    title5 << "Drift Length BC3#" << std::setw(2) << i;

    HB1( 100*i+0, title1.str().c_str(), MaxWireBC3+1, 0., double(MaxWireBC3+1) );
    HB1( 100*i+1, title2.str().c_str(), MaxWireBC3+1, 0., double(MaxWireBC3+1) );
    HB1( 100*i+2, title3.str().c_str(), NbinBcOutTdc, MinBcOutTdc, MaxBcOutTdc );
    HB1( 100*i+3, title4.str().c_str(), NbinBcOutDT, MinBcOutDT, MaxBcOutDT );
    HB1( 100*i+4, title5.str().c_str(), NbinBcOutDL, MinBcOutDL, MaxBcOutDL );

    for (int wire=1; wire<=MaxWireBC3; wire++) {
      std::ostringstream title11, title12, title13, title14;
      title11 << "Tdc BC3#" << std::setw(2) << i << " Wire#" << wire;
      HB1( 10000*i+wire, title11.str().c_str(), NbinBcOutTdc, MinBcOutTdc, MaxBcOutTdc );
      title12 << "Drift Time BC3#" << std::setw(2) << i << " Wire#" << wire;
      HB1( 10000*i+1000+wire, title12.str().c_str(), NbinBcOutDT, MinBcOutDT, MaxBcOutDT );
      title13 << "Drift Length BC3#" << std::setw(2) << i << " Wire#" << wire;
      HB1( 10000*i+2000+wire, title13.str().c_str(), NbinBcOutDL, MinBcOutDL, MaxBcOutDL );
      // Added
      title14 << "Drift Time BC3 after tracking#" << std::setw(2) << i << " Wire#" << wire;
      HB1( 10000*i+5000+wire, title14.str().c_str(), NbinBcOutDT, MinBcOutDT, MaxBcOutDT );
    }
  }

  // BC4
  for( int i=1; i<=NumOfLayersBc; ++i ){
    std::ostringstream title1, title2, title3, title4, title5;
    title1 << "#Hits BC4#" << std::setw(2) << i;
    title2 << "Hitpat BC4#" << std::setw(2) << i;
    title3 << "Tdc BC4#" << std::setw(2) << i;
    title4 << "Drift Time BC4#" << std::setw(2) << i;
    title5 << "Drift Length BC4#" << std::setw(2) << i;
    HB1( 100*(i+6)+0, title1.str().c_str(), MaxWireBC4+1, 0., double(MaxWireBC4+1) );
    HB1( 100*(i+6)+1, title2.str().c_str(), MaxWireBC4+1, 0., double(MaxWireBC4+1) );
    HB1( 100*(i+6)+2, title3.str().c_str(), NbinBcOutTdc, MinBcOutTdc, MaxBcOutTdc );
    HB1( 100*(i+6)+3, title4.str().c_str(), NbinBcOutDT, MinBcOutDT, MaxBcOutDT );
    HB1( 100*(i+6)+4, title5.str().c_str(), NbinBcOutDL, MinBcOutDL, MaxBcOutDL );

    for (int wire=1; wire<=MaxWireBC4; wire++) {
      std::ostringstream title11, title12, title13, title14;
      title11 << "Tdc BC4#" << std::setw(2) << i << " Wire#" << wire;
      HB1( 10000*(i+6)+wire, title11.str().c_str(), NbinBcOutTdc, MinBcOutTdc, MaxBcOutTdc );
      title12 << "Drift Time BC4#" << std::setw(2) << i << " Wire#" << wire;
      HB1( 10000*(i+6)+1000+wire, title12.str().c_str(), NbinBcOutDT, MinBcOutDT, MaxBcOutDT );
      title13 << "Drift Length BC4#" << std::setw(2) << i << " Wire#" << wire;
      HB1( 10000*(i+6)+2000+wire, title13.str().c_str(), NbinBcOutDL, MinBcOutDL, MaxBcOutDL );
      // added
      title14 << "Drift Time BC4 after tracking#" << std::setw(2) << i << " Wire#" << wire;
      HB1( 10000*(i+6)+5000+wire, title14.str().c_str(), NbinBcOutDT, MinBcOutDT, MaxBcOutDT );

    }
  }

  // Tracking Histgrams
  HB1( 10, "#Tracks BcOut", 10, 0., 10. );
  HB1( 11, "#Hits of Track BcOut", 15, 0., 15. );
  HB1( 12, "Chisqr BcOut", 500, 0., 50. ); 
  HB1( 13, "LayerId BcOut", 15, 0., 15. );
  HB1( 14, "X0 BcOut", 400, -100., 100. ); 
  HB1( 15, "Y0 BcOut", 400, -100., 100. );
  HB1( 16, "U0 BcOut", 200, -0.20, 0.20 );
  HB1( 17, "V0 BcOut", 200, -0.20, 0.20 );
  HB2( 18, "U0%X0 BcOut", 100, -100., 100., 100, -0.20, 0.20 );
  HB2( 19, "V0%Y0 BcOut", 100, -100., 100., 100, -0.20, 0.20 );
  HB2( 20, "X0%Y0 BcOut", 100, -100., 100., 100, -100, 100 );

  HB1( 21, "Xtgt BcOut", 400, -100., 100. ); 
  HB1( 22, "Ytgt BcOut", 400, -100., 100. );
  HB1( 23, "Utgt BcOut", 200, -0.20, 0.20 );
  HB1( 24, "Vtgt BcOut", 200, -0.20, 0.20 );
  HB2( 25, "Utgt%Xtgt BcOut", 100, -100., 100., 100, -0.20, 0.20 );
  HB2( 26, "Vtgt%Ytgt BcOut", 100, -100., 100., 100, -0.20, 0.20 );
  HB2( 27, "Xtgt%Ytgt BcOut", 100, -100., 100., 100, -100, 100 );

  HB1( 31, "Xbac BcOut", 400, -100., 100. ); 
  HB1( 32, "Ybac BcOut", 400, -100., 100. );
  HB1( 33, "Ubac BcOut", 200, -0.20, 0.20 );
  HB1( 34, "Vbac BcOut", 200, -0.20, 0.20 );
  HB2( 35, "Ubac%Xbac BcOut", 100, -100., 100., 100, -0.20, 0.20 );
  HB2( 36, "Vbac%Ybac BcOut", 100, -100., 100., 100, -0.20, 0.20 );
  HB2( 37, "Xbac%Ybac BcOut", 100, -100., 100., 100, -100, 100 );

  HB2( 51, "X-X' 245 BcOut", 400, -100., 100., 120, -60, 60); 
  HB2( 52, "X-X' 600 BcOut", 400, -100., 100., 120, -60, 60); 
  HB2( 53, "X-X' 1200 BcOut", 400, -100., 100., 120, -60, 60); 
  HB2( 54, "X-X' 1600 BcOut", 400, -100., 100., 120, -60, 60); 

  HB2( 61, "X-X 600 BcOut", 400, -100., 100., 400, -100, 100); 
  HB2( 62, "X-X 1200 BcOut", 400, -100., 100., 400, -100, 100); 
  HB2( 63, "X-X 1600 BcOut", 400, -100., 100., 400, -100, 100); 

  // Plane eff
  HB1( 38, "Plane Eff", 36, 0, 36);
  HB1( 40, "#Tracks BcOut w/ chisqr cut", 10, 0., 10. );

  for( int i=1; i<=NumOfLayersBcOut; ++i ){
    std::ostringstream title0, title1, title2, title3, title4;
    std::ostringstream title5, title6, title7, title8, title9;
    std::ostringstream title10, title11, title12, title13, title14;
    
    title0 << "wire for LayerId = " << std::setw(2) << i<< " [Track]";
    HB1( 100*i+11, title0.str().c_str(), 64, 0., 64. );
    title1 << "drift time for LayerId = " << std::setw(2) << i<< " [Track]";
    HB1( 100*i+12, title1.str().c_str(), NbinBcOutDT, MinBcOutDT, MaxBcOutDT );
    title2 << "drift length for LayerId = " << std::setw(2) << i<< " [Track]";
    HB1( 100*i+13, title2.str().c_str(), NbinBcOutDL, MinBcOutDL, MaxBcOutDL );
    title3 << "Position BcOut" << std::setw(2) << i;
    HB1( 100*i+14, title3.str().c_str(), 100, -250., 250. );
    title4 << "Residual BcOut" << std::setw(2) << i;
    HB1( 100*i+15, title4.str().c_str(), 200, -1.0, 1.0 );
    title5 << "Resid%Pos BcOut" << std::setw(2) << i;
    HB2( 100*i+16, title5.str().c_str(), 250, -250., 250., 100, -1.0, 1.0 );
    title6 << "Y%Xcal BcOut" << std::setw(2) << i;
    HB2( 100*i+17, title6.str().c_str(), 100, -250., 250., 100, -250., 250. );
    title7 << "Res%dl BcOut" << std::setw(2) << i;
    HB2( 100*i+18, title7.str().c_str(), 160, -1.6, 1.6, 100, -0.5, 0.5 );
    title8 << "Hit Pos%Drift Time BcOut" << std::setw(2) << i;
    HB2( 100*i+19, title8.str().c_str(), NbinBcOutDT, MinBcOutDT, MaxBcOutDT, 100, -2., 2.);

    title9 << "Drift Length%Drift Time BcOut" << std::setw(2) << i;
    HBProf( 100*i+20, title9.str().c_str(), 1000, -5, 50, 0, 12);
    HB2( 100*i+22, title9.str().c_str(), 850, -5, 80, 250, 0, 5);

    title4 << " w/o self plane"; 
    HB1( 100*i+21, title4.str().c_str(), 200, -5.0, 5.0 );

    title10 << "Residual BcOut (0<theta<15)" << std::setw(2) << i;
    HB1( 100*i+71, title10.str().c_str(), 200, -5.0, 5.0 );

    title11 << "Residual BcOut (15<theta<30)" << std::setw(2) << i;
    HB1( 100*i+72, title11.str().c_str(), 200, -5.0, 5.0 );

    title12 << "Residual BcOut (30<theta<45)" << std::setw(2) << i;
    HB1( 100*i+73, title12.str().c_str(), 200, -5.0, 5.0 );

    title13 << "Residual BcOut (45<theta)" << std::setw(2) << i;
    HB1( 100*i+74, title13.str().c_str(), 200, -5.0, 5.0 );

    for (int j=1; j<=64; j++) {
      std::ostringstream title;
      title << "XT of Layer " << std::setw(2) << i << " Wire #"<< std::setw(4) << j;
      HBProf( 100000*i+3000+j, title.str().c_str(), 101, -4.4, 100, -5,40);
      HB2( 100000*i+4000+j, title.str().c_str(), 100, -4, 4, 100, -5,40);
    }
    
    //     title9 << " w/o self plane"; 
    //     HB2( 100*i+23, title9.str().c_str(), 100, -5., 5., 150, -50, 100);
    //     HBProf( 100*i+24, title9.str().c_str(), 100, -50, 150, 0, 5);
    //     HB2( 100*i+25, title9.str().c_str(), 100, -50, 150, 100,0, 3.5);
    //     title3 << " w/o self plane"; 
    //     HB2( 100*i+26, title3.str().c_str(), 50, -100., 100., 50, -1.0, 1.0 );
  }
  
  ////////////////////////////////////////////
  //Tree
  HBTree("tree","tree of Counter");
  TTree *tree = dynamic_cast<TTree *>(gFile->Get("tree"));

  tree->Branch("runnum",    &event.runnum,    "runnum/I");
  tree->Branch("evnum",     &event.evnum,     "evnum/I");
  tree->Branch("trigpat",    event.trigpat,   "trigpat[20]/I");
  tree->Branch("trigflag",   event.trigflag,   "trigflag[20]/I");

  //Hodoscope
  tree->Branch("nhBh2",   &event.nhBh2,   "nhBh2/I");
  tree->Branch("rawnhBh2",&event.rawnhBh2,"rawnhBh2/I");
  tree->Branch("T0Seg",   &event.T0Seg,   "T0Seg/D");
  tree->Branch("Bh2Seg",   event.Bh2Seg,  "Bh2Seg[nhBh2]/D");
  tree->Branch("tBh2",     event.tBh2,    "tBh2[nhBh2]/D");
  tree->Branch("deBh2",    event.deBh2,   "deBh2[nhBh2]/D");
  tree->Branch("bh2ut",    event.bh2ut,   "bh1ut[7]/I");
  tree->Branch("bh2dt",    event.bh2dt,   "bh1dt[7]/I");
  tree->Branch("rawBh2Seg",event.rawBh2Seg,"rawBh2Seg[rawnhBh2]/I");

  tree->Branch("nhBh1",   &event.nhBh1,   "nhBh1/I");
  tree->Branch("Bh1Seg",   event.Bh1Seg,  "Bh1Seg[nhBh1]/D");
  tree->Branch("tBh1",     event.tBh1,    "tBh1[nhBh1]/D");
  tree->Branch("deBh1",    event.deBh1,   "deBh1[nhBh1]/D");
  tree->Branch("btof",     event.btof,    "btof[nhBh1]/D");
  tree->Branch("mp_bh1",  &event.multip_bh1, "mp_bh1/I");
  tree->Branch("mp_bh2",  &event.multip_bh2, "mp_bh2/I");

  tree->Branch("ntrack",   &event.ntrack,   "ntrack/I");
  tree->Branch("chosqr",    event.chisqr,   "chisqr[ntrack]/D");
  tree->Branch("x_bh2",     event.x_Bh2,    "x_bh2[ntrack]/D");
  tree->Branch("x0",        event.x0,       "x0[ntrack]/D");
  tree->Branch("y0",        event.y0,       "y0[ntrack]/D");
  tree->Branch("u0",        event.u0,       "u0[ntrack]/D");
  tree->Branch("v0",        event.v0,       "v0[ntrack]/D");

  tree->Branch("posBC3X0",  event.pos[0],   "posBC3X0[50]/D");
  tree->Branch("posBC3X1",  event.pos[1],   "posBC3X1[50]/D");
  tree->Branch("posBC3V0",  event.pos[2],   "posBC3V0[50]/D");
  tree->Branch("posBC3V1",  event.pos[3],   "posBC3V1[50]/D");
  tree->Branch("posBC3U0",  event.pos[4],   "posBC3U0[50]/D");
  tree->Branch("posBC3U1",  event.pos[5],   "posBC3U1[50]/D");

  tree->Branch("posBC4V0",  event.pos[6],   "posBC4V0[50]/D");
  tree->Branch("posBC4V1",  event.pos[7],   "posBC4V1[50]/D");
  tree->Branch("posBC4U0",  event.pos[8],   "posBC4U0[50]/D");
  tree->Branch("posBC4U1",  event.pos[9],   "posBC4U1[50]/D");
  tree->Branch("posBC4X0",  event.pos[10],  "posBC4X0[50]/D");
  tree->Branch("posBC4X1",  event.pos[11],  "posBC4X1[50]/D");

  tree->Branch("BCOUT_TDC", event.BCOUT_TDC, "BCOUT_TDC[12][100]/I");

  tree->Branch("Hits_wF", event.NHits_wF, "Hits_wF[7][12]/I");

  return true;
}

bool ConfMan::InitializeParameterFiles( void )
{
  bool status = true;
  
  if( HodoParamFileName_!="" )
    HodoParamManager_ = new HodoParamMan(HodoParamFileName_);
  if(HodoParamManager_) HodoParamManager_->Initialize();

  if( HodoPHCFileName_!="" )
    HodoPHCManager_ = new HodoPHCMan(HodoPHCFileName_);
  if( HodoPHCManager_ ) HodoPHCManager_->Initialize();

  DCGeomManager_ = & DCGeomMan::GetInstance();
  if( DCGeomFileName_!="" )
    DCGeomManager_->Initialize(DCGeomFileName_);
  else
    DCGeomManager_->Initialize();

  if( DCTdcCalibFileName_!="" )
    DCTdcCalibManager_ = new DCTdcCalibMan(DCTdcCalibFileName_);
  if(DCTdcCalibManager_) DCTdcCalibManager_->Initialize();

  if( DCDriftParamFileName_!="" )
    DCDriftParamManager_ = new DCDriftParamMan(DCDriftParamFileName_);
  if(DCDriftParamManager_) DCDriftParamManager_->Initialize();

  if(FileManagerName_ != ""){
    h_Utility::FileManager& gFMan = h_Utility::FileManager::GetInstance();
    if(-1 == gFMan.stInitialize(FileManagerName_)){
      status = false;
    }
  }

  return status;
}
