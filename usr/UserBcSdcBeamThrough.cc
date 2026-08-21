/*
  UserBcSdcBeamThrough.cc

  2012/1/24
*/

#include "VEvent.hh"

#include <iostream>
#include <sstream>
#include <cmath>
#include <string>
#include <iomanip>

#include "UnpackerManager.hh"
#include "ConfMan.hh"
#include "HistHelper.hh"
#include "DetectorID.hh"
#include "RawData.hh"
#include "DCRawHit.hh"
#include "SksLib.hh"

const int BEAM_TRIG  = 1;
const int PIK_TRIG   = 2;
const int PIPI_TRIG  = 11;
const int PIP_TRIG   = 12;
const int BH2_TRIG   = 4;

const double Deg2Rad = acos(-1.)/180.;
const double Rad2Deg = 180./acos(-1.);
const double TdcLow  =    0.;
const double TdcHigh = 1100.;

#define HodoCut 1

 VEvent::VEvent()
 {
 }

 VEvent::~VEvent()
 {
 }

class EventBcSdcBeamThrough 
  : public VEvent
{

private:
  RawData *rawData;
  DCAnalyzer *DCAna;
  HodoAnalyzer *hodoAna;

public:
  EventBcSdcBeamThrough();
  ~EventBcSdcBeamThrough();

  bool ProcessingBegin();
  bool ProcessingEnd();
  bool ProcessingNormal();
  bool InitializeHistograms();
  void InitializeEvent();
};

EventBcSdcBeamThrough::EventBcSdcBeamThrough()
  : VEvent(),
    rawData(0), 
    DCAna(new DCAnalyzer()),
    hodoAna(new HodoAnalyzer())
{
}

EventBcSdcBeamThrough::~EventBcSdcBeamThrough()
{
  if (DCAna)   delete DCAna;
  if (hodoAna) delete hodoAna;
  if (rawData) delete rawData;
}

bool EventBcSdcBeamThrough::ProcessingBegin()
{
 return true;
}

#ifndef MaxHits 
#define MaxHits 20
#endif

//For Tree
struct Event{
  int runnum;
  int evnum;
  int trigpat[NumOfMisc];
  int trigflag[NumOfMisc];

  //Beam
  double btof[MaxHits];

  int nhBh2;
  double Bh2Seg[MaxHits];
  double tBh2[MaxHits];
  double deBh2[MaxHits];

  int nhBh1;
  double Bh1Seg[MaxHits];
  double tBh1[MaxHits];
  double deBh1[MaxHits];

  int ntrackb;
  double chisqrb[MaxHits];  
  double x0b[MaxHits];  
  double y0b[MaxHits];  
  double u0b[MaxHits];  
  double v0b[MaxHits];  

  int ntracks;
  double chisqrs[MaxHits];
  double x0s[MaxHits];
  double y0s[MaxHits];
  double u0s[MaxHits];
  double v0s[MaxHits];

  int ntrackc;
  double chisqrc[MaxHits];
};
static Event event;

const double MinDeltaEBH2 = 0.5;
const double MaxDeltaEBH2 = 4.0;

const double MinDeltaEBH1 = 0.5;
const double MaxDeltaEBH1 = 4.0;

const double MinBeamToF  = -2.0;
const double MaxBeamToF  =  2.0;

const double MaxMultiHitSdcIn = 100.0;
const double MaxMultiHitBcOut = 100.0;

bool EventBcSdcBeamThrough::ProcessingNormal()
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
  int trig_type=0;
  int trig=0;
  {
    const HodoRHitContainer &cont=rawData->GetMiscRawHC();
    int nh=cont.size();
    for( int i=0; i<nh; ++i ){
      HodoRawHit *hit=cont[i];
      int seg=hit->SegmentId()+1;
      int T=hit->GetTdc1();

      //pi beam
      if( seg==1 && T>0 ) trig_type = seg;
      //(pi, K)
      //if( seg==2 && T>0 ) trig_type = seg;
      //(pi, pi)
      //if( seg==11 && T>0 ) trig_type = seg;
      event.trigflag[seg-1] = T;

      if( T>0 ) event.trigpat[i] = seg;
    }
  }

  //------------------------Cut
  //  if (trig_type != BEAM_TRIG )  return true;
  //if (trig_type != PIK_TRIG )  return true;
  //if (trig_type != PIPI_TRIG )  return true;
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
    }
    //if( ncOk<1 ) return true;
  }

  HF1( 1, 3. );

  //////////////BH1 Analysis
  hodoAna->DecodeBH1Hits(rawData); 
  int ncBh1=hodoAna->GetNClustersBH1();
  event.nhBh1=ncBh1;
  {
    int ncOk=0;
    for( int i=0; i<ncBh1; ++i ){
      HodoCluster *cl=hodoAna->GetClusterBH1(i);
      double btof= (cl->CMeanTime())-time0;
      double dEbh1 = cl->DeltaE();
      //------------------------Cut
      if( dEbh1<MinDeltaEBH1 || MaxDeltaEBH1<dEbh1 ) continue;
      event.Bh1Seg[i]=cl->MeanSeg()+1;
      event.tBh1[i]=cl->CMeanTime();
      event.deBh1[i]=cl->DeltaE();
      event.btof[i] = btof;
      //------------------------Cut
      if( MinBeamToF<btof && btof<MaxBeamToF ){
	++ncOk;
      }
      else{
	cl->GoodForAnalysis( false );
      }
    }
    //------------------------Cut
    if( ncOk<1 ){ tree->Fill(); return true; }
  }
  HF1( 1, 4. );
#endif

  //tree->Fill();return true;

  HF1( 1, 10. );

  DCAna->DecodeRawHits( rawData );  
  double multi_BcOut=0.;
  double multi_SdcIn=0.;

 //BC3&BC4
  {
    for( int layer=1; layer<=NumOfLayersBcOut; ++layer ){
      const DCHitContainer &contOut =DCAna->GetBcOutHC(layer);
      int nhOut=contOut.size();
      multi_BcOut += double(nhOut);
      HF1( 100*(layer+50), nhOut );
      for( int i=0; i<nhOut; ++i ){
	DCHit *hit=contOut[i];
	double wire=hit->GetWire();
	
	HF1( 100*(layer+50)+1, wire-0.5 );
	int nhtdc = hit->GetTdcSize();
	int tdc1st = -1;
	for( int k=0; k<nhtdc; k++ ){
	  int tdc = hit->GetTdcVal(k);
	  if( (TdcLow< tdc && tdc < TdcHigh) && tdc > tdc1st ) tdc1st=tdc;
	}
	HF1( 100*(layer+50)+2, tdc1st );
	HF1( 10000*(layer+50)+int(wire), tdc1st );

	int nhdt = hit->GetDriftTimeSize();
	for( int k=0; k<nhdt; k++ ){
	  double dt = hit->GetDriftTime(k);
	  HF1( 100*(layer+50)+3, dt );
	  HF1( 10000*(layer+50)+1000+int(wire), dt );
	}
	int nhdl = hit->GetDriftTimeSize();
	for( int k=0; k<nhdl; k++ ){
	  double dl = hit->GetDriftLength(k);
	  HF1( 100*(layer+50)+4, dl );
	}
      }
    }    
  }

  //SDC1&SDC2
  {
    for( int layer=1; layer<=NumOfLayersSdcIn; ++layer ){
      const DCHitContainer &contIn =DCAna->GetSdcInHC(layer);
      int nhIn=contIn.size();
      multi_SdcIn += double(nhIn);
      HF1( 100*layer, nhIn );
      for( int i=0; i<nhIn; ++i ){
	DCHit *hit=contIn[i];
	double wire=hit->GetWire();
	
	HF1( 100*layer+1, wire-0.5 );
	int nhtdc = hit->GetTdcSize();
	int tdc1st = -1;
	for( int k=0; k<nhtdc; k++ ){
	  int tdc = hit->GetTdcVal(k);
	  if( (TdcLow< tdc && tdc < TdcHigh) && tdc > tdc1st ) tdc1st=tdc;
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
    }    
  }
  if( multi_SdcIn/double(NumOfLayersSdcIn) > MaxMultiHitSdcIn ||
      multi_BcOut/double(NumOfLayersBcOut) > MaxMultiHitBcOut)
   { tree->Fill(); return true; }

 HF1( 1, 11. );

#if 1
  // Bc Out
  //  std::cout << "==========TrackSearch BcOut============" << std::endl;
  DCAna->TrackSearchBcOut();
  int nhitsb=0;
  int ntb=DCAna->GetNtracksBcOut();
  HF1( 50, double(ntb) );
  for( int it=0; it<ntb; ++it ){
    DCLocalTrack *tp=DCAna->GetTrackBcOut(it);
    int nh=tp->GetNHit();
    double chisqr=tp->GetChiSquare();
    double x0=tp->GetX0(), y0=tp->GetY0();
    double u0=tp->GetU0(), v0=tp->GetV0();

    double cost = 1./sqrt(1.+u0*u0+v0*v0);
    double theta = acos(cost)*Rad2Deg;

    HF1( 51, double(nh) );
    HF1( 52, chisqr ); 
    HF1( 54, x0 ); HF1( 55, y0 );
    HF1( 56, u0 ); HF1( 57, v0 );
    HF2( 58, x0, u0 ); HF2( 59, y0, v0 );
    HF2( 60, x0, y0 );

    double xtgt=tp->GetX(1624.7), ytgt=tp->GetY(1624.7);
    double utgt=u0, vtgt=v0;
    HF1( 61, xtgt ); HF1( 62, ytgt );
    HF1( 63, utgt ); HF1( 64, vtgt );
    HF2( 65, xtgt, utgt ); HF2( 66, ytgt, vtgt );
    HF2( 67, xtgt, ytgt );

    event.chisqrb[nhitsb]=chisqr;
    event.x0b[nhitsb]=xtgt;
    event.y0b[nhitsb]=ytgt;
    event.u0b[nhitsb]=utgt;
    event.v0b[nhitsb]=vtgt;
    nhitsb++;

    for( int ih=0; ih<nh; ++ih ){
      DCLTrackHit *hit=tp->GetHit(ih);
      int layerId=hit->GetLayer()-62; 
      HF1( 53, layerId );
      double wire=hit->GetWire();
      double dt=hit->GetDriftTime(), dl=hit->GetDriftLength();
      HF1( 100*layerId+11, wire-0.5 );
      HF1( 100*layerId+12, dt );
      HF1( 100*layerId+13, dl );
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
  event.ntrackb = nhitsb;
#endif

#if 1
  // Sdc In
  //std::cout << "==========TrackSearch SdcIn============" << std::endl;
  DCAna->TrackSearchSdcIn();
  int nhitss=0;
  int nts=DCAna->GetNtracksSdcIn();
  HF1( 10, double(nts) );
  for( int it=0; it<nts; ++it ){
    DCLocalTrack *tp=DCAna->GetTrackSdcIn(it);
    int nh=tp->GetNHit();
    double chisqr=tp->GetChiSquare();
    double x0=tp->GetX0(), y0=tp->GetY0();
    double u0=tp->GetU0(), v0=tp->GetV0();

    double cost = 1./sqrt(1.+u0*u0+v0*v0);
    double theta = acos(cost)*Rad2Deg;

    HF1( 11, double(nh) );
    HF1( 12, chisqr ); 
    HF1( 14, x0 ); HF1( 15, y0 );
    HF1( 16, u0 ); HF1( 17, v0 );
    HF2( 18, x0, u0 ); HF2( 19, y0, v0 );
    HF2( 20, x0, y0 );

    double xtgt=tp->GetX(0.), ytgt=tp->GetY(0.);
    double utgt=u0, vtgt=v0;
    HF1( 21, xtgt ); HF1( 22, ytgt );
    HF1( 23, utgt ); HF1( 24, vtgt );
    HF2( 25, xtgt, utgt ); HF2( 26, ytgt, vtgt );
    HF2( 27, xtgt, ytgt );

    event.chisqrs[nhitss]=chisqr;
    event.x0s[nhitss]=xtgt;
    event.y0s[nhitss]=ytgt;
    event.u0s[nhitss]=utgt;
    event.v0s[nhitss]=vtgt;
    nhitss++;

    for( int ih=0; ih<nh; ++ih ){
      DCLTrackHit *hit=tp->GetHit(ih);
      int layerId=hit->GetLayer(); 
      HF1( 13, layerId );
      double wire=hit->GetWire();
      double dt=hit->GetDriftTime(), dl=hit->GetDriftLength();
      HF1( 100*layerId+11, wire-0.5 );
      HF1( 100*layerId+12, dt );
      HF1( 100*layerId+13, dl );
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
  event.ntracks = nhitss;
#endif

#if 1
  //////////////BcOut + SdcIn combined straight tracking
  // (for the DCGeom alignment with the SKS-off beam through)
  DCAna->TrackSearchBcOutSdcIn();
  int ntc=DCAna->GetNtracksBcOutSdcIn();
  HF1( 90, double(ntc) );
  if( ntc>MaxHits ) ntc=MaxHits;
  event.ntrackc = ntc;
  for( int it=0; it<ntc; ++it ){
    DCLocalTrack *tp=DCAna->GetTrackBcOutSdcIn(it);
    if( !tp ) continue;
    int nh=tp->GetNHit();
    double chisqr=tp->GetChiSquare();
    HF1( 91, double(nh) );
    HF1( 92, chisqr );
    event.chisqrc[it]=chisqr;
    if( chisqr>30. ) continue;
    double u0c=tp->GetU0(), v0c=tp->GetV0();
    for( int ih=0; ih<nh; ++ih ){
      DCLTrackHit *hit=tp->GetHit(ih);
      int layer=hit->GetLayer();
      int layerId = (layer>=100) ? layer-62 : layer;
      double res=hit->GetResidual();
      double pos=hit->GetLocalHitPos();
      double wp=hit->GetWirePosition();
      double dl=hit->GetDriftLength();
      double sign=1.;
      if( pos-wp<0. ) sign=-1;
      HF1( 100*layerId+25, res );
      HF2( 100*layerId+26, sign*dl, res );
      // z-consistency check: a wrong plane z shows up as a linear
      // correlation of the residual with the track slope along the
      // plane's measurement direction (slope = effective delta-z)
      double aa=hit->GetTiltAngle()*Deg2Rad;
      HF2( 100*layerId+27, u0c*cos(aa)+v0c*sin(aa), res );
    }
  }
#endif

  HF1( 1, 12. );

  tree->Fill();

  return true;
}

void EventBcSdcBeamThrough::InitializeEvent( void )
{
  event.runnum = -1;
  event.evnum = -1;
  // size branches must be 0 (not -1): TTree::Fill crashes on
  // negative array sizes when the fill section is skipped
  event.nhBh2    = 0;
  event.nhBh1    = 0;
  event.ntrackb  = 0;
  event.ntracks  = 0;
  event.ntrackc  = 0;
    
  for( int it=0; it<MaxHits; it++){
    event.chisqrb[it] = -1.0;
    event.x0b[it] = -9999.0;
    event.y0b[it] = -9999.0;
    event.u0b[it] = -9999.0;
    event.v0b[it] = -9999.0;

    event.chisqrs[it] = -1.0;
    event.x0s[it] = -9999.0;
    event.y0s[it] = -9999.0;
    event.u0s[it] = -9999.0;
    event.v0s[it] = -9999.0;
  }

  for( int it=0; it<NumOfMisc; it++){
    event.trigpat[it] = -1;
    event.trigflag[it] = -1;
  }

  for( int it=0; it<MaxHits; it++){
    event.Bh2Seg[it] = -1;
    event.tBh2[it]  = -9999.0;
    event.deBh2[it] = -9999.0;

    event.Bh1Seg[it] = -1;
    event.tBh1[it]  = -9999.0;
    event.deBh1[it] = -9999.0;
    event.btof[it]  = -1;
  }
}

bool EventBcSdcBeamThrough::ProcessingEnd()
{
  return true;
}

VEvent *ConfMan::EventAllocator()
{
  return new EventBcSdcBeamThrough;
}


const int NbinTdc   = 2000;
const double MinTdc =    0.;
const double MaxTdc = 2000.;

const int NbinDT   =   100;
const double MinDT =  -30.;
const double MaxDT =   80.;

const int NbinDL   = 100;
const double MinDL =  -1.;
const double MaxDL =   3.;

bool ConfMan:: InitializeHistograms()
{  
  HB1( 1, "Status", 20, 0., 20. );

  //***********************Chamber
  // BC3
  for( int i=1+50; i<=NumOfLayersBc+50; ++i ){
    std::ostringstream title1, title2, title3, title4, title5;
    title1 << "#Hits BC3#" << std::setw(2) << i;
    title2 << "Hitpat BC3#" << std::setw(2) << i;
    title3 << "Tdc BC3#" << std::setw(2) << i;
    title4 << "Drift Time BC3#" << std::setw(2) << i;
    title5 << "Drift Length BC3#" << std::setw(2) << i;
    
    HB1( 100*i+0, title1.str().c_str(), MaxWireBC3+1, 0., double(MaxWireBC3+1) );
    HB1( 100*i+1, title2.str().c_str(), MaxWireBC3+1, 0., double(MaxWireBC3+1) );
    HB1( 100*i+2, title3.str().c_str(), NbinTdc, MinTdc, MaxTdc );
    HB1( 100*i+3, title4.str().c_str(), NbinDT, MinDT, MaxDT );
    HB1( 100*i+4, title5.str().c_str(), NbinDL, MinDL, MaxDL );
    
    for (int wire=1; wire<=MaxWireBC3; wire++) {
      std::ostringstream title11, title12, title13;
      title11 << "Tdc BC3#" << std::setw(2) << i << " Wire#" << wire;
      HB1( 10000*i+wire, title11.str().c_str(), NbinTdc, MinTdc, MaxTdc );
      title12 << "Drift Time BC3#" << std::setw(2) << i << " Wire#" << wire;
      HB1( 10000*i+1000+wire, title12.str().c_str(), NbinDT, MinDT, MaxDT );
      title13 << "Drift Length BC3#" << std::setw(2) << i << " Wire#" << wire;
      HB1( 10000*i+2000+wire, title13.str().c_str(), NbinDL, MinDL, MaxDL );
    }
  }

  // BC4
  for( int i=1+50; i<=NumOfLayersBc+50; ++i ){
    std::ostringstream title1, title2, title3, title4, title5;
    title1 << "#Hits BC4#" << std::setw(2) << i;
    title2 << "Hitpat BC4#" << std::setw(2) << i;
    title3 << "Tdc BC4#" << std::setw(2) << i;
    title4 << "Drift Time BC4#" << std::setw(2) << i;
    title5 << "Drift Length BC4#" << std::setw(2) << i;
    HB1( 100*(i+6)+0, title1.str().c_str(), MaxWireBC4+1, 0., double(MaxWireBC4+1) );
    HB1( 100*(i+6)+1, title2.str().c_str(), MaxWireBC4+1, 0., double(MaxWireBC4+1) );
    HB1( 100*(i+6)+2, title3.str().c_str(), NbinTdc, MinTdc, MaxTdc );
    HB1( 100*(i+6)+3, title4.str().c_str(), NbinDT, MinDT, MaxDT );
    HB1( 100*(i+6)+4, title5.str().c_str(), NbinDL, MinDL, MaxDL );

    for (int wire=1; wire<=MaxWireBC4; wire++) {
      std::ostringstream title11, title12, title13;
      title11 << "Tdc BC4#" << std::setw(2) << i << " Wire#" << wire;
      HB1( 10000*(i+6)+wire, title11.str().c_str(), NbinTdc, MinTdc, MaxTdc );
      title12 << "Drift Time BC4#" << std::setw(2) << i << " Wire#" << wire;
      HB1( 10000*(i+6)+1000+wire, title12.str().c_str(), NbinDT, MinDT, MaxDT );
      title13 << "Drift Length BC4#" << std::setw(2) << i << " Wire#" << wire;
      HB1( 10000*(i+6)+2000+wire, title13.str().c_str(), NbinDL, MinDL, MaxDL );

    }
  }

  // SDC1
  for( int i=1; i<=NumOfLayersSdc-2; ++i ){
    std::ostringstream title1, title2, title3, title4, title5;
    title1 << "#Hits SDC1#" << std::setw(2) << i;
    title2 << "Hitpat SDC1#" << std::setw(2) << i;
    title3 << "Tdc SDC1#" << std::setw(2) << i;
    title4 << "Drift Time SDC1#" << std::setw(2) << i;
    title5 << "Drift Length SDC1#" << std::setw(2) << i;

    HB1( 100*i+0, title1.str().c_str(), MaxWireSDC1+1, 0., double(MaxWireSDC1+1) );
    HB1( 100*i+1, title2.str().c_str(), MaxWireSDC1+1, 0., double(MaxWireSDC1+1) );
    HB1( 100*i+2, title3.str().c_str(), NbinTdc, MinTdc, MaxTdc );
    HB1( 100*i+3, title4.str().c_str(), NbinDT, MinDT, MaxDT );
    HB1( 100*i+4, title5.str().c_str(), NbinDL, MinDL, MaxDL );

    for (int wire=1; wire<=MaxWireSDC1; wire++) {
      std::ostringstream title11, title12, title13;
      title11 << "Tdc SDC1#" << std::setw(2) << i << " Wire#" << wire;
      HB1( 10000*i+wire, title11.str().c_str(), NbinTdc, MinTdc, MaxTdc );
      title12 << "Drift Time SDC1#" << std::setw(2) << i << " Wire#" << wire;
      HB1( 10000*i+1000+wire, title12.str().c_str(), NbinDT, MinDT, MaxDT );
      title13 << "Drift Length SDC1#" << std::setw(2) << i << " Wire#" << wire;
      HB1( 10000*i+2000+wire, title13.str().c_str(), NbinDL, MinDL, MaxDL );
    }
  }

  // SDC2
  for( int i=1; i<=NumOfLayersSdc; ++i ){
    std::ostringstream title1, title2, title3, title4, title5;
    title1 << "#Hits SDC2#" << std::setw(2) << i;
    title2 << "Hitpat SDC2#" << std::setw(2) << i;
    title3 << "Tdc SDC2#" << std::setw(2) << i;
    title4 << "Drift Time SDC2#" << std::setw(2) << i;
    title5 << "Drift Length SDC2#" << std::setw(2) << i;
    HB1( 100*(i+4)+0, title1.str().c_str(), MaxWireSDC2+1, 0., double(MaxWireSDC2+1) );
    HB1( 100*(i+4)+1, title2.str().c_str(), MaxWireSDC2+1, 0., double(MaxWireSDC2+1) );
    HB1( 100*(i+4)+2, title3.str().c_str(), NbinTdc, MinTdc, MaxTdc );
    HB1( 100*(i+4)+3, title4.str().c_str(), NbinDT, MinDT, MaxDT );
    HB1( 100*(i+4)+4, title5.str().c_str(), NbinDL, MinDL, MaxDL );

    for (int wire=1; wire<=MaxWireSDC2; wire++) {
      std::ostringstream title11, title12, title13;
      title11 << "Tdc SDC2#" << std::setw(2) << i << " Wire#" << wire;
      HB1( 10000*(i+4)+wire, title11.str().c_str(), NbinTdc, MinTdc, MaxTdc );
      title12 << "Drift Time SDC2#" << std::setw(2) << i << " Wire#" << wire;
      HB1( 10000*(i+4)+1000+wire, title12.str().c_str(), NbinDT, MinDT, MaxDT );
      title13 << "Drift Length SDC2#" << std::setw(2) << i << " Wire#" << wire;
      HB1( 10000*(i+4)+2000+wire, title13.str().c_str(), NbinDL, MinDL, MaxDL );

    }
  }

 // Tracking Histgrams BcOut
  HB1( 50, "#Tracks BcOut", 10, 0., 10. );
  HB1( 51, "#Hits of Track BcOut", 15, 0., 15. );
  HB1( 52, "Chisqr BcOut", 500, 0., 50. ); 
  HB1( 53, "LayerId BcOut", 15, 50., 65. );
  HB1( 54, "X0 BcOut", 400, -100., 100. ); 
  HB1( 55, "Y0 BcOut", 400, -100., 100. );
  HB1( 56, "U0 BcOut", 200, -0.20, 0.20 );
  HB1( 57, "V0 BcOut", 200, -0.20, 0.20 );
  HB2( 58, "U0%X0 BcOut", 100, -100., 100., 100, -0.20, 0.20 );
  HB2( 59, "V0%Y0 BcOut", 100, -100., 100., 100, -0.20, 0.20 );
  HB2( 60, "X0%Y0 BcOut", 100, -100., 100., 100, -100, 100 );

  HB1( 61, "Xtgt BcOut", 400, -100., 100. ); 
  HB1( 62, "Ytgt BcOut", 400, -100., 100. );
  HB1( 63, "Utgt BcOut", 200, -0.20, 0.20 );
  HB1( 64, "Vtgt BcOut", 200, -0.20, 0.20 );
  HB2( 65, "Utgt%Xtgt BcOut", 100, -100., 100., 100, -0.20, 0.20 );
  HB2( 66, "Vtgt%Ytgt BcOut", 100, -100., 100., 100, -0.20, 0.20 );
  HB2( 67, "Xtgt%Ytgt BcOut", 100, -100., 100., 100, -100, 100 );

  for( int i=1+50; i<=NumOfLayersBcOut+50; ++i ){
    std::ostringstream title0, title1, title2, title3, title4;
    std::ostringstream title5, title6, title7, title8, title9;
    std::ostringstream title10, title11, title12, title13, title14;
    
    title0 << "wire for LayerId = " << std::setw(2) << i<< " [Track]";
    HB1( 100*i+11, title0.str().c_str(), 64, 0., 64. );
    title1 << "drift time for LayerId = " << std::setw(2) << i<< " [Track]";
    HB1( 100*i+12, title1.str().c_str(), 200, -10, 150 );
    title2 << "drift length for LayerId = " << std::setw(2) << i<< " [Track]";
    HB1( 100*i+13, title2.str().c_str(), 50, -0.5, 3);
    title3 << "Position BcOut" << std::setw(2) << i;
    HB1( 100*i+14, title3.str().c_str(), 100, -250., 250. ); 
    title4 << "Residual BcOut" << std::setw(2) << i;
    HB1( 100*i+15, title4.str().c_str(), 400, -2.0, 2.0 );
    title5 << "Resid%Pos BcOut" << std::setw(2) << i;
    HB2( 100*i+16, title5.str().c_str(), 250, -250., 250., 100, -1.0, 1.0 );
    title6 << "Y%Xcal BcOut" << std::setw(2) << i;
    HB2( 100*i+17, title6.str().c_str(), 100, -250., 250., 100, -250., 250. );
    title7 << "Res%dl BcOut" << std::setw(2) << i;
    HB2( 100*i+18, title7.str().c_str(), 100, -3., 3., 100, -1.0, 1.0 );
    title8 << "Hit Pos%Drift Time BcOut" << std::setw(2) << i;
    HB2( 100*i+19, title8.str().c_str(), 100, -5., 100., 100, -3, 3);

    title9 << "Drift Length%Drift Time BcOut" << std::setw(2) << i;
    HBProf( 100*i+20, title9.str().c_str(), 100, -5, 50, 0, 12);
    HB2( 100*i+22, title9.str().c_str(), 100, -5, 50, 50, 0, 12);

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
      HBProf( 100000*i+3000+j, title.str().c_str(), 101, -12.12, 12.12, -30,300);
      HB2( 100000*i+4000+j, title.str().c_str(), 100, -12, 12, 100, -30,300);
    }
    
    //     title9 << " w/o self plane"; 
    //     HB2( 100*i+23, title9.str().c_str(), 100, -5., 5., 150, -50, 100);
    //     HBProf( 100*i+24, title9.str().c_str(), 100, -50, 150, 0, 5);
    //     HB2( 100*i+25, title9.str().c_str(), 100, -50, 150, 100,0, 3.5);
    //     title3 << " w/o self plane"; 
    //     HB2( 100*i+26, title3.str().c_str(), 50, -100., 100., 50, -1.0, 1.0 );
  }

  // Tracking Histgrams SdcIn
  HB1( 10, "#Tracks SdcIn", 10, 0., 10. );
  HB1( 11, "#Hits of Track SdcIn", 15, 0., 15. );
  HB1( 12, "Chisqr SdcIn", 500, 0., 50. ); 
  HB1( 13, "LayerId SdcIn", 15, 0., 15. );
  HB1( 14, "X0 SdcIn", 400, -100., 100. ); 
  HB1( 15, "Y0 SdcIn", 400, -100., 100. );
  HB1( 16, "U0 SdcIn", 200, -0.20, 0.20 );
  HB1( 17, "V0 SdcIn", 200, -0.20, 0.20 );
  HB2( 18, "U0%X0 SdcIn", 100, -100., 100., 100, -0.20, 0.20 );
  HB2( 19, "V0%Y0 SdcIn", 100, -100., 100., 100, -0.20, 0.20 );
  HB2( 20, "X0%Y0 SdcIn", 100, -100., 100., 100, -100, 100 );

  HB1( 21, "Xtgt SdcIn", 400, -100., 100. ); 
  HB1( 22, "Ytgt SdcIn", 400, -100., 100. );
  HB1( 23, "Utgt SdcIn", 200, -0.20, 0.20 );
  HB1( 24, "Vtgt SdcIn", 200, -0.20, 0.20 );
  HB2( 25, "Utgt%Xtgt SdcIn", 100, -100., 100., 100, -0.20, 0.20 );
  HB2( 26, "Vtgt%Ytgt SdcIn", 100, -100., 100., 100, -0.20, 0.20 );
  HB2( 27, "Xtgt%Ytgt SdcIn", 100, -100., 100., 100, -100, 100 );

  for( int i=1; i<=NumOfLayersSdcIn; ++i ){
    std::ostringstream title0, title1, title2, title3, title4;
    std::ostringstream title5, title6, title7, title8, title9;
    std::ostringstream title10, title11, title12, title13, title14;
    
    title0 << "wire for LayerId = " << std::setw(2) << i<< " [Track]";
    HB1( 100*i+11, title0.str().c_str(), 100, 0., 100. );
    title1 << "drift time for LayerId = " << std::setw(2) << i<< " [Track]";
    HB1( 100*i+12, title1.str().c_str(), 200, -10, 150 );
    title2 << "drift length for LayerId = " << std::setw(2) << i<< " [Track]";
    HB1( 100*i+13, title2.str().c_str(), 50, -0.5, 3);
    title3 << "Position SdcIn" << std::setw(2) << i;
    HB1( 100*i+14, title3.str().c_str(), 100, -250., 250. ); 
    title4 << "Residual SdcIn" << std::setw(2) << i;
    HB1( 100*i+15, title4.str().c_str(), 400, -2.0, 2.0 );
    title5 << "Resid%Pos SdcIn" << std::setw(2) << i;
    HB2( 100*i+16, title5.str().c_str(), 250, -250., 250., 100, -1.0, 1.0 );
    title6 << "Y%Xcal SdcIn" << std::setw(2) << i;
    HB2( 100*i+17, title6.str().c_str(), 100, -250., 250., 100, -250., 250. );
    title7 << "Res%dl SdcIn" << std::setw(2) << i;
    HB2( 100*i+18, title7.str().c_str(), 100, -3., 3., 100, -1.0, 1.0 );
    title8 << "Hit Pos%Drift Time SdcIn" << std::setw(2) << i;
    HB2( 100*i+19, title8.str().c_str(), 100, -5., 100., 100, -3, 3);

    title9 << "Drift Length%Drift Time SdcIn" << std::setw(2) << i;
    HBProf( 100*i+20, title9.str().c_str(), 100, -5, 50, 0, 12);
    HB2( 100*i+22, title9.str().c_str(), 100, -5, 50, 50, 0, 12);

    title4 << " w/o self plane"; 
    HB1( 100*i+21, title4.str().c_str(), 200, -5.0, 5.0 );

    title10 << "Residual SdcIn (0<theta<15)" << std::setw(2) << i;
    HB1( 100*i+71, title10.str().c_str(), 200, -5.0, 5.0 );

    title11 << "Residual SdcIn (15<theta<30)" << std::setw(2) << i;
    HB1( 100*i+72, title11.str().c_str(), 200, -5.0, 5.0 );

    title12 << "Residual SdcIn (30<theta<45)" << std::setw(2) << i;
    HB1( 100*i+73, title12.str().c_str(), 200, -5.0, 5.0 );

    title13 << "Residual SdcIn (45<theta)" << std::setw(2) << i;
    HB1( 100*i+74, title13.str().c_str(), 200, -5.0, 5.0 );

    for (int j=1; j<=64; j++) {
      std::ostringstream title;
      title << "XT of Layer " << std::setw(2) << i << " Wire #"<< std::setw(4) << j;
      HBProf( 100000*i+3000+j, title.str().c_str(), 101, -12.12, 12.12, -30,300);
      HB2( 100000*i+4000+j, title.str().c_str(), 100, -12, 12, 100, -30,300);
    }
    
    //     title9 << " w/o self plane"; 
    //     HB2( 100*i+23, title9.str().c_str(), 100, -5., 5., 150, -50, 100);
    //     HBProf( 100*i+24, title9.str().c_str(), 100, -50, 150, 0, 5);
    //     HB2( 100*i+25, title9.str().c_str(), 100, -50, 150, 100,0, 3.5);
    //     title3 << " w/o self plane"; 
    //     HB2( 100*i+26, title3.str().c_str(), 50, -100., 100., 50, -1.0, 1.0 );
  }

  
  // BcOut + SdcIn combined straight tracking (alignment)
  HB1( 90, "#Tracks BcOutSdcIn", 10, 0., 10. );
  HB1( 91, "#Hits of Track BcOutSdcIn", 25, 0., 25. );
  HB1( 92, "Chisqr BcOutSdcIn", 500, 0., 50. );
  for( int j=0; j<22; ++j ){
    // SdcIn layers 1-10, BcOut layers 51-62
    int L = (j<10) ? j+1 : j+41;
    std::ostringstream t1, t2;
    // the combined-fit residual in the x view is several mm wide
    // (multiple scattering at the target x defocused beam): the
    // range must not clip it, or the mean becomes meaningless
    t1 << "Residual BcOutSdcIn layer " << std::setw(2) << L;
    HB1( 100*L+25, t1.str().c_str(), 500, -5.0, 5.0 );
    t2 << "Res%dl BcOutSdcIn layer " << std::setw(2) << L;
    HB2( 100*L+26, t2.str().c_str(), 120, -3., 3., 100, -5.0, 5.0 );
    std::ostringstream t3;
    t3 << "Res%slope BcOutSdcIn layer " << std::setw(2) << L;
    HB2( 100*L+27, t3.str().c_str(), 100, -0.1, 0.1, 100, -5.0, 5.0 );
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
  tree->Branch("Bh2Seg",   event.Bh2Seg,  "Bh2Seg[nhBh2]/D");
  tree->Branch("tBh2",     event.tBh2,    "tBh2[nhBh2]/D");
  tree->Branch("deBh2",    event.deBh2,   "deBh2[nhBh2]/D");

  tree->Branch("nhBh1",   &event.nhBh1,   "nhBh1/I");
  tree->Branch("Bh1Seg",   event.Bh1Seg,  "Bh1Seg[nhBh1]/D");
  tree->Branch("tBh1",     event.tBh1,    "tBh1[nhBh1]/D");
  tree->Branch("deBh1",    event.deBh1,   "deBh1[nhBh1]/D");
  tree->Branch("btof",     event.btof,    "btof[nhBh1]/D");

  tree->Branch("ntrackb",   &event.ntrackb,   "ntrackb/I");
  tree->Branch("chosqrb",    event.chisqrb,   "chisqrb[ntrackb]/D");
  tree->Branch("x0b",        event.x0b,       "x0b[ntrackb]/D");
  tree->Branch("y0b",        event.y0b,       "y0b[ntrackb]/D");
  tree->Branch("u0b",        event.u0b,       "u0b[ntrackb]/D");
  tree->Branch("v0b",        event.v0b,       "v0b[ntrackb]/D");

  tree->Branch("ntracks",   &event.ntracks,   "ntracks/I");
  tree->Branch("chosqrs",    event.chisqrs,   "chisqrs[ntracks]/D");
  tree->Branch("x0s",        event.x0s,       "x0s[ntracks]/D");
  tree->Branch("y0s",        event.y0s,       "y0s[ntracks]/D");
  tree->Branch("u0s",        event.u0s,       "u0s[ntracks]/D");
  tree->Branch("v0s",        event.v0s,       "v0s[ntracks]/D");

  tree->Branch("ntrackc",   &event.ntrackc,   "ntrackc/I");
  tree->Branch("chisqrc",    event.chisqrc,   "chisqrc[ntrackc]/D");

  return true;
}

bool ConfMan::InitializeParameterFiles( void )
{
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

  return true;
}
