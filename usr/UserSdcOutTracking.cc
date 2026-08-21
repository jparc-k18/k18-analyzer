/*
  UserSdcOutTracking.cc

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

const double MaxMultiHitSdcOut = 3.0;

const double TdcLow  =   0.;
const double TdcHigh = 2000.;

#define HodoCut 0
#define SdcOutTracking 1
#define TDCrearrange 0

 VEvent::VEvent()
 {
 }

 VEvent::~VEvent()
 {
 }

class EventSdcOutTracking 
  : public VEvent
{

private:
  RawData *rawData;
  DCAnalyzer *DCAna;
  HodoAnalyzer *hodoAna;

public:
  EventSdcOutTracking();
  ~EventSdcOutTracking();

  bool ProcessingBegin();
  bool ProcessingEnd();
  bool ProcessingNormal();
  bool InitializeHistograms();
  void InitializeEvent();
};

EventSdcOutTracking::EventSdcOutTracking()
  : VEvent(),
    rawData(0), 
    DCAna(new DCAnalyzer()),
    hodoAna(new HodoAnalyzer())
{
}

EventSdcOutTracking::~EventSdcOutTracking()
{
  if (DCAna)   delete DCAna;
  if (hodoAna) delete hodoAna;
  if (rawData) delete rawData;
}

bool EventSdcOutTracking::ProcessingBegin()
{
 return true;
}

#ifndef MaxHits 
#define MaxHits 20
#endif
// For Dst
struct Dst{
  int ntrack;
  double chisqr[MaxHits];
  double x0[MaxHits];
  double y0[MaxHits];
  double u0[MaxHits];
  double v0[MaxHits];
};

//TDC rearrange
double diff_track_seg[MaxHits];

//For Tree
struct Event{
  int runnum;
  int evnum;
  int trigpat[NumOfMisc];
  int trigflag[NumOfMisc];

  double btof[MaxHits];

  int nhBh2;
  double Bh2Seg[MaxHits];
  double tBh2[MaxHits];
  double deBh2[MaxHits];

  int nhBh1;
  double Bh1Seg[MaxHits];
  double tBh1[MaxHits];
  double deBh1[MaxHits];

  int nhTof;
  double TofSeg[MaxHits];
  double tTof[MaxHits];
  double dtTof[MaxHits];
  double deTof[MaxHits];

  int nhLc;
  double LcSeg[MaxHits];
  double tLc[MaxHits];
  double dtLc[MaxHits];
  double deLc[MaxHits];

  int nhAc;
  int AcSeg[NumOfSegAC];
  double tAc[NumOfSegAC];
  double deAc[NumOfSegAC];
  double npeAc;

  int ntrack;
  double chisqr[MaxHits];  
  double x0[MaxHits];  
  double y0[MaxHits];  
  double u0[MaxHits];  
  double v0[MaxHits];
  double diff_track_seg[MaxHits];
  double res31[MaxHits];
  double res31080[MaxHits];
  double xtrackingtof[MaxHits];
  double ytrackingtof[MaxHits];

  int nhTofmt;
  int TofmtSeg[NumOfSegTOF];
  double tofmtt[NumOfSegTOF];
  double tofmtT[NumOfSegTOF];
};
static Event event;
static Dst dst;

// BH2 Cut Parameters
const double MinDeltaEBH2 = 0.5; //for Pion
const double MaxDeltaEBH2 = 3.0; //for Pion

// BH1 Cut Parameters
const double MinDeltaEBH1 = 0.5; //for Pion
const double MaxDeltaEBH1 = 3.0; //for Pion

// BH1-BH2
const double MinBeamToF  = -1.0;
const double MaxBeamToF  =  1.0;

// TOF Cut Parameters
const double MinDeltaETof = 0.40; //for Kaon
const double MaxDeltaETof = 4.00; //for Kaon
const double MinTimeTof =  15.0;
const double MaxTimeTof =  35.0;

// LC Cut Parameters
const double MinDeltaELc = 0.35; //for Kaon
const double MaxDeltaELc = 4.00; //for Kaon
const double MinTimeLc =  20.0;
const double MaxTimeLc =  40.0;

// TOF-LC  
const double MinToFT2L =  -2.0;
const double MaxToFT2L =  15.0;
const int MinToF2LSeg  =  -5;
const int MaxToF2LSeg  =  20;

bool EventSdcOutTracking::ProcessingNormal()
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

  HF1( 1, 0. );

  //------------------------Cut
  //if (trig_type != PIK_TRIG )  return true;
  //if (trig_type != BEAM_TRIG )  return true;
  //if (trig_type != PIPI_TRIG )  return true;
  HF1( 1, 1. );

#if HodoCut 
  /*
  //////////////BH2 Analysis
  hodoAna->DecodeBH2Hits(rawData);
  int ncBh2=hodoAna->GetNClustersBH2();
  event.nhBh2=ncBh2;
  //if( ncBh2==0 ) return true;

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
    //if( ncOk<1 ) return true;
  }

  HF1( 1, 2. );
  */
  //////////////Tof Analysis
  hodoAna->DecodeTOFHits(rawData);
  int ncTof=hodoAna->GetNClustersTOF();
  event.nhTof=ncTof;
  {
    int ncOk=0;
    for( int i=0; i<ncTof; ++i ){
      HodoCluster *cl=hodoAna->GetClusterTOF(i);
      double de=cl->DeltaE();
      double t=cl->CMeanTime(); 
      double dt=cl->TimeDif();
      event.TofSeg[i]=cl->MeanSeg()+1;
      event.tTof[i]=t;
      event.dtTof[i]=dt;
      event.deTof[i]=de;
      //------------------------Cut
      if( MinDeltaETof<de && de<MaxDeltaETof &&
	  MinTimeTof  <t  && t< MaxTimeTof ){
	++ncOk;
      }	
      else{
	//cl->GoodForAnalysis( false );
      }
    }
    //HF1( 311, double(ncOk) );
    //------------------------Cut
    if( ncOk<1 ) {
      //return true;
    }
  }

  HF1( 1, 3. );
  /*
  //////////////LC Analysis
  hodoAna->DecodeLCHits(rawData);
  int ncLc=hodoAna->GetNClustersLC();
  event.nhLc=ncLc;
  {
    int ncOk=0;
    for( int i=0; i<ncLc; ++i ){
      HodoCluster *cl=hodoAna->GetClusterLC(i);
      double de=cl->DeltaE();
      double t=cl->CMeanTime()-time0;
      double dt=cl->TimeDif();
      event.LcSeg[i]=cl->MeanSeg()+1;
      event.tLc[i]=t;
      event.dtLc[i]=dt;
      event.deLc[i]=de;
      //------------------------Cut
      if( MinDeltaELc<de && de<MaxDeltaELc &&
	  MinTimeLc < t  && t <MaxTimeLc ){
	++ncOk;
      }	
      else{
	//cl->GoodForAnalysis( false );
      }
    }
    //HF1( 411, double(ncOk) );
    //------------------------Cut
    if( ncOk<1 ) {
      //return true;
    }
  }

  HF1( 1, 4. );

  //////////////Tof-Lc
  /*
  {
    int nc=0;
    for( int iTof=0; iTof<ncTof; ++iTof ){
      HodoCluster *clTof=hodoAna->GetClusterTOF(iTof);
      //if( !clTof || !clTof->GoodForAnalysis() ) continue;
      double ttof=clTof->CMeanTime()-time0,
	segTof=clTof->MeanSeg()+1;
      for( int iLc=0; iLc<ncLc; ++iLc ){
	HodoCluster *clLc=hodoAna->GetClusterLC(iLc);
	//if( !clLc || !clLc->GoodForAnalysis() ) continue;
	double tlc=clLc->CMeanTime()-time0,
	  segLc=clLc->MeanSeg()+1;
	//------------------------Cut
	if( MinToF2LSeg<=(segTof-segLc) && (segTof-segLc)<=MaxToF2LSeg &&
	    MinToFT2L<(tlc-ttof) && (tlc-ttof)<MaxToFT2L ){
	  ++nc;
	}
	else{
	  //clTof->GoodForAnalysis( false );
	  //clLc->GoodForAnalysis( false );
	}
      }
    }
    //------------------------Cut
    if( nc<1 ) {
      //return true;
    }
  }

  HF1( 1, 5. );
  */
  
  // AC
  /*
  hodoAna->DecodeACHits(rawData);
  int nh=hodoAna->GetNHitsAC();
  double npeac;
  for( int i=0; i<nh; ++i ){
    Hodo1Hit *hit=hodoAna->GetHitAC(i);
    if(!hit) continue;
    int seg=hit->SegmentId()+1;
    double a=hit->GetA(), t=hit->GetT(), ct=hit->GetCT();
    
    event.AcSeg[i] = seg;
    event.deAc[i] = a;
    event.tAc[i] = t;
    if( t >0 ) npeac += a;
  }
  event.npeAc = npeac;

  */  
  HF1( 1, 6. );
#endif
  
  HF1( 1, 10. );

  //////////////SDC3&4 number of hit layer
  DCAna->DecodeSdcOutHits( rawData ); 

  int IdTof = DCGeomMan::GetInstance().GetTofId();
  double zTof = DCGeomMan::GetInstance().GetLocalZ( IdTof );

  //SDC3&SDC4
  double multi_SdcOut=0.;
  {
    for( int layer=1; layer<=NumOfLayersSdcOut; ++layer ){
      const DCHitContainer &contOut =DCAna->GetSdcOutHC(layer);
      int nhOut=contOut.size();
      multi_SdcOut += double(nhOut);
      HF1( 100*layer, nhOut );
      for( int i=0; i<nhOut; ++i ){
	DCHit *hit=contOut[i];
	double wire=hit->GetWire();

	HF1( 100*layer+1, wire-0.5 );
	int nhtdc = hit->GetTdcSize();

	//std::cout << "nhtdc=" << nhtdc << std::endl;
	for( int k=0; k<nhtdc; k++ ){
	  int tdc = hit->GetTdcVal(k);
	  HF1( 100*layer+2, tdc );
	  HF1( 10000*layer+int(wire), tdc );
	}
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
	  HF1( 10000*layer+2000+int(wire), dl );
	}
      }
    }    
  } 
  //if( multi_SdcOut/double(NumOfLayersSdcOut) > MaxMultiHitSdcOut )
  //return true;

  HF1( 1, 11. );

#if SdcOutTracking
  // Sdc Out
  //std::cout << "==========TrackSearch SdcOut============" << std::endl;
  DCAna->TrackSearchSdcOut();
#endif 

  int tdc_rearrange_flag = 0;
#if TDCrearrange
  tdc_rearrange_flag = 1;
  {
    const HodoRHitContainer &cont = rawData->GetTOFMTRawHC();
    int nh = cont.size();
    event.nhTofmt = nh;
    HF1(1300, double(nh));
    int nh1=0, nh2=0;
    for( int i=0; i<nh; ++i){
      HodoRawHit *hit=cont[i];
      int seg=hit->SegmentId()+1;
      int T= hit->GetTdcUp();
      //std::cout<<"tofmt tdc:   "<<T<<std::endl;
      //std::cout<<"tofmt seg:   "<<seg<<std::endl;
      HF1(130000+seg, double(T));
      event.TofmtSeg[i] = seg;
      HF1(130000+seg, double(T));
      HF1(1302, double(T));
      HF1(1301, seg-0.5);
      event.tofmtt[i] = T;
    }
  }

  hodoAna->DecodeTOFMTHits( rawData );
  {
    int nh=hodoAna->GetNHitsTOFMT();
    //event.nhTOFMT=nhTOFMT;
    {
      for( int i=0; i<nh; ++i ){
	Hodo1Hit *cl=hodoAna->GetHitTOFMT(i);
	int seg=cl->SegmentId()+1;
	double t=cl->GetT();
	//std::cout<<"tofmt time:   "<<t<<std::endl;
	HF1(1303, t);
	HF1(131000+seg, t);
        event.tofmtT[i]=t;
	}	
    }
  }


  bool rearrange_status = DCAna->TDCrearrangeSDCout(DCAna, hodoAna);
  multi_SdcOut=0.;
  {
    for( int layer=1; layer<=NumOfLayersSdcOut; ++layer ){
      const DCHitContainer &contOut =DCAna->GetSdcOutHC(layer);
      int nhOut=contOut.size();
      multi_SdcOut += double(nhOut);
      for( int i=0; i<nhOut; ++i ){
	DCHit *hit=contOut[i];
	double wire=hit->GetWire();
	int nhdt = hit->GetDriftTimeSize();
	for( int k=0; k<nhdt; k++ ){
	  double dt = hit->GetDriftTime(k);
	  HF1( 100*layer+41, dt );  
	  HF1( 10000*layer+6000+int(wire), dt );
	}
	int nhdl = hit->GetDriftTimeSize();
	for( int k=0; k<nhdl; k++ ){
	  double dl = hit->GetDriftLength(k);
	  HF1( 100*layer+42, dl );	 
	  HF1( 10000*layer+7000+int(wire), dl );
	}
      }
    }    
  }
  if(rearrange_status) 
    {
      DCAna->TrackSearchSdcOut();
      std::cout<<"###########TDC is Rearranged#############"<<std::endl;
    }
  else std::cout<<"#########No rearrange###############"<<std::endl;

#endif

#if SdcOutTracking
  int nt=DCAna->GetNtracksSdcOut();
  if( nt>MaxHits ) nt=MaxHits;
  event.ntrack=nt;
  HF1( 10, double(nt) );
  for( int it=0; it<nt; ++it ){
    DCLocalTrack *tp=DCAna->GetTrackSdcOut(it);
    int nh=tp->GetNHit();
    double chisqr=tp->GetChiSquare();
    double x0=tp->GetX0(), y0=tp->GetY0();
    double u0=tp->GetU0(), v0=tp->GetV0();

    double cost = 1./sqrt(1.+u0*u0+v0*v0);
    double theta = acos(cost)*Rad2Deg;

    event.chisqr[it]=chisqr;
    event.x0[it]=x0;
    event.y0[it]=y0;
    event.u0[it]=u0;
    event.v0[it]=v0;
    if(tdc_rearrange_flag==1){
      event.diff_track_seg[it] = DCAna->GetDiffTrackSeg(it);}
    HF1( 11, double(nh) );
    HF1( 12, chisqr ); 
    HF1( 14, x0 ); HF1( 15, y0 );
    HF1( 16, u0 ); HF1( 17, v0 );
    HF2( 18, x0, u0 ); HF2( 19, y0, v0 );
    HF2( 20, x0, y0 );
    
    double xtof=tp->GetX( zTof ), ytof=tp->GetY( zTof );
    double utof=u0, vtof=v0;
    HF1( 21, xtof ); HF1( 22, ytof );
    HF1( 23, utof ); HF1( 24, vtof );
    HF2( 25, xtof, utof ); HF2( 26, ytof, vtof );
    HF2( 27, xtof, ytof );
    event.xtrackingtof[it] = xtof;
    event.ytrackingtof[it] = ytof;


    for( int ih=0; ih<nh; ++ih ){
      DCLTrackHit *hit=tp->GetHit(ih);
      int layerId=hit->GetLayer()-30;  

      HF1( 13, layerId );
      double wire=hit->GetWire();
	//if(layerId == 3&&wire==1.)
	//std::cout<<"go"<<std::endl; 
      double dt=hit->GetDriftTime(), dl=hit->GetDriftLength();
      HF1( 100*layerId+11, wire-0.5 );
      HF1( 100*layerId+12, dt );
      HF1( 100*layerId+13, dl );


      if(nh>=9){
	HF1 (100*layerId + 31, wire-0.5);
	HF1 (100*layerId + 32, dt);
	HF1 (100*layerId + 33, dl);
	HF1 (10000*layerId + 5000 + int(wire), dt);
	}	



      double xcal=hit->GetXcal(), ycal=hit->GetYcal();
      double pos=hit->GetLocalHitPos(), res=hit->GetResidual();
      HF1( 100*layerId+14, pos );
      HF1( 100*layerId+15, res );
      if(layerId==1)
	{
	  event.res31[it]=res;
	  //std::cout<<"I found residula"<<std::endl;
	  if(wire==80.)
	    {
	      event.res31080[it]=res;
	    }
	}
      HF2( 100*layerId+16, pos, res );
      HF2( 100*layerId+17, xcal, ycal);
      HF2( 100*layerId+35, wire-0.5, res );
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
#endif

  HF1( 1, 12. );

  tree->Fill();

  return true;
}

void EventSdcOutTracking::InitializeEvent( void )
{
  event.runnum = -1;
  event.evnum = -1;
  // size branches must be 0 (not -1): TTree::Fill crashes on
  // negative array sizes when the fill section is skipped
  event.ntrack   = 0;
  event.nhBh2  = 0;
  event.nhBh1  = 0;
  event.nhTof  = 0;
  event.nhLc   = 0;
  event.nhAc   = 0;
  event.npeAc  = -9999.0;
  event.nhTofmt = 0;

  for( int it=0; it<MaxHits; it++){
    event.Bh2Seg[it] = -1;
    event.tBh2[it] = -9999.0;
    event.deBh2[it] = -9999.0;

    event.Bh1Seg[it] = -1;
    event.tBh1[it] = -9999.0;
    event.deBh1[it] = -9999.0;
    event.btof[it]  = -1;

    event.TofSeg[it] = -1;
    event.tTof[it] = -9999.0;
    event.dtTof[it] = -9999.0;
    event.deTof[it] = -9999.0;

    event.LcSeg[it] = -1;
    event.tLc[it] = -9999.0;
    event.dtLc[it] = -9999.0;
    event.deLc[it] = -9999.0;

    event.AcSeg[it] = -1;
    event.tAc[it] = -9999.0;
    event.deAc[it] = -9999.0;
    event.TofmtSeg[it] = -1;
    event.tofmtt[it] = -9999.;
    event.tofmtT[it] = -9999.;
    event.xtrackingtof[it] = -9999.;
    event.ytrackingtof[it] = -9999.;
  }
    
  for( int it=0; it<MaxHits; it++){
    event.chisqr[it] = -1.0;
    event.x0[it] = -9999.0;
    event.y0[it] = -9999.0;
    event.u0[it] = -9999.0;
    event.v0[it] = -9999.0;
    event.diff_track_seg[it] = -9999.0;
  }

  for( int it=0; it<NumOfMisc; it++){
    event.trigpat[it] = -1;
    event.trigflag[it] = -1;
  }
}

bool EventSdcOutTracking::ProcessingEnd()
{
  return true;
}

VEvent *ConfMan::EventAllocator()
{
  return new EventSdcOutTracking;
}


const int NbinSdcOutTdc   = 2000;
const double MinSdcOutTdc =    0.;
const double MaxSdcOutTdc = 2000.;

// 1 ns/bin; the SDC3/4 TDC LSB is ~0.6 ns with wire-by-wire
// offsets, so the layer-summed drift time has no comb structure
const int NbinSdcOutDT   = 440;
const double MinSdcOutDT =  -40.;
const double MaxSdcOutDT =  400.;

// drift length is clamped to [0, 10] mm (DriftLength6)
const int NbinSdcOutDL   = 130;
const double MinSdcOutDL =  -1.;
const double MaxSdcOutDL =  12.;

bool ConfMan:: InitializeHistograms()
{  
  HB1( 1, "Status", 20, 0., 20. );

  //***********************Chamber
  // SDC3
  for( int i=1; i<=NumOfLayersSdc; ++i ){
    std::ostringstream title1, title2, title3, title4, title5, title41, title42;
    title1 << "#Hits SDC3#" << std::setw(2) << i;
    title2 << "Hitpat SDC3#" << std::setw(2) << i;
    title3 << "Tdc SDC3#" << std::setw(2) << i;
    title4 << "Drift Time SDC3#" << std::setw(2) << i;
    title5 << "Drift Length SDC3#" << std::setw(2) << i;
    title41 << "Drift Time after rearrange SDC3#" << std::setw(2) << i;
    title42 << "Drift Length after rearrange SDC3#" << std::setw(2) << i;

    HB1( 100*i+0, title1.str().c_str(), MaxWireSDC3+1, 0., double(MaxWireSDC3+1) );
    HB1( 100*i+1, title2.str().c_str(), MaxWireSDC3+1, 0., double(MaxWireSDC3+1) );
    HB1( 100*i+2, title3.str().c_str(), NbinSdcOutTdc, MinSdcOutTdc, MaxSdcOutTdc );
    HB1( 100*i+3, title4.str().c_str(), NbinSdcOutDT, MinSdcOutDT, MaxSdcOutDT );
    HB1( 100*i+4, title5.str().c_str(), NbinSdcOutDL, MinSdcOutDL, MaxSdcOutDL );
    HB1( 100*i+41, title41.str().c_str(), NbinSdcOutDT, MinSdcOutDT, MaxSdcOutDT );
    HB1( 100*i+42, title42.str().c_str(), NbinSdcOutDL, MinSdcOutDL, MaxSdcOutDL );
    
    int NumOfWires;
    if (i==2 || i==5)
      NumOfWires = MaxWireSDC3X;
    else 
      NumOfWires = MaxWireSDC3V;
    for (int wire=1; wire<=NumOfWires; wire++) {
      std::ostringstream title11, title12, title13, title100, title14, title15;
      title11 << "Tdc SDC3#" << std::setw(2) << i << " Wire#" << wire;
      HB1( 10000*i+wire, title11.str().c_str(), NbinSdcOutTdc, MinSdcOutTdc, MaxSdcOutTdc );
      title12 << "Drift Time SDC3#" << std::setw(2) << i << " Wire#" << wire;
      HB1( 10000*i+1000+wire, title12.str().c_str(), NbinSdcOutDT, MinSdcOutDT, MaxSdcOutDT );
      title13 << "Drift Length SDC3#" << std::setw(2) << i << " Wire#" << wire;
      HB1( 10000*i+2000+wire, title13.str().c_str(), NbinSdcOutDL, MinSdcOutDL, MaxSdcOutDL );
      title14 << "Drift Time after rearrange SDC3#" << std::setw(2) << i << " Wire#" << wire;
      HB1( 10000*i+6000+wire, title14.str().c_str(), NbinSdcOutDT, MinSdcOutDT, MaxSdcOutDT );
      title15 << "Drift Length after rearrange SDC3#" << std::setw(2) << i << " Wire#" << wire;
      HB1( 10000*i+7000+wire, title15.str().c_str(), NbinSdcOutDL, MinSdcOutDL, MaxSdcOutDL );
      title100 << "Drift Time (Tracking) SDC3#" << std::setw(2) << i << " Wire#" << wire;
      HB1( 10000*i+5000+wire, title100.str().c_str(), NbinSdcOutDT, MinSdcOutDT, MaxSdcOutDT );

    }

  }

  // SDC4
  for( int i=1; i<=NumOfLayersSdc; ++i ){
    std::ostringstream title1, title2, title3, title4, title5, title41, title42;
    title1 << "#Hits SDC4#" << std::setw(2) << i;
    title2 << "Hitpat SDC4#" << std::setw(2) << i;
    title3 << "Tdc SDC4#" << std::setw(2) << i;
    title4 << "Drift Time SDC4#" << std::setw(2) << i;
    title5 << "Drift Length SDC4#" << std::setw(2) << i;
    title41 << "Drift Time after rearrage SDC4#" << std::setw(2) << i;
    title42 << "Drift Length after rearrage SDC4#" << std::setw(2) << i;
    HB1( 100*(i+6)+0, title1.str().c_str(), MaxWireSDC4+1, 0., double(MaxWireSDC4+1) );
    HB1( 100*(i+6)+1, title2.str().c_str(), MaxWireSDC4+1, 0., double(MaxWireSDC4+1) );
    HB1( 100*(i+6)+2, title3.str().c_str(), NbinSdcOutTdc, MinSdcOutTdc, MaxSdcOutTdc );
    HB1( 100*(i+6)+3, title4.str().c_str(), NbinSdcOutDT, MinSdcOutDT, MaxSdcOutDT );
    HB1( 100*(i+6)+4, title5.str().c_str(), NbinSdcOutDL, MinSdcOutDL, MaxSdcOutDL );
    HB1( 100*(i+6)+41, title41.str().c_str(), NbinSdcOutDT, MinSdcOutDT, MaxSdcOutDT );
    HB1( 100*(i+6)+42, title42.str().c_str(), NbinSdcOutDL, MinSdcOutDL, MaxSdcOutDL );
    

    int NumOfWires;
    if (i==2 || i==5)
      NumOfWires = MaxWireSDC4X;
    else 
      NumOfWires = MaxWireSDC4V;
    for (int wire=1; wire<=NumOfWires; wire++) {
      std::ostringstream title11, title12, title13, title101, title14, title15;
      title11 << "Tdc SDC4#" << std::setw(2) << i << " Wire#" << wire;
      HB1( 10000*(i+6)+wire, title11.str().c_str(), NbinSdcOutTdc, MinSdcOutTdc, MaxSdcOutTdc );
      title12 << "Drift Time SDC4#" << std::setw(2) << i << " Wire#" << wire;
      HB1( 10000*(i+6)+1000+wire, title12.str().c_str(), NbinSdcOutDT, MinSdcOutDT, MaxSdcOutDT );
      title13 << "Drift Length SDC4#" << std::setw(2) << i << " Wire#" << wire;
      HB1( 10000*(i+6)+2000+wire, title13.str().c_str(), NbinSdcOutDL, MinSdcOutDL, MaxSdcOutDL );
      title14 << "Drift Time after rearrange SDC4#" << std::setw(2) << i << " Wire#" << wire;
      HB1( 10000*(i+6)+6000+wire, title14.str().c_str(), NbinSdcOutDT, MinSdcOutDT, MaxSdcOutDT );
      title15 << "Drift Length after rearrange SDC4#" << std::setw(2) << i << " Wire#" << wire;
      HB1( 10000*(i+6)+7000+wire, title15.str().c_str(), NbinSdcOutDL, MinSdcOutDL, MaxSdcOutDL );
      title101 << "Drift Time (Tracking) SDC4#" << std::setw(2) << i << " Wire#" << wire;
      HB1( 10000*(i+6)+5000+wire, title101.str().c_str(), NbinSdcOutDT, MinSdcOutDT, MaxSdcOutDT );
    }
  }

  //TOFMT Histograms
  {
    HB1( 1300, "#Hits TOFMT",32,0.,32.);
    HB1( 1301, "Hit Pattern TOFMT",32,0.,32.);
    HB1( 1302, "TDC TOFMT", NbinSdcOutTdc, MinSdcOutTdc, MaxSdcOutTdc);
    HB1( 1303, "Time TOFMT", 2000,-1000., 1000.);

    for( int i =1 ; i<33 ; i++)
      {
	std::ostringstream title1, title2;
	title1 << "Tdc TOFMT Seg#" << std::setw(2)<<i;
	HB1 (130000+i, title1.str().c_str(),NbinSdcOutTdc, MinSdcOutTdc, MaxSdcOutTdc);
	title2 << "Time TOFMT seg#" << std::setw(2)<<i;
	HB1 (130000+1000+i, title2.str().c_str(),2000,-1000.,1000.);
      }
  }

   // Tracking Histgrams
  HB1( 10, "#Tracks SdcOut", 10, 0., 10. );
  HB1( 11, "#Hits of Track SdcOut", 20, 0., 20. );
  HB1( 12, "Chisqr SdcOut", 500, 0., 50. ); 
  HB1( 13, "LayerId SdcOut", 20, 0., 13. );
  HB1( 14, "X0 SdcOut", 1400, -1200., 1200. ); 
  HB1( 15, "Y0 SdcOut", 1000, -500., 500. );
  HB1( 16, "U0 SdcOut", 200, -0.35, 0.35 );
  HB1( 17, "V0 SdcOut", 200, -0.20, 0.20 );
  HB2( 18, "U0%X0 SdcOut", 120, -1200., 1200., 100, -0.40, 0.40 );
  HB2( 19, "V0%Y0 SdcOut", 100, -500., 500., 100, -0.20, 0.20 );
  HB2( 20, "X0%Y0 SdcOut", 100, -1200., 1200., 100, -500, 500 );

  HB1( 21, "Xtof SdcOut", 1400, -1200., 1200. ); 
  HB1( 22, "Ytof SdcOut", 1000, -500., 500. );
  HB1( 23, "Utof SdcOut", 200, -0.35, 0.35 );
  HB1( 24, "Vtof SdcOut", 200, -0.20, 0.20 );
  HB2( 25, "Utof%Xtof SdcOut", 100, -1200., 1200., 100, -0.40, 0.40 );
  HB2( 26, "Vtof%Ytof SdcOut", 100, -500., 500., 100, -0.20, 0.20 );
  HB2( 27, "Xtof%Ytof SdcOut", 100, -1200., 1200., 100, -500, 500 );

  for( int i=1; i<=NumOfLayersSdcOut; ++i ){
    std::ostringstream title0, title1, title2, title3, title4, title41, title42;
    std::ostringstream title5, title6, title7, title8, title9;
    std::ostringstream title10, title11, title12, title13, title14, title102, title103, title104, title105, title106;
    
    title0 << "wire for LayerId = " << std::setw(2) << i<< " [Track]";
    HB1( 100*i+11, title0.str().c_str(), 120, 0., 120. );
    title1 << "drift time for LayerId = " << std::setw(2) << i<< " [Track]";
    HB1( 100*i+12, title1.str().c_str(), NbinSdcOutDT, MinSdcOutDT, MaxSdcOutDT );
    title2 << "drift length for LayerId = " << std::setw(2) << i<< " [Track]";
    HB1( 100*i+13, title2.str().c_str(), NbinSdcOutDL, MinSdcOutDL, MaxSdcOutDL );
    title3 << "Position SdcOut" << std::setw(2) << i;
    HB1( 100*i+14, title3.str().c_str(), 100, -1000., 1000. );
    title4 << "Residual SdcOut" << std::setw(2) << i;
    HB1( 100*i+15, title4.str().c_str(), 400, -2.0, 2.0 );
    title5 << "Resid%Pos SdcOut" << std::setw(2) << i;
    HB2( 100*i+16, title5.str().c_str(), 400, -1000., 1000., 100, -1.0, 1.0 );
    title6 << "Y%Xcal SdcOut" << std::setw(2) << i;
    HB2( 100*i+17, title6.str().c_str(), 100, -1000., 1000., 100, -1000., 1000. );
    title7 << "Res%dl SdcOut" << std::setw(2) << i;
    HB2( 100*i+18, title7.str().c_str(), 220, -11., 11., 120, -1.5, 1.5 );
    title8 << "Hit Pos%Drift Time SdcOut" << std::setw(2) << i;
    HB2( 100*i+19, title8.str().c_str(), NbinSdcOutDT, MinSdcOutDT, MaxSdcOutDT, 100, -20, 20);
    title102 << "wire for LayerId = " << std::setw(2) << i << "[Track & nhit>9]";
    HB1( 100*i+31, title102.str().c_str(), 120, 0., 120.);
    title103 << "drift time for LayerId = " << std::setw(2) << i << "[Track & nhit>9]";
    HB1( 100*i+32, title103.str().c_str(), 500, -100, 550.);
    title104 << "drift length for LayerId = " << std::setw(2) << i << "[Track & nhit>9]";
    HB1( 100*i+33, title104.str().c_str(), 100, -5, 15);
    title105 << "Residual SdcOut = " << std::setw(2) << i << "[Track & nhit>9]";
    HB1( 100*i+34, title105.str().c_str(), 1000, -5., 5.);
    title106 << "Residual SdcOut each wire = " << std::setw(2) << i;
    HB2( 100*i+35, title106.str().c_str(), 120, 0., 120., 1000, -5., 5.);



    title9 << "Drift Length%Drift Time SdcOut" << std::setw(2) << i;
    HBProf( 100*i+20, title9.str().c_str(), 100, -50, 300, 0, 12);
    HB2( 100*i+22, title9.str().c_str(), 100, -50, 300, 100,0, 12);

    title4 << " w/o self plane"; 
    HB1( 100*i+21, title4.str().c_str(), 200, -5.0, 5.0 );

    title10 << "Residual SdcOut (0<theta<15)" << std::setw(2) << i;
    HB1( 100*i+71, title10.str().c_str(), 200, -5.0, 5.0 );

    title11 << "Residual SdcOut (15<theta<30)" << std::setw(2) << i;
    HB1( 100*i+72, title11.str().c_str(), 200, -5.0, 5.0 );

    title12 << "Residual SdcOut (30<theta<45)" << std::setw(2) << i;
    HB1( 100*i+73, title12.str().c_str(), 200, -5.0, 5.0 );

    title13 << "Residual SdcOut (45<theta)" << std::setw(2) << i;
    HB1( 100*i+74, title13.str().c_str(), 200, -5.0, 5.0 );

    for (int j=1; j<=120; j++) {
      std::ostringstream title;
      title << "XT of Layer " << std::setw(2) << i << " Wire #"<< std::setw(4) << j;
      HBProf( 100000*i+3000+j, title.str().c_str(), 101, -12.12, 12.12, -30,300);
      HB2( 100000*i+4000+j, title.str().c_str(), 100, -12, 12, 100, -30,300);
    }
    /*
    title9 << " w/o self plane"; 
    HB2( 100*i+23, title9.str().c_str(), 100, -5., 5., 150, -50, 100);
    HBProf( 100*i+24, title9.str().c_str(), 100, -50, 150, 0, 5);
    HB2( 100*i+25, title9.str().c_str(), 100, -50, 150, 100,0, 3.5);
    title3 << " w/o self plane"; 
    HB2( 100*i+26, title3.str().c_str(), 50, -100., 100., 50, -1.0, 1.0 );
    */
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

  tree->Branch("nhTof",   &event.nhTof,   "nhTof/I");
  tree->Branch("TofSeg",   event.TofSeg,  "TofSeg[nhTof]/D");
  tree->Branch("tTof",     event.tTof,    "tTof[nhTof]/D");
  tree->Branch("dtTof",    event.dtTof,   "dtTof[nhTof]/D");
  tree->Branch("deTof",    event.deTof,   "deTof[nhTof]/D");

  tree->Branch("nhLc",   &event.nhLc,   "nhLc/I");
  tree->Branch("LcSeg",   event.LcSeg,  "LcSeg[nhLc]/D");
  tree->Branch("tLc",     event.tLc,    "tLc[nhLc]/D");
  tree->Branch("dtLc",    event.dtLc,   "dtLc[nhLc]/D");
  tree->Branch("deLc",    event.deLc,   "deLc[nhLc]/D");

  tree->Branch("nhAc",   &event.nhAc,   "nhAc/I");
  tree->Branch("AcSeg",   event.AcSeg,  "AcSeg[nhAc]/I");
  tree->Branch("tAc",     event.tAc,    "tAc[nhAc]/D");
  tree->Branch("deAc",    event.deAc,   "deAc[nhAc]/D");
  tree->Branch("npeAc",  &event.npeAc,  "npeAc/D");

  tree->Branch("ntrack",   &event.ntrack,   "ntrack/I");
  tree->Branch("chisqr",    event.chisqr,   "chisqr[ntrack]/D");
  tree->Branch("x0",        event.x0,       "x0[ntrack]/D");
  tree->Branch("y0",        event.y0,       "y0[ntrack]/D");
  tree->Branch("u0",        event.u0,       "u0[ntrack]/D");
  tree->Branch("v0",        event.v0,       "v0[ntrack]/D");
  tree->Branch("diff_track_seg", event.diff_track_seg, "diff_track_seg[ntrack]/D");
  tree->Branch("res31", event.res31, "res31[ntrack]/D");
  tree->Branch("res31080", event.res31080, "res31080[ntrack]/D");
  tree->Branch("xtrackingtof", event.xtrackingtof, "xtrackingtof[ntrack]/D");
  tree->Branch("ytrackingtof", event.ytrackingtof, "ytrackingtof[ntrack]/D");

  tree->Branch("nhTofmt", &event.nhTofmt, "nhTofmt/I");
  tree->Branch("TofmtSeg", event.TofmtSeg, "TofmtSeg[nhTofmt]/I");
  tree->Branch("tofmtt", event.tofmtt, "tofmtt[nhTofmt]/D");
  tree->Branch("tofmtT", event.tofmtT, "tofmtT[nhTofmt]/D");
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
