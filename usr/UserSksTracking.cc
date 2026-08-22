/*
  UserSksTracking.cc

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
#include "Parameter.hh"

const int BEAM_TRIG  = 1;
const int PIK_TRIG   = 3;
const int PIPI_TRIG  = 11;
const int PIP_TRIG   = 12;
const int BH2_TRIG   = 4;

const double Deg2Rad = acos(-1.)/180.;
const double Rad2Deg = 180./acos(-1.);

#define In  1
#define Out 1
#define SKS 1

 VEvent::VEvent()
 {
 }

 VEvent::~VEvent()
 {
 }

class EventSksTracking 
  : public VEvent
{

private:
  RawData *rawData;
  DCAnalyzer *DCAna;
  HodoAnalyzer *hodoAna;

public:
  EventSksTracking();
  ~EventSksTracking();

  bool ProcessingBegin();
  bool ProcessingEnd();
  bool ProcessingNormal();
  bool InitializeHistograms();
  void InitializeEvent();
};

EventSksTracking::EventSksTracking()
  : VEvent(),
    rawData(0), 
    DCAna(new DCAnalyzer()),
    hodoAna(new HodoAnalyzer())
{
}

EventSksTracking::~EventSksTracking()
{
  if (DCAna)   delete DCAna;
  if (hodoAna) delete hodoAna;
  if (rawData) delete rawData;
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
  double deTof[MaxHits];

  int nhLc;
  double LcSeg[MaxHits];
  double tLc[MaxHits];
  double deLc[MaxHits];

  int ntIn;
  int nhIn[MaxHits]; 
  double chisqrIn[MaxHits];
  double u0In[MaxHits];
  double v0In[MaxHits];
  double x0In[MaxHits];
  double y0In[MaxHits];

  int ntOut;
  int nhOut[MaxHits]; 
  double chisqrOut[MaxHits];
  double u0Out[MaxHits];
  double v0Out[MaxHits];
  double x0Out[MaxHits];
  double y0Out[MaxHits];

  int ntSks;
  int nhSks[MaxHits]; 
  double chisqrSks[MaxHits];
  double path[MaxHits];
  double p[MaxHits];
  double m2[MaxHits];

  double xtgt[MaxHits];  
  double ytgt[MaxHits];  
  double utgt[MaxHits];  
  double vtgt[MaxHits]; 
  double theta[MaxHits]; 
  double phi[MaxHits]; 

  double layer1ResL[MaxHits];
  double layer2ResL[MaxHits];
  double layer3ResL[MaxHits];
  double layer4ResL[MaxHits];

  double layer5ResL[MaxHits];
  double layer6ResL[MaxHits];
  double layer7ResL[MaxHits];
  double layer8ResL[MaxHits];
  double layer9ResL[MaxHits];
  double layer10ResL[MaxHits];

  double layer31ResL[MaxHits];
  double layer32ResL[MaxHits];
  double layer33ResL[MaxHits];
  double layer34ResL[MaxHits];
  double layer35ResL[MaxHits];
  double layer36ResL[MaxHits];

  double layer37ResL[MaxHits];
  double layer38ResL[MaxHits];
  double layer39ResL[MaxHits];
  double layer40ResL[MaxHits];
  double layer41ResL[MaxHits];
  double layer42ResL[MaxHits];

  double layer1ResG[MaxHits];
  double layer2ResG[MaxHits];
  double layer3ResG[MaxHits];
  double layer4ResG[MaxHits];

  double layer5ResG[MaxHits];
  double layer6ResG[MaxHits];
  double layer7ResG[MaxHits];
  double layer8ResG[MaxHits];
  double layer9ResG[MaxHits];
  double layer10ResG[MaxHits];

  double layer31ResG[MaxHits];
  double layer32ResG[MaxHits];
  double layer33ResG[MaxHits];
  double layer34ResG[MaxHits];
  double layer35ResG[MaxHits];
  double layer36ResG[MaxHits];

  double layer37ResG[MaxHits];
  double layer38ResG[MaxHits];
  double layer39ResG[MaxHits];
  double layer40ResG[MaxHits];
  double layer41ResG[MaxHits];
  double layer42ResG[MaxHits];
};
static Event event;


bool EventSksTracking::ProcessingBegin()

{
  //  ConfMan::InitializeHistograms(); 
 return true;
}

bool EventSksTracking::ProcessingNormal()
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

  //////////////BH2 Analysis
  hodoAna->DecodeBH2Hits(rawData);
  int ncBh2=hodoAna->GetNClustersBH2();
  event.nhBh2=ncBh2;
  if( ncBh2==0 ){ tree->Fill(); return true; }

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
    if( ncOk<1 ){ tree->Fill(); return true; }
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

  //////////////Tof Analysis
  hodoAna->DecodeTOFHits(rawData);
  int ncTof=hodoAna->GetNClustersTOF();


    event.nhTof=ncTof;
  {
    int ncOk=0;
    for( int i=0; i<ncTof; ++i ){
      HodoCluster *cl=hodoAna->GetClusterTOF(i);
      double de=cl->DeltaE();
      double t=cl->CMeanTime()-time0;
      event.TofSeg[i]=cl->MeanSeg()+1;
      event.tTof[i]=t;
      event.deTof[i]=de;
      //------------------------Cut
      if( MinDeltaETof<de && de<MaxDeltaETof &&
	  MinTimeTof  <t  && t< MaxTimeTof ){
	++ncOk;
      }	
      else{
	//	cl->GoodForAnalysis( false );
      }
    }
    HF1( 311, double(ncOk) );
    //------------------------Cut
    if( ncOk<1 ) {
      //      return true;
    }
  }

  HF1( 1, 3. );

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
      event.LcSeg[i]=cl->MeanSeg()+1;
      event.tLc[i]=t;
      event.deLc[i]=de;
      //------------------------Cut
      if( MinDeltaELc<de && de<MaxDeltaELc &&
	  MinTimeLc < t  && t <MaxTimeLc ){
	++ncOk;
      }	
      else{
	//	cl->GoodForAnalysis( false );
      }
    }
    HF1( 411, double(ncOk) );
    //------------------------Cut
    if( ncOk<1 ) {
      //      return true;
    }
  }
  HF1( 1, 4. );

  //////////////Tof-Lc
  HF1( 1, 5. );

  HF1( 1, 10. );


  DCAna->DecodeRawHits( rawData );
  
  int IdTof = DCGeomMan::GetInstance().GetTofId();
  int IdLc  = DCGeomMan::GetInstance().GetLcId();
  double zTof = DCGeomMan::GetInstance().GetLocalZ( IdTof ); 
  double zLc = DCGeomMan::GetInstance().GetLocalZ( IdLc ); 

  double multi_SdcIn=0.;
  double multi_SdcOut=0.;
  
  //////////////SDC1&2 number of hit layer
  {
    for( int layer=1; layer<=NumOfLayersSdcIn; ++layer ){
      const DCHitContainer &contIn =DCAna->GetSdcInHC(layer);
      int nhIn=contIn.size();
      multi_SdcIn += double(nhIn);
    }
  }
  //------------------------Cut
  if( multi_SdcIn/double(NumOfLayersSdcIn) > MaxMultiHitSdcIn )
   { tree->Fill(); return true; }
  
  //////////////SDC3&4 number of hit layer
  {
    for( int layer=1; layer<=NumOfLayersSdcOut; ++layer ){
      const DCHitContainer &contOut =DCAna->GetSdcOutHC(layer);
      int nhOut=contOut.size();
      multi_SdcOut += double(nhOut);
    }
  }
  //------------------------Cut
  if( multi_SdcOut/double(NumOfLayersSdcOut) > MaxMultiHitSdcOut )
    { tree->Fill(); return true; }

  HF1( 1, 11. );

  //////////////SdcIn tracking
  DCAna->TrackSearchSdcIn();

  int ntIn=DCAna->GetNtracksSdcIn();
  if( ntIn>MaxHits ) ntIn=MaxHits;
  event.ntIn=ntIn;
  HF1( 10, double(ntIn) );
  for( int it=0; it<ntIn; ++it ){
    DCLocalTrack *tp=DCAna->GetTrackSdcIn(it);
    int nh=tp->GetNHit();
    double chisqr=tp->GetChiSquare();
    double x0=tp->GetX0(), y0=tp->GetY0();
    double u0=tp->GetU0(), v0=tp->GetV0();

    HF1( 11, double(nh) );
    HF1( 12, chisqr ); 
    HF1( 14, x0 ); HF1( 15, y0 );
    HF1( 16, u0 ); HF1( 17, v0 );
    HF2( 18, x0, u0 ); HF2( 19, y0, v0 );
    HF2( 20, x0, y0 );
    
    event.nhIn[it] = nh;
    event.chisqrIn[it] = chisqr;
    event.x0In[it] = x0;
    event.y0In[it] = y0;
    event.u0In[it] = u0;
    event.v0In[it] = v0;
   
    for( int ih=0; ih<nh; ++ih ){
      DCLTrackHit *hit=tp->GetHit(ih);
      int layerId_f=hit->GetLayer();
      int layerId=layerId_f;
      HF1( 13, layerId );
      double wire=hit->GetWire();
      double wire_position=hit->GetWirePosition();
      double dt=hit->GetDriftTime(), dl=hit->GetDriftLength();
      double tdc_value = hit->GetTdcVal();
      if( layerId<=4) 
	{
	//HF1( 100*layerId+1, wire);
	HF1( 100*layerId+1, wire_position);
	HF1( 100*layerId+2, tdc_value);
        HF1( 100*layerId+3, dl );
	}
      else 
       {
       HF1( 100*layerId+1, wire-0.5 );
       HF1( 100*layerId+2, dt );
       HF1( 100*layerId+3, dl );
       }
      double xcal=hit->GetXcal(), ycal=hit->GetYcal();
      double pos=hit->GetLocalHitPos(), res=hit->GetResidual();
      HF1( 100*layerId+4, pos );
      HF1( 100*layerId+5, res );
      HF2( 100*layerId+6, pos, res );
      HF2( 100*layerId+7, xcal, ycal);

      if( layerId==1 ) event.layer1ResL[it] = res;
      if( layerId==2 ) event.layer2ResL[it] = res;
      if( layerId==3 ) event.layer3ResL[it] = res;
      if( layerId==4 ) event.layer4ResL[it] = res;

      if( layerId==5 ) event.layer5ResL[it] = res;
      if( layerId==6 ) event.layer6ResL[it] = res;
      if( layerId==7 ) event.layer7ResL[it] = res;
      if( layerId==8 ) event.layer8ResL[it] = res;
      if( layerId==9 ) event.layer9ResL[it] = res;
      if( layerId==10) event.layer10ResL[it] = res;
    }
  } 
  //if( ntIn<1 ) return true;
  //if( !(ntIn==1) ) return true;

  HF1( 1, 12. );  
  
  //////////////SdcIn vs Tof&LC cut for Proton event
  {
    int ntOk=0;
    for( int it=0; it<ntIn; ++it ){
      DCLocalTrack *tp=DCAna->GetTrackSdcIn(it);
      if(!tp) continue;
    
      int nh=tp->GetNHit();
      double chisqr=tp->GetChiSquare();
      double u0=tp->GetU0(), v0=tp->GetV0();
      double x0=tp->GetX0(), y0=tp->GetY0();

      bool condTof=false, condLc=false;
      for( int j=0; j<ncTof; ++j ){
	HodoCluster *clTof=hodoAna->GetClusterTOF(j);
	if( !clTof || !clTof->GoodForAnalysis() ) continue;
	double ttof=clTof->CMeanTime()-time0;
	//------------------------Cut
	if( MinModTimeTof< ttof+14.0*u0 && ttof+14.0*u0 <MaxModTimeTof )
	  condTof=true;
      }
      for( int j=0; j<ncLc; ++j ){
	HodoCluster *clLc=hodoAna->GetClusterLC(j);
	if( !clLc || !clLc->GoodForAnalysis() ) continue;
	double tlc=clLc->CMeanTime()-time0;
	//------------------------Cut
	if( MinModTimeLc< tlc+14.0*u0 && tlc+14.0*u0 <MaxModTimeLc )
	  condLc=true;
      }
      if( condTof && condLc ){
	++ntOk;
	for( int j=0; j<ncTof; ++j ){
	  HodoCluster *clTof=hodoAna->GetClusterTOF(j);
	  if( !clTof || !clTof->GoodForAnalysis() ) continue;
	  double ttof=clTof->CMeanTime()-time0;
	}
	for( int j=0; j<ncLc; ++j ){
	  HodoCluster *clLc=hodoAna->GetClusterLC(j);
	  if( !clLc || !clLc->GoodForAnalysis() ) continue;
	  double tlc=clLc->CMeanTime()-time0;
	}
	//if(ntOk>0) tp->GoodForTracking( false );
      }	
      else {
	//	tp->GoodForTracking( false );
      }
    }
    //if( ntOk<1 ) return true;
  }

  HF1( 1, 13. );

  //////////////SDCOut tracking
  DCAna->TrackSearchSdcOut();

  int ntOut=DCAna->GetNtracksSdcOut();
  if( ntOut>MaxHits ) ntOut=MaxHits;
  event.ntOut=ntOut;
  HF1( 30, double(ntOut) );
  for( int it=0; it<ntOut; ++it ){
    DCLocalTrack *tp=DCAna->GetTrackSdcOut(it);
    int nh=tp->GetNHit();
    double chisqr=tp->GetChiSquare();
    double u0=tp->GetU0(), v0=tp->GetV0();
    double x0=tp->GetX(zTof), y0=tp->GetY(zTof);

    HF1( 31, double(nh) );
    HF1( 32, chisqr ); 
    HF1( 34, x0 ); HF1( 35, y0 );
    HF1( 36, u0 ); HF1( 37, v0 );
    HF2( 38, x0, u0 ); HF2( 39, y0, v0 );
    HF2( 40, x0, y0 );

    event.nhOut[it] = nh;
    event.chisqrOut[it] = chisqr;
    event.x0Out[it] = x0;
    event.y0Out[it] = y0;
    event.u0Out[it] = u0;
    event.v0Out[it] = v0;
   
    for( int ih=0; ih<nh; ++ih ){
      DCLTrackHit *hit=tp->GetHit(ih);
      int layerId=hit->GetLayer();
      HF1( 33, layerId );
      double wire=hit->GetWire();
      double dt=hit->GetDriftTime(), dl=hit->GetDriftLength();
      HF1( 100*layerId+1, wire-0.5 );
      HF1( 100*layerId+2, dt );
      HF1( 100*layerId+3, dl );
      double xcal=hit->GetXcal(), ycal=hit->GetYcal();
      double pos=hit->GetLocalHitPos(), res=hit->GetResidual();
      HF1( 100*layerId+4, pos );
      HF1( 100*layerId+5, res );
      HF2( 100*layerId+6, pos, res );
      HF2( 100*layerId+7, xcal, ycal);

      if( layerId==31 ) event.layer31ResL[it] = res;
      if( layerId==32 ) event.layer32ResL[it] = res;
      if( layerId==33 ) event.layer33ResL[it] = res;
      if( layerId==34 ) event.layer34ResL[it] = res;
      if( layerId==35 ) event.layer35ResL[it] = res;
      if( layerId==36 ) event.layer36ResL[it] = res;

      if( layerId==37 ) event.layer37ResL[it] = res;
      if( layerId==38 ) event.layer38ResL[it] = res;
      if( layerId==39 ) event.layer39ResL[it] = res;
      if( layerId==40 ) event.layer40ResL[it] = res;
      if( layerId==41 ) event.layer41ResL[it] = res;
      if( layerId==42 ) event.layer42ResL[it] = res;
    }
  }

//   for( int i1=0; i1<ntIn; ++i1 ){
//     DCLocalTrack *trIn=DCAna->GetTrackSdcIn(i1);
//     double yin=trIn->GetY(500.), vin=trIn->GetV0();
//     for( int i2=0; i2<ntOut; ++i2 ){
//       DCLocalTrack *trOut=DCAna->GetTrackSdcOut(i2);
//       double yout=trOut->GetY(3800.), vout=trOut->GetV0();
//       HF2( 20001, yin, yout ); HF2( 20002, vin, vout );
//       HF2( 20003, vin, yout ); HF2( 20004, vout, yin );
//     }
//   }
  //if( ntOut<1 ) return true;

  HF1( 1, 14. );

  //tree->Fill(); return true;

  HF1( 1, 20. );

  //if( !(ntIn==1 && ntOut ==1) ) return true;

  HF1( 1, 21. );

  //////////////SKS tracking
  DCAna->TrackSearchSks();

  int ntSks=DCAna->GetNTracksSks();
  HF1( 50, double(ntSks) );
  if( ntSks>MaxHits ) ntSks=MaxHits;
  event.ntSks = ntSks;
  for( int i=0; i<ntSks; ++i ){
    SksTrack *tp=DCAna->GetSksTrack(i);
    if(!tp) continue;
    int nh=tp->GetNHits();
    double chisqr=tp->chisqr();
    ThreeVector Ppos=tp->PrimaryPosition();
    ThreeVector Pmom=tp->PrimaryMomentum();
    double pathL=tp->PathLengthToTOF();
    double xt=Ppos.x(), yt=Ppos.y();
    double p=Pmom.Mag();
    double ut=Pmom.x()/p, vt=Pmom.y()/p;
    double cost = 1./sqrt(1.+ut*ut+vt*vt);
    double theta = acos(cost)*Rad2Deg;
    double phi   = atan2( ut, vt );

    HF1( 51, double(nh) );
    HF1( 52, chisqr );
    HF1( 54, xt ); HF1( 55, yt ); HF1( 56, ut ); HF1( 57,vt );
    HF2( 58, xt, ut ); HF2( 59, yt, vt ); HF2( 60, xt, yt );
    HF1( 61, p ); HF1( 62, pathL );

    event.nhSks[i] = nh;
    event.chisqrSks[i] = chisqr;
    event.path[i] = pathL;
    event.p[i] = p;

    event.xtgt[i] = xt;
    event.ytgt[i] = yt;
    event.utgt[i] = ut;
    event.vtgt[i] = vt; 

    event.theta[i] = theta;
    event.phi[i] = phi;

    // scattered-particle flight time target -> TOF: the RK path
    // length starts at the target, so the beam leg BH2 -> target
    // is subtracted using the SURVEYED positions in DCGEO
    // (K18Target id 130 and BH2 id 202, same beam-line frame) and
    // the BEAM velocity from PK18 (kaon). The old hardcoded
    // Parameter.hh offset (E13 era) is not used.
    const DCGeomMan & geomBt = DCGeomMan::GetInstance();
    double lbt = geomBt.GetLocalZ( 130 ) - geomBt.GetLocalZ( 202 );
    double pk18 = ConfMan::GetConfManager()->K18Momentum();
    double bBeam = pk18/sqrt(pk18*pk18 + KaonMass*KaonMass);
    double tBeam = lbt/(bBeam*299.792458);  // [mm]/[mm/ns] = ns
    double m2 = -9999.;
    for( int j=0; j<ncTof; ++j ){
      HodoCluster *clTof=hodoAna->GetClusterTOF(j);
      if( !clTof ) continue;
      double time = clTof->CMeanTime()-time0-tBeam;

      if( time>0 ){
	m2=MassSquare( p, pathL, time );
	HF1( 63, m2 );
      }
    }
    event.m2[i] = m2;
   
    for( int j=0; j<nh; ++j ){
      TrackHit *hit=tp->GetHit(j);
      if(!hit) continue;
      int layerId=hit->GetLayer();
      //int layerId=layerId_f-1;
      HF1( 53, layerId );
      double wire=hit->GetHit()->GetWire();
      double dt=hit->GetHit()->GetDriftTime();
      double dl=hit->GetHit()->GetDriftLength();
      double pos=hit->GetCalLPos(), res=hit->GetResidual();
      DCLTrackHit *lhit=hit->GetHit();
      double xcal=lhit->GetXcal(), ycal=lhit->GetYcal();
      HF1( 100*layerId+11, double(wire)-0.5 ); 
      HF1( 100*layerId+12, dt );  HF1( 100*layerId+13, dl ); 
      HF1( 100*layerId+14, pos ); HF1( 100*layerId+15, res );
      HF2( 100*layerId+16, pos, res );
      HF2( 100*layerId+17, xcal, ycal );

      if( layerId==1 ) event.layer1ResG[i] = res;
      if( layerId==2 ) event.layer2ResG[i] = res;
      if( layerId==3 ) event.layer3ResG[i] = res;
      if( layerId==4 ) event.layer4ResG[i] = res;
			     				   
      if( layerId==5 ) event.layer5ResG[i] = res;
      if( layerId==6 ) event.layer6ResG[i] = res;
      if( layerId==7 ) event.layer7ResG[i] = res;
      if( layerId==8 ) event.layer8ResG[i] = res;
      if( layerId==9 ) event.layer9ResG[i] = res;
      if( layerId==10) event.layer10ResG[i] = res;

      if( layerId==31 ) event.layer31ResG[i] = res;
      if( layerId==32 ) event.layer32ResG[i] = res;
      if( layerId==33 ) event.layer33ResG[i] = res;
      if( layerId==34 ) event.layer34ResG[i] = res;
      if( layerId==35 ) event.layer35ResG[i] = res;
      if( layerId==36 ) event.layer36ResG[i] = res;
			      			         
      if( layerId==37 ) event.layer37ResG[i] = res;
      if( layerId==38 ) event.layer38ResG[i] = res;
      if( layerId==39 ) event.layer39ResG[i] = res;
      if( layerId==40 ) event.layer40ResG[i] = res;
      if( layerId==41 ) event.layer41ResG[i] = res;
      if( layerId==42 ) event.layer42ResG[i] = res;

      //      std::cout<<"Layer="     << layerId
      // 	       <<" Residual=" << res << std::endl;
    }

    DCLocalTrack *trIn =tp->GetLocalTrackIn();
    DCLocalTrack *trOut=tp->GetLocalTrackOut();
    if(trIn){
      int nhIn=trIn->GetNHit();
      double x0in=trIn->GetX0(), y0in=trIn->GetY0();
      double u0in=trIn->GetU0(), v0in=trIn->GetV0();
      double chiin=trIn->GetChiSquare();
      HF1( 21, double(nhIn) ); HF1( 22, chiin );
      HF1( 24, x0in ); HF1( 25, y0in ); HF1( 26, u0in ); HF1( 27, v0in );
      HF2( 28, x0in, u0in ); HF2( 29, y0in, v0in );
      for( int jin=0; jin<nhIn; ++jin ){
	int layer=trIn->GetHit(jin)->GetLayer();
	HF1( 23, layer );
      }
    }
    if(trOut){
      int nhOut=trOut->GetNHit();
      double x0out=trOut->GetX(zTof), y0out=trOut->GetY(zTof);
      double u0out=trOut->GetU0(), v0out=trOut->GetV0();
      double chiout=trOut->GetChiSquare();
      HF1( 41, double(nhOut) ); HF1( 42, chiout );
      HF1( 44, x0out ); HF1( 45, y0out ); HF1( 46, u0out ); HF1( 47, v0out );
      HF2( 48, x0out, u0out ); HF2( 49, y0out, v0out );
      for( int jout=0; jout<nhOut; ++jout ){
	int layer=trOut->GetHit(jout)->GetLayer();
	HF1( 43, layer );
      }
    }
  }

  for( int i=0; i<ntSks; ++i ){
    SksTrack *tp=DCAna->GetSksTrack(i);
    if(!tp) continue;
    DCLocalTrack *trIn =tp->GetLocalTrackIn();
    DCLocalTrack *trOut=tp->GetLocalTrackOut();
    if( !trIn || !trOut ) continue;

    double yin=trIn->GetY(500.), vin=trIn->GetV0();
    double yout=trOut->GetY(3800.), vout=trOut->GetV0();
    HF2( 20021, yin, yout ); HF2( 20022, vin, vout );
    HF2( 20023, vin, yout ); HF2( 20024, vout, yin );
  }

  HF1( 1, 22. );
  
  tree->Fill();

  return true;
}

void EventSksTracking::InitializeEvent( void )
{
  event.runnum = -1;
  event.evnum = -1;
  // size branches must be 0 (not -1): TTree::Fill crashes on
  // negative array sizes when the fill section is skipped
  event.ntIn     = 0;
  event.ntOut    = 0;
  event.ntSks    = 0;
  event.nhBh2  = 0;
  event.nhBh1  = 0;
  event.nhTof  = 0;
  event.nhLc   = 0;

  for( int it=0; it<NumOfMisc; it++){  
    event.trigpat[it] = -1;
    event.trigflag[it] = -1;
  }

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
    event.deTof[it] = -9999.0;

    event.LcSeg[it] = -1;
    event.tLc[it] = -9999.0;
    event.deLc[it] = -9999.0;
  }

  for( int it=0; it<MaxHits; it++){
    event.nhIn[it]     = -1;
    event.chisqrIn[it] = -1.0;
    event.x0In[it] = -9999.0;
    event.y0In[it] = -9999.0;
    event.u0In[it] = -9999.0;
    event.v0In[it] = -9999.0;
  }

  for( int it=0; it<MaxHits; it++){
    event.nhOut[it]     = -1; 
    event.chisqrOut[it] = -1.0;
    event.x0Out[it] = -9999.0;
    event.y0Out[it] = -9999.0;
    event.u0Out[it] = -9999.0;
    event.v0Out[it] = -9999.0;
  }

  for( int it=0; it<MaxHits; it++){
    event.nhSks[it]     = -1; 
    event.chisqrSks[it] = -1.0;
    event.path[it]  = -9999.0;
    event.p[it]     = -9999.0;
    event.m2[it]    = -9999.0;

    event.xtgt[it]  = -9999.0;
    event.ytgt[it]  = -9999.0;
    event.utgt[it]  = -9999.0;
    event.vtgt[it]  = -9999.0;
    event.theta[it] = -9999.0;
    event.phi[it]   = -9999.0;
  }

  for( int it=0; it<MaxHits; it++){
    event.layer1ResL[it]  = -9999.0;
    event.layer2ResL[it]  = -9999.0;
    event.layer3ResL[it]  = -9999.0;
    event.layer4ResL[it]  = -9999.0;
	  	
    event.layer5ResL[it]  = -9999.0;
    event.layer6ResL[it]  = -9999.0;
    event.layer7ResL[it]  = -9999.0;
    event.layer8ResL[it]  = -9999.0;
    event.layer9ResL[it]  = -9999.0;
    event.layer10ResL[it]  = -9999.0;

    event.layer31ResL[it]  = -9999.0;
    event.layer32ResL[it]  = -9999.0;
    event.layer33ResL[it]  = -9999.0;
    event.layer34ResL[it]  = -9999.0;
    event.layer35ResL[it]  = -9999.0;
    event.layer36ResL[it]  = -9999.0;
	  	 			     
    event.layer37ResL[it]  = -9999.0;
    event.layer38ResL[it]  = -9999.0;
    event.layer39ResL[it]  = -9999.0;
    event.layer40ResL[it]  = -9999.0;
    event.layer41ResL[it]  = -9999.0;
    event.layer42ResL[it]  = -9999.0;
  }

  for( int it=0; it<MaxHits; it++){
    event.layer1ResG[it]  = -9999.0;
    event.layer2ResG[it]  = -9999.0;
    event.layer3ResG[it]  = -9999.0;
    event.layer4ResG[it]  = -9999.0;
	  	
    event.layer5ResG[it]  = -9999.0;
    event.layer6ResG[it]  = -9999.0;
    event.layer7ResG[it]  = -9999.0;
    event.layer8ResG[it]  = -9999.0;
    event.layer9ResG[it]  = -9999.0;
    event.layer10ResG[it]  = -9999.0;

    event.layer31ResG[it]  = -9999.0;
    event.layer32ResG[it]  = -9999.0;
    event.layer33ResG[it]  = -9999.0;
    event.layer34ResG[it]  = -9999.0;
    event.layer35ResG[it]  = -9999.0;
    event.layer36ResG[it]  = -9999.0;
	  	 			     
    event.layer37ResG[it]  = -9999.0;
    event.layer38ResG[it]  = -9999.0;
    event.layer39ResG[it]  = -9999.0;
    event.layer40ResG[it]  = -9999.0;
    event.layer41ResG[it]  = -9999.0;
    event.layer42ResG[it]  = -9999.0;
  }
}

bool EventSksTracking::ProcessingEnd()
{
  return true;
}

VEvent *ConfMan::EventAllocator()
{
  return new EventSksTracking;
}

const int NBin1DTSdcIn  =  90;
const double MinDTSdcIn = -10.;
const double MaxDTSdcIn =  80.;
const int NBin1DLSdcIn  =  100;
const double MinDLSdcIn = -1.0;
const double MaxDLSdcIn =  3.0;

const int NBin1DTSdcOut =   400;
const double MinDTSdcOut = -100.;
const double MaxDTSdcOut =  300.;
const int NBin1DLSdcOut =   100;
const double MinDLSdcOut = -5.0;
const double MaxDLSdcOut = 15.0;

bool ConfMan:: InitializeHistograms()
{  
  HB1( 1, "Status", 30, 0., 30. );

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
  HB2( 20, "X0%Y0 SdcIn", 100, -100., 100., 100, -100., 100. );

  HB1( 21, "#Hits of Track SdcIn [SksTrack]", 15, 0., 15. );
  HB1( 22, "Chisqr SdcIn [SksTrack]", 500, 0., 50. ); 
  HB1( 23, "LayerId SdcIn [SksTrack]", 15, 0., 15. );
  HB1( 24, "X0 SdcIn [SksTrack]", 400, -100., 100. ); 
  HB1( 25, "Y0 SdcIn [SksTrack]", 400, -100., 100. );
  HB1( 26, "U0 SdcIn [SksTrack]", 200, -0.20, 0.20 );
  HB1( 27, "V0 SdcIn [SksTrack]", 200, -0.20, 0.20 );
  HB2( 28, "U0%X0 SdcIn [SksTrack]", 100, -100., 100., 100, -0.20, 0.20 );
  HB2( 29, "V0%Y0 SdcIn [SksTrack]", 100, -100., 100., 100, -0.20, 0.20 );
  //HB2( 30, "X0%Y0 SdcIn [SksTrack]", 100, -100., 100., 100, -100., 100. );

  HB1( 30, "#Tracks SdcOut", 10, 0., 10. );
  HB1( 31, "#Hits of Track SdcOut", 20, 0., 20. );
  HB1( 32, "Chisqr SdcOut", 500, 0., 50. ); 
  HB1( 33, "LayerId SdcOut", 20, 30., 50. );
  HB1( 34, "X0 SdcOut", 1400, -1200., 1200. ); 
  HB1( 35, "Y0 SdcOut", 1000, -500., 500. );
  HB1( 36, "U0 SdcOut",  700, -0.35, 0.35 );
  HB1( 37, "V0 SdcOut",  200, -0.20, 0.20 );
  HB2( 38, "U0%X0 SdcOut", 120, -1200., 1200., 100, -0.40, 0.40 );
  HB2( 39, "V0%Y0 SdcOut", 100, -500., 500., 100, -0.20, 0.20 );
  HB2( 40, "X0%Y0 SdcOut", 100, -1200., 1200., 100, -500., 500. );

  HB1( 41, "#Hits of Track SdcOut [SksTrack]", 20, 0., 20. );
  HB1( 42, "Chisqr SdcOut [SksTrack]", 500, 0., 50. ); 
  HB1( 43, "LayerId SdcOut [SksTrack]", 20, 30., 50. );
  HB1( 44, "X0 SdcOut [SksTrack]", 1400, -1200., 1200. ); 
  HB1( 45, "Y0 SdcOut [SksTrack]", 1000, -500., 500. );
  HB1( 46, "U0 SdcOut [SksTrack]",  700, -0.35, 0.35 );
  HB1( 47, "V0 SdcOut [SksTrack]",  200, -0.10, 0.10 );
  HB2( 48, "U0%X0 SdcOut [SksTrack]", 120, -600., 600., 100, -0.40, 0.40 );
  HB2( 49, "V0%Y0 SdcOut [SksTrack]", 100, -500., 500., 100, -0.10, 0.10 );
  //HB2( 50, "X0%Y0 SdcOut [SksTrack]", 100, -700., 700., 100, -500., 500. );

  HB1( 50, "#Tracks SKS", 10, 0., 10. );
  HB1( 51, "#Hits of SksTrack", 30, 0., 30. );
  HB1( 52, "Chisqr SksTrack", 500, 0., 100. );
  HB1( 53, "LayerId SksTrack", 50, 0., 50. );
  HB1( 54, "Xtgt SksTrack", 200, -100., 100. );
  HB1( 55, "Ytgt SksTrack", 200, -100., 100. );
  HB1( 56, "Utgt SksTrack", 300, -0.30, 0.30 );
  HB1( 57, "Vtgt SksTrack", 300, -0.20, 0.20 );
  HB2( 58, "U%Xtgt SksTrack", 100, -100., 100., 100, -0.25, 0.25 );
  HB2( 59, "V%Ytgt SksTrack", 100, -100., 100., 100, -0.10, 0.10 );
  HB2( 60, "Y%Xtgt SksTrack", 100, -100., 100., 100, -100., 100. );
  HB1( 61, "P SksTrack", 1400, 0.80, 2.2 );
  HB1( 62, "PathLength SksTrack", 600, 3000., 6000. );
  HB1( 63, "MassSqr", 600, -0.2, 1.2 );

  // SDC1
  for( int i=1; i<=4; ++i ){
    std::ostringstream title1, title2, title3;
    title1 << "HitPat Sdc" << std::setw(2) << i;
    HB1( 100*i+1, title1.str().c_str(), 300., -150., 150. );
    title2 << "DriftTime Sdc" << std::setw(2) << i;
    HB1( 100*i+2, title2.str().c_str(), 2000, -20, 30 );
    title3 << "DriftLength Sdc" << std::setw(2) << i;
    HB1( 100*i+3, title3.str().c_str(), NBin1DLSdcIn, MinDLSdcIn, MaxDLSdcIn );

    std::ostringstream title4, title5, title6, title7;
    title4 << "Position Sdc" << std::setw(2) << i;
    HB1( 100*i+4, title4.str().c_str(), 500, -100., 100. );
    title5 << "Residual Sdc" << std::setw(2) << i;
    HB1( 100*i+5, title5.str().c_str(), 200, -2.0, 2.0 );
    title6 << "Resid%Pos Sdc" << std::setw(2) << i;
    HB2( 100*i+6, title6.str().c_str(), 50, -100., 100., 50, -1.0, 1.0 );
    title7 << "Y%Xcal Sdc" << std::setw(2) << i;
    HB2( 100*i+7, title7.str().c_str(), 100, -100., 100., 50, -50., 50. );

    title1 << " [SksTrack]";
    HB1( 100*i+11, title1.str().c_str(), 300., -150., 150. );
    title2 << " [SksTrack]";
    HB1( 100*i+12, title2.str().c_str(), NBin1DTSdcIn, MinDTSdcIn, MaxDTSdcIn );
    title3 << " [SksTrack]";
    HB1( 100*i+13, title3.str().c_str(), NBin1DLSdcIn, MinDLSdcIn, MaxDLSdcIn );
    title4 << " [SksTrack]";
    HB1( 100*i+14, title4.str().c_str(), 500, -100., 100. );
    title5 << " [SksTrack]";
    HB1( 100*i+15, title5.str().c_str(), 200, -2.0, 2.0 );
    title6 << " [SksTrack]";
    HB2( 100*i+16, title6.str().c_str(), 50, -100., 100., 50, -2.0, 2.0 );
    title7 << " [SksTrack]";
    HB2( 100*i+17, title7.str().c_str(), 100, -100., 100., 50, -50., 50. );
  }

  // SDC2
  for( int i=5; i<=10; ++i ){
    std::ostringstream title1, title2, title3;
    title1 << "HitPat Sdc" << std::setw(2) << i;
    HB1( 100*i+1, title1.str().c_str(), 96, 0., 96. );
    title2 << "DriftTime Sdc" << std::setw(2) << i;
    HB1( 100*i+2, title2.str().c_str(), NBin1DTSdcIn, MinDTSdcIn, MaxDTSdcIn );
    title3 << "DriftLength Sdc" << std::setw(2) << i;
    HB1( 100*i+3, title3.str().c_str(), NBin1DLSdcIn, MinDLSdcIn, MaxDLSdcIn );

    std::ostringstream title4, title5, title6, title7;
    title4 << "Position Sdc" << std::setw(2) << i;
    HB1( 100*i+4, title4.str().c_str(), 1000, -200., 200. );
    title5 << "Residual Sdc" << std::setw(2) << i;
    HB1( 100*i+5, title5.str().c_str(), 200, -2.0, 2.0 );
    title6 << "Resid%Pos Sdc" << std::setw(2) << i;
    HB2( 100*i+6, title6.str().c_str(), 100, -200., 200., 50, -1.0, 1.0 );
    title7 << "Y%Xcal Sdc" << std::setw(2) << i;
    HB2( 100*i+7, title7.str().c_str(), 100, -200., 200., 50, -100., 100. );

    title1 << " [SksTrack]";
    HB1( 100*i+11, title1.str().c_str(), 96, 0., 96. );
    title2 << " [SksTrack]";
    HB1( 100*i+12, title2.str().c_str(), NBin1DTSdcIn, MinDTSdcIn, MaxDTSdcIn );
    title3 << " [SksTrack]";
    HB1( 100*i+13, title3.str().c_str(), NBin1DLSdcIn, MinDLSdcIn, MaxDLSdcIn );
    title4 << " [SksTrack]";
    HB1( 100*i+14, title4.str().c_str(), 1000, -250., 250. );
    title5 << " [SksTrack]";
    HB1( 100*i+15, title5.str().c_str(), 200, -2.0, 2.0 );
    title6 << " [SksTrack]";
    HB2( 100*i+16, title6.str().c_str(), 100, -200., 200., 50, -2.0, 2.0 );
    title7 << " [SksTrack]";
    HB2( 100*i+17, title7.str().c_str(), 100, -250., 250., 50, -100., 100. );
  }

  // SDC34
  for( int i=31; i<=42; ++i ){
    std::ostringstream title1, title2, title3;
    title1 << "HitPat Sdc" << std::setw(2) << i;
    HB1( 100*i+1, title1.str().c_str(), 120, 0., 120. );
    title2 << "DriftTime Sdc" << std::setw(2) << i;
    HB1( 100*i+2, title2.str().c_str(), NBin1DTSdcOut, MinDTSdcOut, MaxDTSdcOut );
    title3 << "DriftLength Sdc" << std::setw(2) << i;
    HB1( 100*i+3, title3.str().c_str(), NBin1DLSdcOut, MinDLSdcOut, MaxDLSdcOut );

    std::ostringstream title4, title5, title6, title7;
    title4 << "Position Sdc" << std::setw(2) << i;
    HB1( 100*i+4, title4.str().c_str(), 1000, -600., 600. );
    title5 << "Residual Sdc" << std::setw(2) << i;
    HB1( 100*i+5, title5.str().c_str(), 200, -2.0, 2.0 );
    title6 << "Resid%Pos Sdc" << std::setw(2) << i;
    HB2( 100*i+6, title6.str().c_str(), 100, -600., 600., 50, -1.0, 1.0 );
    title7 << "Y%Xcal Sdc" << std::setw(2) << i;
    HB2( 100*i+7, title7.str().c_str(), 100, -600., 600., 100, -600., 600. );

    title1 << " [SksTrack]";
    HB1( 100*i+11, title1.str().c_str(), 120, 0., 120. );
    title2 << " [SksTrack]";
    HB1( 100*i+12, title2.str().c_str(), NBin1DTSdcOut, MinDTSdcOut, MaxDTSdcOut );
    title3 << " [SksTrack]";
    HB1( 100*i+13, title3.str().c_str(), NBin1DLSdcOut, MinDLSdcOut, MaxDLSdcOut );
    title4 << " [SksTrack]";
    HB1( 100*i+14, title4.str().c_str(), 1000, -600., 600. );
    title5 << " [SksTrack]";
    HB1( 100*i+15, title5.str().c_str(), 200, -2.0, 2.0 );
    title6 << " [SksTrack]";
    HB2( 100*i+16, title6.str().c_str(), 100, -600., 600., 50, -2.0, 2.0 );
    title7 << " [SksTrack]";
    HB2( 100*i+17, title7.str().c_str(), 100, -600., 600., 100, -600., 600. );
  }

  HB2( 20001, "Yout%Yin", 100, -150., 150., 120, -300., 300. );
  HB2( 20002, "Vout%Vin", 100, -0.05, 0.05, 100, -0.1, 0.1 );
  HB2( 20003, "Yout%Vin", 100, -0.05, 0.05, 100, -300., 300. );
  HB2( 20004, "Yin%Vout", 100, -0.10, 0.10, 100, -150., 150. );

  HB2( 20021, "Yout%Yin [SksTrack]", 100, -150., 150., 120, -300., 300. );
  HB2( 20022, "Vout%Vin [SksTrack]", 100, -0.05, 0.05, 100, -0.1, 0.1 );
  HB2( 20023, "Yout%Vin [SksTrack]", 100, -0.05, 0.05, 100, -300., 300. );
  HB2( 20024, "Yin%Vout [SksTrack]", 100, -0.10, 0.10, 100, -150., 150. );
  
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
  tree->Branch("deTof",    event.deTof,   "deTof[nhTof]/D");

  tree->Branch("nhLc",   &event.nhLc,   "nhLc/I");
  tree->Branch("LcSeg",   event.LcSeg,  "LcSeg[nhLc]/D");
  tree->Branch("tLc",     event.tLc,    "tLc[nhLc]/D");
  tree->Branch("deLc",    event.deLc,   "deLc[nhLc]/D");

  //Tracking  
  tree->Branch("ntIn",    &event.ntIn,     "ntIn/I");
  tree->Branch("nhIn",     event.nhIn,     "nhIn[ntIn]/I");
  tree->Branch("chisqrIn", event.chisqrIn, "chisqrIn[ntIn]/D");
  tree->Branch("x0In",     event.x0In,     "x0In[ntIn]/D");
  tree->Branch("y0In",     event.y0In,     "y0In[ntIn]/D");
  tree->Branch("u0In",     event.u0In,     "u0In[ntIn]/D");
  tree->Branch("v0In",     event.v0In,     "v0In[ntIn]/D");

  tree->Branch("layer1ResL",  event.layer1ResL,  "layer1ResL[ntIn]/D");
  tree->Branch("layer2ResL",  event.layer2ResL,  "layer2ResL[ntIn]/D");
  tree->Branch("layer3ResL",  event.layer3ResL,  "layer3ResL[ntIn]/D");
  tree->Branch("layer4ResL",  event.layer4ResL,  "layer4ResL[ntIn]/D");
		      		  	        	
  tree->Branch("layer5ResL",  event.layer5ResL,  "layer5ResL[ntIn]/D");
  tree->Branch("layer6ResL",  event.layer6ResL,  "layer6ResL[ntIn]/D");
  tree->Branch("layer7ResL",  event.layer7ResL,  "layer7ResL[ntIn]/D");
  tree->Branch("layer8ResL",  event.layer8ResL,  "layer8ResL[ntIn]/D");
  tree->Branch("layer9ResL",  event.layer9ResL,  "layer9ResL[ntIn]/D");
  tree->Branch("layer10ResL", event.layer10ResL, "layer10ResL[ntIn]/D");

  tree->Branch("ntOut",   &event.ntOut,     "ntOut/I");
  tree->Branch("nhOut",    event.nhOut,     "nhOut[ntOut]/I");
  tree->Branch("chisqrOut",event.chisqrOut, "chisqrIn[ntOut]/D");
  tree->Branch("x0Out",    event.x0Out,     "x0Out[ntOut]/D");
  tree->Branch("y0Out",    event.y0Out,     "y0Out[ntOut]/D");
  tree->Branch("u0Out",    event.u0Out,     "u0Out[ntOut]/D");
  tree->Branch("v0Out",    event.v0Out,     "v0Out[ntOut]/D");

  tree->Branch("layer31ResL",  event.layer31ResL,  "layer31ResL[ntOut]/D");
  tree->Branch("layer32ResL",  event.layer32ResL,  "layer32ResL[ntOut]/D");
  tree->Branch("layer33ResL",  event.layer33ResL,  "layer33ResL[ntOut]/D");
  tree->Branch("layer34ResL",  event.layer34ResL,  "layer34ResL[ntOut]/D");
  tree->Branch("layer35ResL",  event.layer35ResL,  "layer35ResL[ntOut]/D");
  tree->Branch("layer36ResL",  event.layer36ResL,  "layer36ResL[ntOut]/D");
		       		   	    		 
  tree->Branch("layer37ResL",  event.layer37ResL,  "layer37ResL[ntOut]/D");
  tree->Branch("layer38ResL",  event.layer38ResL,  "layer38ResL[ntOut]/D");
  tree->Branch("layer39ResL",  event.layer39ResL,  "layer39ResL[ntOut]/D");
  tree->Branch("layer40ResL",  event.layer40ResL,  "layer40ResL[ntOut]/D");
  tree->Branch("layer41ResL",  event.layer41ResL,  "layer41ResL[ntOut]/D");
  tree->Branch("layer42ResL",  event.layer42ResL,  "layer42ResL[ntOut]/D");

  tree->Branch("ntSks",      &event.ntSks,     "ntSks/I");
  tree->Branch("nhSks",       event.nhSks,     "nhSks[ntSks]/I");
  tree->Branch("chisqrSks",   event.chisqrSks, "chisqrSks[ntSks]/D");
  tree->Branch("path",        event.path,      "path[ntSks]/D");
  tree->Branch("p",           event.p,         "p[ntSks]/D");
  tree->Branch("m2",          event.m2,        "m2[ntSks]/D");

  tree->Branch("xtgt",    event.xtgt,   "xtgt[ntSks]/D");
  tree->Branch("ytgt",    event.ytgt,   "ytgt[ntSks]/D");
  tree->Branch("utgt",    event.utgt,   "utgt[ntSks]/D");
  tree->Branch("vtgt",    event.vtgt,   "vtgt[ntSks]/D");

  tree->Branch("theta",   event.theta,  "theta[ntSks]/D");
  tree->Branch("phi",     event.phi,    "phi[ntSks]/D");

  tree->Branch("layer1ResG",  event.layer1ResG,  "layer1ResG[ntSks]/D");
  tree->Branch("layer2ResG",  event.layer2ResG,  "layer2ResG[ntSks]/D");
  tree->Branch("layer3ResG",  event.layer3ResG,  "layer3ResG[ntSks]/D");
  tree->Branch("layer4ResG",  event.layer4ResG,  "layer4ResG[ntSks]/D");
				    		  	
  tree->Branch("layer5ResG",  event.layer5ResG,  "layer5ResG[ntSks]/D");
  tree->Branch("layer6ResG",  event.layer6ResG,  "layer6ResG[ntSks]/D");
  tree->Branch("layer7ResG",  event.layer7ResG,  "layer7ResG[ntSks]/D");
  tree->Branch("layer8ResG",  event.layer8ResG,  "layer8ResG[ntSks]/D");
  tree->Branch("layer9ResG",  event.layer9ResG,  "layer9ResG[ntSks]/D");
  tree->Branch("layer10ResG", event.layer10ResG, "layer10ResG[ntSks]/D");

  tree->Branch("layer31ResG",  event.layer31ResG,  "layer31ResG[ntSks]/D");
  tree->Branch("layer32ResG",  event.layer32ResG,  "layer32ResG[ntSks]/D");
  tree->Branch("layer33ResG",  event.layer33ResG,  "layer33ResG[ntSks]/D");
  tree->Branch("layer34ResG",  event.layer34ResG,  "layer34ResG[ntSks]/D");
  tree->Branch("layer35ResG",  event.layer35ResG,  "layer35ResG[ntSks]/D");
  tree->Branch("layer36ResG",  event.layer36ResG,  "layer36ResG[ntSks]/D");
		       		            	           
  tree->Branch("layer37ResG",  event.layer37ResG,  "layer37ResG[ntSks]/D");
  tree->Branch("layer38ResG",  event.layer38ResG,  "layer38ResG[ntSks]/D");
  tree->Branch("layer39ResG",  event.layer39ResG,  "layer39ResG[ntSks]/D");
  tree->Branch("layer40ResG",  event.layer40ResG,  "layer40ResG[ntSks]/D");
  tree->Branch("layer41ResG",  event.layer41ResG,  "layer41ResG[ntSks]/D");
  tree->Branch("layer42ResG",  event.layer42ResG,  "layer42ResG[ntSks]/D");

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

  if( FieldMapFileName_!="" )
    FieldMan::GetInstance().Initialize( FieldMapFileName_ );

  return true;
}
