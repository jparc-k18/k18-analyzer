/*
  UserKPiAna.cc

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
#include "K18Track_BFT.hh"
#include "FiberCluster.hh"

#include "FileManager.hh"
#include "ParameterManager.hh"
#include "Parameter.hh"

#define BFT 1

const int BEAM_TRIG  = 1;
const int PIK_TRIG   = 2;
const int PIPI_TRIG  = 11;
const int PIP_TRIG   = 12;
const int WOMT_TRIG  = 13;
const int BH2_TRIG   = 4;

const double Deg2Rad = acos(-1.)/180.;
const double Rad2Deg = 180./acos(-1.);

VEvent::VEvent()
 {
 }

 VEvent::~VEvent()
 {
 }

class EventKPiAna 
  : public VEvent
{

private:
  RawData *rawData;
  DCAnalyzer *DCAna;
  HodoAnalyzer *hodoAna;

public:
  EventKPiAna();
  ~EventKPiAna();

  bool ProcessingBegin();
  bool ProcessingEnd();
  bool ProcessingNormal();
  bool InitializeHistograms();
  void InitializeEvent();
};

EventKPiAna::EventKPiAna()
  : VEvent(),
    rawData(0), 
    DCAna(new DCAnalyzer()),
    hodoAna(new HodoAnalyzer())
{
}

EventKPiAna::~EventKPiAna()
{
  if (DCAna)   delete DCAna;
  if (hodoAna) delete hodoAna;
  if (rawData) delete rawData;
}

#ifndef MaxHits 
#define MaxHits 50
#endif

//For Tree
struct Event{
  int runnum;
  int evnum;
  //Trigger
  int trigpat[NumOfMisc];
  int trigflag[NumOfMisc];

  int matrixpat[NumOfMatrix];
  int matrixflag[NumOfMatrix];

  //Hodoscope
  int nhitsBh2;
  int nhBH2;
  double BH2Seg[MaxHits];
  double tBH2[MaxHits];
  double deBH2[MaxHits];

  int nhBH1;
  double BH1Seg[MaxHits];
  double tBH1[MaxHits];
  double deBH1[MaxHits];
  double btof[MaxHits];

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

  int Ac1Tdc[NumOfSegAC1];
  int multip_bh1, multip_bh2;


  int nhBAC;
  int BACSeg[MaxHits];
  double tBAC[NumOfSegBAC];
  double aBAC[NumOfSegBAC];

  //DC Beam
  int ntBcIn;
  int nhBcIn[MaxHits]; 
  double chisqrBcIn[MaxHits];
  double x0BcIn[MaxHits];
  double y0BcIn[MaxHits];
  double u0BcIn[MaxHits];
  double v0BcIn[MaxHits];
  double xbh1BcIn[MaxHits];
  double ybh1BcIn[MaxHits];

  int ntBcOut;
  int nhBcOut[MaxHits]; 
  double chisqrBcOut[MaxHits];
  double x0BcOut[MaxHits];
  double y0BcOut[MaxHits];
  double u0BcOut[MaxHits];
  double v0BcOut[MaxHits];
  double xtgtBcOut[MaxHits];
  double ytgtBcOut[MaxHits];
  double xbh2BcOut[MaxHits];
  double ybh2BcOut[MaxHits];

  int    ntK18;
  int    nhK18[MaxHits]; 
  double chisqrK18[MaxHits]; 
  double pK18[MaxHits]; 
  double xtb[MaxHits];
  double ytb[MaxHits];
  double utb[MaxHits];
  double vtb[MaxHits];
  double thetab[MaxHits];

  //DC SKS
  int ntSdcIn;
  int nhSdcIn[MaxHits]; 
  double chisqrSdcIn[MaxHits];
  double x0SdcIn[MaxHits];
  double y0SdcIn[MaxHits];
  double u0SdcIn[MaxHits];
  double v0SdcIn[MaxHits];

  int ntSdcOut;
  int nhSdcOut[MaxHits]; 
  double chisqrSdcOut[MaxHits];
  double u0SdcOut[MaxHits];
  double v0SdcOut[MaxHits];
  double xtofSdcOut[MaxHits];
  double ytofSdcOut[MaxHits];

  int    ntSks;
  int    nhSks[MaxHits]; 
  double chisqrSks[MaxHits]; 
  double path[MaxHits]; 
  double pSks[MaxHits]; 
  double m2[MaxHits]; 
  double xts[MaxHits];
  double yts[MaxHits];
  double uts[MaxHits];
  double vts[MaxHits];
  double thetas[MaxHits];
  
  int    ncBFT;
  double cxpos[MaxHits];
  double ctime_BFT[MaxHits];
  //Reaction
  int    nPi;
  int    nK;
  int    nKPi;
  double vtx[MaxHits];
  double vty[MaxHits];
  double vtz[MaxHits];
  double closeDist[MaxHits];
  double theta[MaxHits];
  double MissMass1[MaxHits];
  double MissMass2[MaxHits];  
  double MissMass3[MaxHits];
  double MissMass4[MaxHits];
  double MissMassMuon[MaxHits];
  double MissMassPion[MaxHits];
  double MissMassCorr[MaxHits];
  double MissMassCorrDE[MaxHits];
  double BE1[MaxHits];
  double BE2[MaxHits];
  double BE3[MaxHits];
  double BE4[MaxHits];
  double theta_CM[MaxHits];
  double cost_CM[MaxHits];

  double xk[MaxHits];
  double yk[MaxHits];
  double uk[MaxHits];
  double vk[MaxHits];
  double xpi[MaxHits];
  double ypi[MaxHits];
  double upi[MaxHits];
  double vpi[MaxHits];

  double pOrg[MaxHits];
  double pCalc[MaxHits];
  double pCorr[MaxHits];
  double pCorrDE[MaxHits];

};
static Event event;


bool EventKPiAna::ProcessingBegin()
{
  using namespace h_Utility;
  return true;
}

bool EventKPiAna::ProcessingNormal()
{

  //using namespace h_Utility;
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
  bool fl_trigger = false;
  {
    const HodoRHitContainer &cont=rawData->GetMiscRawHC();
    int nh=cont.size();
    for( int i=0; i<nh; ++i ){
      HodoRawHit *hit=cont[i];
      int seg=hit->SegmentId()+1;
      int T=hit->GetTdc1();
      event.trigflag[seg-1] = T;
      if( T>0 ){
	event.trigpat[i] = seg;
	//if(TargetTrig == seg){ fl_trigger = true;}
      }
    }
  }

  //Matrix
  int matrix_flag=0;
  {
    const HodoRHitContainer &cont=rawData->GetMatrixRawHC();
    int nh=cont.size();
    for( int i=0; i<nh; ++i ){
      HodoRawHit *hit=cont[i];
      int seg=hit->SegmentId()+1;
      int T=hit->GetTdc1();

      if( T>0 ) event.matrixpat[i] = seg;
      event.matrixflag[seg-1] = T;
      if( seg==1 && T>0 ){ matrix_flag = seg; }
    }
  }

  HF1( 1, 0. );

  //------------------------Cut
  //if(!fl_trigger){return true;}
  HF1( 1, 1. );

  //if ( matrix_flag==0 ) return true;
  //   std::cout << "KPi==" << trig_type <<std::endl;
  //   std::cout << "Flag==" << matrix_flag <<std::endl;

  ////////////////////////////Hodoscope
  //////////////BH2 Analysis
  hodoAna->DecodeBH2Hits(rawData);
  int ncBh2=hodoAna->GetNClustersBH2();
  event.nhBH2=ncBh2;
  HF1( 101, double(ncBh2) );
  //------------------------Cut
  if( ncBh2==0 ) { return true; }

  BH2Cluster *clBH2Time0=hodoAna->GetClusterBH2(0);
  double time0=clBH2Time0->CTime0();
  std::vector<BH2Cluster*> valid_bh2;
  {
    int ncOk = 0;
    double mint=clBH2Time0->CMeanTime();
    HF1( 102, clBH2Time0->ClusterSize() );
    HF1( 103, clBH2Time0->MeanSeg()+1-0.5 );
    HF1( 104, mint );
    HF1( 105, clBH2Time0->DeltaE() );
    for( int i=0; i<ncBh2; ++i ){
      BH2Cluster *cl=hodoAna->GetClusterBH2(i);
      double t=cl->CMeanTime();
      double dEbh2 = cl->DeltaE();
      //------------------------Cut
      //if( dEbh2<MinDeltaEBH2 || MaxDeltaEBH2<dEbh2 ) continue;
      ++ncOk;
      HF1( 112, cl->ClusterSize() );
      HF1( 113, cl->MeanSeg()+1-0.5 );
      HF1( 114, t );
      HF1( 115, cl->DeltaE() );

      if(i<MaxHits){
	event.BH2Seg[i]=cl->MeanSeg()+1;
	event.tBH2[i]=t;
	event.deBH2[i]=cl->DeltaE();
      }
      if( fabs(t)<fabs(mint) ){
	clBH2Time0=cl;
	mint=t; time0=clBH2Time0->CTime0();
      }
      // Add my cut ---------------------------------------------------------
      //if(-1.<t && t < 1.){ valid_bh2.push_back(cl); }
      // Add my cut ---------------------------------------------------------

    }

    // Add my cut ---------------------------------------------------------
    //if(false
    //   || valid_bh2.size() > 1
    //   || clBH2Time0->ClusterSize() > 1
    //   ){
    //  event.multip_bh2 = 1;
    //}
    
    double trig = clBH2Time0->CMeanTime();
    static const double min_t0 = -0.5;
    static const double max_t0 =  0.5;
    //if(min_t0 > trig || trig > max_t0 || ncOk<1){ return true;}
    // Add my cut ---------------------------------------------------------
  }

  HF1( 122, clBH2Time0->ClusterSize() );
  HF1( 123, clBH2Time0->MeanSeg()+1-0.5 );
  HF1( 124, clBH2Time0->CMeanTime() );
  HF1( 125, clBH2Time0->DeltaE() );

  HF1( 1, 2. );

  //////////////BH1 Analysis
  hodoAna->DecodeBH1Hits(rawData); 
  int ncBh1=hodoAna->GetNClustersBH1();
  event.nhBH1=ncBh1;
  if( ncBh1==0 ) { return true; }

    HodoCluster *clBH1Time0=hodoAna->GetClusterBH1(0);
  std::vector<HodoCluster*> valid_bh1;
  {
    int ncOk=0;
    double min_tof = clBH1Time0->CMeanTime() - time0;
    HF1( 201, double(ncBh1) );
    for( int i=0; i<ncBh1; ++i ){
      HodoCluster *cl=hodoAna->GetClusterBH1(i);
      double btof= (cl->CMeanTime())-time0;
      double dEbh1 = cl->DeltaE();
      //------------------------Cut
      //if( !(MinDeltaEBH1<dEbh1 && dEbh1<MaxDeltaEBH1) ) continue;

      HF1( 202, cl->ClusterSize() );
      HF1( 203, cl->MeanSeg()+1-0.5 );
      HF1( 204, cl->CMeanTime() );
      HF1( 205, cl->DeltaE() );
      HF1( 206, btof );

      if(i<MaxHits){
	event.BH1Seg[i]=cl->MeanSeg()+1;
	event.tBH1[i]=cl->CMeanTime();
	event.deBH1[i]=cl->DeltaE();
	event.btof[i] = btof;
      }
      //------------------------Cut
      if( btof>MinBeamToF && btof<MaxBeamToF ){
	valid_bh1.push_back(cl);
	++ncOk;
	HF1( 212, cl->ClusterSize() );
	HF1( 213, cl->MeanSeg()+1-0.5 );
	HF1( 214, cl->CMeanTime() );
	HF1( 215, cl->DeltaE() );
	if( fabs(btof)<fabs(min_tof) ){
	  clBH1Time0 = cl;
	  min_tof = btof;
	}
      }
      else{
	//cl->GoodForAnalysis( false );
      }
    }
    HF1( 211, double(ncOk) );
    //------------------------Cut
    if( ncOk<1 ) { return true; }

  }

  //if(1 == event.multip_bh1 || 1 == event.multip_bh2){return true;}
  HF1( 1, 3. );

  //std::cout<<"here!!!"<<std::endl;


  //////////////BAC
  hodoAna->DecodeBACHits(rawData);
  int nhBAC = hodoAna->GetNHitsBAC();
  event.nhBAC = nhBAC;
  {
    for(int i=0; i<nhBAC; i++)
      {
	Hodo1Hit *hit = hodoAna->GetHitBAC(i);
	int seg = hit->SegmentId();
	double A= hit->GetA();
	double T= hit->GetT();
	event.BACSeg[i] = seg;
	event.aBAC[i]=A;
	event.tBAC[i]=T;
      }
  }


  //////////////Tof Analysis
  hodoAna->DecodeTOFHits(rawData);
  int ncTof=hodoAna->GetNClustersTOF();
  event.nhTof=ncTof;
  HF1( 301, double(ncTof) );
  {
    int ncOk=0;
    for( int i=0; i<ncTof; ++i ){
      HodoCluster *cl=hodoAna->GetClusterTOF(i);
      double de=cl->DeltaE();
      double t=cl->CMeanTime()-time0, dt=cl->TimeDif();
      HF1( 302, cl->ClusterSize() );
      HF1( 303, cl->MeanSeg()+1-0.5 );
      HF1( 304, t );
      HF1( 305, de );

      if(i<MaxHits){
	event.TofSeg[i]=cl->MeanSeg()+1;
	event.tTof[i]=t;
	event.dtTof[i]=dt;
	event.deTof[i]=de;
      }
      //------------------------Cut
      if( true 
	  //	  && MinDeltaETof<de && de<MaxDeltaETof
	  && MinTimeTof  <t  && t<MaxTimeTof
	  //	  && (cl->MeanSeg()+1)>8
	  ){
	++ncOk;
	HF1( 312, cl->ClusterSize() );
	HF1( 313, cl->MeanSeg()+1-0.5 );
	HF1( 314, t );
	HF1( 315, de );
      }	
      else{
	//cl->GoodForAnalysis( false );
      }
    }
    HF1( 311, double(ncOk) );
    //------------------------Cut
    if( ncOk<1 ) {
      //return true;
    }
  }
  HF1( 1, 4. );

  //////////////LC Analysis
  hodoAna->DecodeLCHits(rawData);
  int ncLc=hodoAna->GetNClustersLC();
  event.nhLc=ncLc;
  HF1( 401, double(ncLc) );
  {
    int ncOk=0;
    for( int i=0; i<ncLc; ++i ){
      HodoCluster *cl=hodoAna->GetClusterLC(i);
      double de=cl->DeltaE();
      double t=cl->CMeanTime()-time0, dt=cl->TimeDif();
      HF1( 402, cl->ClusterSize() );
      HF1( 403, cl->MeanSeg()+1-0.5 );
      HF1( 404, t );
      HF1( 405, de );

      if(i<MaxHits){
	event.LcSeg[i]=cl->MeanSeg()+1;
	event.tLc[i]=t;	
	event.dtLc[i]=dt;
	event.deLc[i]=de;
      }
      //------------------------Cut
      if( true 
	  && MinDeltaELc<de && de<MaxDeltaELc
	  && MinTimeLc  <t  && t<MaxTimeLc 
	  //	  && (cl->MeanSeg()+1)>7
	  ){
	++ncOk;
	HF1( 412, cl->ClusterSize() );
	HF1( 413, cl->MeanSeg()+1-0.5 );
	HF1( 414, t );
	HF1( 415, de );
      }	
      else{
	//cl->GoodForAnalysis( false );
      }
    }
    HF1( 411, double(ncOk) );
    //------------------------Cut
    if( ncOk<1 ) {
      //return true;
    }
  }
  HF1( 1, 5. );


 /////////////////////////DC
  double multi_BcIn=0., multi_BcOut=0.;
  double multi_SdcIn=0., multi_SdcOut=0.;

  DCAna->DecodeRawHits( rawData ); 
  int IdTof = DCGeomMan::GetInstance().GetTofId(); 
  double zTof = DCGeomMan::GetInstance().GetLocalZ( IdTof );
  double zK18Tgt = DCGeomMan::GetInstance().GetLocalZ( IdK18Target );

  //////////////SDC3&4 number of hit in one layer not 0
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
  HF1(1, 8.);

  //////////////BC3&4 number of hit in one layer
  {
    for( int layer=1; layer<=NumOfLayersBcOut; ++layer ){
      const DCHitContainer &contOut =DCAna->GetBcOutHC(layer);
      int nhOut=contOut.size();
      multi_BcOut += double(nhOut);
    }
  }
  //------------------------Cut
  if( multi_BcOut/double(NumOfLayersBcOut) > MaxMultiHitBcOut )
   { tree->Fill(); return true; }
  HF1(1, 9.);
  //////////////SDC1&2 number of hit in one layer not 0
  {
    for( int layer=1; layer<=NumOfLayersSdcIn; ++layer ){
      const DCHitContainer &contIn =DCAna->GetSdcInHC(layer);
      int nhIn=contIn.size();
      multi_SdcIn  += double(nhIn);
    }
  }
  //------------------------Cut
  if( multi_SdcIn/double(NumOfLayersSdcIn) > MaxMultiHitSdcIn )
    { tree->Fill(); return true; }
  
  HF1( 1, 11. );

  //////////////SdcIn Analysis
  DCAna->TrackSearchSdcIn();
  int ntSdcIn=DCAna->GetNtracksSdcIn();
  event.ntSdcIn=ntSdcIn;
  HF1( 1001, double(ntSdcIn) );
  //------------------------Cut
  if( ntSdcIn<1 || ntSdcIn>MaxNTrackSdcIn ) {
    { tree->Fill(); return true; }
  }

  HF1( 1, 12. );
  {
    int ntOk=0;
    for( int i=0; i<ntSdcIn; ++i ){
      DCLocalTrack *tp=DCAna->GetTrackSdcIn(i);
      if(!tp) continue;
      int nh=tp->GetNHit();
      double chisqr=tp->GetChiSquare();
      double u0=tp->GetU0(), v0=tp->GetV0();
      double x0=tp->GetX(0.), y0=tp->GetY(0.);
      
      if(i<MaxHits){
	event.nhSdcIn[i]=nh;
	event.chisqrSdcIn[i]=chisqr;
	event.u0SdcIn[i]=u0;
	event.v0SdcIn[i]=v0;
	event.x0SdcIn[i]=x0;
	event.y0SdcIn[i]=y0;
      }
      
      HF1( 1002, double(nh) );
      HF1( 1003, chisqr );
      HF1( 1004, x0 ); HF1( 1005, y0 );
      HF1( 1006, u0 ); HF1( 1007, v0 );
      HF2( 1008, x0, u0 ); HF2( 1009, y0, v0 );
      HF2( 1010, x0, y0 );
    
    }
    if( ntOk<1 ) {
      //return true;
    }
    
  }      

  HF1( 1, 13. );

  //tree->Fill(); return true;
    //////////////SdcOut Analysis
  DCAna->TrackSearchSdcOut();
  int ntSdcOut=DCAna->GetNtracksSdcOut();
  event.ntSdcOut=ntSdcOut;
  HF1( 1201, double(ntSdcOut) );
  //------------------------Cut
  if( ntSdcOut<1 || ntSdcOut>MaxNTrackSdcOut ) {
    { tree->Fill(); return true; }
  }
  HF1( 1, 14. );

  {
    int ntOk=0;
    for( int i=0; i<ntSdcOut; ++i ){
      DCLocalTrack *tp=DCAna->GetTrackSdcOut(i);
      if(!tp) continue;
      int nh=tp->GetNHit();
      double chisqr=tp->GetChiSquare();
      double u0=tp->GetU0(), v0=tp->GetV0();
      double x0=tp->GetX( zTof ), y0=tp->GetY( zTof );
    
      if( i<MaxHits ){
	event.nhSdcOut[i]=nh;
	event.chisqrSdcOut[i]=chisqr;
	event.u0SdcOut[i]=u0;
	event.v0SdcOut[i]=v0;
	event.xtofSdcOut[i]=x0;
	event.ytofSdcOut[i]=y0;
      }
      
      HF1( 1202, double(nh) );
      HF1( 1203, chisqr );
      HF1( 1204, x0 ); HF1( 1205, y0 );
      HF1( 1206, u0 ); HF1( 1207, v0 );
      HF2( 1208, x0, u0 ); HF2( 1209, y0, v0 );
      HF2( 1210, x0, y0 );

      for( int ih=0; ih<nh; ++ih ){
	DCLTrackHit *hit=tp->GetHit(ih);
	int layerId=hit->GetLayer();
	// 	  HF1( 33, layerId );
	}
      if( chisqr<MaxChisqrSdcOut ) ++ntOk;
      else 
	{
	  //tp->GoodForTracking( false );
	}
    }
  //    HF1( 30, double(ntOk) );
   
    if( ntOk<1 ) {
      { tree->Fill(); return true; }
    }
  }
  HF1( 1, 15. );

  HF1( 1, 20. );
  
  
  //std::cout<<"***************************************"<<std::endl;
  //////////////SKSTracking Analysis
  DCAna->TrackSearchSks();
  int ntSks=DCAna->GetNTracksSks();
  event.ntSks=ntSks;
  // std::cout<<ntSks<<std::endl;
  HF1( 3001, double(ntSks) );
  //------------------------Cut
  if( ntSks<1 || ntSks>MaxNTrackSks ) {
    { tree->Fill(); return true; }
  }
  HF1( 1, 21. );

  {
    int ntOk=0;
    for( int i=0; i<ntSks; ++i ){
      SksTrack *track=DCAna->GetSksTrack(i);
      if(!track) continue;
      int nh=track->GetNHits();
      double chisqr=track->chisqr();
      ThreeVector Ppos=track->PrimaryPosition();
      ThreeVector Pmom=track->PrimaryMomentum();
      double path=track->PathLengthToTOF();
      double xt=Ppos.x(), yt=Ppos.y();
      double p=Pmom.Mag();
      double ut=Pmom.x()/p, vt=Pmom.y()/p;

      if(i<MaxHits){
	event.nhSks[i] = nh;
	event.chisqrSks[i] = chisqr;
	event.path[i] = path;
	event.pSks[i] = p;
	event.xts[i] = xt;
	event.yts[i] = yt;
	event.uts[i] = ut;
	event.vts[i] = vt;
      }

      HF1( 3002, double(nh) );
      HF1( 3003, chisqr );
      HF1( 3004, xt ); HF1( 3005, yt ); 
      HF1( 3006, ut ); HF1( 3007, vt );
      HF2( 3008, xt, ut ); HF2( 3009, yt, vt );
      HF2( 3010, xt, yt );
      HF1( 3011, p ); 
      HF1( 3012, path );

      double m2;
      for ( int j=0; j<ncTof; ++j){
	HodoCluster *clTof=hodoAna->GetClusterTOF(j);
	if (!clTof || !clTof->GoodForAnalysis()) continue;
	m2=MassSquare( p, path, clTof->CMeanTime()-time0+offset);
	if (i<MaxHits){
	  event.m2[i] = m2;
	}
	HF1 (3013, m2);
      }

      //------------------------Cut
      //if( chisqr<MaxChisqrSks ){
      if( MinMassSquare<m2 && m2<MaxMassSquare && chisqr<MaxChisqrSks ){
	++ntOk;
	HF1( 3102, double(nh) );
	HF1( 3103, chisqr );
	HF1( 3104, xt ); HF1( 3105, yt ); 
	HF1( 3106, ut ); HF1( 3107, vt );
	HF2( 3108, xt, ut ); HF2( 3109, yt, vt );
	HF2( 3110, xt, yt );
	HF1( 3111, p ); 
	HF1( 3112, path );

	for( int j=0; j<ncTof; ++j ){
	  HodoCluster *clTof=hodoAna->GetClusterTOF(j);
	  if( !clTof || !clTof->GoodForAnalysis() ) continue;
	  double m2=MassSquare( p, path, clTof->CMeanTime()-time0+offset );
	  HF1( 3113, m2 );
	  //std::cout<<"Time="<< clTof->CMeanTime()-time0+offset <<std::endl;
	  //std::cout<<"******* m2= "<< m2 <<std::endl;
	}
      }
      else {
	track->GoodForAnalysis( false );
      }
    }
    HF1( 3101, double(ntOk) );
    //------------------------Cut
    if( ntOk<1 ) {
      { tree->Fill(); return true; }
    }
  }
  
  HF1( 1, 22. );

  HF1( 1, 30. );

  //////////////BcOut Analysis
  /*
  {
    std::vector<std::vector<DCHitContainer> > bcOutCandidates;
    BH2Filter& gFilter_bh2 = BH2Filter::GetInstance();
    //    gFilter_bh2.Apply(*hodoAna, *DCAna, bcOutCandidates);
    int T0Seg = int(clBH2Time0->MeanSeg());
    gFilter_bh2.Apply_withT0Seg(T0Seg, *DCAna, bcOutCandidates);
    DCAna->TrackSearchBcOut(bcOutCandidates);
  }
  */
  DCAna->TrackSearchBcOut();
  
  int ntBcOut=DCAna->GetNtracksBcOut();
  event.ntBcOut=ntBcOut;
  HF1( 1101, double(ntBcOut) );
  //std::cout<<"BcOut Num:   "<<ntBcOut<<std::endl;
  //------------------------Cut
  if( ntBcOut<1 || ntBcOut>MaxNTrackBcOut ) {
    { tree->Fill(); return true; }
  }
  HF1( 1, 31. );

  {
    int ntOk=0;
    for( int i=0; i<ntBcOut; ++i ){
      DCLocalTrack *tp=DCAna->GetTrackBcOut(i);
      if(!tp) continue;
      int nh=tp->GetNHit();
      double chisqr=tp->GetChiSquare();
      double u0=tp->GetU0(), v0=tp->GetV0();
      double x0=tp->GetX(0.), y0=tp->GetY(0.);
      double xtgt=tp->GetX( zK18Tgt ), ytgt=tp->GetY( zK18Tgt );
      
      if( i<MaxHits ){
	event.nhBcOut[i]=nh;
	event.chisqrBcOut[i]=chisqr;
	event.u0BcOut[i]=u0;
	event.v0BcOut[i]=v0;
	event.x0BcOut[i]=x0;
	event.y0BcOut[i]=y0;
	event.xtgtBcOut[i]=xtgt;
	event.ytgtBcOut[i]=ytgt;
      }
      
      HF1( 1102, double(nh) );
      HF1( 1103, chisqr );
      HF1( 1104, x0 ); HF1( 1105, y0 );
      HF1( 1106, u0 ); HF1( 1107, v0 );
      HF2( 1108, x0, u0 ); HF2( 1109, y0, v0 );
      HF2( 1110, x0, y0 );
      HF1( 1111, xtgt ); HF1( 1112, ytgt );
      HF2( 1113, xtgt, ytgt );
      //------------------------Cut
      if( chisqr<MaxChisqrBcOut ){
	++ntOk;
      }
    }
    if( ntOk<1 ) {
      { tree->Fill(); return true; }
    }
  }
  HF1( 1, 32. );
  
#if BFT
  // BFT analysis
 hodoAna->DecodeBFTHits(rawData);

  //BFT TimeCut
 hodoAna->TimeCutBFT( -8.,8.);
 HF1 ( 1, 33.);
 std::vector<double> cont_x;
 {
    //ParamMan& BH1Cut = FileManager::GetInstance().GetParamMan("CutBh12");
    //BH1Cut.stChangeParamSet("PositionCutBFT");
    //int clpos_bh1 = static_cast<int>(2*clBH1Time0->MeanSeg());
    //double min    = BH1Cut.GetParameter(clpos_bh1, "min");
    //double max    = BH1Cut.GetParameter(clpos_bh1, "max");
    double min = -160. , max = 160.; 
    int nc = hodoAna->GetNClustersBFT();
    for(int i = 0; i<nc; ++i){
      FiberCluster* cl = hodoAna->GetClusterBFT(i);
      double ctime = cl->CMeanTime();
      double x = cl->MeanPosition();
      if(min < x && x < max){
	event.cxpos[i] = x;
	cont_x.push_back(x);
	event.ctime_BFT[i]=ctime;
      }
    }
    
    if(0 == cont_x.size()){ return true;}
    int ncBFT = cont_x.size();
    event.ncBFT = ncBFT;
  }
  if (cont_x.size()>5)
    {return true;}
  HF1( 1, 35. );
  
  // K18 with D2U
  DCAna->TrackSearchK18_BFT(cont_x);
  int ntK18=DCAna->GetNTracksK18_BFT();
  event.ntK18 = ntK18;
  HF1( 2201, double(ntK18) );
  if(ntK18 > 20){
    std::cout << "#W : " 
	      << " Too many K18 tracks : ntK18 = "
	      << ntK18
	      << std::endl;
    ntK18 = 20;
  }
  //------------------------Cut  
  if( ntK18<1 || ntK18>MaxNTrackK18) {
    { tree->Fill(); return true; }
  }
  
  HF1( 1, 36. );

  for( int i=0; i<ntK18; ++i ){
    K18Track_BFT *tp=DCAna->GetK18Track_BFT(i);
    if(!tp) continue;
    double pk18=tp->PD2U_3rd();
    event.pK18[i] = pk18;
    HF1( 2204, pk18 );
  }
#endif

  HF1( 1, 37. );

  HF1( 1, 40. );

  std::vector <ThreeVector> PionPCont, PionXCont;
  std::vector <ThreeVector> KaonPCont, KaonXCont;
    
  //////////////Pion cut
  for( int itSks=0; itSks<ntSks; ++itSks ){
    SksTrack *track=DCAna->GetSksTrack(itSks);
    if( !track || !track->GoodForAnalysis() ) continue;
    DCLocalTrack *trIn =track->GetLocalTrackIn();
    DCLocalTrack *trOut=track->GetLocalTrackOut();
    HodoCluster *clTof=0;  
    double xtof=trOut->GetX(zTof), ytof=trOut->GetY(zTof);
    double mindifTof=MinXDifTof;
    for( int j=0; j<ncTof; ++j ){
      HodoCluster *cl=hodoAna->GetClusterTOF(j);
      //if( !cl || !cl->GoodForAnalysis() ) continue;
      double x=(cl->MeanSeg()-15.5)*70.;
      //if( fabs(x-xtof)<mindifTof ){
      clTof=cl; //mindifTof=fabs(x-xtof);
      //}
      //}
    
    if( !clTof ) continue;
    int nh=track->GetNHits();
    double chisqr=track->chisqr();
    ThreeVector Ppos=track->PrimaryPosition();
    ThreeVector Pmom=track->PrimaryMomentum();
    //Calibration
    double p0 = Pmom.Mag()+pSKS_offset+pS_offset;
 
    double u0 = Pmom.x()/p0, v0 = Pmom.y()/p0;
    double pt = p0/sqrt(1.+u0*u0+v0*v0);

    ThreeVector PposCorr( Ppos.x()+x_off, Ppos.y()+y_off, Ppos.z() );
    ThreeVector PmomCorr( pt*(u0+u_off), pt*(v0+v_off), pt);

    double pathL=track->PathLengthToTOF();
    double xt=PposCorr.x(), yt=PposCorr.y();
    double p=PmomCorr.Mag();
    double ut=PmomCorr.x()/p, vt=PmomCorr.y()/p;
    double m2=MassSquare( p0, pathL, clTof->CMeanTime()-time0+offset );
    HF1( 3202, double(nh) );
    HF1( 3203, chisqr );
    HF1( 3204, xt ); HF1( 3205, yt ); 
    HF1( 3206, ut ); HF1( 3207, vt );
    HF2( 3208, xt, ut ); HF2( 3209, yt, vt );
    HF2( 3210, xt, yt );
    HF1( 3211, p ); 
    HF1( 3212, pathL );
    HF1( 3213, m2 );

    //------------------------Cut    
    if( MinMassSquare<m2 && m2<MaxMassSquare ){
      HF1( 4202, double(nh) );
      HF1( 4203, chisqr );
      HF1( 4204, xt ); HF1( 4205, yt ); 
      HF1( 4206, ut ); HF1( 4207, vt );
      HF2( 4208, xt, ut ); HF2( 4209, yt, vt );
      HF2( 4210, xt, yt );
      HF1( 4211, p ); 
      HF1( 4212, pathL );
      PionPCont.push_back(PmomCorr); PionXCont.push_back(PposCorr);
    }
  }
  if( PionPCont.empty() ) {
    { tree->Fill(); return true; }
  }

  HF1( 1, 41. );

  //Kaon cut  
  for( int itK18=0; itK18<ntK18; ++itK18 ){
#if BFT
    K18Track_BFT *track=DCAna->GetK18Track_BFT(itK18);
    if( !track || !track->GoodForAnalysis() ) continue;
    DCLocalTrack *trOut=track->TrackOut();

    //Calibration
    double p_offset = pK18_offset;

    double p=track->PD2U()+p_offset, x=track->Xtgt(), y=track->Ytgt();
    //double p=track->PD2U_3rd()*pB_offset+p_offset, x=track->Xtgt(), y=track->Ytgt();
    double u=track->Utgt(), v=track->Vtgt();
    double loss_bh2 = 1.09392e-3;
    p = p - loss_bh2;
    
    double pt=p/sqrt(1.+u*u+v*v);
    ThreeVector Pos( x, y, 0. );
    ThreeVector Mom( pt*u, pt*v, pt );
    double xo=trOut->GetX0(), yo=trOut->GetY0();

    HodoCluster *clBh1=0;
    for( int j=0; j<ncBh1; ++j ){
      HodoCluster *cl=hodoAna->GetClusterBH1(j);
      clBh1=cl;
    }
    
    if(!clBh1) continue;

    
    HF1( 4104, p ); 
    HF1( 4105, x ); HF1( 4106, y ); 
    HF1( 4107, xo ); HF1( 4108, yo ); HF1( 4109, u ); HF1( 4110, v );

    KaonPCont.push_back(Mom); KaonXCont.push_back(Pos);
#endif
  }
  //------------------------Cut
  if( KaonPCont.empty() ) {
    { tree->Fill(); return true; }
  }
  HF1( 1, 42. );

  //MissingMass
  int nK=KaonPCont.size(), nPi=PionPCont.size();

  event.nK=nK;
  event.nPi=nPi;
  event.nKPi=nK*nPi;
  HF1( 4101, double(nPi));
  HF1( 4201, double(nK) );

  int nkpi=0;
  for( int ipi=0; ipi<nPi; ++ipi ){
    ThreeVector ppi=PionPCont[ipi], xpi=PionXCont[ipi];
    for( int ik=0; ik<nK; ++ik ){
      
      ThreeVector pk=KaonPCont[ik], xk=KaonXCont[ik];

      ThreeVector vert=VertexPoint( xk, xpi, pk, ppi );
      double closedist=closeDist( xk, xpi, pk, ppi );

      double us=ppi.x()/ppi.z(), vs=ppi.y()/ppi.z();
      double ub=pk.x()/pk.z(), vb=pk.y()/pk.z();
      double cost=pk*ppi/(pk.Mag()*ppi.Mag());

      //For correction @ QM02 for narrow region, 4th pol for V
      double pu0 =  0.00416686;
      double pu1 =  0.0142352;
      double pu2 = -0.00556256;
      double pu3 = -0.18308;
      double pu4 =  0.472155;
      double pu5 = -1.79675;

      //      double pv0 =  0.00369582;
      //      double pv1 =  0.0137256;
      //      double pv2 =  0.429082;
      //      double pv3 =  0.238147;
      //      double pv4 =-37.2991;
      double pv0 =  0.00382395;
      double pv1 =  0.0144088;
      double pv2 =  0.23244;

      double ppi0= ppi.Mag();
      double pCorr = ppi0 + 
	(( pu1*us + pu2*us*us + pu3*us*us*us + pu4*us*us*us*us +pu5*us*us*us*us*us) +
	 //	 (  pv1*vs + pv2*vs*vs + pv3*vs*vs*vs + pv4*vs*vs*vs*vs ))
	 (  pv1*vs + pv2*vs*vs))
	+ (pu0+pv0);
      
      ThreeVector ppiCorr(pCorr*ppi.x()/ppi.Mag(),
			 pCorr*ppi.y()/ppi.Mag(),
			 pCorr*ppi.z()/ppi.Mag());
      
      ThreeVector pkCorrDE = CorrElossIn( pk, xk, vert, KaonMass);
      
      ThreeVector ppiCorrDE  = CorrElossOut( ppiCorr, xpi, vert, PionMass);
      //ThreeVector ppiCorrDE  = CorrElossIn( pk, xk, vert, KaonMass);
      
      LorentzVector LvK( pk, sqrt(KaonMass*KaonMass+pk.Mag2()) );
      LorentzVector LvKCorrDE( pkCorrDE, sqrt(KaonMass*KaonMass+pkCorrDE.Mag2()) );

      LorentzVector LvPi( ppi, sqrt(PionMass*PionMass+ppi.Mag2()) );
      LorentzVector LvPiCorr( ppiCorr, sqrt(PionMass*PionMass+ppiCorr.Mag2()) );
      LorentzVector LvPiCorrDE( ppiCorrDE, sqrt(PionMass*PionMass+ppiCorrDE.Mag2()) );
      //added by me
      LorentzVector LvMuon( ppi, sqrt(MuonMass*MuonMass+ppi.Mag2()));


      LorentzVector LvCore1( 0., 0., 0., CoreMass1 );
      LorentzVector LvCore2( 0., 0., 0., CoreMass2 );
      LorentzVector LvCore3( 0., 0., 0., CoreMass3 );
      LorentzVector LvCore4( 0., 0., 0., CoreMass4 );

      LorentzVector LvT1( 0., 0., 0., TargetMass1 );
      LorentzVector LvT2( 0., 0., 0., TargetMass2 );
      LorentzVector LvT3( 0., 0., 0., TargetMass3 );
      LorentzVector LvT4( 0., 0., 0., TargetMass4 );
      LorentzVector LvT0( 0., 0., 0., 0.);

      LorentzVector LvRc1 = LvK+LvT1-LvPi;
      LorentzVector LvRc2 = LvK+LvT2-LvPi;
      LorentzVector LvRc3 = LvK+LvT3-LvPi;
      LorentzVector LvRc4 = LvK+LvT4-LvPi;
      LorentzVector LvRc_pion = LvK+LvT0-LvPi;
      LorentzVector LvRc_muon = LvK+LvT0-LvMuon; 

      LorentzVector LvRcCorr = LvK+LvT1-LvPiCorr;
      LorentzVector LvRcCorrDE = LvKCorrDE+LvT1-LvPiCorrDE;
      double MisMass1 = LvRc1.Mag();//-LvC.Mag();
      double MisMass2 = LvRc2.Mag();
      double MisMass3 = LvRc3.Mag();
      double MisMass4 = LvRc4.Mag();

      double MisMass_pion = LvRc_pion.Mag();
      double MisMass_muon = LvRc_muon.Mag();
      double MisMassCorr = LvRcCorr.Mag();//-LvC.Mag();
      double MisMassCorrDE = LvRcCorrDE.Mag();//-LvC.Mag();

      double BE1=LvRc1.Mag()-( LvCore1.Mag()+LambdaMass );
      double BE2=LvRc2.Mag()-( LvCore2.Mag()+LambdaMass );
      double BE3=LvRc3.Mag()-( LvCore3.Mag()+LambdaMass );
      double BE4=LvRc4.Mag()-( LvCore4.Mag()+LambdaMass );

      double BECorr=LvRcCorr.Mag()-( LvCore1.Mag()+PionMass );
      double BECorrDE=LvRcCorrDE.Mag()-( LvCore1.Mag()+PionMass );
      
      std::cout<<"******* Missing Mass1= "<< MisMass1 <<std::endl; 

      //Primary frame
      LorentzVector PrimaryLv =  LvK+LvT1;
      double TotalEnergyCM = PrimaryLv.Mag();
      ThreeVector beta( PrimaryLv.Vect()/PrimaryLv.E() );
      
      //CM
      double TotalMomCM 
	= 0.5*sqrt(( TotalEnergyCM*TotalEnergyCM
		     -( PionMass+SigmaPMass )*( PionMass+SigmaPMass ))
		   *( TotalEnergyCM*TotalEnergyCM
		      -( PionMass-SigmaPMass )*( PionMass-SigmaPMass )))/TotalEnergyCM;
      
      double costLab = cost;
      double cottLab = costLab/sqrt(1.-costLab*costLab);
      double bt=beta.Mag(), gamma=1./sqrt(1.-bt*bt);
      double gbep=gamma*bt*sqrt(TotalMomCM*TotalMomCM+PionMass*PionMass)/TotalMomCM;
      double a  = gamma*gamma+cottLab*cottLab;
      double bp = gamma*gbep;
      double c  = gbep*gbep-cottLab*cottLab;
      double dd=bp*bp-a*c;
      
      if( dd<0. ){
	std::cerr << "dd<0." << std::endl;
	dd=0.;
	//exit(-1);
      }
      
      double costCM=(sqrt(dd)-bp)/a;
      if( costCM>1. || costCM<-1. ){
	std::cerr << "costCM>1. || costCM<-1." << std::endl;
	//exit(-1);
	costCM=-1.;
      }
      double sintCM=sqrt(1.-costCM*costCM);
      double PionMom = TotalMomCM*sintCM/sqrt(1.-costLab*costLab);

      if (nkpi<MaxHits) {
	event.vtx[nkpi]=vert.x();
	event.vty[nkpi]=vert.y();
	event.vtz[nkpi]=vert.z();
	event.closeDist[nkpi]=closedist;
	event.theta[nkpi]=acos(cost)*Rad2Deg;
	event.theta_CM[nkpi] =acos(costCM)*Rad2Deg;
	event.cost_CM[nkpi] =costCM;

	event.MissMass1[nkpi]=MisMass1;
	event.MissMass2[nkpi]=MisMass2;
	event.MissMass3[nkpi]=MisMass3;
	event.MissMass4[nkpi]=MisMass4;
	event.MissMassCorr[nkpi] = MisMassCorr;
	event.MissMassCorrDE[nkpi] = MisMassCorrDE;
	event.MissMassMuon[nkpi] = MisMass_muon;
	event.MissMassPion[nkpi] = MisMass_pion;
	event.BE1[nkpi]=BE1;
	event.BE2[nkpi]=BE2;
	event.BE3[nkpi]=BE3;
	event.BE4[nkpi]=BE4;

	event.xk[nkpi] = xk.x();
	event.yk[nkpi] = xk.y();
	event.uk[nkpi] = us;
	event.vk[nkpi] = vs;

	event.xpi[nkpi] = xpi.x();
	event.ypi[nkpi] = xpi.y();
	event.upi[nkpi] = ub;
	event.vpi[nkpi] = vb;
	event.pOrg[nkpi] = ppi0;
	event.pCalc[nkpi] = PionMom;
	event.pCorr[nkpi] = pCorr;
	event.pCorrDE[nkpi] = pkCorrDE.Mag();
	nkpi++;
      }
      
      HF1( 5001, vert.z() );

      HF1( 5002, MisMass1 );
      HF2( 5011, MisMass1, us );
      HF2( 5012, MisMass1, vs );
      HF2( 5013, MisMass1, ub );
      HF2( 5014, MisMass1, vb );
    }
  }  

  HF1( 1, 50. );

  tree->Fill();
  
  //Final Hodoscope histograms
  for( int i=1; i<ncBh2; ++i ){
    BH2Cluster *cl=hodoAna->GetClusterBH2(i);
    double t=cl->CMeanTime();
    double dEbh2 = cl->DeltaE();
    if( !(MinDeltaEBH2<dEbh2 && dEbh2<MaxDeltaEBH2) ) continue;

    HF1( 152, cl->ClusterSize() );
    HF1( 153, cl->MeanSeg()+1-0.5 );
    HF1( 154, t );
    HF1( 155, cl->DeltaE() );
  }

  for( int i=0; i<ncBh1; ++i ){
    HodoCluster *cl=hodoAna->GetClusterBH1(i);
    double btof= (cl->CMeanTime())-time0;
    if( !cl || !cl->GoodForAnalysis() ) continue;

    HF1( 252, cl->ClusterSize() );
    HF1( 253, cl->MeanSeg()+1-0.5 );
    HF1( 254, cl->CMeanTime() );
    HF1( 255, cl->DeltaE() );
    HF1( 256, btof );
  }

  for( int i=0; i<ncTof; ++i ){
    HodoCluster *cl=hodoAna->GetClusterTOF(i);
    double de=cl->DeltaE();
    double t=cl->CMeanTime()-time0;
    if( !cl || !cl->GoodForAnalysis() ) continue;

    HF1( 352, cl->ClusterSize() );
    HF1( 353, cl->MeanSeg()+1-0.5 );
    HF1( 354, t );
    HF1( 355, de );
  }

  }
  
  
  { tree->Fill(); return true; }
}

void EventKPiAna::InitializeEvent( void )
{
  event.runnum = -1;
  event.evnum = -1;
  //Trigger
  for( int it=0; it<NumOfMisc; it++){
    event.trigpat[it] = -1;
    event.trigflag[it] = -1;
  }

  for( int it=0; it<NumOfMatrix; it++){
    event.matrixpat[it] = -1;
    event.matrixflag[it] = -1;
  }

  //Hodoscope  
  event.nhitsBh2  = -1;
  event.nhBH2  = 0;
  event.nhBH1  = 0;
  event.nhTof  = 0;
  event.nhLc   = 0;
  event.multip_bh1 = 0;
  event.multip_bh2 = 0;

  for( int it=0; it<MaxHits; it++){
    event.BH2Seg[it] = -1;
    event.tBH2[it] = -9999.0;
    event.deBH2[it] = -9999.0;

    event.BH1Seg[it] = -1;
    event.tBH1[it] = -9999.0;
    event.deBH1[it] = -9999.0;
    event.btof[it]  = -1;

    event.TofSeg[it] = -1;
    event.tTof[it] = -9999.0;
    event.dtTof[it] = -9999.0;
    event.deTof[it] = -9999.0;

    event.LcSeg[it] = -1;
    event.tLc[it] = -9999.0;
    event.dtLc[it] = -9999.0;
    event.deLc[it] = -9999.0;
  }

  for(int i = 0; i<NumOfSegAC1; ++i){
    event.Ac1Tdc[i] = -1;
  }

  //DC
  event.ntBcIn   = 0;
  event.ntBcOut  = 0;
  event.ntSdcIn  = 0;
  event.ntSdcOut = 0;
  event.ntK18    = 0;
  event.ntSks    = 0;
  event.ncBFT    = 0;

  //Beam DC    
  for( int it=0; it<MaxHits; it++){
    event.nhBcIn[it]     = -1;
    event.chisqrBcIn[it] = -1.0;
    event.x0BcIn[it] = -9999.0;
    event.y0BcIn[it] = -9999.0;
    event.u0BcIn[it] = -9999.0;
    event.v0BcIn[it] = -9999.0;
    
    event.ctime_BFT[it] = -9999.0;
    event.xbh1BcIn[it] = -9999.0;
    event.ybh1BcIn[it] = -9999.0;

    event.cxpos[it]    = -999.0;
  }

  for( int it=0; it<MaxHits; it++){
    event.nhBcOut[it]     = -1;
    event.chisqrBcOut[it] = -1.0;
    event.x0BcOut[it] = -9999.0;
    event.y0BcOut[it] = -9999.0;
    event.u0BcOut[it] = -9999.0;
    event.v0BcOut[it] = -9999.0;

    event.xtgtBcOut[it] = -9999.0;
    event.ytgtBcOut[it] = -9999.0;
    event.xbh2BcOut[it] = -9999.0;
    event.ybh2BcOut[it] = -9999.0;
  }

  for( int it=0; it<MaxHits; it++){
    event.nhK18[it]     = -1;
    event.chisqrK18[it] = -1.0;
    event.xtb[it] = -9999.0;
    event.ytb[it] = -9999.0;
    event.utb[it] = -9999.0;
    event.vtb[it] = -9999.0;

    event.pK18[it]   = -9999.0;
    event.thetab[it] = -9999.0;
  }

  //SKS DC
  for( int it=0; it<MaxHits; it++){
    event.nhSdcIn[it]     = -1;
    event.chisqrSdcIn[it] = -1.0;
    event.x0SdcIn[it] = -9999.0;
    event.y0SdcIn[it] = -9999.0;
    event.u0SdcIn[it] = -9999.0;
    event.v0SdcIn[it] = -9999.0;
  }


  event.nhBAC = 0;
  for( int it=0; it<MaxHits; it++)
    {
      event.BACSeg[it] = -1;
      event.tBAC[it] = -9999.;
      event.aBAC[it] = -9999.;
    }

  for( int it=0; it<MaxHits; it++){
    event.nhSdcOut[it]     = -1;
    event.chisqrSdcOut[it] = -1.0;
    event.u0SdcOut[it] = -9999.0;
    event.v0SdcOut[it] = -9999.0;

    event.xtofSdcOut[it] = -9999.0;
    event.ytofSdcOut[it] = -9999.0;
  }

  for( int it=0; it<MaxHits; it++){
    event.nhSks[it]     = -1;
    event.chisqrSks[it] = -1.0;
    event.xts[it] = -9999.0;
    event.yts[it] = -9999.0;
    event.uts[it] = -9999.0;
    event.vts[it] = -9999.0;

    event.pSks[it] = -9999.0;
    event.path[it] = -9999.0;
    event.m2[it]   = -9999.0;
    event.thetas[it] = -9999.0;
  }

  //Reaction
  event.nPi      = -1;
  event.nK       = -1;
  event.nKPi     = 0;

  for( int it=0; it<MaxHits; it++){
    event.vtx[it]       = -9999.0; 
    event.vty[it]       = -9999.0; 
    event.vtz[it]       = -9999.0; 
    event.closeDist[it] = -9999.0; 
    event.theta[it]     = -9999.0; 
    event.theta_CM[it]  = -9999.0; 
    event.cost_CM[it]   = -9999.0; 
    event.MissMass1[it]  = -9999.0;
    event.MissMass2[it]  = -9999.0;
    event.MissMass3[it]  = -9999.0;
    event.MissMass4[it]  = -9999.0;
    event.MissMassMuon[it] = -9999.0;
    event.MissMassPion[it] = -9999.0;    
    event.MissMassCorr[it]  = -9999.0; 
    event.MissMassCorrDE[it]  = -9999.0;
    event.BE1[it]        = -9999.0; 
    event.BE2[it]        = -9999.0; 
    event.BE3[it]        = -9999.0; 
    event.BE4[it]        = -9999.0; 

    event.xk[it] = -9999.0;
    event.yk[it] = -9999.0;
    event.uk[it] = -9999.0;
    event.vk[it] = -9999.0;
    event.xpi[it] = -9999.0;
    event.ypi[it] = -9999.0;
    event.upi[it] = -9999.0;
    event.vpi[it] = -9999.0;

    event.pOrg[it] = -9999.0;
    event.pCalc[it] = -9999.0;
    event.pCorr[it] = -9999.0;
    event.pCorrDE[it] = -9999.0;
  }

}

bool EventKPiAna::ProcessingEnd()
{
  return true;
}

VEvent *ConfMan::EventAllocator()
{
  return new EventKPiAna;
}

bool ConfMan:: InitializeHistograms()
{  
  HB1( 1, "Status", 60, 0., 60. );

  HB1( 101, "#Clusters BH2",   7, 0., 7. );
  HB1( 102, "ClusterSize BH2", 7, 0., 7. );
  HB1( 103, "HitPat BH2", 10, 0., 10. );
  HB1( 104, "MeanTime BH2", 200, -10., 10. );
  HB1( 105, "Delta-E BH2", 200, -0.5, 4.5 );

  HB1( 112, "ClusterSize BH2", 7, 0., 7. );
  HB1( 113, "HitPat BH2", 10, 0., 10. );
  HB1( 114, "MeanTime BH2", 200, -10., 10. );
  HB1( 115, "Delta-E BH2", 200, -0.5, 4.5 );

  HB1( 122, "ClusterSize BH2 [T0]", 7, 0., 7. );
  HB1( 123, "HitPat BH2 [T0]", 10, 0., 10. );
  HB1( 124, "MeanTime BH2 [T0]", 200, -10., 10. );
  HB1( 125, "Delta-E BH2 [T0]", 200, -0.5, 4.5 );

  HB1( 152, "ClusterSize BH2 [piK]", 7, 0., 7. );
  HB1( 153, "HitPat BH2 [piK]", 10, 0., 10. );
  HB1( 154, "MeanTime BH2 [piK]", 200, -10., 10. );
  HB1( 155, "Delta-E BH2 [piK]", 200, -0.5, 4.5 );

  HB1( 201, "#Clusters BH1",  11, 0., 11. );
  HB1( 202, "ClusterSize BH1",11, 0., 11. );
  HB1( 203, "HitPat BH1", 11, 0., 11. );
  HB1( 204, "MeanTime BH1", 200, -10., 10. );
  HB1( 205, "Delta-E BH1", 200, -0.5, 4.5 );
  HB1( 206, "Beam ToF", 200, -10., 10. );

  HB1( 211, "#Clusters BH1 [pi]",  11, 0., 11. );
  HB1( 212, "ClusterSize BH1 [pi]",11, 0., 11. );
  HB1( 213, "HitPat BH1 [pi]", 11, 0., 11. );
  HB1( 214, "MeanTime BH1 [pi]", 200, -10., 10. );
  HB1( 215, "Delta-E BH1 [pi]", 200, -0.5, 4.5 );

  HB1( 252, "ClusterSize BH1 [piK]",11, 0., 11. );
  HB1( 253, "HitPat BH1 [piK]", 11, 0., 11. );
  HB1( 254, "MeanTime BH1 [piK]", 200, -10., 10. );
  HB1( 255, "Delta-E BH1 [piK]", 200, -0.5, 4.5 );
  HB1( 256, "Beam ToF [piK]", 200, -10., 10. );

  HB1( 301, "#Clusters Tof",  32, 0., 32. );
  HB1( 302, "ClusterSize Tof",32, 0., 32. );
  HB1( 303, "HitPat Tof", 32, 0., 32. );
  HB1( 304, "TimeOfFlight Tof", 500, -50., 100. );
  HB1( 305, "Delta-E Tof", 200, -0.5, 4.5 );

  HB1( 311, "#Clusters Tof [Good]",  32, 0., 32. );
  HB1( 312, "ClusterSize Tof [Good]",32, 0., 32. );
  HB1( 313, "HitPat Tof [Good]", 32, 0., 32. );
  HB1( 314, "TimeOfFlight Tof [Good]", 500, -50., 100. );
  HB1( 315, "Delta-E Tof [Good]", 200, -0.5, 4.5 );

  HB1( 352, "ClusterSize Tof [piK]",32, 0., 32. );
  HB1( 353, "HitPat Tof [piK]", 32, 0., 32. );
  HB1( 354, "TimeOfFlight Tof [piK]", 500, -50., 100. );
  HB1( 355, "Delta-E Tof [piK]", 200, -0.5, 4.5 );

  HB1( 401, "#Clusters Lc",  28, 0., 28. );
  HB1( 402, "ClusterSize Lc",28, 0., 28. );
  HB1( 403, "HitPat Lc", 28, 0., 28. );
  HB1( 404, "TimeOfFlight Lc", 500, -50., 200. );
  HB1( 405, "Delta-E Lc", 200, -0.5, 4.5 );

  HB1( 411, "#Clusters Lc [Good]",  28, 0., 28. );
  HB1( 412, "ClusterSize Lc [Good]",28, 0., 28. );
  HB1( 413, "HitPat Lc [Good]", 28, 0., 28. );
  HB1( 414, "TimeOfFlight Lc [Good]", 500, -50., 200. );
  HB1( 415, "Delta-E Lc [Good]", 200, -0.5, 4.5 );

  HB1( 452, "ClusterSize Lc [piK]",28, 0., 28. );
  HB1( 453, "HitPat Lc [piK]", 28, 0., 28. );
  HB1( 454, "TimeOfFlight Lc [piK]", 500, -50., 200. );
  HB1( 455, "Delta-E Lc [piK]", 200, -0.5, 4.5 );

  HB2( 501, "SegLC%SegTOF", 32, 0., 32., 28, 0., 28. );
  HB2( 502, "TLc%TTof", 100,-50., 100., 100, -50., 200. );
  HB1( 503, "TLc-TTof", 500, -30., 100. );
  HB1( 504, "TofSeg-LCSeg", 25, -10., 15. );

  HB1( 510, "#GoodCombi Tof-LC", 20, 0., 20. ); 
  HB2( 511, "SegLC%SegTOF [Good]", 32, 0., 32., 28, 0., 28. );
  HB2( 512, "TLc%TTof [Good]", 100,-20., 100., 100, -20., 100. );
  HB1( 513, "TLc-TTof [Good]", 500, -5., 100. );
  HB1( 514, "TofSeg-LCSeg [Good]", 25, -10., 15. );

  HB2( 551, "SegLC%SegTOF [piK]", 32, 0., 32., 28, 0., 28. );
  HB2( 552, "TLc%TTof [piK]", 100,-20., 100., 100, -20., 100. );
  HB1( 553, "TLc-TTof [piK]", 500, -5., 100. );
  HB1( 554, "TofSeg-LCSeg [piK]", 25, -10., 15. );

  HB1( 1001, "#Tracks SdcIn", 10, 0., 10. );
  HB1( 1002, "#Hits SdcIn", 20, 0., 20. );
  HB1( 1003, "Chisqr SdcIn", 200, 0., 100. );
  HB1( 1004, "X0 SdcIn", 500, -200., 200. ); 
  HB1( 1005, "Y0 SdcIn", 500, -200., 200. );
  HB1( 1006, "U0 SdcIn",  700, -0.35, 0.35 );
  HB1( 1007, "V0 SdcIn",  400, -0.20, 0.20 );
  HB2( 1008, "X0%U0 SdcIn", 120, -200., 200., 100, -0.35, 0.35 );
  HB2( 1009, "Y0%V0 SdcIn", 100, -200., 200., 100, -0.20, 0.20 );
  HB2( 1010, "X0%Y0 SdcIn", 100, -200., 200., 100, -200., 200. );

  HB1( 1101, "#Tracks BcOut", 10, 0., 10. );
  HB1( 1102, "#Hits BcOut", 20, 0., 20. );
  HB1( 1103, "Chisqr BcOut", 200, 0., 100. );
  HB1( 1104, "X0 BcOut", 500, -200., 200. ); 
  HB1( 1105, "Y0 BcOut", 500, -200., 200. );
  HB1( 1106, "U0 BcOut",  700, -0.35, 0.35 );
  HB1( 1107, "V0 BcOut",  400, -0.20, 0.20 );
  HB2( 1108, "X0%U0 BcOut", 120, -200., 200., 100, -0.35, 0.35 );
  HB2( 1109, "Y0%V0 BcOut", 100, -200., 200., 100, -0.20, 0.20 );
  HB2( 1110, "X0%Y0 BcOut", 100, -200., 200., 100, -200., 200. );
  HB1( 1111, "Xtgt BcOut", 500, -200., 200. ); 
  HB1( 1112, "Ytgt BcOut", 500, -200., 200. );
  HB2( 1113, "Xtgt%Ytgt BcOut", 100, -200., 200., 100, -200., 200. );

  HB1( 1201, "#Tracks SdcOut", 10, 0., 10. );
  HB1( 1202, "#Hits SdcOut", 20, 0., 20. );
  HB1( 1203, "Chisqr SdcOut", 200, 0., 100. );
  HB1( 1204, "X0 SdcOut", 600, -1200., 1200. ); 
  HB1( 1205, "Y0 SdcOut", 600, -600., 600. );
  HB1( 1206, "U0 SdcOut",  700, -0.35, 0.35 );
  HB1( 1207, "V0 SdcOut",  400, -0.20, 0.20 );
  HB2( 1208, "X0%U0 SdcOut", 120, -1200., 1200., 100, -0.35, 0.35 );
  HB2( 1209, "Y0%V0 SdcOut", 100,  -600.,  600., 100, -0.20, 0.20 );
  HB2( 1210, "X0%Y0 SdcOut", 100, -1200., 1200., 100, -600., 600. );

  HB1( 1301, "#Tracks BcIn", 10, 0., 10. );
  HB1( 1302, "#Hits BcIn", 20, 0., 20. );
  HB1( 1303, "Chisqr BcIn", 200, 0., 100. );
  HB1( 1304, "X0 BcIn", 500, -200., 200. ); 
  HB1( 1305, "Y0 BcIn", 500, -200., 200. );
  HB1( 1306, "U0 BcIn",  400, -0.35, 0.35 );
  HB1( 1307, "V0 BcIn",  200, -0.20, 0.20 );
  HB2( 1308, "X0%U0 BcIn", 120, -200., 200., 100, -0.35, 0.35 );
  HB2( 1309, "Y0%V0 BcIn", 100, -200., 200., 100, -0.20, 0.20 );
  HB2( 1310, "X0%Y0 BcIn", 100, -200., 200., 100, -200., 200. );

  HB1( 2201, "#Tracks K18", 10, 0., 10. );
  HB1( 2202, "#Hits K18", 30, 0., 30. );
  HB1( 2203, "Chisqr K18", 500, 0., 50. );
  HB1( 2204, "P K18", 1000, 0.5, 2.0 );
  HB1( 2251, "#Tracks K18 [Good]", 10, 0., 10. );

  HB1( 3001, "#Tracks Sks", 10, 0., 10. );
  HB1( 3002, "#Hits Sks", 30, 0., 30. );
  HB1( 3003, "Chisqr Sks", 500, 0., 500. );
  HB1( 3004, "Xtgt Sks", 500, -200., 200. );
  HB1( 3005, "Ytgt Sks", 500, -100., 100. );
  HB1( 3006, "Utgt Sks", 400, -0.35, 0.35 );
  HB1( 3007, "Vtgt Sks", 200, -0.20, 0.20 );
  HB2( 3008, "Xtgt%U Sks", 100, -200., 200., 100, -0.35, 0.35 );
  HB2( 3009, "Ytgt%V Sks", 100, -100., 100., 100, -0.20, 0.20 );
  HB2( 3010, "Xtgt%Ytgt Sks", 100, -200., 200., 100, -100., 100. );
  HB1( 3011, "P Sks", 200, 0.0, 1.0 );
  HB1( 3012, "PathLength Sks", 600, 3000., 6000. );
  HB1( 3013, "MassSqr", 600, -1.2, 1.2 );

  HB1( 3101, "#Tracks Sks [Good]", 10, 0., 10. );
  HB1( 3102, "#Hits Sks [Good]", 30, 0., 30. );
  HB1( 3103, "Chisqr Sks [Good]", 500, 0., 500. );
  HB1( 3104, "Xtgt Sks [Good]", 500, -200., 200. );
  HB1( 3105, "Ytgt Sks [Good]", 500, -100., 100. );
  HB1( 3106, "Utgt Sks [Good]", 700, -0.35, 0.35 );
  HB1( 3107, "Vtgt Sks [Good]", 400, -0.20, 0.20 );
  HB2( 3108, "Xtgt%U Sks [Good]", 100, -200., 200., 100, -0.35, 0.35 );
  HB2( 3109, "Ytgt%V Sks [Good]", 100, -100., 100., 100, -0.20, 0.20 );
  HB2( 3110, "Xtgt%Ytgt Sks [Good]", 100, -200., 200., 100, -100., 100. );
  HB1( 3111, "P Sks [Good]", 200, 0.0, 1.0 );
  HB1( 3112, "PathLength Sks [Good]", 600, 3000., 6000. );
  HB1( 3113, "MassSqr", 600, -1.2, 1.2 );

  HB1( 3201, "#Tracks Sks [Good2]", 10, 0., 10. );
  HB1( 3202, "#Hits Sks [Good2]", 30, 0., 30. );
  HB1( 3203, "Chisqr Sks [Good2]", 500, 0., 500. );
  HB1( 3204, "Xtgt Sks [Good2]", 500, -200., 200. );
  HB1( 3205, "Ytgt Sks [Good2]", 500, -100., 100. );
  HB1( 3206, "Utgt Sks [Good2]", 700, -0.35, 0.35 );
  HB1( 3207, "Vtgt Sks [Good2]", 400, -0.20, 0.20 );
  HB2( 3208, "Xtgt%U Sks [Good2]", 100, -200., 200., 100, -0.35, 0.35 );
  HB2( 3209, "Ytgt%V Sks [Good2]", 100, -100., 100., 100, -0.20, 0.20 );
  HB2( 3210, "Xtgt%Ytgt Sks [Good2]", 100, -200., 200., 100, -100., 100. );
  HB1( 3211, "P Sks [Good2]", 200, 0.0, 1.0 );
  HB1( 3212, "PathLength Sks [Good2]", 600, 3000., 6000. );
  HB1( 3213, "MassSqr", 600, -1.2, 1.2 );

  HB1( 4101, "#Tracks K18 [SksP]", 10, 0., 10. );
  HB1( 4102, "#Hits K18 [SksP]", 30, 0., 30. );
  HB1( 4103, "Chisqr K18 [SksP]", 500, 0., 50. );
  HB1( 4104, "P K18 [SksP]", 500, 0.5, 2.0 );
  HB1( 4105, "Xtgt K18 [SksP]", 500, -200., 200. );
  HB1( 4106, "Ytgt K18 [SksP]", 500, -100., 100. );
  HB1( 4107, "Xout K18 [SksP]", 500, -200., 200. ); 
  HB1( 4108, "Yout K18 [SksP]", 500, -100., 100. );
  HB1( 4109, "Uout K18 [SksP]", 400, -0.35, 0.35 );
  HB1( 4110, "Vout K18 [SksP]", 200, -0.20, 0.20 );
  HB1( 4111, "Xin  K18 [SksP]", 500, -100., 100. ); 
  HB1( 4112, "Yin  K18 [SksP]", 500, -100., 100. );
  HB1( 4113, "Uin  K18 [SksP]", 400, -0.35, 0.35 );
  HB1( 4114, "Vin  K18 [SksP]", 200, -0.20, 0.20 );

  HB1( 4201, "#Tracks Sks [Proton]", 10, 0., 10. );
  HB1( 4202, "#Hits Sks [Proton]", 30, 0., 30. );
  HB1( 4203, "Chisqr Sks [Proton]", 500, 0., 500. );
  HB1( 4204, "Xtgt Sks [Proton]", 500, -200., 200. );
  HB1( 4205, "Ytgt Sks [Proton]", 500, -100., 100. );
  HB1( 4206, "Utgt Sks [Proton]", 700, -0.35, 0.35 );
  HB1( 4207, "Vtgt Sks [Proton]", 400, -0.20, 0.20 );
  HB2( 4208, "Xtgt%U Sks [Proton]", 100, -200., 200., 100, -0.35, 0.35 );
  HB2( 4209, "Ytgt%V Sks [Proton]", 100, -100., 100., 100, -0.20, 0.20 );
  HB2( 4210, "Xtgt%Ytgt Sks [Proton]", 100, -200., 200., 100, -100., 100. );
  HB1( 4211, "P Sks [Proton]", 200, 0.0, 1.0 );
  HB1( 4212, "PathLength Sks [Proton]", 600, 3000., 6000. );
  HB1( 4213, "MassSqr", 600, -1.2, 1.2 );

  HB1( 5001, "Zvert [KPi]", 1000, -1000., 1000. ); 
  HB1( 5002, "MissingMass1 [KPi]", 1000, 0.0, 2.0 );

  HB2( 5011, "MissingMass%Us", 200, 0.0, 2.50, 100, -0.40, 0.40 );
  HB2( 5012, "MissingMass%Vs", 200, 0.0, 2.50, 100, -0.20, 0.20 );
  HB2( 5013, "MissingMass%Ub", 200, 0.0, 2.50, 100, -0.30, 0.30 );
  HB2( 5014, "MissingMass%Vb", 200, 0.0, 2.50, 100, -0.10, 0.10 );
  
  ////////////////////////////////////////////
  //Tree
  HBTree("tree","tree of Counter");
  TTree *tree = dynamic_cast<TTree *>(gFile->Get("tree"));

  //Trigger  
  tree->Branch("runnum",    &event.runnum,    "runnum/I");
  tree->Branch("evnum",     &event.evnum,     "evnum/I");
  tree->Branch("trigpat",    event.trigpat,   "trigpat[20]/I");
  tree->Branch("trigflag",   event.trigflag,   "trigflag[20]/I");

  //Misc
  tree->Branch("matrixpat",  event.matrixpat,   "matrixpat[10]/I");
  tree->Branch("matrixflag", event.matrixflag,  "matrixflag[10]/I");

  //Hodoscope
  tree->Branch("nhitsBh2",&event.nhitsBh2, "nhitsBh2/I");
  tree->Branch("nhBH2",   &event.nhBH2,   "nhBH2/I");
  tree->Branch("BH2Seg",   event.BH2Seg,  "BH2Seg[nhBH2]/D");
  tree->Branch("tBH2",     event.tBH2,    "tBH2[nhBH2]/D");
  tree->Branch("deBH2",    event.deBH2,   "deBH2[nhBH2]/D");

  tree->Branch("nhBH1",   &event.nhBH1,   "nhBH1/I");
  tree->Branch("BH1Seg",   event.BH1Seg,  "BH1Seg[nhBH1]/D");
  tree->Branch("tBH1",     event.tBH1,    "tBH1[nhBH1]/D");
  tree->Branch("deBH1",    event.deBH1,   "deBH1[nhBH1]/D");
  tree->Branch("btof",     event.btof,    "btof[nhBH1]/D");

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

  tree->Branch("Ac1Tdc",  event.Ac1Tdc, "Ac1Tdc[18]/I");


  tree->Branch("nhBAC", &event.nhBAC, "nhBAC/I");
  tree->Branch("BACSeg", event.BACSeg, "BACSeg[nhBAC]/I");
  tree->Branch("tBAC", event.tBAC, "tBAC[nhBAC]/D");
  tree->Branch("aBAC", event.aBAC, "aBAC[nhBAC]/D");


  //Beam DC 
  tree->Branch("ntBcIn",    &event.ntBcIn,     "ntBcIn/I");
  tree->Branch("nhBcIn",     event.nhBcIn,     "nhBcIn[ntBcIn]/I");
  tree->Branch("chisqrBcIn", event.chisqrBcIn, "chisqrBcIn[ntBcIn]/D");
  tree->Branch("x0BcIn",     event.x0BcIn,     "x0BcIn[ntBcIn]/D");
  tree->Branch("y0BcIn",     event.y0BcIn,     "y0BcIn[ntBcIn]/D");
  tree->Branch("u0BcIn",     event.u0BcIn,     "u0BcIn[ntBcIn]/D");
  tree->Branch("v0BcIn",     event.v0BcIn,     "v0BcIn[ntBcIn]/D");
  tree->Branch("xbh1BcIn",   event.xbh1BcIn,   "xbh1BcIn[ntBcIn]/D");
  tree->Branch("ybh1BcIn",   event.ybh1BcIn,   "ybh1BcIn[ntBcIn]/D");

  tree->Branch("ntBcOut",   &event.ntBcOut,     "ntBcOut/I");
  tree->Branch("nhBcOut",    event.nhBcOut,     "nhBcOut[ntBcOut]/I");
  tree->Branch("chisqrBcOut",event.chisqrBcOut, "chisqrBcOut[ntBcOut]/D");
  tree->Branch("x0BcOut",    event.x0BcOut,     "x0BcOut[ntBcOut]/D");
  tree->Branch("y0BcOut",    event.y0BcOut,     "y0BcOut[ntBcOut]/D");
  tree->Branch("u0BcOut",    event.u0BcOut,     "u0BcOut[ntBcOut]/D");
  tree->Branch("v0BcOut",    event.v0BcOut,     "v0BcOut[ntBcOut]/D");
  tree->Branch("xtgtBcOut",  event.xtgtBcOut,   "xtgtBcOut[ntBcOut]/D");
  tree->Branch("ytgtBcOut",  event.ytgtBcOut,   "ytgtBcOut[ntBcOut]/D");
  tree->Branch("xbh2BcOut",  event.xbh2BcOut,   "xbh2BcOut[ntBcOut]/D");
  tree->Branch("ybh2BcOut",  event.ybh2BcOut,   "ybh2BcOut[ntBcOut]/D");

  tree->Branch("ntK18",      &event.ntK18,     "ntK18/I");
  tree->Branch("nhK18",       event.nhK18,     "nhK18[ntK18]/I");
  tree->Branch("chisqrK18",   event.chisqrK18, "chisqrK18[ntK18]/D");
  tree->Branch("pK18",        event.pK18,      "pK18[ntK18]/D");
  tree->Branch("xtb",         event.xtb,       "xtb[ntK18]/D");
  tree->Branch("ytb",         event.ytb,       "ytb[ntK18]/D");
  tree->Branch("utb",         event.utb,       "utb[ntK18]/D");
  tree->Branch("vtb",         event.vtb,       "vtb[ntK18]/D");
  tree->Branch("thetab",      event.thetab,    "thetab[ntK18]/D");

  tree->Branch("ncBFT",      &event.ncBFT,     "ncBFT/I");
  tree->Branch("cxpos_BFT",   event.cxpos,     "cxpos_BFT[ncBFT]/D");
  tree->Branch("ctime_BFT", event.ctime_BFT, "ctime_BFT[ncBFT]/D");
  //SKS
  tree->Branch("ntSdcIn",    &event.ntSdcIn,     "ntSdcIn/I");
  tree->Branch("nhSdcIn",     event.nhSdcIn,     "nhSdcIn[ntSdcIn]/I");
  tree->Branch("chisqrSdcIn", event.chisqrSdcIn, "chisqrSdcIn[ntSdcIn]/D");
  tree->Branch("x0SdcIn",     event.x0SdcIn,     "x0SdcIn[ntSdcIn]/D");
  tree->Branch("y0SdcIn",     event.y0SdcIn,     "y0SdcIn[ntSdcIn]/D");
  tree->Branch("u0SdcIn",     event.u0SdcIn,     "u0SdcIn[ntSdcIn]/D");
  tree->Branch("v0SdcIn",     event.v0SdcIn,     "v0SdcIn[ntSdcIn]/D");

  tree->Branch("ntSdcOut",   &event.ntSdcOut,     "ntSdcOut/I");
  tree->Branch("nhSdcOut",    event.nhSdcOut,     "nhSdcOut[ntSdcOut]/I");
  tree->Branch("chisqrSdcOut",event.chisqrSdcOut, "chisqrOut[ntSdcOut]/D");
  tree->Branch("xtofSdcOut",  event.xtofSdcOut,   "xtofSdcOut[ntSdcOut]/D");
  tree->Branch("ytofSdcOut",  event.ytofSdcOut,   "ytofSdcOut[ntSdcOut]/D");
  tree->Branch("u0SdcOut",    event.u0SdcOut,     "u0SdcOut[ntSdcOut]/D");
  tree->Branch("v0SdcOut",    event.v0SdcOut,     "v0SdcOut[ntSdcOut]/D");

  tree->Branch("ntSks",      &event.ntSks,     "ntSks/I");
  tree->Branch("nhSks",       event.nhSks,     "nhSks[ntSks]/I");
  tree->Branch("chisqrSks",   event.chisqrSks, "chisqrSks[ntSks]/D");
  tree->Branch("path",        event.path,      "path[ntSks]/D");
  tree->Branch("pSks",        event.pSks,      "pSks[ntSks]/D");
  tree->Branch("m2",          event.m2,        "m2[ntSks]/D");
  tree->Branch("xts",         event.xts,       "xts[ntSks]/D");
  tree->Branch("yts",         event.yts,       "yts[ntSks]/D");
  tree->Branch("uts",         event.uts,       "uts[ntSks]/D");
  tree->Branch("vts",         event.vts,       "vts[ntSks]/D");
  tree->Branch("thetas",      event.thetas,    "thetas[ntSks]/D");

  //Reaction
  tree->Branch("nPi",    &event.nPi,     "nPi/I");
  tree->Branch("nK",     &event.nK,      "nK/I");
  tree->Branch("nKPi",   &event.nKPi,    "nKPi/I");
  tree->Branch("vtx",     event.vtx,     "vtx[nKPi]/D");
  tree->Branch("vty",     event.vty,     "vty[nKPi]/D");
  tree->Branch("vtz",     event.vtz,     "vtz[nKPi]/D");
  tree->Branch("closeDist",event.closeDist,"closeDist[nKPi]/D");
  tree->Branch("theta",   event.theta,    "theta[nKPi]/D");
  tree->Branch("MissMass1",event.MissMass1, "MissMass1[nKPi]/D");
  tree->Branch("MissMass2",event.MissMass2, "MissMass2[nKPi]/D");
  tree->Branch("MissMass3",event.MissMass3, "MissMass3[nKPi]/D");
  tree->Branch("MissMass4",event.MissMass4, "MissMass4[nKPi]/D");
  tree->Branch("MissMassPion",event.MissMassPion, "MissMassPion[nKPi]/D");
  tree->Branch("MissMassMuon",event.MissMassMuon, "MissMassMuon[nKPi]/D");
  tree->Branch("MissMassCorr",event.MissMassCorr, "MissMassCorr[nKPi]/D");
  tree->Branch("MissMassCorrDE",event.MissMassCorrDE, "MissMassCorrDE[nKPi]/D");
  tree->Branch("BE1",      event.BE1,       "BE1[nKPi]/D");
  tree->Branch("BE2",      event.BE2,       "BE2[nKPi]/D");
  tree->Branch("BE3",      event.BE3,       "BE3[nKPi]/D");
  tree->Branch("BE4",      event.BE4,       "BE4[nKPi]/D");
  tree->Branch("theta_CM", event.theta_CM, "theta_CM[nKPi]/D");
  tree->Branch("cost_CM",  event.cost_CM,  "cost_CM[nKPi]/D");

  tree->Branch("xk",         event.xk,       "xk[nKPi]/D");
  tree->Branch("yk",         event.yk,       "yk[nKPi]/D");
  tree->Branch("uk",         event.uk,       "uk[nKPi]/D");
  tree->Branch("vk",         event.vk,       "vk[nKPi]/D");
  tree->Branch("xpi",        event.xpi,      "xpi[nKPi]/D");
  tree->Branch("ypi",        event.ypi,      "ypi[nKPi]/D");
  tree->Branch("upi",        event.upi,      "upi[nKPi]/D");
  tree->Branch("vpi",        event.vpi,      "vpi[nKPi]/D");

  tree->Branch("pOrg",       event.pOrg,      "pOrg[nKPi]/D");
  tree->Branch("pCalc",      event.pCalc,     "pCalc[nKPi]/D");
  tree->Branch("pCorr",      event.pCorr,     "pCorr[nKPi]/D");
  tree->Branch("pCorrDE",    event.pCorrDE,   "pCorrDE[nKPi]/D");
  

  return true;
}

bool ConfMan::InitializeParameterFiles( void )
{
  using namespace h_Utility;
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

  if( FieldMapFileName_!="" )
    FieldMan::GetInstance().Initialize( FieldMapFileName_ );

  if( K18MatrixFileName_!="" )
    K18Matrix_ = new K18TransMatrix(K18MatrixFileName_);
  if(K18Matrix_) K18Matrix_->Initialize();

  if(FileManagerName_ != ""){
    if(-1 == FileManager::GetInstance().stInitialize(FileManagerName_)){
      status = false;
    }
  }

  return status;
}
