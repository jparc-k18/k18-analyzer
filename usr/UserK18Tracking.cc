/*
  UserK18Tracking.cc

  2012/1/24
*/

#include "VEvent.hh"

#include <iostream>
#include <sstream>
#include <cmath>

#include "UnpackerManager.hh"
#include "ConfMan.hh"
#include "DCGeomMan.hh"
#include "HistHelper.hh"
#include "DetectorID.hh"
#include "RawData.hh"
#include "HodoRawHit.hh"
#include "SksLib.hh"

#include "Hodo2Hit.hh"
#include "Hodo1Hit.hh"
#include "BH2Hit.hh"
#include "BH2Cluster.hh"
#include "FiberHit.hh"
#include "HodoCluster.hh"
#include "FiberCluster.hh"
#include "HodoAnalyzer.hh"
#include "DCAnalyzer.hh"
#include "K18Track_BFT.hh"
#include "HodoParamMan.hh"
#include "HodoPHCMan.hh"

#include "DCRawHit.hh"

#include "FileManager.hh"
#include "ParameterManager.hh"

const double Deg2Rad = acos(-1.)/180.;
const double Rad2Deg = 180./acos(-1.);

const double MaxMultiHitBcIn = 8.0;
const double MaxMultiHitBcOut = 40.0;

const double TdcLow  = -10.;
const double TdcHigh = 300.;

// BH2 Cut Parameters
const double MinDeltaEBH2 = 0.5; //for Pion
const double MaxDeltaEBH2 = 3.0; //for Pion

// BH1 Cut Parameters
const double MinDeltaEBH1 = 0.5; //for Pion
const double MaxDeltaEBH1 = 3.0; //for Pion

// BH1-BH2
const double MinBeamToF  = -2.0;
const double MaxBeamToF  =  2.0;

 VEvent::VEvent()
 {
 }

 VEvent::~VEvent()
 {
 }

class EventBFT 
  : public VEvent
{

private:
  RawData *rawData;
  DCAnalyzer *DCAna;
  HodoAnalyzer *hodoAna;

public:
  EventBFT();
  ~EventBFT();

  bool ProcessingBegin();
  bool ProcessingEnd();
  bool ProcessingNormal();
  bool InitializeHistograms();
  void InitializeEvent();
};

EventBFT::EventBFT()
  : VEvent(),
    rawData(0),
    DCAna(new DCAnalyzer()),
    hodoAna(new HodoAnalyzer)
{
}

EventBFT::~EventBFT()
{
  if (hodoAna) delete hodoAna;
  if (DCAna)   delete DCAna;
  if (rawData) delete rawData;

  hodoAna = NULL;
  DCAna   = NULL;
  rawData = NULL;
}

bool EventBFT::ProcessingBegin()
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
  int anaFlag;

  int matrixpat[NumOfMatrix];
  int matrixflag[NumOfMatrix];

  int ncOkBh1;
  int multip_bh1;
  int multip_bh2;
  int nhBH1;
  double btof[MaxHits];
  int MaxSizeCopper;
  int sizebft[NumOfPlaneBFT][NumOfSegBFT];
  int bft[NumOfPlaneBFT][NumOfSegBFT][8];

  int sizeTime[NumOfPlaneBFT];
  int HitSeg[NumOfPlaneBFT][NumOfSegBFT];
  double xpos[NumOfPlaneBFT][NumOfSegBFT];
  double time[NumOfPlaneBFT][NumOfSegBFT];

  int sizeCl;
  double ctime[NumOfSegBFT];
  double cxpos[NumOfSegBFT];

  // BcOut
  int ntOut;
  int nhOut[MaxHits];
  double chisqrOut[MaxHits];
  double x0Out[MaxHits];
  double y0Out[MaxHits];
  double u0Out[MaxHits];
  double v0Out[MaxHits];

  // K18
  int ntK18; 
  double p_2nd[MaxHits];
  double p_3rd[MaxHits];
  double delta_2nd[MaxHits];
  double delta_3rd[MaxHits];

  double xtgt[MaxHits];
  double ytgt[MaxHits];
  double utgt[MaxHits];
  double vtgt[MaxHits];

  double theta[MaxHits];
  double phi[MaxHits];
};
static Event event;

bool EventBFT::ProcessingNormal()
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

  //**************************************************************************
  //******************RawData

  //Misc
  
  bool fl_trigger = false;
  {
    const HodoRHitContainer &cont=rawData->GetMiscRawHC();
    //ParamMan& Trig = FileManager::GetInstance().GetParamMan("TriggerFlag");
    //Trig.stChangeParamSet("TriggerFlag");
    //static int TargetTrig = (int)Trig.GetParameter(0, "PIKLV2");
    int nh=cont.size();
    for( int i=0; i<nh; ++i ){
      HodoRawHit *hit=cont[i];
      int seg=hit->SegmentId()+1;
      int T=hit->GetTdc1();

      if( T>0 ){
	event.trigpat[i] = seg;
	//if(TargetTrig == seg){ fl_trigger = true;}
      }
      event.trigflag[seg-1] = T;
    }
  }

  //Matrix
  {
    const HodoRHitContainer &cont=rawData->GetMatrixRawHC();
    int nh=cont.size();
    for( int i=0; i<nh; ++i ){
      HodoRawHit *hit=cont[i];
      int seg=hit->SegmentId()+1;
      int T=hit->GetTdc1();

      if( T>0 ) event.matrixpat[i] = seg;
      event.matrixflag[seg-1] = T;
    }
  }

  HF1( 1, 0. );

  //if(!fl_trigger){return true;}
  
  HF1( 1, 1. );

 //////////////BH2 time 0
  
  hodoAna->DecodeBH2Hits(rawData);
  int ncBh2=hodoAna->GetNClustersBH2();
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
      //if( dEbh2<MinDeltaEBH2 || MaxDeltaEBH2<dEbh2 ) continue;
      ++ncOk;
      if( fabs(t)<fabs(mint) ){
	clBH2Time0 = cl;
	mint=t; time0=clBH2Time0->CTime0();
      }

      if(-1.<t && t < 1.){
	valid_bh2.push_back(cl);
      }
    }
    if(false
       || valid_bh2.size() > 1
       || clBH2Time0->ClusterSize() > 1
       ){
      event.multip_bh2 = 1;
    }

    double trig = clBH2Time0->CMeanTime();
    double min_t0 = -0.5;
    double max_t0 =  0.5;
    //if(min_t0 > trig || trig > max_t0 || ncOk<1){return true;}
  }

  HF1( 1, 3. );

  //////////////BH1 Analysis
  hodoAna->DecodeBH1Hits(rawData); 
  int ncBh1=hodoAna->GetNClustersBH1();
  event.nhBH1 = ncBh1;
  if(0 == ncBh1){ tree->Fill(); return true; };
  
  HF1( 1, 4. );

  HodoCluster* clBH1Time0 = hodoAna->GetClusterBH1(0);
  std::vector<HodoCluster*> valid_bh1;
  {
    int ncOk=0;
    double min_tof = clBH1Time0->CMeanTime() - time0;
    for( int i=0; i<ncBh1; ++i ){
      HodoCluster *cl=hodoAna->GetClusterBH1(i);
      double btof= (cl->CMeanTime())-time0;
      event.btof[i]=btof;
      double dEbh1 = cl->DeltaE();
      //------------------------Cut
      //if( dEbh1<MinDeltaEBH1 || MaxDeltaEBH1<dEbh1 ) continue;

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
	//cl->GoodForAnalysis( false );
      }
    }
    //------------------------Cut
    if( ncOk<1 ){
      // return true;
    }

    if(false
       || valid_bh1.size() > 1
       || clBH1Time0->ClusterSize() > 2
       ){
      event.multip_bh1 = 1;
    }
  }
  HF1( 1, 5. );
  //if(1 == event.multip_bh1 || 1 == event.multip_bh2){return true;}
  HF1( 1, 6. );
  
  // BFT (Copper)
  {
    for(int p = 0; p<NumOfPlaneBFT; ++p){
      const HodoRHitContainer &cont = rawData->GetBFTRawHC(p);
      int nhit = cont.size();
      HF1 (100*(p+1)+1, nhit);
      //      std::cout << "Nhits " << nhit << std::endl;
      for(int i = 0; i<nhit; ++i){
	HodoRawHit *hit = cont.at(i);
	int MHit = hit->GetSizeTdcUp();
	int seg = hit->SegmentId();

	if(MHit > event.MaxSizeCopper){
	  event.MaxSizeCopper = MHit;
	}
	event.sizebft[p][seg] = MHit;
	if(MHit > 8){
	  std::cout << "#W : Mhit is over 8" << std::endl;
	  MHit = 8;
	}

	bool fl_tdc_ok = false;
	for(int m = 0; m<MHit; ++m){
	  int buf = hit->GetTdcUp(m);
	  event.bft[p][seg][m] = buf;
	}
	
	HF1 (100*(p+1)+2, seg);
      }
    }
  }

  // BFT FiberHit Copper
  hodoAna->DecodeBFTHits(rawData);
  hodoAna->TimeCutBFT(-10., 10.);
  int MaxClSize = 0;
  std::vector<double> cont_x;
  {
    /*
    ParamMan& BH1Cut = FileManager::GetInstance().GetParamMan("CutBh12");
    BH1Cut.stChangeParamSet("PositionCutBFT");
    int clpos_bh1 = static_cast<int>(2*clBH1Time0->MeanSeg());
    double min    = BH1Cut.GetParameter(clpos_bh1, "min");
    double max    = BH1Cut.GetParameter(clpos_bh1, "max");
    */
    int nc = hodoAna->GetNClustersBFT();
    for(int i = 0; i<nc; ++i){
      FiberCluster* cl = hodoAna->GetClusterBFT(i);
      double x = cl->MeanPosition();

      //if(min < x && x < max){
	cont_x.push_back(x);
	event.cxpos[i] = x;
	//}
    }
    
    //    if(0 == cont_x.size()){ return true;}
    event.sizeCl = cont_x.size();
  }
  
  HF1( 1, 7.);

  // BcOut tracking
  HF1( 1, 10. );

  DCAna->DecodeRawHits( rawData );  

  //////////////BC3&4 number of hit in one layer not 0
  double multi_BcOut=0.;
  {
    for( int layer=1; layer<=NumOfLayersBcOut; ++layer ){
      const DCHitContainer &contOut =DCAna->GetBcOutHC(layer);
      int nhOut=contOut.size();
      multi_BcOut += double(nhOut);
    }
  }
  if( multi_BcOut/double(NumOfLayersBcOut) > MaxMultiHitBcOut )
    //    return true;

  HF1( 1, 11. );


  //////////////BCOut tracking
#if 1
  //DCAna->ChiSqrCutBcOut(10);
  DCAna->TrackSearchBcOut();
  int ntOut=DCAna->GetNtracksBcOut();
  if( ntOut>MaxHits ) ntOut=MaxHits;
  event.ntOut=ntOut;
  int ntOutOk=0;
  HF1( 30, double(ntOut) );
  for( int it=0; it<ntOut; ++it ){
    DCLocalTrack *tp=DCAna->GetTrackBcOut(it);
    int nh=tp->GetNHit();
    double chisqr=tp->GetChiSquare();
    double u0=tp->GetU0(), v0=tp->GetV0();
    double x0=tp->GetX(0.), y0=tp->GetY(0.);

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

    if( chisqr < 20 ) ntOutOk++;

    for( int ih=0; ih<nh; ++ih ){
      DCLTrackHit *hit=tp->GetHit(ih);
      int layerId=hit->GetLayer()-100;
      HF1( 33, layerId );
    }
  }
  //  if( ntOutOk<1 ) return true;
#endif
  HF1( 1, 12. );

  // K18 D2U

  HF1( 1, 20. );

#if 1
  DCAna->TrackSearchK18_BFT(cont_x);

  int ntK18=DCAna->GetNTracksK18_BFT();
  HF1( 50, double(ntK18) );
  if(ntK18 > MaxHits){
    std::cout << "#W : "
	      << " Too many K18 tracks : ntK18 = "
	      << ntK18
	      << std::endl;
    ntK18 = MaxHits;
  }
  event.ntK18 = ntK18;
  for( int i=0; i<ntK18; ++i ){
    K18Track_BFT *tp=DCAna->GetK18Track_BFT(i);
    if(!tp) continue;
    double xt=tp->Xtgt(), yt=tp->Ytgt();
    double ut=tp->Utgt(), vt=tp->Vtgt();
    double p_2nd=tp->PD2U();
    double p_3rd=tp->PD2U_3rd();
    double delta_2nd=tp->DeltaD2U();
    double delta_3rd=tp->DeltaD2U_3rd();
    double cost = 1./sqrt(1.+ut*ut+vt*vt);
    double theta = acos(cost)*Rad2Deg;
    double phi   = atan2( ut, vt );

    HF1( 54, xt ); HF1( 55, yt ); HF1( 56, ut ); HF1( 57,vt );
    HF2( 58, xt, ut ); HF2( 59, yt, vt ); HF2( 60, xt, yt );
    HF1( 61, p_3rd ); HF1( 62, delta_3rd );


    event.p_2nd[i] = p_2nd;
    event.p_3rd[i] = p_3rd;
    event.delta_2nd[i] = delta_2nd;
    event.delta_3rd[i] = delta_3rd;

    event.xtgt[i] = xt;
    event.ytgt[i] = yt;
    event.utgt[i] = ut;
    event.vtgt[i] = vt; 

    event.theta[i] = theta;
    event.phi[i] = phi;

  }

#endif

  HF1( 1, 22. );

  if(0 != ntK18){event.anaFlag = 1;};

  tree->Fill();

  return true;
}

void EventBFT::InitializeEvent( void )
{
  event.runnum = -1;
  // size branches must be 0 (not -1): TTree::Fill crashes on
  // negative array sizes when the fill section is skipped
  event.evnum = -1;
  event.MaxSizeCopper  = -1;
  event.sizeCl = 0;
  event.ncOkBh1= -1;
  event.ntK18 = 0;
  event.ntOut = 0;
  event.nhBH1 = 0;
  event.multip_bh1 = 0;
  event.multip_bh2 = 0;
  event.anaFlag  = 0;

  for(int ip = 0; ip<NumOfPlaneBFT; ++ip){
    event.sizeTime[ip] = 0;
    for( int it=0; it<NumOfSegBFT; it++){
      event.sizebft[ip][it] = -1;
      event.HitSeg[ip][it] = -1;
      event.time[ip][it] = -999;
      event.ctime[it] = -999;
      event.xpos[ip][it] = -999;
      event.cxpos[it] = -999;
      for(int m = 0; m<8; ++m){
	event.bft[ip][it][m] = -999;
      }
    }
  }

  for( int it=0; it<NumOfMisc; it++){
    event.trigpat[it] = -1;
    event.trigflag[it] = -1;
  }

  for( int it=0; it<NumOfMatrix; it++){
    event.matrixpat[it] = -1;
    event.matrixflag[it] = -1;
  }

  for(int i = 0; i<MaxHits; ++i){
    event.nhOut[i] = -1;
    event.chisqrOut[i] = -999.;
    event.x0Out[i] = -999.;
    event.y0Out[i] = -999.;
    event.u0Out[i] = -999.;
    event.v0Out[i] = -999.;

    event.p_2nd[i] = -999.;
    event.p_3rd[i] = -999.;
    event.delta_2nd[i] = -999.;
    event.delta_3rd[i] = -999.;

    event.xtgt[i] = -999.;
    event.ytgt[i] = -999.;
    event.utgt[i] = -999.;
    event.vtgt[i] = -999.;

    event.theta[i] = -999.;
    event.phi[i] = -999.;
  }
}

bool EventBFT::ProcessingEnd()
{
  return true;
}

VEvent *ConfMan::EventAllocator()
{
  return new EventBFT;
}

// const int NbinAdc = 1024;
// const double MinAdc  =    0.;
// const double MaxAdc  = 4096.;

// const int NbinTdc = 1024;
// const double MinTdc  =    0.;
// const double MaxTdc  = 4096.;

bool ConfMan:: InitializeHistograms()
{
  HB1( 1, "Status", 30, 0., 30. );

  //BFT
  HB1( 100+1, "N hits U", 160, 0, 160 );
  HB1( 100+2, "Hit pat U", 160, 0, 160 );
  HB1( 200+1, "N hits D", 160, 0, 160 );
  HB1( 200+2, "Hit pat D", 160, 0, 160 );

  // FiberHit
  HB1( 300+1, "N hits U", 160, 0, 160 );
  HB1( 300+2, "Hit pat U", 160, 0, 160 );
  HB1( 400+1, "N hits D", 160, 0, 160 );
  HB1( 400+2, "Hit pat D", 160, 0, 160 );
  HB1(1000+2, "X pos U",  180, -90, 90);
  HB1(2000+2, "X pos D",  180, -90, 90);
  HB1(1000+15, "X res U",  200, -10, 10);
  HB1(2000+15, "X res D",  200, -10, 10);

  // Cluster
  HB1( 10,    "ClusterSize", 5,0,5);
  HB1( 14,    "X pos Cl", 180, -90, 90);
  HB1( 15,    "Residual Cl", 200,-10,10);

  // BcIn
  HB2(20, "XY",  250, -125, 125, 200, -100, 100);

  // BcOut
  HB1( 30, "#Tracks BcOut", 10, 0., 10. );
  HB1( 31, "#Hits of Track BcOut", 20, 0., 20. );
  HB1( 32, "Chisqr BcOut", 500, 0., 50. ); 
  HB1( 33, "LayerId BcOut", 15, 12., 27. );
  HB1( 34, "X0 BcOut", 400, -100., 100. ); 
  HB1( 35, "Y0 BcOut", 400, -100., 100. );
  HB1( 36, "U0 BcOut",  200, -0.20, 0.20 );
  HB1( 37, "V0 BcOut",  200, -0.20, 0.20 );
  HB2( 38, "U0%X0 BcOut", 100, -100., 100., 100, -0.20, 0.20 );
  HB2( 39, "V0%Y0 BcOut", 100, -100., 100., 100, -0.20, 0.20 );
  HB2( 40, "X0%Y0 BcOut", 100, -100., 100., 100, -100., 100. );

  // K18
  HB1( 50, "#Tracks K18", 20, 0., 20. );
  HB1( 51, "#Hits of K18Track", 30, 0., 30. );
  HB1( 52, "Chisqr K18Track", 500, 0., 100. );
  HB1( 53, "LayerId K18Track", 50, 0., 50. );
  HB1( 54, "Xtgt K18Track", 200, -100., 100. );
  HB1( 55, "Ytgt K18Track", 200, -100., 100. );
  HB1( 56, "Utgt K18Track", 300, -0.30, 0.30 );
  HB1( 57, "Vtgt K18Track", 300, -0.20, 0.20 );
  HB2( 58, "U%Xtgt K18Track", 100, -100., 100., 100, -0.25, 0.25 );
  HB2( 59, "V%Ytgt K18Track", 100, -100., 100., 100, -0.10, 0.10 );
  HB2( 60, "Y%Xtgt K18Track", 100, -100., 100., 100, -100., 100. );
  HB1( 61, "P K18Track", 600, 1.50, 2.1 );
  HB1( 62, "dP K18Track", 400, -0.04, 0.04 );

  //Tree
  HBTree("tree","tree of Counter");
  TTree *tree = dynamic_cast<TTree *>(gFile->Get("tree"));

  //Misc
  tree->Branch("runnum",    &event.runnum,    "runnum/I");
  tree->Branch("evnum",     &event.evnum,     "evnum/I");
  tree->Branch("trigpat",    event.trigpat,   "trigpat[20]/I");
  tree->Branch("trigflag",   event.trigflag,  "trigflag[20]/I");
  tree->Branch("anaFlag",   &event.anaFlag,   "anaFlag/I");

  //Misc
  tree->Branch("matrixpat",  event.matrixpat,   "matrixpat[10]/I");
  tree->Branch("matrixflag", event.matrixflag,  "matrixflag[10]/I");

  //BH1
  tree->Branch("ncOkBh1",   &event.ncOkBh1,     "ncOkBh1/I");
  tree->Branch("nhBH1", &event.nhBH1, "nhBH1/I");
  tree->Branch("btof", event.btof, "btof[nhBH1]/D");

  //BFT (Copper)
  tree->Branch("MSC",     &event.MaxSizeCopper,  "MSC/I");
  tree->Branch("sizeBft_u",  event.sizebft[0],   "sizebft_u[160]/I");
  tree->Branch("sizeBft_d",  event.sizebft[1],   "sizebft_d[160]/I");
  tree->Branch("bft_u",      event.bft[0],       "bft_u[160][8]/I");
  tree->Branch("bft_d",      event.bft[1],       "bft_d[160][8]/I");  

  tree->Branch("sizeTime_u",&event.sizeTime[0],  "sizeTime_u/I");
  tree->Branch("sizeTime_d",&event.sizeTime[1],  "sizeTime_d/I");
  tree->Branch("seg_u",      event.HitSeg[0],    "seg_u[sizeTime_u]/I");
  tree->Branch("seg_d",      event.HitSeg[1],    "seg_d[sizeTime_d]/I");
  tree->Branch("time_u",     event.time[0],      "time_u[sizeTime_u]/D");
  tree->Branch("time_d",     event.time[1],      "time_d[sizeTime_d]/D");

  tree->Branch("sizeCl",    &event.sizeCl,       "sizeCl/I");
  tree->Branch("ctime",      event.ctime,        "ctime[sizeCl]/D");
  tree->Branch("cxpos",      event.cxpos,        "cxpos[sizeCl]/D");

  // BcOut
  tree->Branch("ntOut",   &event.ntOut,     "ntOut/I");
  tree->Branch("nhOut",    event.nhOut,     "nhOut[ntOut]/I");
  tree->Branch("chisqrOut",event.chisqrOut, "chisqrOut[ntOut]/D");
  tree->Branch("x0Out",    event.x0Out,     "x0Out[ntOut]/D");
  tree->Branch("y0Out",    event.y0Out,     "y0Out[ntOut]/D");
  tree->Branch("u0Out",    event.u0Out,     "u0Out[ntOut]/D");
  tree->Branch("v0Out",    event.v0Out,     "v0Out[ntOut]/D");

  // K18
  tree->Branch("ntK18",      &event.ntK18,     "ntK18/I");
  tree->Branch("p_2nd",       event.p_2nd,     "p_2nd[ntK18]/D");
  tree->Branch("p_3rd",       event.p_3rd,     "p_3rd[ntK18]/D");
  tree->Branch("delta_2nd",   event.delta_2nd, "delta_2nd[ntK18]/D");
  tree->Branch("delta_3rd",   event.delta_3rd, "delta_3rd[ntK18]/D");

  tree->Branch("xtgt",    event.xtgt,   "xtgt[ntK18]/D");
  tree->Branch("ytgt",    event.ytgt,   "ytgt[ntK18]/D");
  tree->Branch("utgt",    event.utgt,   "utgt[ntK18]/D");
  tree->Branch("vtgt",    event.vtgt,   "vtgt[ntK18]/D");

  tree->Branch("theta",   event.theta,  "theta[ntK18]/D");
  tree->Branch("phi",     event.phi,    "phi[ntK18]/D");

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

  if( DCTdcCalibFileName_!="" )
    DCTdcCalibManager_ = new DCTdcCalibMan(DCTdcCalibFileName_);
  if(DCTdcCalibManager_) DCTdcCalibManager_->Initialize();

  if( DCDriftParamFileName_!="" )
    DCDriftParamManager_ = new DCDriftParamMan(DCDriftParamFileName_);
  if(DCDriftParamManager_) DCDriftParamManager_->Initialize();

  DCGeomManager_ = & DCGeomMan::GetInstance();
  if( DCGeomFileName_!="" )
    DCGeomManager_->Initialize(DCGeomFileName_);
  else
    DCGeomManager_->Initialize();

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
