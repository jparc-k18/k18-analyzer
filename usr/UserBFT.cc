/*
  UserHodoscope.cc

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
#include "HodoParamMan.hh"
#include "HodoPHCMan.hh"

#include "DCRawHit.hh"

#include "FileManager.hh"
#include "ParameterManager.hh"

const double Deg2Rad = acos(-1.)/180.;
const double Rad2Deg = 180./acos(-1.);

const double MaxMultiHitBcIn = 8.0;

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
#define MaxHits 512
#endif

//For Tree
struct Event{
  int runnum;
  int evnum;
  int trigpat[NumOfMisc];
  int trigflag[NumOfMisc];

  int matrixpat[NumOfMatrix];
  int matrixflag[NumOfMatrix];

  int ncOkBh1;
  int multip_bh1;
  int multip_bh2;

  int MaxSizeCopper;
  int sizebft[NumOfPlaneBFT][NumOfSegBFT];
  int bft[NumOfPlaneBFT][NumOfSegBFT][8];

  int sizeTime[NumOfPlaneBFT];
  int HitSeg[NumOfPlaneBFT][NumOfSegBFT];
  double xpos[NumOfPlaneBFT][NumOfSegBFT];
  double time[NumOfPlaneBFT][NumOfSegBFT];

  int sizeCl;
  double ctime[MaxHits];
  double cxpos[MaxHits];
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
    // ParamMan& Trig = FileManager::GetInstance().GetParamMan("TriggerFlag");
    //Trig.stChangeParamSet("TriggerFlag");
    //static int TargetTrig = (int)Trig.GetParameter(0, "TARGET");
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

  HF1(1, 0);
  //if(!fl_trigger){return true;}
  HF1(1, 1);

 //////////////BH2 time 0
  hodoAna->DecodeBH2Hits(rawData);
  int ncBh2=hodoAna->GetNClustersBH2();
  HF1(20003, ncBh2);
  if( ncBh2==0 ){ tree->Fill(); return true; }

  HF1(1, 2);

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
      //      if( dEbh2<MinDeltaEBH2 || MaxDeltaEBH2<dEbh2 ) continue;
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
    double min_t0 = -1.0;
    double max_t0 =  1.0;
    if(min_t0 > trig || trig > max_t0 || ncOk<1){
      HF1(20002, trig);

      if(min_t0 > trig || trig > max_t0){
	if(ncOk<1){
	  HF1(2, 3);
	}else{
	  HF1(2, 1);
	}
      }else if(ncOk<1){
      	HF1(2, 2);
      }
      
      { tree->Fill(); return true; }
    }else{
      HF1(20000, ncOk);
      HF1(20001, trig);
    }
  }

  HF1(1, 3);

 //////////////BH1 Analysis
  hodoAna->DecodeBH1Hits(rawData); 
  int ncBh1=hodoAna->GetNClustersBH1();
  HF1(10003, ncBh1);
  if(0 == ncBh1){ tree->Fill(); return true; };

  HF1(1, 4);

  HodoCluster* clBH1Time0 = hodoAna->GetClusterBH1(0);
  std::vector<HodoCluster*> valid_bh1;
  {
    int ncOk=0;
    double min_tof = clBH1Time0->CMeanTime() - time0;
    for( int i=0; i<ncBh1; ++i ){
      HodoCluster *cl=hodoAna->GetClusterBH1(i);
      double btof= cl->CMeanTime()-time0;
      double dEbh1 = cl->DeltaE();
      //------------------------Cut
      if( dEbh1<MinDeltaEBH1 || MaxDeltaEBH1<dEbh1 ) continue;

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
    
    event.ncOkBh1 = ncOk;
    HF1(10000, ncOk);
    //------------------------Cut
    if( ncOk<1 ){
      { tree->Fill(); return true; }
    }

    if(false
       || valid_bh1.size() > 1
       || clBH1Time0->ClusterSize() > 2
       ){
      event.multip_bh1 = 1;
    }
  }

  HF1(1, 5);

  if(1 == event.multip_bh1 || 1 == event.multip_bh2){ tree->Fill(); return true; }

  HF1(1, 6);

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
	  HF1( 10006+p, double(buf) );
	  HF1( 11000+1000*p+seg+1, double(buf) );
	}
	
	HF1 (100*(p+1)+2, seg);
      }
    }
  }

  // BFT FiberHit Copper
  hodoAna->DecodeBFTHits(rawData);
  int NofP = 0;

  std::vector<FiberCluster*> valid_bft;
  {
    int nh = hodoAna->GetNHitsBFT();
    int nh_u = 0, nh_d = 0;
    for(int i = 0; i<nh; ++i){
      FiberHit* hit = hodoAna->GetHitBFT(i);
      int pairid = hit->PairId();
      
      int p = 0, seg = 0;
      if(0 == pairid%2){
	// U layer
	p   = 0;
	seg = pairid/2;
	++nh_u;
      }else{
	// D layer
	p   = 1;
	seg = (pairid-1)/2;
	++nh_d;
      }
      
      HF1(100*(p+3) +2, seg);
      event.HitSeg[p][i];
      double time = hit->GetCTime();
      event.time[p][i] = time;
      {
	// PHC calibration: tot vs time relative to the BH2 time zero
	double w = hit->GetWidth();
	double t = hit->GetTime() - time0;
	double ct = time - time0;
	HF1( 10021+p, ct );
	HF2( 10023+p, w, t );
	HF2( 10025+p, w, ct );
	HF2( 15001+1000*p+seg, w, t );
	HF2( 17001+1000*p+seg, w, ct );
      }
      double x = hit->GetPosition();
      HF1(1000*(p+1)+2, x);
      event.xpos[p][i] = x;
    }
    
    HF1(100*(0+3) +1, nh_u);
    HF1(100*(1+3) +1, nh_d);
    event.sizeTime[0] = nh_u;
    event.sizeTime[1] = nh_d;

    // Copper Cluster
    int nc = hodoAna->GetNClustersBFT();
    if(nc >= MaxHits){
      std::cout << "nc is too big : " << nc<<std::endl;
      nc = MaxHits;
      //      getchar();
    }
    NofP = nc;
    event.sizeCl = nc;
    HF1(10, nc);
    for(int i = 0; i<nc; ++i){
      FiberCluster* cl = hodoAna->GetClusterBFT(i);
      double ctime = cl->CMeanTime();
      HF1(10030, ctime - time0);
      if(-9 < ctime && ctime <6){
	event.ctime[i] = ctime;
	int clsize = cl->ClusterSize();
	HF1(11, clsize);
	double x = cl->MeanPosition();
	HF1(14, x);
	HF2(16, clBH1Time0->MeanSeg(), x);
	event.cxpos[i] = x;
	
	valid_bft.push_back(cl);
      }
    }

    HF1(20, (double)valid_bft.size());

    // correlation BH1
    // the CutBh12 parameter file (FILEMAN) does not exist in E05;
    // skip the BH1 position matching when it is not initialized
    if(0 == event.multip_bh1 && ConfMan::GetConfManager() &&
       !ConfMan::GetConfManager()->GetFileManagerName().empty()){
      int nc_wBh1 = 0;

      ParamMan& BH1Cut = FileManager::GetInstance().GetParamMan("CutBh12");
      BH1Cut.stChangeParamSet("PositionCutBFT");
      for(int i = 0; i<valid_bft.size(); ++i){
	int clpos_bh1 = static_cast<int>(2*clBH1Time0->MeanSeg());
	double min = BH1Cut.GetParameter(clpos_bh1, "min");
	double max = BH1Cut.GetParameter(clpos_bh1, "max");

	double x = valid_bft.at(i)->MeanPosition();
	if(min < x && x < max){
	  ++nc_wBh1;
	  HF2(17, clBH1Time0->MeanSeg(), x);
	}
      }

      HF1(30, nc_wBh1);
    }
    
  }
  
  //////////////BC1&2 number of hit layer
  DCAna->DecodeRawHits( rawData );  

  //BC1&BC2
  double multi_BcIn=0.;
  {
    for( int layer=1; layer<=NumOfLayersBcIn; ++layer ){
      const DCHitContainer &contIn =DCAna->GetBcInHC(layer);
      const DCRHitContainer &raw_contIn =rawData->GetBcInRawHC(layer);
      int nhIn=contIn.size();
      multi_BcIn += double(nhIn);

      int plane_eff = (layer-1)*3;
      bool fl_valid_sig = false; 
      int Current = 0;
      for( int i=0; i<nhIn; ++i ){
	DCHit *hit     = contIn[i];
	DCRawHit *raw_hit = raw_contIn[i];
	double wire=hit->GetWire();
	int csize=hit->GetClusterSize();

	int nhtdc = hit->GetTdcSize();
	int tdc1st = -1;
	for( int k=0; k<nhtdc; k++ ){
	  int tdc = raw_hit->GetTdc(k);
	  if( (TdcLow< tdc && tdc < TdcHigh) && tdc > tdc1st ){
	    tdc1st=tdc;
	    fl_valid_sig = true;
	  }
	}

	int nhdt = hit->GetDriftTimeSize();
	for( int k=0; k<nhdt; k++ ){
	  double dt = hit->GetDriftTime(k);
	}
	int nhdl = hit->GetDriftTimeSize();
	for( int k=0; k<nhdl; k++ ){
	  double dl = hit->GetDriftLength(k);
	}
      }

      if(fl_valid_sig){
	++plane_eff;
      }
    }    
  }
  if( multi_BcIn/double(NumOfLayersBcIn) > MaxMultiHitBcIn )
    { tree->Fill(); return true; }

#if 0
  DCAna->TrackSearchBcIn();
  
  int nt=DCAna->GetNtracksBcIn();
  for( int it=0; it<nt; ++it ){
    DCLocalTrack *tp=DCAna->GetTrackBcIn(it);
    int nh=tp->GetNHit();
    double chisqr=tp->GetChiSquare();
    double x0=tp->GetX0(), y0=tp->GetY0();
    double u0=tp->GetU0(), v0=tp->GetV0();

    double cost = 1./sqrt(1.+u0*u0+v0*v0);
    double theta = acos(cost)*Rad2Deg;
    
    double xtgt=tp->GetX(-555.), ytgt=tp->GetY(555.);
    double utgt=u0, vtgt=v0;
    double xbac=tp->GetX(885.), ybac=tp->GetY(885.);
    double ubac=u0, vbac=v0;

    double xbc2x1  = tp->GetX(-20.);
    double xbc1bc2 = tp->GetX(-235.);
    double x_bh1   = tp->GetX(-527.25);
    double x_bft   = tp->GetX(-199.0);
    double y_bft   = tp->GetY(-199.0);
    
    HF2(20, x_bft, y_bft);
    
    //    if(NofP == 1 && nt == 1 && event.multip_bh1 == 0 && event.multip_bh2 == 0){
    {
      for(int ic = 0; ic<NofP; ++ic){
	double diff = x_bft-event.cxpos[ic];
	if(fabs(diff)< 3.){
	  HF1(15, diff);
	  HF2(25, event.cxpos[ic], diff);
	}
      }
    }
  }

#endif

  tree->Fill();
  return true;
}

void EventBFT::InitializeEvent( void )
{
  event.runnum = -1;
  event.evnum = -1;
  event.MaxSizeCopper  = -1;
  event.sizeCl = 0;
  event.ncOkBh1 = -1;
  event.multip_bh1 = 0;
  event.multip_bh2 = 0;
  
  for(int ip = 0; ip<NumOfPlaneBFT; ++ip){
    event.sizeTime[ip] = 0;
    for( int it=0; it<NumOfSegBFT; it++){
      event.sizebft[ip][it] = -1;
      event.HitSeg[ip][it] = -1;
      event.time[ip][it] = -999;
      event.xpos[ip][it] = -999;
      for(int m = 0; m<8; ++m){
	event.bft[ip][it][m] = -999;
      }
    }
  }

  for(int i = 0; i<MaxHits; ++i){
    event.cxpos[i] = -999;
    event.ctime[i] = -999;
  }

  for( int it=0; it<NumOfMisc; it++){
    event.trigpat[it] = -1;
    event.trigflag[it] = -1;
  }

  for( int it=0; it<NumOfMatrix; it++){
    event.matrixpat[it] = -1;
    event.matrixflag[it] = -1;
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
  HB1( 1, "Status", 20, 0., 20. );
  HB1( 2, "analys flag", 20, 0., 20. );

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

  // calibration histograms
  HB1( 10006, "BFT Tdc U all", 1000, 0., 1000. );
  HB1( 10007, "BFT Tdc D all", 1000, 0., 1000. );
  HB1( 10021, "BFT CTime-Time0 U all", 400, -20., 20. );
  HB1( 10022, "BFT CTime-Time0 D all", 400, -20., 20. );
  HB1( 10030, "BFT cluster CTime-Time0", 800, -20., 20. );
  HB2( 10023, "BFT Time-Time0%Tot U all",  100, 0., 100., 200, -20., 20. );
  HB2( 10024, "BFT Time-Time0%Tot D all",  100, 0., 100., 200, -20., 20. );
  HB2( 10025, "BFT CTime-Time0%Tot U all", 100, 0., 100., 200, -20., 20. );
  HB2( 10026, "BFT CTime-Time0%Tot D all", 100, 0., 100., 200, -20., 20. );
  for( int p=0; p<NumOfPlaneBFT; ++p ){
    for( int i=0; i<NumOfSegBFT; ++i ){
      std::stringstream t1, t2, t3;
      const char* ud = (p==0) ? "U" : "D";
      t1 << "BFT Tdc " << ud << "-" << i;
      t2 << "BFT Time-Time0%Tot " << ud << "-" << i;
      t3 << "BFT CTime-Time0%Tot " << ud << "-" << i;
      HB1( 11000+1000*p+i+1, t1.str().c_str(), 1000, 0., 1000. );
      HB2( 15001+1000*p+i, t2.str().c_str(), 50, 0., 100., 100, -20., 20. );
      HB2( 17001+1000*p+i, t3.str().c_str(), 50, 0., 100., 100, -20., 20. );
    }
  }
  HB1(1000+15, "X res U",  200, -10, 10);
  HB1(2000+15, "X res D",  200, -10, 10);

  // Cluster
  HB1( 10,    "Ncluster raw", 10,0,10);
  HB1( 11,    "ClusterSize", 5,0,5);
  HB1( 14,    "X pos Cl", 180, -90, 90);
  HB1( 15,    "Residual Cl", 200,-10,10);
  HB2( 16,    "Time0 cor", 24, 0, 12, 160, -80, 80);
  HB2( 17,    "Time0 cor w/ Bh1", 24, 0, 12, 160, -80, 80);
  HB1( 20,    "Ncluster valid", 10,0,10);
  HB2( 25,    "resi % pos", 320, -160, 160, 100, -5, 5);
  HB1( 30,    "Ncluster valid w/ Bh1", 10,0,10);

  // BH1
  HB1( 10000, "nhvalid", 10,0,10);
  
  // BH2
  HB1( 20000, "nhvalid", 10,0,10);
  HB1( 20001, "time0 true", 3000,-100,100);
  HB1( 20002, "time0 false", 3000,-100,100);

  //Tree
  HBTree("tree","tree of Counter");
  TTree *tree = dynamic_cast<TTree *>(gFile->Get("tree"));

  //Misc
  tree->Branch("runnum",    &event.runnum,    "runnum/I");
  tree->Branch("evnum",     &event.evnum,     "evnum/I");
  tree->Branch("trigpat",    event.trigpat,   "trigpat[20]/I");
  tree->Branch("trigflag",   event.trigflag,  "trigflag[20]/I");

  //Misc
  tree->Branch("matrixpat",  event.matrixpat,   "matrixpat[10]/I");
  tree->Branch("matrixflag", event.matrixflag,  "matrixflag[10]/I");

  //BH1
  tree->Branch("ncOkBh1",   &event.ncOkBh1,     "ncOkBh1/I");
  tree->Branch("mp_bh1",  &event.multip_bh1, "mp_bh1/I");
  tree->Branch("mp_bh2",  &event.multip_bh2, "mp_bh2/I");

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

  if(FileManagerName_ != ""){
    if(-1 == FileManager::GetInstance().stInitialize(FileManagerName_)){
      status = false;
    }
  }

  return status;
}
