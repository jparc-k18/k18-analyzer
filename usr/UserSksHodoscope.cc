/*
  UserSksHodoscope.cc

  2012/1/24
*/

#include "VEvent.hh"

#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <cmath>

#include "UnpackerManager.hh"
#include "ConfMan.hh"
#include "HistHelper.hh"
#include "DetectorID.hh"
#include "RawData.hh"
#include "DCRawHit.hh"
#include "SksLib.hh"

#include "FileManager.hh"
#include "ParameterManager.hh"

const double Deg2Rad = acos(-1.)/180.;
const double Rad2Deg = 180./acos(-1.);

static const int BH2Hid = 20000;
static const int TOFHid = 50000;
static const int ACHid  = 60000;
static const int LCHid  = 70000;

static const double LightVel = 299.792458;  /* mm/nsec */
static const double PionMass   = 0.139570;
static const double KaonMass   = 0.493677;
static const double ProtonMass = 0.9382720;
// TrackChiSqr
const double MaxChisqrSks    = 100.;
const double MaxChisqrBdcIn  = 100.;
const double MaxChisqrBdcOut = 100.;
const double MaxChisqrSdcIn  = 100.;
const double MaxChisqrSdcOut = 100.;

//Time difference between BH2 and Target [ns]  
double offset = 1.3108917166; // for kaon @ 1.8 GeV/c
// TOF Cut Parameters
double MinDeltaETof = 0.; //for Kaon
double MaxDeltaETof = 0.; //for Kaon
double MinTimeTof =  0.;
double MaxTimeTof =  0.;

// LC Cut Parameters
double MinDeltaELc = 0.; //for Kaon
double MaxDeltaELc = 0.; //for Kaon
double MinTimeLc =  0.;
double MaxTimeLc =  0.;

// TOF-LC  
double MinToFT2L =  0.;
double MaxToFT2L =  0.;
int MinToF2LSeg  =  0;
int MaxToF2LSeg  =  0;

// DC multiplisity cut
double MaxMultiHitSdcIn  = 30.;
double MaxMultiHitSdcOut = 30.;

//SdcIn vs TOF&LC cut
double ModFactorTof  = 0.;
double ModFactorLc   = 0.;
double MinModTimeTof = 0.;
double MaxModTimeTof = 0.;
double MinModTimeLc  = 0.;
double MaxModTimeLc  = 0.;

// TOF-SdcOut
double MinXDifTof  = 0.;
double MaxXDifTof  = 0.;
double MinYDifTof  = 0.;
double MaxYDifTof  = 0.;

// LC-SdcOut
double MinXDifLc  = 0.;
double MaxXDifLc  = 0.;
double MinYDifLc  = 0.;
double MaxYDifLc  = 0.;

VEvent::VEvent()
{
}

VEvent::~VEvent()
{
}

class EventSksHodoscope : public VEvent
{

private:
  RawData *rawData;
  DCAnalyzer *DCAna;
  HodoAnalyzer *hodoAna;

public:
  EventSksHodoscope();
  ~EventSksHodoscope();

  bool ProcessingBegin();
  bool ProcessingEnd();
  bool ProcessingNormal();
  bool InitializeHistograms();
  void InitializeEvent();
};

EventSksHodoscope::EventSksHodoscope()
  : VEvent(), 
    rawData(0), 
    DCAna(new DCAnalyzer()), 
    hodoAna(new HodoAnalyzer())
{
}

EventSksHodoscope::~EventSksHodoscope()
{
  if (DCAna)   delete DCAna;
  if (hodoAna) delete hodoAna;
  if (rawData) delete rawData;
}

#ifndef MaxHits 
#define MaxHits 10
#endif

//For Tree
struct Event{
  int runnum;
  int trigpat[NumOfMisc];
  int trigflag[NumOfMisc];

  int matrixpat[NumOfMatrix];
  int matrixflag[NumOfMatrix];

  int ntIn;
  int nhIn[MaxHits]; 
  double chisqrIn[MaxHits];
  double u0In[MaxHits];
  double v0In[MaxHits];
  double x0In[MaxHits];
  double y0In[MaxHits];

  int nhx1;
  int x1wire[MaxHits];
  int nhx2;
  int x2wire[MaxHits];

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

  //Hodoscope
  int nhTof;
  int TofSeg[MaxHits];
  double tTof[MaxHits];
  double deTof[MaxHits];
  int evnum;
  double time0;
  int tofua[NumOfSegTOF];
  int tofda[NumOfSegTOF];
  int tofut[NumOfSegTOF];
  int tofdt[NumOfSegTOF];
  double tofctu[NumOfSegTOF];
  double tofctd[NumOfSegTOF];

  int nhLc;
  int LcSeg[MaxHits];
  double tLc[MaxHits];
  double deLc[MaxHits];
  int lcua[NumOfSegLC];
  int lcda[NumOfSegLC];

  int nhAc;
  int AcSeg[NumOfSegAC];
  double tAc[NumOfSegAC];
  double deAc[NumOfSegAC];
  double npeAc;

  //SKS calc hodo
  enum _Particle{Pi, Kaon, Proton, sizeParticle};
  double tTofCalc[sizeParticle];
  double tLcCalc[sizeParticle];
  double ctTofSeg[NumOfSegTOF];
  double ctLcSeg[NumOfSegLC];
};
static Event event;

bool EventSksHodoscope::ProcessingBegin()
{
  using namespace h_Utility;
  static bool fl_init = false;
  if(fl_init){
    ParamMan& CutSKS = FileManager::GetInstance().GetParamMan("CutSKS");
    int select = (int)CutSKS.GetParameter(0, "select", "select");
    CutSKS.stChangeParamSet("CutSKS");
    
    // TOF offset
    offset = CutSKS.GetParameter(select, "offset");

    // TOF Cut Parameters
    MinDeltaETof = CutSKS.GetParameter(select, "MinDeltaETof");
    MaxDeltaETof = CutSKS.GetParameter(select, "MaxDeltaETof");
    MinTimeTof   = CutSKS.GetParameter(select, "MinTimeTof");
    MaxTimeTof   = CutSKS.GetParameter(select, "MaxTimeTof");

    // LC Cut Parameters
    MinDeltaELc  = CutSKS.GetParameter(select, "MinDeltaELc");
    MaxDeltaELc  = CutSKS.GetParameter(select, "MaxDeltaELc");
    MinTimeLc    = CutSKS.GetParameter(select, "MinTimeLc");
    MaxTimeLc    = CutSKS.GetParameter(select, "MaxTimeLc");

    // TOF-LC  
    MinToFT2L    = CutSKS.GetParameter(select, "MinToFT2L");
    MaxToFT2L    = CutSKS.GetParameter(select, "MaxToFT2L");
    MinToF2LSeg  = (int)CutSKS.GetParameter(select, "MinToF2LSeg");
    MaxToF2LSeg  = (int)CutSKS.GetParameter(select, "MaxToF2LSeg");

    // DC multiplisity cut
    MaxMultiHitSdcIn  = CutSKS.GetParameter(select, "MaxHitSdcIn");
    MaxMultiHitSdcOut = CutSKS.GetParameter(select, "MaxHitSdcOut");

    //SdcIn vs TOF&LC cut
    ModFactorTof      = CutSKS.GetParameter(select, "ModFactorTof");
    ModFactorLc       = CutSKS.GetParameter(select, "ModFactorLc");
    MinModTimeTof     = CutSKS.GetParameter(select, "MinModTimeTof");
    MaxModTimeTof     = CutSKS.GetParameter(select, "MaxModTimeTof");
    MinModTimeLc      = CutSKS.GetParameter(select, "MinModTimeLc");
    MaxModTimeLc      = CutSKS.GetParameter(select, "MaxModTimeLc");

    // TOF-SdcOut
    MinXDifTof  = CutSKS.GetParameter(select, "MinXDifTof");
    MaxXDifTof  = CutSKS.GetParameter(select, "MaxXDifTof");
    MinYDifTof  = CutSKS.GetParameter(select, "MinYDifTof");
    MaxYDifTof  = CutSKS.GetParameter(select, "MaxYDifTof");

    // LC-SdcOut
    MinXDifLc   = CutSKS.GetParameter(select, "MinXDifLc");
    MaxXDifLc   = CutSKS.GetParameter(select, "MaxXDifLc");
    MinYDifLc   = CutSKS.GetParameter(select, "MinYDifLc");
    MaxYDifLc   = CutSKS.GetParameter(select, "MaxYDifLc");
    
    fl_init = false;
  }

 return true;
}

bool EventSksHodoscope::ProcessingNormal()
{
  using namespace h_Utility;
  static const std::string funcname = "[EventSksHodoscope::ProcessingNormal]";

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
    //ParamMan& Trig = FileManager::GetInstance().GetParamMan("TriggerFlag");
    //Trig.stChangeParamSet("TriggerFlag");
    //static int TargetTrig = (int)Trig.GetParameter(0, "TARGET");
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

  //------------------------Cut
  //if(!fl_trigger){return true;}

  HF1( 1, 1. );

  //////////////BH2=>time0
  hodoAna->DecodeBH2Hits(rawData);
  int nhBh2=hodoAna->GetNHitsBH2();
  HF1( BH2Hid+10, double(nhBh2) );
  if( nhBh2==0 ){
  //  if( nhBh2!=1 ){
    //std::cout<<"#D NHitsBH2 != 1 [return]"<<std::endl; 
    { tree->Fill(); return true; }
  }
  double time0=hodoAna->GetHitBH2(0)->CTime0();
  event.time0 = time0;
  HF1( BH2Hid+50, time0 );

  HF1( 1, 2. );

  //////////////pre Event Selection by TOF
  hodoAna->DecodeTOFHits(rawData);
  int ncTof=hodoAna->GetNClustersTOF();
  event.nhTof=ncTof;
  {
    int ncOk=0;
    for( int i=0; i<ncTof; ++i ){
      HodoCluster *cl=hodoAna->GetClusterTOF(i);
      int cs=cl->ClusterSize();
      double ms=cl->MeanSeg()+1, cmt=cl->CMeanTime(),
	de=cl->DeltaE(), dt=cl->TimeDif(),
	x=(ms-16.5)*70., y=dt*50.;
      double t=cmt-time0;	

      event.TofSeg[i]=int(ms);
      event.tTof[i]=t;
      event.deTof[i]=de;

      //------------------------Cut
      if( true 
	  //	  && MinDeltaETof<de 
	  //	  && de<MaxDeltaETof 
	  && MinTimeTof  < t  
	  && t < MaxTimeTof
	  ){
        ++ncOk;
	HF1( TOFHid+31, double(cs) );
	HF1( TOFHid+32, ms-0.5 );
	HF1( TOFHid+33, cmt ); HF1( TOFHid+34, de ); HF1( TOFHid+35, dt );
	HF1( TOFHid+36, x ); HF1( TOFHid+37, y ); HF2( TOFHid+38, x, y );
	HF1( TOFHid+39, t );
      }
      else{
	//cl->GoodForAnalysis( false );
      }
    }
    if( ncOk<1 ){
      //std::cout<<"#D Good TofCluster was not found. [return]"<<std::endl;
      //return true;
    }
    HF1( TOFHid+30, double(ncTof) );
  }

  HF1( 1, 3. );

  //////////////pre Event Selection by LC
  hodoAna->DecodeLCHits(rawData);
  int ncLc=hodoAna->GetNClustersLC();
  event.nhLc=ncLc;
  {
    int ncOk=0;
    for( int i=0; i<ncLc; ++i ){
      HodoCluster *cl=hodoAna->GetClusterLC(i);
      int cs=cl->ClusterSize();
      double ms=cl->MeanSeg()+1, cmt=cl->CMeanTime(),
	de=cl->DeltaE(), dt=cl->TimeDif(),
	x=(ms-14.5)*100., y=dt*50.;
      double t=cl->CMeanTime()-time0;

      event.LcSeg[i]=int(ms);
      event.tLc[i]=t;
      event.deLc[i]=de;

      //------------------------Cut
      if( true 
	  //	  && MinDeltaELc<de 
	  //	  && de<MaxDeltaELc 
	  && MinTimeLc <t  
	  && t<MaxTimeLc 
	  ){
        ++ncOk;
	HF1( LCHid+31, double(cs) );
	HF1( LCHid+32, ms-0.5 );
	HF1( LCHid+33, cmt ); HF1( LCHid+34, de ); HF1( LCHid+35, dt );
	HF1( LCHid+36, x ); HF1( LCHid+37, y ); HF2( LCHid+38, x, y );
	HF1( LCHid+39, t );
      }
      else{
	//cl->GoodForAnalysis( false );
      }
    }
    if( ncOk<1 ){
      //std::cout<<"#D Good LcCluster was not found. [return]"<<std::endl;
      //return true;
    }
    HF1( LCHid+30, double(ncLc) );
  }

  HF1( 1, 4. );
  
  //////////////pre Event Selection by TOF-LC
  /*
  {
    int nc=0;
    for( int iTof=0; iTof<ncTof; ++iTof ){
      HodoCluster *clTof=hodoAna->GetClusterTOF(iTof);
      if( !clTof || !clTof->GoodForAnalysis() ) continue;
      double ttof=clTof->CMeanTime()-time0,
        segTof=clTof->MeanSeg()+1;
      for( int iLc=0; iLc<ncLc; ++iLc ){
        HodoCluster *clLc=hodoAna->GetClusterLC(iLc);
        if( !clLc || !clLc->GoodForAnalysis() ) continue;
        double tlc=clLc->CMeanTime()-time0,
          segLc=clLc->MeanSeg()+1;
	//------------------------Cut
	if( true 
	    && MinToF2LSeg<=segTof-segLc && segTof-segLc<=MaxToF2LSeg
	    && MinToFT2L<tlc-ttof && tlc-ttof<MaxToFT2L
	    ){
          ++nc;
	  HF2( 111, segTof-0.5, segLc-0.5 );
	  HF1( 113, segTof-segLc );
	  HF1( 115, tlc-ttof );
        }
        else{
	  clTof->GoodForAnalysis( false );
	  clLc->GoodForAnalysis( false );
        }
      }
    }
    if( nc<1 ){
      //std::cout<<"#D Good TOF-LC-Combination was not found. [return]"<<std::endl;
      { tree->Fill(); return true; }
    }
  }
  HF1( 1, 5. );
  */
  //tree->Fill();return true;
  
  DCAna->DecodeRawHits( rawData );
  
  int IdTof = DCGeomMan::GetInstance().GetTofId();
  int IdLc  = DCGeomMan::GetInstance().GetLcId();
  double zTof = DCGeomMan::GetInstance().GetLocalZ( IdTof ); 
  double zLc = DCGeomMan::GetInstance().GetLocalZ( IdLc ); 

  double multi_SdcIn=0.;
  double multi_SdcOut=0.;

  //////////////SDC1&2 number of hit in one layer not 0
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

  //////////////SDC3&4 number of hit in one layer not 0
  //////////////SDC3-X1 vs SDC4-X2 chekc
  {
    int nhX1=0, nhX2=0;
    for( int layer=1; layer<=NumOfLayersSdcOut; ++layer ){
      const DCHitContainer &contOut =DCAna->GetSdcOutHC(layer);
      int nhOut=contOut.size();
      multi_SdcOut += double(nhOut);
      for( int i=0; i<nhOut; ++i ){
	DCHit *hit=contOut[i];
	double wire=hit->GetWire();
	if( layer == 2 ){ 
	  event.x1wire[nhX1] = int(wire);
	  nhX1++;
	}
	if( layer == 11 ){ 
	  event.x2wire[nhX2] = int(wire);
	  nhX2++;
	}
      }
    }
    event.nhx1=nhX1;
    event.nhx2=nhX2;
  }
  //------------------------Cut
  if( multi_SdcOut/double(NumOfLayersSdcOut) > MaxMultiHitSdcOut )
    { tree->Fill(); return true; }

  HF1( 1, 6. );

  //tree->Fill();return true;

  //////////////SdcIn tracking
  //std::cout<<"#D Begin TrackSearchSdcIn()"<<std::endl;
  DCAna->TrackSearchSdcIn();
  //std::cout<<"#D End TrackSearchSdcIn()"<<std::endl;

  int ntIn=DCAna->GetNtracksSdcIn();
  event.ntIn=ntIn;
  {
    int ntOk=0;
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
	int layerId=hit->GetLayer();
	HF1( 13, layerId );
      }
      //------------------------Cut
      if( chisqr<MaxChisqrSdcIn ) ++ntOk;
      else
        tp->GoodForTracking( false );
    }
    HF1( 10, double(ntOk) );
    if( ntOk<1 ){ tree->Fill(); return true; }
  }

  HF1( 1, 7. );
  /*
  //////////////SdcIn vs Tof&LC cut
  {
    int ntOk=0;
    for( int it=0; it<ntIn; ++it ){
      DCLocalTrack *tp=DCAna->GetTrackSdcIn(it);
      if(!tp) continue;
    
      int nh=tp->GetNHit();
      double chisqr=tp->GetChiSquare();
      double u0=tp->GetU0(), v0=tp->GetV0();
      double x0=tp->GetX0(), y0=tp->GetY0();
      
      //evDisp
      //       if (FlagEvDisp) {
      // 	const EvDisp & evDisp = EvDisp::GetInstance();
      // 	evDisp.DrawSdcInLocalTrack(tp);
      //       }
      
      bool condTof=false, condLc=false;
      for( int j=0; j<ncTof; ++j ){
	HodoCluster *clTof=hodoAna->GetClusterTOF(j);
	if( !clTof || !clTof->GoodForAnalysis() ) continue;
	double ttof=clTof->CMeanTime()-time0;
	HF2( 1111, u0, ttof ); HF2( 1113, u0, ttof+ModFactorTof*u0 );
	//------------------------Cut
	if(true
	   && MinModTimeTof< ttof+ModFactorTof*u0 
	   && ttof+ModFactorTof*u0 <MaxModTimeTof 
	   ){
	  condTof=true;
	}
      }
      for( int j=0; j<ncLc; ++j ){
	HodoCluster *clLc=hodoAna->GetClusterLC(j);
	if( !clLc || !clLc->GoodForAnalysis() ) continue;
	double tlc=clLc->CMeanTime()-time0;
	HF2( 1112, u0, tlc ); HF2( 1114, u0, tlc+ModFactorLc*u0 );
	//------------------------Cut
	if(true
	   && MinModTimeLc< tlc+ModFactorLc*u0 
	   && tlc+ModFactorLc*u0 <MaxModTimeLc
	   ){
	  condLc=true;
	}
      }
      if( condTof && condLc ){
	++ntOk;
	HF1( 1152, double(nh) );
	HF1( 1153, chisqr );
	HF1( 1154, x0 ); HF1( 1155, y0 );
	HF1( 1156, u0 ); HF1( 1157, v0 );
	HF2( 1158, x0, u0 ); HF2( 1159, y0, v0 );

	for( int j=0; j<ncTof; ++j ){
	  HodoCluster *clTof=hodoAna->GetClusterTOF(j);
	  if( !clTof || !clTof->GoodForAnalysis() ) continue;
	  double ttof=clTof->CMeanTime()-time0;
	  HF2( 1161, u0, ttof ); HF2( 1163, u0, ttof+12.5*u0 );
	}
	for( int j=0; j<ncLc; ++j ){
	  HodoCluster *clLc=hodoAna->GetClusterLC(j);
	  if( !clLc || !clLc->GoodForAnalysis() ) continue;
	  double tlc=clLc->CMeanTime()-time0;
	  HF2( 1162, u0, tlc ); HF2( 1164, u0, tlc+12.5*u0 );
	}
	// For quick tracking
	//	if(ntOk>0) tp->GoodForTracking( false );
      }	
      else {
	tp->GoodForTracking( false );
      }
    }
    if( ntOk<1 ){ tree->Fill(); return true; }
  }
  */
  HF1( 1, 8. );

  //////////////SdcOut tracking
  //std::cout<<"#D Begin TrackSearchSdcOut()"<<std::endl;
  DCAna->TrackSearchSdcOut();
  //std::cout<<"#D End TrackSearchSdcOut()"<<std::endl;
  
  int ntOut=DCAna->GetNtracksSdcOut();
  event.ntOut=ntOut;
  {
    int ntOk=0;
    for( int it=0; it<ntOut; ++it ){
      DCLocalTrack *tp=DCAna->GetTrackSdcOut(it);
      int nh=tp->GetNHit();
      double chisqr=tp->GetChiSquare();
      double u0=tp->GetU0(), v0=tp->GetV0();
      double x0=tp->GetX(zTof), y0=tp->GetY(zTof);
      double xlc=tp->GetX(zLc), ylc=tp->GetY(zLc);

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
   
      bool condTof=false, condLc=false;
      for( int j=0; j<ncTof; ++j ){
	HodoCluster *clTof=hodoAna->GetClusterTOF(j);
        //if( !clTof || !clTof->GoodForAnalysis() ) continue;
        double x=(clTof->MeanSeg()-15.5)*70.;
        double y=(clTof->TimeDif())*50.;
	HF1( TOFHid+41, x0-x ); HF1( TOFHid+42, y0-y );
	//------------------------Cut
	//	if( MinXDifTof<x0-x && x0-x<MaxXDifTof &&
	//            MinYDifTof<y0-y && y0-y<MaxYDifTof )
          condTof=true;
      }
      for( int j=0; j<ncLc; ++j ){
        HodoCluster *clLc=hodoAna->GetClusterLC(j);
        //if( !clLc || !clLc->GoodForAnalysis() ) continue;
        double x=(clLc->MeanSeg()-13.5)*100.;
        double y=(clLc->TimeDif())*50.;
	HF1( LCHid+41, xlc-x ); HF1( LCHid+42, ylc-y );
	//------------------------Cut
	//	if( MinXDifLc<xlc-x && xlc-x<MaxXDifLc &&
	//            MinYDifLc<ylc-y && ylc-y<MaxYDifLc )
          condLc=true;
      }
      if( condTof && condLc ){
	HF1( 31, double(nh) ); HF1( 32, chisqr ); 
	HF1( 34, x0 ); HF1( 35, y0 );
	HF1( 36, u0 ); HF1( 37, v0 );
	HF2( 38, x0, u0 ); HF2( 39, y0, v0 );
	for( int ih=0; ih<nh; ++ih ){
	  DCLTrackHit *hit=tp->GetHit(ih);
	  int layerId=hit->GetLayer();
	  HF1( 33, layerId );
	}
        if( chisqr<MaxChisqrSdcOut ) ++ntOk;
	else 
	  tp->GoodForTracking( false );
      }
      else{
        tp->GoodForTracking( false );
      }
    }
    HF1( 30, double(ntOk) );
    //if( ntOk<1 ) return true;
  }

  HF1( 1, 9. );

  //tree->Fill(); return true;

  //////////////SKS tracking
  //std::cout<<"#D Begin TrackSearchSks()"<<std::endl;
  DCAna->TrackSearchSks();
  //std::cout<<"#D End TrackSearchSks()"<<std::endl;

  int ntSks=DCAna->GetNTracksSks();
  event.ntSks = ntSks; 
  HF1( 50, double(ntSks) );
  {
    int ntOk=0;
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

      //std::cout<<"******************************"<<std::endl;
      double m2;
      for( int j=0; j<ncTof; ++j ){
	HodoCluster *clTof=hodoAna->GetClusterTOF(j);
	if( !clTof ) continue;
	double time = clTof->CMeanTime()-time0+offset;
	
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
	HF1(53, layerId );
      }
      //------------------------Cut
      if( chisqr<MaxChisqrSks ) ++ntOk;
      else continue;
      
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
	//HF2( 64, x0out, p );
      }
      //HF1( 50, double(ntOk) );
      if( ntOk<1 ){ tree->Fill(); return true; }
    }
  }

  HF1( 1, 10. );

  //------------------------Cut
  if( ntSks!=1 ){ tree->Fill(); return true; }

  HF1( 1, 11. );

  SksTrack *track=DCAna->GetSksTrack(0);
  double plTof=track->PathLengthToTOF();
  double plLc=track->PathLengthTotal();
  double p=track->PrimaryMomentum().Mag();
  double tpiTof=plTof*sqrt(p*p+PionMass*PionMass)/p/LightVel;
  double tpiLc=plLc*sqrt(p*p+PionMass*PionMass)/p/LightVel;
  double tkTof=plTof*sqrt(p*p+KaonMass*KaonMass)/p/LightVel;
  double tkLc=plLc*sqrt(p*p+KaonMass*KaonMass)/p/LightVel;
  double tpTof=plTof*sqrt(p*p+ProtonMass*ProtonMass)/p/LightVel;
  double tpLc=plLc*sqrt(p*p+ProtonMass*ProtonMass)/p/LightVel;

  HF1( 71, plTof ); HF1( 72, plLc ); HF1( 73, p );
  HF1( 74, tpiTof ); HF1( 75, tpiLc );
  HF1( 76, tkTof ); HF1( 77, tkLc );
  HF1( 78, tpTof ); HF1( 79, tpLc );

  event.tTofCalc[Event::Pi]     = tpiTof;
  event.tTofCalc[Event::Kaon]   = tkTof;
  event.tTofCalc[Event::Proton] = tpTof;

  event.tLcCalc[Event::Pi]      = tpiLc;
  event.tLcCalc[Event::Kaon]    = tkLc;
  event.tLcCalc[Event::Proton]  = tpLc;

  // TOF
  {
    int nh=hodoAna->GetNHitsTOF();
    HF1( TOFHid+50, double(nh) );
     for( int i=0; i<nh; ++i ){
      Hodo2Hit *hit=hodoAna->GetHitTOF(i);
      if(!hit) continue;
      int seg=hit->SegmentId()+1;
      HF1( TOFHid+51, seg-0.5 );
      double ctu=hit->GetCTUp(), ctd=hit->GetCTDown();
      double cmt=hit->CMeanTime(), t= cmt-time0+offset;//cmt-time0;
      double de=hit->DeltaE();
      double m2=MassSquare( p, plTof, t );
      double beta=plTof/t/LightVel;
      HF1( TOFHid+100*seg+31, cmt );
      HF1( TOFHid+100*seg+32, ctd-ctu );
      HF1( TOFHid+100*seg+33, t );
      HF1( TOFHid+100*seg+34, t-tpiTof );
      HF1( TOFHid+100*seg+35, t-tkTof );
      HF1( TOFHid+100*seg+36, t-tpTof );
      HF1( TOFHid+100*seg+37, m2 );
      HF1( TOFHid+100*seg+38, beta );
      HF1( TOFHid+52, cmt ); HF1( TOFHid+53, de );
      HF1( TOFHid+54, ctd-ctu ); HF1( TOFHid+55, t );
      HF1( TOFHid+56, t-tpiTof ); HF1( TOFHid+57, t-tkTof );
      HF1( TOFHid+58, t-tpTof ); HF1( TOFHid+59, m2 );
      HF1( TOFHid+60, beta );
      event.ctTofSeg[seg-1] = t;
      event.tofctu[seg-1] = ctu;
      event.tofctd[seg-1] = ctd;

    }

     const HodoRHitContainer &cont=rawData->GetTOFRawHC();
     int NofHit = cont.size();
     for(int i = 0; i<NofHit; ++i){
       HodoRawHit *hit = cont[i];
       int         seg = hit->SegmentId();
       event.tofua[seg] = hit->GetAdcUp();
       event.tofda[seg] = hit->GetAdcDown();
       event.tofut[seg] = hit->GetTdcUp();
       event.tofdt[seg] = hit->GetTdcDown();
 
     }
  }

  // LC
  {
    int nh=hodoAna->GetNHitsLC();
    HF1( LCHid+50, double(nh) );
    int nh2=0;
    for( int i=0; i<nh; ++i ){
      Hodo2Hit *hit=hodoAna->GetHitLC(i);
      if(!hit) continue;
      int seg=hit->SegmentId()+1;
      HF1( LCHid+51, seg-0.5 );
      double ctu=hit->GetCTUp(), ctd=hit->GetCTDown();
      double cmt=hit->CMeanTime(), t= cmt-time0+offset;//cmt-time0;
      double de=hit->DeltaE();
      double m2=MassSquare( p, plLc, t );
      double beta=plLc/t/LightVel;
      HF1( LCHid+100*seg+31, cmt );
      HF1( LCHid+100*seg+32, ctd-ctu );
      HF1( LCHid+100*seg+33, t );
      HF1( LCHid+100*seg+34, t-tpiLc );
      HF1( LCHid+100*seg+35, t-tkLc );
      HF1( LCHid+100*seg+36, t-tpLc );
      HF1( LCHid+100*seg+37, m2 );
      HF1( LCHid+100*seg+38, beta );
      HF1( LCHid+52, cmt ); HF1( LCHid+53, de );
      HF1( LCHid+54, ctd-ctu ); HF1( LCHid+55, t );
      HF1( LCHid+56, t-tpiLc ); HF1( LCHid+57, t-tkLc );
      HF1( LCHid+58, t-tpLc ); HF1( LCHid+59, m2 );
      HF1( LCHid+60, beta );
      event.ctLcSeg[seg-1] = t;
    }

     const HodoRHitContainer &cont=rawData->GetLCRawHC();
     int NofHit = cont.size();
     for(int i = 0; i<NofHit; ++i){
       HodoRawHit *hit = cont[i];
       int         seg = hit->SegmentId();
       event.lcua[seg] = hit->GetAdcUp();
       event.lcda[seg] = hit->GetAdcDown();
     }
  }

  // AC
  /*
  {
    hodoAna->DecodeACHits(rawData);
    int nh=hodoAna->GetNHitsAC();
    HF1( ACHid+10, double(nh) );
    event.nhAc = nh;
    int nh2=0;
    double npeac;
    for( int i=0; i<nh; ++i ){
      Hodo1Hit *hit=hodoAna->GetHitAC(i);
      if(!hit) continue;
      int seg=hit->SegmentId()+1;
      HF1( ACHid+11, seg-0.5 );
      double a=hit->GetA(), t=hit->GetT(), ct=hit->GetCT();
      HF1( ACHid+12, a );
      HF1( ACHid+13, t );
    
      event.AcSeg[i] = seg;
      event.deAc[i] = a;
      event.tAc[i] = t;
      if( t >0 ) npeac += a;
      if(a>0.1)	++nh2;
      event.npeAc = npeac;
    
      HF1( ACHid+20, double(nh2) );
    }
  }
  
  // TOF-LC-AC
  {
    int nhTof=hodoAna->GetNHitsTOF();
    int nhLc =hodoAna->GetNHitsLC();
    int nhAc =hodoAna->GetNHitsAC();
    for( int iTof=0; iTof<nhTof; ++iTof ){
      Hodo2Hit *hitTof=hodoAna->GetHitTOF(iTof);
      if(!hitTof) continue;
      int segTof=hitTof->SegmentId()+1; 
      double cmtTof=hitTof->CMeanTime();
      double deTof=hitTof->DeltaE();
      double m2=MassSquare( p, plTof, cmtTof-time0+offset );
      HF2( 301, m2, deTof );
      for( int iLc=0; iLc<nhLc; ++iLc ){
        Hodo2Hit *hitLc =hodoAna->GetHitLC(iLc);
        if(!hitLc) continue;
        int segLc=hitLc->SegmentId()+1;
        double cmtLc=hitLc->CMeanTime();
        double deLc=hitLc->DeltaE();
        HF2( 201, segTof-0.5, segLc-0.5 );
        HF1( 203, cmtLc-cmtTof );
	HF2( 302, m2, deLc );
	HF2( 303, m2, cmtLc-cmtTof );
      }
      for( int iAc=0; iAc<nhAc; ++iAc ){
	Hodo1Hit *hit=hodoAna->GetHitAC(iAc);
	if(!hit) continue;
	double a=hit->GetA();
	HF2( 304, m2, a);
      }
    }
  }
  */
  HF1( 1, 12. );

  tree->Fill();

  return true;
}

void EventSksHodoscope::InitializeEvent( void )
{
  event.runnum = -1;
  event.evnum = -1;
  event.time0 = -9999.;
  event.ntIn     = 0;
  event.ntOut    = 0;
  event.ntSks    = 0;
  event.nhx1     = 0;
  event.nhx2     = 0;

  for( int it=0; it<NumOfMisc; it++){
    event.trigpat[it] = -1;
    event.trigflag[it] = -1;
  }

  for( int it=0; it<NumOfMatrix; it++){
    event.matrixpat[it] = -1;
    event.matrixflag[it] = -1;
  }
    
  for( int it=0; it<MaxHits; it++){
    event.x1wire[it] = -1;
    event.x2wire[it] = -1;
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

  //Hodoscope
  event.nhTof  = 0;
  event.nhLc   = 0;
  event.nhAc   = 0;
  event.npeAc  = -9999.0;

  for( int it=0; it<MaxHits; it++){
    event.TofSeg[it] = -1;
    event.tTof[it] = -9999.0;
    event.deTof[it] = -9999.0;

    event.LcSeg[it] = -1;
    event.tLc[it] = -9999.0;
    event.deLc[it] = -9999.0;

    event.AcSeg[it] = -1;
    event.tAc[it] = -9999.0;
    event.deAc[it] = -9999.0;
  }

  // SKS hodo calc
  for(int i = 0; i<Event::sizeParticle; ++i){
    event.tTofCalc[i] = -9999.;
    event.tLcCalc[i]  = -9999.;
  }

  for(int i = 0; i<NumOfSegTOF; ++i){
    event.ctTofSeg[i] = -9999.;
    event.tofua[i]    = -1;
    event.tofda[i]    = -1;
    event.tofut[i]    = -1;
    event.tofdt[i]    = -1;
    event.tofctu[i]    = -1;
    event.tofctd[i]    = -1;

  }

  for(int i = 0; i<NumOfSegLC; ++i){
    event.ctLcSeg[i] = -9999.;
    event.lcua[i]    = -1;
    event.lcda[i]    = -1;
  }
}

bool EventSksHodoscope::ProcessingEnd()
{
  
  return true;
}

VEvent *ConfMan::EventAllocator()
{
  return new EventSksHodoscope;
}

bool ConfMan::InitializeHistograms( void )
{
  HB1( 1, "Status", 20, 0., 20. );

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
  HB1( 61, "P SksTrack", 400, 0.50, 0.90 );
  HB1( 62, "PathLength SksTrack", 600, 3000., 6000. );
  HB1( 63, "MassSqr", 600, -0.2, 1.2 );
  //   HB1( 62, "PathLength SksTrack", 400, 3000., 7000. );
  //   HB2( 63, "P%PathLength SksTrack", 100, 3000., 7000., 100, 0.60, 0.90 );
  //   HB2( 64, "P%X(zTof) SksTrack", 200, -1200., 1200., 100, 0.60, 0.90 );

  HB1( 71, "PathLength to TOF", 400, 3000., 7000. );
  HB1( 72, "PathLength to LC", 400, 3000., 7000. );
  HB1( 73, "P", 300, 0.60, 0.90 );
  HB1( 74, "Tpi Tof", 500, -5., 45. );
  HB1( 75, "Tpi Lc", 500, -5., 45. );
  HB1( 76, "Tk Tof", 500, -5., 45. );
  HB1( 77, "Tk Lc", 500, -5., 45. );
  HB1( 78, "Tp Tof", 500, -5., 45. );
  HB1( 79, "Tp Lc", 500, -5., 45. );

  // BH2
  HB1( BH2Hid+10, "#Hits BH2[Hodo]", NumOfSegBH2+1, 0., double(NumOfSegBH2+1) );
  HB1( BH2Hid+50, "CTime0 BH2", 200, -10., 10.);

  // Tof
  HB1( TOFHid+30, "#Clusters Tof", NumOfSegTOF+1, 0., double(NumOfSegTOF+1) );
  HB1( TOFHid+31, "ClusterSize Tof", 5, 0., 5. ); 
  HB1( TOFHid+32, "HitPat Cluster Tof", 2*NumOfSegTOF, 0., double(NumOfSegTOF) );
  HB1( TOFHid+33, "CMeamTime Cluster Tof", 500, -5., 45. ); 
  HB1( TOFHid+34, "DeltaE Cluster Tof", 200, -0.5, 4.5 );
  HB1( TOFHid+35, "TimeDif Cluster Tof", 200, -10., 10. );
  HB1( TOFHid+36, "PosX Cluster Tof", 2*NumOfSegTOF, -double(NumOfSegTOF)*35., double(NumOfSegTOF)*35. );
  HB1( TOFHid+37, "PosY Cluster Tof", 200, -600., 600. );
  HB2( TOFHid+38, "PosY%PosX Cluster Tof", 2*NumOfSegTOF, -double(NumOfSegTOF)*35., double(NumOfSegTOF)*35., 200, -600., 600. );
  HB1( TOFHid+39, "CMeamTime-Time0 Cluster Tof", 500, -5., 45. ); 

  HB1( TOFHid+41, "X(zTof)-PosXTof", 500, -500., 500. );
  HB1( TOFHid+42, "Y(zTof)-PosYTof", 500, -500., 500. );

  HB1( TOFHid+50, "#Hits Tof[SksHodo]", NumOfSegTOF+1, 0., double(NumOfSegTOF+1) );
  HB1( TOFHid+51, "Hitpat Tof[SksHodo]", NumOfSegTOF, 0., double(NumOfSegTOF) );
  HB1( TOFHid+52, "CMeanTime Tof[SksHodo]", 500, -5., 45. );
  HB1( TOFHid+53, "dE Tof[SksHodo]", 200, -0.5, 4.5 );
  HB1( TOFHid+54, "Tdown-Tup Tof[SksHodo]", 200, -10.0, 10.0 );
  HB1( TOFHid+55, "TofTime[SksHodo]", 500, -5., 45. );
  HB1( TOFHid+56, "TofTime-Tpi", 500, -25., 25. );
  HB1( TOFHid+57, "TofTime-Tk", 500, -25., 25. );
  HB1( TOFHid+58, "TofTime-Tp", 500, -25., 25. );
  HB1( TOFHid+59, "MassSquare Tof", 600, -0.1, 1.1 );
  HB1( TOFHid+60, "beta Tof", 500, 0., 1.5 );

  for( int i=1; i<=NumOfSegTOF; ++i ){
    std::stringstream title1,title2;
    title1 << "Tof-" << i << " CMeanTime[SksHodo]";
    HB1( TOFHid+100*i+31, title1.str().c_str(), 500, -5., 45. );
    title2 << "TOF-" << i << " Tdown-Tup[SksHodo]";
    HB1( TOFHid+100*i+32, title2.str().c_str(), 200, -10.0, 10.0 );
    std::stringstream title3,title4,title5,title6;
    title3 << "Tof-" << i << " TofTime[SksHodo]";
    HB1( TOFHid+100*i+33, title3.str().c_str(), 500, -5., 45. );
    title4 << "Tof-" << i << " TofTime-Tpi";
    HB1( TOFHid+100*i+34, title4.str().c_str(), 500, -25., 25. );
    title5 << "Tof-" << i << " TofTime-Tk";
    HB1( TOFHid+100*i+35, title5.str().c_str(), 500, -25., 25. );
    title6 << "Tof-" << i << " TofTime-Tp";
    HB1( TOFHid+100*i+36, title6.str().c_str(), 500, -25., 25. );
    std::stringstream title7,title8;
    title7 << "Tof-" << i << " MassSquare";
    HB1( TOFHid+100*i+37, title7.str().c_str(), 600, -0.1, 1.1 );
    title8 << "Tof-" << i << " beta";
    HB1( TOFHid+100*i+38, title8.str().c_str(), 500, 0., 1.5 );
  }

  // Lc
  HB1( LCHid+30, "#Clusters Lc", NumOfSegLC+1, 0., double(NumOfSegLC+1) );
  HB1( LCHid+31, "ClusterSize Lc", 5, 0., 5. ); 
  HB1( LCHid+32, "HitPat Cluster Lc", 2*NumOfSegLC, 0., double(NumOfSegLC) );
  HB1( LCHid+33, "CMeamTime Cluster Lc", 500, -5., 45. ); 
  HB1( LCHid+34, "DeltaE Cluster Lc", 200, -0.5, 4.5 );
  HB1( LCHid+35, "TimeDif Cluster Lc", 200, -15., 15. );
  HB1( LCHid+36, "PosX Cluster Lc", 2*NumOfSegLC, -double(NumOfSegLC)*50., double(NumOfSegLC)*50. );
  HB1( LCHid+37, "PosY Cluster Lc", 200, -800., 800. );
  HB2( LCHid+38, "PosY%PosX Cluster Lc", 2*NumOfSegLC, -double(NumOfSegLC)*50., double(NumOfSegLC)*50., 200, -800., 800. );
  HB1( LCHid+39, "CMeamTime-Time0 Cluster Lc", 500, -5., 45. ); 

  HB1( LCHid+41, "X(zLc)-PosXLc", 500, -500., 500. );
  HB1( LCHid+42, "Y(zLc)-PosYLc", 500, -500., 500. );

  HB1( LCHid+50, "#Hits Lc[SksHodo]", NumOfSegLC+1, 0., double(NumOfSegLC+1) );
  HB1( LCHid+51, "Hitpat Lc[SksHodo]", NumOfSegLC, 0., double(NumOfSegLC) );
  HB1( LCHid+52, "CMeanTime Lc[SksHodo]", 500, -5., 45. );
  HB1( LCHid+53, "dE Lc[SksHodo]", 200, -0.5, 4.5 );
  HB1( LCHid+54, "Tdown-Tup Lc[SksHodo]", 200, -10.0, 10.0 );
  HB1( LCHid+55, "LcTime[SksHodo]", 500, -5., 45. );
  HB1( LCHid+56, "LcTime-Tpi", 500, -25., 25. );
  HB1( LCHid+57, "LcTime-Tk", 500, -25., 25. );
  HB1( LCHid+58, "LcTime-Tp", 500, -25., 25. );
  HB1( LCHid+59, "MassSquare Lc", 600, -0.1, 1.1 );
  HB1( LCHid+60, "beta Lc", 500, 0., 1.5 );

  for( int i=1; i<=NumOfSegLC; ++i ){
    std::stringstream title1,title2;
    title1 << "Lc-" << i << " CMeanTime[SksHodo]";
    HB1( LCHid+100*i+31, title1.str().c_str(), 500, -5., 45. );
    title2 << "Lc-" << i << " Tdown-Tup[SksHodo]";
    HB1( LCHid+100*i+32, title2.str().c_str(), 200, -10.0, 10.0 );
    std::stringstream title3,title4,title5,title6;
    title3 << "Lc-" << i << " LcTime[SksHodo]";
    HB1( LCHid+100*i+33, title3.str().c_str(), 500, -5., 45. );
    title4 << "Lc-" << i << " LcTime-Tpi";
    HB1( LCHid+100*i+34, title4.str().c_str(), 500, -25., 25. );
    title5 << "Lc-" << i << " LcTime-Tk";
    HB1( LCHid+100*i+35, title5.str().c_str(), 500, -25., 25. );
    title6 << "Lc-" << i << " LcTime-Tp";
    HB1( LCHid+100*i+36, title6.str().c_str(), 500, -25., 25. );
    std::stringstream title7,title8;
    title7 << "Lc-" << i << " MassSquare";
    HB1( LCHid+100*i+37, title7.str().c_str(), 600, -0.1, 1.1 );
    title8 << "Lc-" << i << " beta";
    HB1( LCHid+100*i+38, title8.str().c_str(), 500, 0., 1.5 );
  }

  // Ac
  HB1( ACHid+10, "#Hits Ac[Hodo]",  NumOfSegAC+1, 0., double(NumOfSegAC+1) );
  HB1( ACHid+11, "Hitpat Ac[Hodo]", NumOfSegAC,   0., double(NumOfSegAC)   );
  HB1( ACHid+12, "dE Ac", 200, -0.5, 4.5 );
  HB1( ACHid+13, "T Ac", 200,  0.0, 100 );
  HB1( ACHid+20, "#Hits Ac[HodoGood]",  NumOfSegAC+1, 0., double(NumOfSegAC+1) );

  // TOF-LC-AC
  HB2( 111, "SegLc%SegTof Cluster", NumOfSegTOF, 0., double(NumOfSegTOF),
       NumOfSegLC, 0., double(NumOfSegLC) );
  HB1( 113, "SegTof-SegLc Cluster", 30, -15., 15. );
  HB1( 115, "LcTime-TofTime Cluster", 300, -5., 25. );

  HB2( 201, "SegLc%SegTof[SksHodo]", NumOfSegTOF, 0., double(NumOfSegTOF),
       NumOfSegLC, 0., double(NumOfSegLC) );
  HB1( 203, "LcTime-TofTime[SksHodo]", 300, -5., 25. );

  HB2( 301, "TofDe%MSqrTof", 120, -0.1, 1.1, 100, -0.5, 4.5 );
  HB2( 302, "LcDe%MSqrTof", 120, -0.1, 1.1, 100, -0.5, 4.5 );
  HB2( 303, "LcTime-TofTime%MSqrTof", 120, -0.1, 1.1, 100, -5., 25 );
  HB2( 304, "AcDe%MSqrTof", 120, -0.1, 1.1, 100, -0.5, 4.5 );

 ////////////////////////////////////////////
  //Tree
  HBTree("tree","tree of Counter");
  TTree *tree = dynamic_cast<TTree *>(gFile->Get("tree"));

  tree->Branch("runnum",    &event.runnum,    "runnum/I");
  tree->Branch("evnum",     &event.evnum,     "evnum/I");
  tree->Branch("time0",     &event.time0,     "time0/D");
  tree->Branch("trigpat",    event.trigpat,   "trigpat[20]/I");
  tree->Branch("trigflag",   event.trigflag,   "trigflag[20]/I");

  tree->Branch("matrixpat",  event.matrixpat,   "matrixpat[10]/I");
  tree->Branch("matrixflag", event.matrixflag,  "matrixflag[10]/I");

  tree->Branch("nhx1",    &event.nhx1,     "nhx1/I");
  tree->Branch("x1wire",   event.x1wire,   "x1wire[nhx1]/I");
  tree->Branch("nhx2",    &event.nhx2,     "nhx2/I");
  tree->Branch("x2wire",   event.x2wire,   "x2wire[nhx2]/I");

  tree->Branch("ntIn",    &event.ntIn,     "ntIn/I");
  tree->Branch("nhIn",     event.nhIn,     "nhIn[ntIn]/I");
  tree->Branch("chisqrIn", event.chisqrIn, "chisqrIn[ntIn]/D");
  tree->Branch("x0In",     event.x0In,     "x0In[ntIn]/D");
  tree->Branch("y0In",     event.y0In,     "y0In[ntIn]/D");
  tree->Branch("u0In",     event.u0In,     "u0In[ntIn]/D");
  tree->Branch("v0In",     event.v0In,     "v0In[ntIn]/D");

  tree->Branch("ntOut",   &event.ntOut,     "ntOut/I");
  tree->Branch("nhOut",    event.nhOut,     "nhOut[ntOut]/I");
  tree->Branch("chisqrOut",event.chisqrOut, "chisqrIn[ntOut]/D");
  tree->Branch("x0Out",    event.x0Out,     "x0Out[ntOut]/D");
  tree->Branch("y0Out",    event.y0Out,     "y0Out[ntOut]/D");
  tree->Branch("u0Out",    event.u0Out,     "u0Out[ntOut]/D");
  tree->Branch("v0Out",    event.v0Out,     "v0Out[ntOut]/D");

  tree->Branch("tTofCalc",    event.tTofCalc,  "tTofCalc[3]/D");
  tree->Branch("tLcCalc",     event.tLcCalc,   "tLcCalc[3]/D");
  tree->Branch("ctTofSeg",    event.ctTofSeg,  "ctTofSeg[32]/D");
  tree->Branch("ctLcSeg",     event.ctLcSeg,   "ctLcSeg[28]/D");
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

  tree->Branch("nhTof",   &event.nhTof,   "nhTof/I");
  tree->Branch("TofSeg",   event.TofSeg,  "TofSeg[nhTof]/I");
  tree->Branch("tTof",     event.tTof,    "tTof[nhTof]/D");
  tree->Branch("deTof",    event.deTof,   "deTof[nhTof]/D");
  tree->Branch("tofua",    event.tofua,   "tofua[32]/I");
  tree->Branch("tofda",    event.tofda,   "tofda[32]/I");
  tree->Branch("tofut",     event.tofut,    "tofut[32]/I");
  tree->Branch("tofdt",     event.tofdt,    "tofdt[32]/I");
  tree->Branch("tofctu",    event.tofctu,   "tofctu[32]/D");
  tree->Branch("tofctd",    event.tofctd,   "tofctd[32]/D");

  tree->Branch("nhLc",   &event.nhLc,   "nhLc/I");
  tree->Branch("LcSeg",   event.LcSeg,  "LcSeg[nhLc]/I");
  tree->Branch("tLc",     event.tLc,    "tLc[nhLc]/D");
  tree->Branch("deLc",    event.deLc,   "deLc[nhLc]/D");
  tree->Branch("lcua",    event.lcua,   "lcua[28]/I");
  tree->Branch("lcda",    event.lcda,   "lcda[28]/I");

  tree->Branch("nhAc",   &event.nhAc,   "nhAc/I");
  tree->Branch("AcSeg",   event.AcSeg,  "AcSeg[nhAc]/I");
  tree->Branch("tAc",     event.tAc,    "tAc[nhAc]/D");
  tree->Branch("deAc",    event.deAc,   "deAc[nhAc]/D");
  tree->Branch("npeAc",  &event.npeAc,  "npeAc/D");

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

  if(FileManagerName_ != ""){
    if(-1 == FileManager::GetInstance().stInitialize(FileManagerName_)){
      status = false;
    }
  }

  return status;
}



