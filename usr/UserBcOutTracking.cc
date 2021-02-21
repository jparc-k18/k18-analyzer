// -*- C++ -*-

#include <cmath>
#include <iomanip>
#include <iostream>
#include <set>
#include <sstream>
#include <string>

#include "ConfMan.hh"
#include "DCRawHit.hh"
#include "DetectorID.hh"
#include "RMAnalyzer.hh"
#include "MathTools.hh"
#include "RootHelper.hh"
#include "KuramaLib.hh"
#include "RawData.hh"
#include "UnpackerManager.hh"
#include "VEvent.hh"

#define HodoCut 0
#define TotCut  1

namespace
{
  using namespace root;
  const std::string& classname("BcOutTracking");
  RMAnalyzer&         gRM     = RMAnalyzer::GetInstance();
  const UserParamMan& gUser   = UserParamMan::GetInstance();
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
class EventBcOutTracking : public VEvent
{
private:
  RawData      *rawData;
  DCAnalyzer   *DCAna;
  HodoAnalyzer *hodoAna;

public:
          EventBcOutTracking( void );
         ~EventBcOutTracking( void );
  Bool_t  ProcessingBegin( void );
  Bool_t  ProcessingEnd( void );
  Bool_t  ProcessingNormal( void );
  Bool_t  InitializeHistograms( void );
  void    InitializeEvent( void );
};

//______________________________________________________________________________
EventBcOutTracking::EventBcOutTracking( void )
  : VEvent(),
    rawData( nullptr ),
    DCAna( new DCAnalyzer ),
    hodoAna( new HodoAnalyzer )
{
}

//______________________________________________________________________________
EventBcOutTracking::~EventBcOutTracking( void )
{
  if( DCAna )   delete DCAna;
  if( hodoAna ) delete hodoAna;
  if( rawData ) delete rawData;
}

//______________________________________________________________________________
Bool_t
EventBcOutTracking::ProcessingBegin( void )
{
  InitializeEvent();
  return true;
}

//______________________________________________________________________________
struct Event
{
  Int_t evnum;
  Int_t trignhits;
  Int_t trigpat[NumOfSegTrig];
  Int_t trigflag[NumOfSegTrig];

  //BH1
  Int_t nhBh1;
  Double_t tBh1[MaxHits];
  Double_t deBh1[MaxHits];
  Double_t Bh1Seg[MaxHits];

  //BH2
  Int_t nhBh2;
  Double_t tBh2[MaxHits];
  Double_t deBh2[MaxHits];
  Double_t Bh2Seg[MaxHits];

  // Time0
  Double_t Time0Seg;
  Double_t deTime0;
  Double_t Time0;
  Double_t CTime0;

  //Beam
  Int_t pid;
  Double_t btof;

  // BcOut
  Int_t nhit[NumOfLayersBcOut];
  Int_t nlayer;
  Double_t pos[NumOfLayersBcOut][MaxHits];

  // BcOutTracking
  Int_t ntrack;
  Double_t chisqr[MaxHits];
  Double_t x_Bh2[MaxHits];
  Double_t x0[MaxHits];
  Double_t y0[MaxHits];
  Double_t u0[MaxHits];
  Double_t v0[MaxHits];

};

//______________________________________________________________________________
namespace root
{
  Event event;
  TH1   *h[MaxHist];
  TTree *tree;
  enum eParticle
    {
      kKaon, kPion, nParticle
    };
}

//______________________________________________________________________________
Bool_t
EventBcOutTracking::ProcessingNormal( void )
{
  static const std::string funcname("["+classname+"::"+__func__+"]");

#if HodoCut
  static const Double_t MinDeBH2   = gUser.GetParameter("DeBH2", 0);
  static const Double_t MaxDeBH2   = gUser.GetParameter("DeBH2", 1);
  static const Double_t MinDeBH1   = gUser.GetParameter("DeBH1", 0);
  static const Double_t MaxDeBH1   = gUser.GetParameter("DeBH1", 1);
  static const Double_t MinBeamToF = gUser.GetParameter("BTOF",  1);
  static const Double_t MaxBeamToF = gUser.GetParameter("BTOF",  1);
#endif
  static const Double_t MaxMultiHitBcOut = gUser.GetParameter("MaxMultiHitBcOut");
#if TotCut
  static const Double_t MinTotBcOut = gUser.GetParameter("MinTotBcOut", 0);
#endif

  rawData = new RawData;
  rawData->DecodeHits();

  gRM.Decode();

  event.evnum = gRM.EventNumber();

  // Trigger Flag
  std::bitset<NumOfSegTrig> trigger_flag;
  {
    const HodoRHitContainer &cont = rawData->GetTrigRawHC();
    Int_t trignhits = 0;
    Int_t nh = cont.size();
    for( Int_t i=0; i<nh; ++i ){
      HodoRawHit *hit = cont[i];
      Int_t seg = hit->SegmentId()+1;
      Int_t tdc = hit->GetTdc1(); // only first hit
      if( tdc>0 ){
	trigger_flag.set( seg-1 );
	event.trigpat[trignhits++] = seg;
	event.trigflag[seg-1]      = tdc;
      }
    }
    event.trignhits = trignhits;
  }
  HF1( 1, 0. );

  if( trigger_flag[trigger::kSpillEnd] ||
      trigger_flag[trigger::kL1SpillOff] )
    return true;

  HF1( 1, 1. );

 //////////////BH2 time 0
  hodoAna->DecodeBH2Hits( rawData );
  Int_t nhBh2 = hodoAna->GetNHitsBH2();
  event.nhBh2 = nhBh2;
#if HodoCut
  if( nhBh2==0 ) return true;
#endif
  HF1( 1, 2 );

  Double_t time0 = -9999.;
  //////////////BH2 Analysis
  for( Int_t i=0; i<nhBh2; ++i ){
    BH2Hit* hit = hodoAna->GetHitBH2( i );
    if( !hit ) continue;
    Double_t seg = hit->SegmentId()+1;
    // Double_t mt  = hit->MeanTime();
    Double_t cmt = hit->CMeanTime();
    // Double_t ct0 = hit->CTime0();
    Double_t dE  = hit->DeltaE();

#if HodoCut
    if( dE<MinDeBH2 || MaxDeBH2<dE ) continue;
#endif
    event.tBh2[i]   = cmt;
    event.deBh2[i]  = dE;
    event.Bh2Seg[i] = seg;
  }

  BH2Cluster *cl_time0 = hodoAna->GetTime0BH2Cluster();
  if( cl_time0 ){
    event.Time0Seg = cl_time0->MeanSeg()+1;
    event.deTime0  = cl_time0->DeltaE();
    event.Time0    = cl_time0->Time0();
    event.CTime0   = cl_time0->CTime0();
    time0          = cl_time0->CTime0();
  } else {
#if HodoCut
    return true;
#endif
  }

  HF1( 1, 3. );

  //////////////BH1 Analysis
  hodoAna->DecodeBH1Hits( rawData );
  Int_t nhBh1 = hodoAna->GetNHitsBH1();
  event.nhBh1 = nhBh1;
#if HodoCut
  if( nhBh1==0 ) return true;
#endif
  HF1( 1, 4 );

  for( Int_t i=0; i<nhBh1; ++i ){
    Hodo2Hit *hit = hodoAna->GetHitBH1( i );
    if( !hit ) continue;
    Double_t cmt = hit->CMeanTime();
    Double_t dE  = hit->DeltaE();
#if HodoCut
    if( dE<MinDeBH1 || MaxDeBH1<dE ) continue;
    if( btof<MinBeamToF || MaxBeamToF<btof ) continue;
#endif
    event.tBh1[i]  = cmt;
    event.deBh1[i] = dE;
  }

  Double_t btof0 = -999.;
  HodoCluster* cl_btof0 = event.Time0Seg > 0 ?
    hodoAna->GetBtof0BH1Cluster( event.CTime0 ) : nullptr;
  if( cl_btof0 ) btof0 = cl_btof0->CMeanTime() - time0;
  event.btof = btof0;

  HF1( 1, 5. );



  //////////////BC3&4 number of hit layer
  DCAna->DecodeRawHits( rawData );
#if TotCut
  DCAna->TotCutBCOut( MinTotBcOut );
#endif
  // DCAna->DriftTimeCutBC34(-10, 50);
  // DCAna->MakeBH2DCHit(event.Time0Seg-1);

  Double_t multi_BcOut = 0.;
  Double_t multiplicity[] = {0, 0};
  struct hit_info
  {
    Double_t dt;
    Double_t pos;
  };
  {
    for( Int_t layer=1; layer<=NumOfLayersBcOut; ++layer ){
      const DCHitContainer &contOut = DCAna->GetBcOutHC( layer );
      Int_t nhOut = contOut.size();
      event.nhit[layer-1] = nhOut;
      if( nhOut>0 ) event.nlayer++;
      if( layer ==  1 ) multiplicity[0] = nhOut;
      if( layer == 12 ) multiplicity[1] = nhOut;
      multi_BcOut += nhOut;
      HF1( 100*layer, nhOut );
      Int_t plane_eff = (layer-1)*3;
      Bool_t fl_valid_sig = false;
      std::vector<hit_info> hit_cont;

      for( Int_t i=0; i<nhOut; ++i ){
	DCHit *hit = contOut[i];
	Double_t wire = hit->GetWire();
	HF1( 100*layer+1, wire-0.5 );
	Int_t nhtdc = hit->GetTdcSize();
	Int_t tdc1st = -1;
	for( Int_t k=0; k<nhtdc; k++ ){
	  Int_t tdc = hit->GetTdcVal( k );
	  HF1( 100*layer+2, tdc );
	  HF1( 10000*layer+int(wire), tdc );
	  HF2( 1000*layer, tdc, wire-0.5 );
	  if( tdc > tdc1st ){
	    tdc1st = tdc;
	    fl_valid_sig = true;
	  }
	}
	HF1( 100*layer+6, tdc1st );

	if( i<MaxHits )
	  event.pos[layer-1][i] = hit->GetWirePosition();

	Int_t nhdt = hit->GetDriftTimeSize();
	int tot1st = -1;
	for( Int_t k=0; k<nhdt; k++ ){
	  Double_t dt = hit->GetDriftTime( k );
	  HF1( 100*layer+3, dt );
	  HF1( 10000*layer+1000+(Int_t)(wire), dt );

	  Double_t tot = hit->GetTot( k );
	  HF1( 100*layer+5, tot);
	  if( tot > tot1st ){
	    tot1st = tot;
	  }

	  hit_info one_hit;
	  one_hit.dt  = dt;
	  one_hit.pos = wire;
	  hit_cont.push_back( one_hit );
	}
	HF1( 100*layer+7, tot1st );

	Int_t nhdl = hit->GetDriftTimeSize();
	for( Int_t k=0; k<nhdl; k++ ){
	  Double_t dl = hit->GetDriftLength( k );
	  HF1( 100*layer+4, dl );
	}
      }
      if( fl_valid_sig ){
	++plane_eff;
      } else {
      }
      HF1(38, plane_eff);

      std::sort( hit_cont.begin(), hit_cont.end(),
		 []( const hit_info& a_info, const hit_info& b_info )->Bool_t
		 { return a_info.dt < b_info.dt; } );

      for( Int_t i=1; i<hit_cont.size(); ++i ){
	HF1( 100*layer+8, hit_cont.at(i).dt -hit_cont.at(0).dt );
	HF1( 100*layer+9, hit_cont.at(i).pos-hit_cont.at(0).pos );
      }
    }// for(layer)
  }
  HF2( 41, multiplicity[0], multiplicity[1] );

  if( multi_BcOut/Double_t(NumOfLayersBcOut) > MaxMultiHitBcOut ){
    std::cout << "#W " << __FILE__ << " L" << __LINE__ << std::endl
	      << "multi_BcOut is too many: " << multi_BcOut << std::endl;
    return true;
  }

  HF1( 1, 11. );

#if 1
  // Bc Out
  //  std::cout << "==========TrackSearch BcOut============" << std::endl;
  Bool_t status_tracking = DCAna->TrackSearchBcOut();
  DCAna->ChiSqrCutBcOut(10);

  Int_t nt=DCAna->GetNtracksBcOut();
  event.ntrack=nt;
  HF1( 10, Double_t(nt) );
  HF1( 40, status_tracking? Double_t(nt) : -1);
  for( Int_t it=0; it<nt; ++it ){
    DCLocalTrack *tp=DCAna->GetTrackBcOut(it);
    Int_t nh=tp->GetNHit();
    Double_t chisqr=tp->GetChiSquare();
    Double_t x0=tp->GetX0(), y0=tp->GetY0();
    Double_t u0=tp->GetU0(), v0=tp->GetV0();
    Double_t cost = 1./std::sqrt(1.+u0*u0+v0*v0);
    Double_t theta = std::acos(cost)*math::Rad2Deg();
    event.chisqr[it]=chisqr;
    event.x0[it]=x0;
    event.y0[it]=y0;
    event.u0[it]=u0;
    event.v0[it]=v0;

    HF1( 11, Double_t(nh) );
    HF1( 12, chisqr );
    HF1( 14, x0 ); HF1( 15, y0 );
    HF1( 16, u0 ); HF1( 17, v0 );
    HF2( 18, x0, u0 ); HF2( 19, y0, v0 );
    HF2( 20, x0, y0 );

    Double_t xtgt=tp->GetX(1800.), ytgt=tp->GetY(1800.);
    Double_t utgt=u0, vtgt=v0;
    HF1( 21, xtgt ); HF1( 22, ytgt );
    HF1( 23, utgt ); HF1( 24, vtgt );
    HF2( 25, xtgt, utgt ); HF2( 26, ytgt, vtgt );
    HF2( 27, xtgt, ytgt );

    Double_t xbac=tp->GetX(603.), ybac=tp->GetY(885.);
    Double_t ubac=u0, vbac=v0;
    HF1( 31, xbac ); HF1( 32, ybac );
    HF1( 33, ubac ); HF1( 34, vbac );
    HF2( 35, xbac, ubac ); HF2( 36, ybac, vbac );
    HF2( 37, xbac, ybac );

    Double_t Xangle = -1000*std::atan(u0);

    HF2(51, -tp->GetX(245.), Xangle);
    HF2(52, -tp->GetX(600.), Xangle);
    HF2(53, -tp->GetX(1200.), Xangle);
    HF2(54, -tp->GetX(1600.), Xangle);

    HF2(61, tp->GetX(600.), tp->GetX(245.));
    HF2(62, tp->GetX(1200.), tp->GetX(245.));
    HF2(63, tp->GetX(1600.), tp->GetX(245.));

    for( Int_t ih=0; ih<nh; ++ih ){
      DCLTrackHit *hit=tp->GetHit(ih);
      Int_t layerId=hit->GetLayer()-112;

      HF1( 13, layerId );
      Double_t wire=hit->GetWire();
      Double_t dt=hit->GetDriftTime(), dl=hit->GetDriftLength();
      HF1( 100*layerId+11, wire-0.5 );
      HF1( 100*layerId+12, dt );
      HF1( 100*layerId+13, dl );
      HF1( 10000*layerId+ 5000 +(Int_t)wire, dt);
      Double_t xcal=hit->GetXcal(), ycal=hit->GetYcal();
      Double_t pos=hit->GetLocalHitPos(), res=hit->GetResidual();
      HF1( 100*layerId+14, pos );
      HF1( 100*layerId+15, res );
      HF2( 100*layerId+16, pos, res );
      HF2( 100*layerId+17, xcal, ycal);
      //      HF1( 100000*layerId+50000+wire, res);
      Double_t wp=hit->GetWirePosition();
      Double_t sign=1.;
      if( pos-wp<0. ) sign=-1;
      HF2( 100*layerId+18, sign*dl, res );
      Double_t xlcal=hit->GetLocalCalPos();
      HF2( 100*layerId+19, dt, xlcal-wp);

      Double_t tot = hit->GetTot();
      HF1( 100*layerId+40, tot);

      if (theta>=0 && theta<15)
	HF1( 100*layerId+71, res );
      else if (theta>=15 && theta<30)
	HF1( 100*layerId+72, res );
      else if (theta>=30 && theta<45)
	HF1( 100*layerId+73, res );
      else if (theta>=45)
	HF1( 100*layerId+74, res );

      if (std::abs(dl-std::abs(xlcal-wp))<2.0) {
	HFProf( 100*layerId+20, dt, std::abs(xlcal-wp));
	HF2( 100*layerId+22, dt, std::abs(xlcal-wp));
	HFProf( 100000*layerId+3000+(Int_t)wire, xlcal-wp,dt);
	HF2( 100000*layerId+4000+(Int_t)wire, xlcal-wp,dt);
      }
    }
  }

#endif

  HF1( 1, 12. );

  return true;
}

//______________________________________________________________________________
void
EventBcOutTracking::InitializeEvent( void )
{
  event.evnum     =  0;
  event.trignhits =  0;
  event.nlayer    = -1;
  event.ntrack    = -1;
  event.nhBh2     = -1;
  event.nhBh1     = -1;

  event.pid       = -1;
  event.btof      = -999.;

  event.Time0Seg  = -1;
  event.deTime0   = -1;
  event.Time0     = -999;
  event.CTime0    = -999;

  for( Int_t it=0; it<MaxHits; it++){

    event.Bh2Seg[it] = -1;
    event.tBh2[it] = -9999.;
    event.deBh2[it] = -9999.;

    event.Bh1Seg[it] = -1;
    event.tBh1[it] = -9999.;
    event.deBh1[it] = -9999.;

    event.chisqr[it] = -1.0;
    event.x0[it] = -9999.0;
    event.y0[it] = -9999.0;
    event.u0[it] = -9999.0;
    event.v0[it] = -9999.0;
  }

  for( Int_t it=0; it<NumOfSegTrig; it++){
    event.trigpat[it] = -1;
    event.trigflag[it] = -1;
  }

  for ( Int_t it=0; it<NumOfLayersBcOut; ++it ){
    event.nhit[it] = -1;
    for( Int_t that=0; that<MaxHits; ++that ){
      event.pos[it][that] = -9999.;
    }
  }
}

//______________________________________________________________________________
Bool_t
EventBcOutTracking::ProcessingEnd( void )
{
  tree->Fill();
  return true;
}

//______________________________________________________________________________
VEvent*
ConfMan::EventAllocator( void )
{
  return new EventBcOutTracking;
}

//______________________________________________________________________________
namespace
{
  const Int_t NbinBcOutTdc   = 1000;
  const Double_t MinBcOutTdc =  200.;
  const Double_t MaxBcOutTdc = 1200.;

  const Int_t NbinBcOutDT   =  120;
  const Double_t MinBcOutDT =  -10.;
  const Double_t MaxBcOutDT =   90.;

  const Int_t NbinBcOutDL   = 80;
  const Double_t MinBcOutDL =  -0.5;
  const Double_t MaxBcOutDL =   3.5;
}
//______________________________________________________________________________
Bool_t
ConfMan:: InitializeHistograms( void )
{
  HB1( 1, "Status", 20, 0., 20. );

  //***********************Chamber
  // BC3
  for( Int_t i=1; i<=NumOfLayersBc; ++i ){
    TString title0 = Form("#Hits BC3#%2d", i);
    TString title1 = Form("Hitpat BC3#%2d", i);
    TString title2 = Form("Tdc BC3#%2d", i);
    TString title3 = Form("Drift Time BC3#%2d", i);
    TString title4 = Form("Drift Length BC3#%2d", i);
    TString title5 = Form("TOT BC3#%2d", i);
    TString title6 = Form("Tdc 1st BC3#%2d",  i);
    TString title7 = Form("TOT 1st BC3#%2d",  i);
    TString title8 = Form("Time interval from 1st hit BC3#%2d", i);
    TString title9 = Form("Position interval from 1st hit BC3#%2d", i);
    HB1( 100*i+0, title0, MaxWireBC3+1, 0., Double_t(MaxWireBC3+1) );
    HB1( 100*i+1, title1, MaxWireBC3+1, 0., Double_t(MaxWireBC3+1) );
    HB1( 100*i+2, title2, NbinBcOutTdc, MinBcOutTdc, MaxBcOutTdc );
    HB1( 100*i+3, title3, NbinBcOutDT, MinBcOutDT, MaxBcOutDT );
    HB1( 100*i+4, title4, NbinBcOutDL, MinBcOutDL, MaxBcOutDL );
    HB1( 100*i+5, title5, 500,    0, 500 );
    HB1( 100*i+6, title6, NbinBcOutTdc, MinBcOutTdc, MaxBcOutTdc );
    HB1( 100*i+7, title7, 500,  0, 500 );
    HB1( 100*i+8, title8, 72,     0, 60 );
    HB1( 100*i+9, title9, 64,   -32, 32 );
    for (Int_t wire=1; wire<=MaxWireBC3; wire++) {
      TString title10 = Form("Tdc BC3#%2d Wire#%d", i, wire);
      TString title11 = Form("Drift Time BC3#%2d Wire#%d", i, wire);
      TString title12 = Form("Drift Length BC3#%2d Wire#%d", i, wire);
      TString title15 = Form("Drift Time BC3#%2d Wire#%d [Track]", i, wire);
      HB1( 10000*i+wire, title10, NbinBcOutTdc, MinBcOutTdc, MaxBcOutTdc );
      HB1( 10000*i+1000+wire, title11, NbinBcOutDT, MinBcOutDT, MaxBcOutDT );
      HB1( 10000*i+2000+wire, title12, NbinBcOutDL, MinBcOutDL, MaxBcOutDL );
      HB1( 10000*i+5000+wire, title15, NbinBcOutDT, MinBcOutDT, MaxBcOutDT );
    }
  }

  // BC4
  for( Int_t i=1; i<=NumOfLayersBc+1; ++i ){
    TString title0 = Form("#Hits BC4#%2d", i);
    TString title1 = Form("Hitpat BC4#%2d", i);
    TString title2 = Form("Tdc BC4#%2d", i);
    TString title3 = Form("Drift Time BC4#%2d", i);
    TString title4 = Form("Drift Length BC4#%2d", i);
    TString title5 = Form("TOT BC4#%2d", i);
    TString title6 = Form("Tdc 1st BC4#%2d",  i);
    TString title7 = Form("TOT 1st BC4#%2d",  i);
    TString title8 = Form("Time interval from 1st hit BC4#%2d", i);
    TString title9 = Form("Position interval from 1st hit BC4#%2d", i);
    HB1( 100*(i+6)+0, title0, MaxWireBC4+1, 0., Double_t(MaxWireBC4+1) );
    HB1( 100*(i+6)+1, title1, MaxWireBC4+1, 0., Double_t(MaxWireBC4+1) );
    HB1( 100*(i+6)+2, title2, NbinBcOutTdc, MinBcOutTdc, MaxBcOutTdc );
    HB1( 100*(i+6)+3, title3, NbinBcOutDT, MinBcOutDT, MaxBcOutDT );
    HB1( 100*(i+6)+4, title4, NbinBcOutDL, MinBcOutDL, MaxBcOutDL );
    HB1( 100*(i+6)+5, title5, 500,    0, 500 );
    HB1( 100*(i+6)+6, title6, NbinBcOutTdc, MinBcOutTdc, MaxBcOutTdc );
    HB1( 100*(i+6)+7, title7, 500,  0, 500 );
    HB1( 100*(i+6)+8, title8, 72,     0, 60 );
    HB1( 100*(i+6)+9, title9, 64,   -32, 32 );
    for (Int_t wire=1; wire<=MaxWireBC4; wire++) {
      TString title10 = Form("Tdc BC4#%2d Wire#%d", i, wire);
      TString title11 = Form("Drift Time BC4#%2d Wire#%d", i, wire);
      TString title12 = Form("Drift Length BC4#%2d Wire#%d", i, wire);
      TString title15 = Form("Drift Time BC4#%2d Wire#%d [Track]", i, wire);
      HB1( 10000*(i+6)+wire, title10, NbinBcOutTdc, MinBcOutTdc, MaxBcOutTdc );
      HB1( 10000*(i+6)+1000+wire, title11, NbinBcOutDT, MinBcOutDT, MaxBcOutDT );
      HB1( 10000*(i+6)+2000+wire, title12, NbinBcOutDL, MinBcOutDL, MaxBcOutDL );
      HB1( 10000*(i+6)+5000+wire, title15, NbinBcOutDT, MinBcOutDT, MaxBcOutDT );
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
  HB1( 38, "Plane Eff", 36, 0, 36);

  HB2( 51, "X-X' 245 BcOut", 400, -100., 100., 120, -60, 60);
  HB2( 52, "X-X' 600 BcOut", 400, -100., 100., 120, -60, 60);
  HB2( 53, "X-X' 1200 BcOut", 400, -100., 100., 120, -60, 60);
  HB2( 54, "X-X' 1600 BcOut", 400, -100., 100., 120, -60, 60);

  HB2( 61, "X-X 600 BcOut", 400, -100., 100., 400, -100, 100);
  HB2( 62, "X-X 1200 BcOut", 400, -100., 100., 400, -100, 100);
  HB2( 63, "X-X 1600 BcOut", 400, -100., 100., 400, -100, 100);

  // Analysis status
  HB1( 40, "Tacking status", 11, -1., 10. );
  HB2( 41, "BC3X0/BC4X1", 20, 0, 20, 20, 0, 20);

  for( Int_t i=1; i<=NumOfLayersBcOut; ++i ){
    TString title11 = Form("HitPat BcOut%2d [Track]", i);
    TString title12 = Form("DriftTime BcOut%2d [Track]", i);
    TString title13 = Form("DriftLength BcOut%2d [Track]", i);
    TString title14 = Form("Position BcOut%2d", i);
    TString title15 = Form("Residual BcOut%2d", i);
    TString title16 = Form("Resid%%Pos BcOut%2d", i);
    TString title17 = Form("Y%%Xcal BcOut%2d", i);
    TString title18 = Form("Res%%dl BcOut%2d", i);
    TString title19 = Form("HitPos%%DriftTime BcOut%2d", i);
    TString title20 = Form("DriftLength%%DriftTime BcOut%2d", i);
    TString title21 = title15 + " [w/o Self]";
    TString title22 = title20;
    TString title40 = Form("TOT BcOut%2d [Track]", i);
    TString title71 = Form("Residual BcOut%2d (0<theta<15)", i);
    TString title72 = Form("Residual BcOut%2d (15<theta<30)", i);
    TString title73 = Form("Residual BcOut%2d (30<theta<45)", i);
    TString title74 = Form("Residual BcOut%2d (45<theta)", i);
    HB1( 100*i+11, title11, 64, 0., 64. );
    HB1( 100*i+12, title12, NbinBcOutDT, MinBcOutDT, MaxBcOutDT );
    HB1( 100*i+13, title13, NbinBcOutDL, MinBcOutDL, MaxBcOutDL );
    HB1( 100*i+14, title14, 100, -250., 250. );
    HB1( 100*i+15, title15, 400, -2.0, 2.0 );
    HB2( 100*i+16, title16, 250, -250., 250., 100, -1.0, 1.0 );
    HB2( 100*i+17, title17, 100, -250., 250., 100, -250., 250. );
    HB2( 100*i+18, title18, 50, 3., 3., 100, -1.0, 1.0 );
    HB2( 100*i+19, title19, 125, -5., 100., 100, -3, 3);
    HBProf( 100*i+20, title20, 100, -5., 50., 0., 12. );
    HB2( 100*i+22, title22,
	 NbinBcOutDT, MinBcOutDT, MaxBcOutDT,
	 NbinBcOutDL, MinBcOutDL, MaxBcOutDL );
    HB1( 100*i+21, title21, 200, -5.0, 5.0 );
    HB1( 100*i+40, title40, 360,    0, 300 );
    HB1( 100*i+71, title71, 200, -5.0, 5.0 );
    HB1( 100*i+72, title72, 200, -5.0, 5.0 );
    HB1( 100*i+73, title73, 200, -5.0, 5.0 );
    HB1( 100*i+74, title74, 200, -5.0, 5.0 );
    for (Int_t j=1; j<=64; j++) {
      TString title = Form("XT of Layer %2d Wire #%4d", i, j);
      HBProf( 100000*i+3000+j, title, 100, -4., 4., -5., 40. );
      HB2( 100000*i+4000+j, title, 100, -4., 4., 100, -5., 40. );
    }
  }

  ////////////////////////////////////////////
  //Tree
  HBTree( "bcout","tree of BcOutTracking" );
  tree->Branch("evnum",     &event.evnum,     "evnum/I");

  tree->Branch("trignhits", &event.trignhits, "trignhits/I");
  tree->Branch("trigpat",    event.trigpat,   Form("trigpat[%d]/I", NumOfSegTrig));
  tree->Branch("trigflag",   event.trigflag,  Form("trigflag[%d]/I", NumOfSegTrig));

  //Hodoscope
  tree->Branch("nhBh1",    &event.nhBh1,   "nhBh1/I");
  tree->Branch("tBh1",      event.tBh1,    Form("tBh1[%d]/I",   MaxHits));
  tree->Branch("deBh1",     event.deBh1,   Form("deBh1[%d]/I",  MaxHits));
  tree->Branch("Bh1Seg",    event.Bh1Seg,  Form("Bh1Seg[%d]/I", MaxHits));

  tree->Branch("nhBh2",    &event.nhBh2,   "nhBh2/I");
  tree->Branch("tBh2",      event.tBh2,    Form("tBh2[%d]/I",   MaxHits));
  tree->Branch("deBh2",     event.deBh2,   Form("deBh2[%d]/I",  MaxHits));
  tree->Branch("Bh2Seg",    event.Bh2Seg,  Form("Bh2Seg[%d]/I", MaxHits));

  tree->Branch("Time0Seg", &event.Time0Seg,  "Time0Seg/D");
  tree->Branch("deTime0",  &event.deTime0,   "deTime0/D");
  tree->Branch("Time0",    &event.Time0,     "Time0/D");
  tree->Branch("CTime0",   &event.CTime0,    "CTime0/D");

  tree->Branch("pid",      &event.pid,      "pid/I");
  tree->Branch("btof",     &event.btof,     "btof/D");

  tree->Branch("nhit",     &event.nhit,     Form("nhit[%d]/I", NumOfLayersBcOut ) );

  tree->Branch("ntrack",   &event.ntrack,   "ntrack/I");
  tree->Branch("chisqr",    event.chisqr,   "chisqr[ntrack]/D");
  tree->Branch("x0",        event.x0,       "x0[ntrack]/D");
  tree->Branch("y0",        event.y0,       "y0[ntrack]/D");
  tree->Branch("u0",        event.u0,       "u0[ntrack]/D");
  tree->Branch("v0",        event.v0,       "v0[ntrack]/D");

  TString layer_name[NumOfLayersBcOut] =
    { "bc3x0", "bc3x1", "bc3v0", "bc3v1", "bc3u0", "bc3u1",
      "bc4u0", "bc4u1", "bc4v0", "bc4v1", "bc4x0", "bc4x1" };
  for( Int_t i=0; i<NumOfLayersBcOut; ++i ){
    tree->Branch( Form("%s_pos", layer_name[i].Data() ), event.pos[i],
		   Form("%s_pos[%d]/D", layer_name[i].Data(), MaxHits ) );
  }

  HPrint();

  return true;
}

//______________________________________________________________________________
Bool_t
ConfMan::InitializeParameterFiles( void )
{
  return
    ( InitializeParameter<DCGeomMan>("DCGEO")        &&
      InitializeParameter<DCDriftParamMan>("DCDRFT") &&
      InitializeParameter<DCTdcCalibMan>("DCTDC")    &&
      InitializeParameter<HodoParamMan>("HDPRM")     &&
      InitializeParameter<HodoPHCMan>("HDPHC")       &&
      InitializeParameter<UserParamMan>("USER")      );
}

//______________________________________________________________________________
Bool_t
ConfMan::FinalizeProcess( void )
{
  return true;
}
