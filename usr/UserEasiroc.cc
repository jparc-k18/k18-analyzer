// -*- C++ -*-

#include "VEvent.hh"

#include <iostream>
#include <sstream>
#include <cmath>
#include <vector>
#include <utility>
#include <algorithm>

#include "BH2Cluster.hh"
#include "BH2Hit.hh"
#include "ConfMan.hh"
#include "DCGeomMan.hh"
#include "DetectorID.hh"
#include "RMAnalyzer.hh"
#include "FiberCluster.hh"
#include "FiberHit.hh"
#include "RootHelper.hh"
#include "HodoHit.hh"
#include "HodoAnalyzer.hh"
#include "HodoCluster.hh"
#include "HodoParamMan.hh"
#include "HodoPHCMan.hh"
#include "HodoRawHit.hh"
#include "S2sLib.hh"
#include "RawData.hh"
#include "UnpackerManager.hh"

#define HodoCut    0 // with BH1/BH2
#define TIME_CUT   1 // in cluster analysis
#define DE_CUT     1 // in cluster analysis for rc
#define FHitBranch 0 // make FiberHit branches (becomes heavy)
#define RawHitRCBranch 1 //make RC RawHit branches (becomes heavy)
#define ClusterHitRCBranch 0 //make RCCluster branches

namespace
{
enum EUorD { kU, kD, kUorD };
using namespace root;
using hddaq::unpacker::GUnpacker;
const auto qnan = TMath::QuietNaN();
const auto& gUnpacker = GUnpacker::get_instance();
const auto& gUser = UserParamMan::GetInstance();
}

//_____________________________________________________________________________
struct Event
{
  Int_t evnum;
  Int_t trigpat[NumOfSegTrig];
  Int_t trigflag[NumOfSegTrig];

  // BH1
  Int_t    bh1_nhits;
  Int_t	bh1seg[MaxHits];
  Int_t    bh1_ncl;
  Int_t    bh1_clsize[NumOfSegBH1];
  Double_t bh1_clseg[NumOfSegBH1];

  // BFT
  Int_t    bft_nhits;
  Int_t    bft_unhits;
  Int_t    bft_dnhits;
  Int_t    bft_uhitpat[NumOfSegBFT];
  Int_t    bft_dhitpat[NumOfSegBFT];
  Double_t bft_utdc[NumOfSegBFT][MaxDepth];
  Double_t bft_dtdc[NumOfSegBFT][MaxDepth];
  Double_t bft_utrailing[NumOfSegBFT][MaxDepth];
  Double_t bft_dtrailing[NumOfSegBFT][MaxDepth];
  Double_t bft_utot[NumOfSegBFT][MaxDepth];
  Double_t bft_dtot[NumOfSegBFT][MaxDepth];
  Int_t    bft_udepth[NumOfSegBFT];
  Int_t    bft_ddepth[NumOfSegBFT];
  Int_t    bft_ncl;
  Int_t    bft_clsize[NumOfSegBFT];
  Double_t bft_ctime[NumOfSegBFT];
  Double_t bft_ctot[NumOfSegBFT];
  Double_t bft_clpos[NumOfSegBFT];
  Double_t bft_clseg[NumOfSegBFT];

  // RC raw
  Int_t    rc_nhits[NumOfPlaneRC];
  Int_t    rc_hitpat[NumOfPlaneRC][NumOfSegRC];
  Double_t rc_tdc[NumOfPlaneRC][NumOfSegRC][kUorD][MaxDepth];
  Double_t rc_adc_high[NumOfPlaneRC][NumOfSegRC][kUorD];
  Double_t rc_adc_low[NumOfPlaneRC][NumOfSegRC][kUorD];
  Double_t rc_ltime[NumOfPlaneRC][NumOfSegRC][kUorD][MaxDepth];
  Double_t rc_ttime[NumOfPlaneRC][NumOfSegRC][kUorD][MaxDepth];
  Double_t rc_tot[NumOfPlaneRC][NumOfSegRC][kUorD][MaxDepth];

  // RC normalized
  Double_t rc_mt[NumOfPlaneRC][NumOfSegRC][MaxDepth];
  Double_t rc_cmt[NumOfPlaneRC][NumOfSegRC][MaxDepth];
  Double_t rc_mtot[NumOfPlaneRC][NumOfSegRC][MaxDepth];
  Double_t rc_de_high[NumOfPlaneRC][NumOfSegRC];
  Double_t rc_de_low[NumOfPlaneRC][NumOfSegRC];

  //RC cluster
  Int_t    rc_ncl;
  Double_t rc_desum;
  Int_t    rc_clsize[MaxCluster];
  Double_t rc_clseg[MaxCluster];
  Double_t rc_cltot[MaxCluster];
  Double_t rc_clde[MaxCluster];
  Double_t rc_cltime[MaxCluster];
  Int_t    rc_clplane[MaxCluster];
  Double_t rc_clpos[MaxCluster];

  void clear();
};

//_____________________________________________________________________________
void
Event::clear()
{
  evnum      = 0;
  bft_nhits  = 0;
  bft_unhits = 0;
  bft_dnhits = 0;
  bft_ncl    = 0;

  for(Int_t it=0; it<NumOfSegTrig; it++){
    trigpat[it]  = -1;
    trigflag[it] = -1;
  }

  bh1_nhits = 0;
  for(Int_t it=0; it<MaxHits; it++){
    bh1seg[it]  = -1;
  }
  bh1_ncl    = 0;
  for(Int_t it=0; it<NumOfSegBH1; it++){
    bh1_clsize[it] = 0;
    bh1_clseg[it] = qnan;
  }

  for(Int_t it=0; it<NumOfSegBFT; it++){
    bft_uhitpat[it] = -1;
    bft_dhitpat[it] = -1;
    bft_udepth[it]  = 0;
    bft_ddepth[it]  = 0;
    for(Int_t that=0; that<MaxDepth; that++){
      bft_utdc[it][that] = qnan;
      bft_dtdc[it][that] = qnan;
      bft_utrailing[it][that] = qnan;
      bft_dtrailing[it][that] = qnan;
      bft_utot[it][that] = qnan;
      bft_dtot[it][that] = qnan;
    }
    bft_clsize[it] = 0;
    bft_ctime[it]  = qnan;
    bft_ctot[it]   = qnan;
    bft_clpos[it]  = qnan;
    bft_clseg[it]  = qnan;
  }

  //rc raw
  for(Int_t p=0; p<NumOfPlaneRC; p++){
    rc_nhits[p] = 0;
    for(Int_t seg=0; seg<NumOfSegRC; seg++){
      rc_hitpat[p][seg] = -1;
      for(Int_t ud=0; ud<kUorD; ud++){
        rc_adc_high[p][seg][ud] = qnan;
        rc_adc_low[p][seg][ud] = qnan;
        for(Int_t i=0; i<MaxDepth; i++){
          rc_tdc[p][seg][ud][i] = qnan;
          rc_ltime[p][seg][ud][i] = qnan;
          rc_ttime[p][seg][ud][i] = qnan;
          rc_tot[p][seg][ud][i] = qnan;
        }
      }
      rc_de_high[p][seg] = qnan;
      rc_de_low[p][seg] = qnan;
      for(Int_t i=0; i<MaxDepth; i++){
        rc_mt[p][seg][i] = qnan;
        rc_cmt[p][seg][i] = qnan;
        rc_mtot[p][seg][i] = qnan;
      }
    }
  }

  //rc cluster
  rc_ncl   = 0;
  rc_desum = qnan;
  for(Int_t it=0; it<MaxCluster; it++){
	rc_clsize[it]  = 0;
	rc_clseg[it]   = qnan;
	rc_cltot[it]   = qnan;
	rc_clde[it]    = qnan;
	rc_cltime[it]  = qnan;
	rc_clplane[it] = qnan;
	rc_clpos[it]   = qnan;
  }
}

//_____________________________________________________________________________
namespace root
{
Event  event;
TH1   *h[MaxHist];
TTree *tree;
enum eDetHid
{
  BFTHid  =  10000,
  RCHid  = 100000,
};
}

//_____________________________________________________________________________
Bool_t
ProcessingBegin()
{
  event.clear();
  return true;
}

//_____________________________________________________________________________
Bool_t
ProcessingNormal()
{
#if HodoCut
  static const Double_t MinDeBH2   = gUser.GetParameter("DeBH2", 0);
  static const Double_t MaxDeBH2   = gUser.GetParameter("DeBH2", 1);
  static const Double_t MinDeBH1   = gUser.GetParameter("DeBH1", 0);
  static const Double_t MaxDeBH1   = gUser.GetParameter("DeBH1", 1);
  static const Double_t MinBeamToF = gUser.GetParameter("BTOF",  0);
  static const Double_t MaxBeamToF = gUser.GetParameter("BTOF",  1);
#endif
  static const Double_t MinTdcBFT  = gUser.GetParameter("TdcBFT",  0);
  static const Double_t MaxTdcBFT  = gUser.GetParameter("TdcBFT",  1);
#if TIME_CUT
  static const Double_t MinTimeBFT = gUser.GetParameter("TimeBFT", 0);
  static const Double_t MaxTimeBFT = gUser.GetParameter("TimeBFT", 1);
  static const Double_t MinTimeRC = gUser.GetParameter("TimeRC", 0);
  static const Double_t MaxTimeRC = gUser.GetParameter("TimeRC", 1);
#endif
#if DE_CUT
  static const Double_t MinDeRC   = gUser.GetParameter("DeRC", 0);
  static const Double_t MaxDeRC   = gUser.GetParameter("DeRC", 1);
#endif

  RawData rawData;
  rawData.DecodeHits("TFlag");
  rawData.DecodeHits("BH1");
  rawData.DecodeHits("BH2");
  rawData.DecodeHits("BFT");
  rawData.DecodeHits("RC");
  HodoAnalyzer hodoAna(rawData);

  event.evnum = gUnpacker.get_event_number();

  HF1(1, 0);

  ///// Trigger Flag
  std::bitset<NumOfSegTrig> trigger_flag;
  {
    for(const auto& hit: rawData.GetHodoRawHC("TFlag")){
      Int_t seg = hit->SegmentId();
      Int_t tdc = hit->GetTdc();
      if(tdc > 0){
	event.trigpat[trigger_flag.count()] = seg;
	event.trigflag[seg] = tdc;
        trigger_flag.set(seg);
	HF1(10, seg);
	HF1(10+seg, tdc);
      }
    }
  }

  if(trigger_flag[trigger::kSpillOnEnd] || trigger_flag[trigger::kSpillOffEnd])
    return true;

  HF1(1, 1);

  ////////// BH2 time 0
  hodoAna.DecodeHits<BH2Hit>("BH2");
  Int_t nhBh2 = hodoAna.GetNHits("BH2");
#if HodoCut
  if(nhBh2==0) return true;
#endif
  HF1(1, 2);
  Double_t time0 = qnan;
  ////////// BH2 Analysis
  for(Int_t i=0; i<nhBh2; ++i){
    const auto& hit = hodoAna.GetHit("BH2", i);
    if(!hit) continue;
    Double_t cmt = hit->CMeanTime();
    Double_t ct0 = hit->CTime0();
    Double_t min_time = qnan;
#if HodoCut
    Double_t dE  = hit->DeltaE();
    if(dE<MinDeBH2 || MaxDeBH2<dE)
      continue;
#endif
    if(std::abs(cmt)<std::abs(min_time)){
      min_time = cmt;
      time0    = ct0;
    }
  }

  HF1(1, 3);

  ////////// BH1 Analysis
  hodoAna.DecodeHits("BH1");
  Int_t nhBh1 = hodoAna.GetNHits("BH1");
#if HodoCut
  if(nhBh1==0) return true;
#endif
  HF1(1, 4);
  Double_t btof0 = qnan;
  Int_t bh1nhits = 0;
  for(Int_t i=0; i<nhBh1; ++i){
    const auto& hit = hodoAna.GetHit("BH1", i);
    if(!hit) continue;
    Double_t cmt  = hit->CMeanTime();
    Double_t btof = cmt - time0;
    Int_t seg = hit->SegmentId()+1;
    event.bh1seg[bh1nhits++] = seg;
#if HodoCut
    Double_t dE   = hit->DeltaE();
    if(dE<MinDeBH1 || MaxDeBH1<dE) continue;
    if(btof<MinBeamToF || MaxBeamToF<btof) continue;
#endif
    if(std::abs(btof)<std::abs(btof0)){
      btof0 = btof;
    }
  }
  event.bh1_nhits  = bh1nhits;

  Int_t nclbh1 = hodoAna.GetNClusters("BH1");
  if(nclbh1 > NumOfSegBH1){
    // std::cout << "#W BFT too much number of clusters" << std::endl;
    nclbh1 = NumOfSegBH1;
  }
  event.bh1_ncl = nclbh1;
  for(Int_t i=0; i<nclbh1; ++i){
    const auto& cl = hodoAna.GetCluster("BH1", i);
    if(!cl) continue;
    Double_t clsize = cl->ClusterSize();
    Double_t seg  = cl->MeanSeg();
    event.bh1_clseg[i] = seg;
    event.bh1_clsize[i] = clsize;
  }

  HF1(1, 5);

  HF1(1, 6);

  ////////// BFT
  hodoAna.DecodeHits<FiberHit>("BFT");
  {
    const auto& U = HodoRawHit::kUp;
    Int_t nh = hodoAna.GetNHits("BFT");
    Int_t unhits = 0;
    Int_t dnhits = 0;
    for(Int_t i=0; i<nh; ++i){
      const auto& hit = hodoAna.GetHit<FiberHit>("BFT", i);
      const auto& rhit = hit->GetRawHit();
      Int_t plane = hit->PlaneId();
      Int_t seg = hit->SegmentId();
      // Raw
      Int_t mh_l = rhit->GetSizeTdcLeading();
      if(plane == 0) event.bft_udepth[seg] = mh_l;
      for(Int_t j=0; j<mh_l; ++j){
        Double_t leading = rhit->GetTdcLeading(U, j);
        HF1(BFTHid +plane + 6, leading);
        HF2(BFTHid +plane + 10, seg, leading);
        HF1(BFTHid +1000*(1+plane)+seg+1, leading);
        if(plane == 0) event.bft_utdc[seg][j] = leading;
        if(plane == 1) event.bft_dtdc[seg][j] = leading;
        if(MinTdcBFT<leading && leading<MaxTdcBFT){
	  HF1(BFTHid +plane+4, seg+0.5);
	  if(plane==0) event.bft_uhitpat[unhits++] = seg;
	  if(plane==1) event.bft_dhitpat[dnhits++] = seg;
        }
      }
      Int_t mh_t = rhit->GetSizeTdcTrailing();
      if(plane == 1) event.bft_ddepth[seg] = mh_l;
      for(Int_t j=0; j<mh_t; ++j){
        Double_t trailing = rhit->GetTdcTrailing(U, j);
        if(plane == 0) event.bft_utrailing[seg][j] = trailing;
        if(plane == 1) event.bft_dtrailing[seg][j] = trailing;
      }
      // Normalized
      Int_t mh = hit->GetEntries();
      for(Int_t j=0; j<mh; ++j){
        Double_t time  = hit->MeanTime(j);
        Double_t ctime = hit->CMeanTime(j);
        Double_t tot   = hit->TOT(U, j);
        HF1(BFTHid +plane+8, tot);
        HF2(BFTHid +plane+12, seg, tot);
        HF1(BFTHid +plane+21, time);
        HF1(BFTHid +plane+31, ctime);
        HF1(BFTHid +1000*(plane+3)+seg+1, tot);
        if(plane == 0) event.bft_utot[seg][j] = tot;
        if(plane == 1) event.bft_dtot[seg][j] = tot;
        if(-10.<time && time<10.){
          HF2(BFTHid +plane+23, tot, time);
          HF2(BFTHid +plane+33, tot, ctime);
          HF2(BFTHid +1000*(plane+5)+seg+1, tot, time);
          HF2(BFTHid +1000*(plane+7)+seg+1, tot, ctime);
        }
      }
    }
    HF1(BFTHid +1, unhits);
    HF1(BFTHid +2, dnhits);
    HF1(BFTHid +3, unhits + dnhits);
    event.bft_unhits = unhits;
    event.bft_dnhits = dnhits;
    event.bft_nhits  = unhits + dnhits;

    // Fiber Cluster
#if TIME_CUT
    hodoAna.TimeCut("BFT", MinTimeBFT, MaxTimeBFT);
#endif
    Int_t ncl = hodoAna.GetNClusters("BFT");
    if(ncl > NumOfSegBFT){
      // std::cout << "#W BFT too much number of clusters" << std::endl;
      ncl = NumOfSegBFT;
    }
    event.bft_ncl = ncl;
    HF1(BFTHid +101, ncl);
    for(Int_t i=0; i<ncl; ++i){
      const auto& cl = hodoAna.GetCluster("BFT", i);
      if(!cl) continue;
      Double_t clsize = cl->ClusterSize();
      Double_t ctime  = cl->CMeanTime();
      Double_t ctot   = cl->TOT();
      Double_t pos    = cl->MeanPosition();
      Double_t seg    = cl->MeanSeg();
      event.bft_clsize[i] = clsize;
      event.bft_ctime[i]  = ctime;
      event.bft_ctot[i]   = ctot;
      event.bft_clpos[i]  = pos;
      event.bft_clseg[i]  = seg;
      HF1(BFTHid +102, clsize);
      HF1(BFTHid +103, ctime);
      HF1(BFTHid +104, ctot);
      HF2(BFTHid +105, ctot, ctime);
      HF1(BFTHid +106, pos);
    }
  }

  ////////// RC
  for(const auto& hit: rawData.GetHodoRawHC("RC")){
    // hit->Print();
    Int_t plane = hit->PlaneId();
    Int_t seg = hit->SegmentId();
    Bool_t isTDC = kFALSE;
    for(Int_t ud=0; ud<kUorD; ++ud){
      auto adc_high = hit->GetAdcHigh(ud);
      auto adc_low = hit->GetAdcLow(ud);
      event.rc_adc_high[plane][seg][ud] = adc_high;
      event.rc_adc_low[plane][seg][ud] = adc_low;
      HF1(RCHid+plane*1000+7+ud, adc_high);
      HF1(RCHid+plane*1000+9+ud, adc_low);
      HF2(RCHid+plane*1000+15+ud, seg, adc_high);
      HF2(RCHid+plane*1000+17+ud, seg, adc_low);
      for(Int_t i=0, n=hit->GetSizeTdcLeading(ud); i<n; ++i){
        auto tdc = hit->GetTdc(ud, i);
	isTDC = kTRUE;
	auto tra = hit->GetTdcTrailing(ud, i);
        event.rc_tdc[plane][seg][ud][i] = tdc;
        HF1(RCHid+plane*1000+3+ud, tdc);
        HF2(RCHid+plane*1000+11+ud, seg, tdc);
	HF2(RCHid+plane*1000+5+ud, seg, n);

        HF2(RCHid+plane*1000+100+13+ud, seg, tra);
	HF2(RCHid+plane*1000+100+15+ud, seg, tdc-tra);	
        // HF1(RCHid+plane*1000+seg+100+ud*100, tdc);
      }
      if (isTDC){
	HF2(RCHid+plane*1000+81+ud, seg, adc_high);
	HF2(RCHid+plane*1000+83+ud, seg, adc_low);
      }
    }
  }
  hodoAna.DecodeHits<FiberHit>("RC");
  for(Int_t i=0, n=hodoAna.GetNHits("RC"); i<n; ++i){
    const auto& hit = hodoAna.GetHit<FiberHit>("RC", i);
    // hit->Print();
    // const auto& rhit = hit->GetRawHit();
    // rhit->Print();
    Int_t plane = hit->PlaneId();
    Int_t seg = hit->SegmentId();
    event.rc_hitpat[plane][event.rc_nhits[plane]++] = seg;
    //HF1(RCHid+plane*1000+2, seg);
    Int_t m = hit->GetEntries();
    for(Int_t j=0; j<m; ++j){
      auto mt = hit->MeanTime(j);
      auto cmt = hit->CMeanTime(j);
      auto mtot = hit->MeanTOT(j);
      event.rc_mt[plane][seg][j] = mt;
      event.rc_cmt[plane][seg][j] = cmt;
      event.rc_mtot[plane][seg][j] = mtot;
      HF1(RCHid+plane*1000+21, mt);
      HF1(RCHid+plane*1000+22, cmt);
      HF1(RCHid+plane*1000+23, mtot);
      HF2(RCHid+plane*1000+31, seg, mt);
      HF2(RCHid+plane*1000+32, seg, cmt);
      HF2(RCHid+plane*1000+33, seg, mtot);
      for(Int_t ud=0; ud<kUorD; ++ud){
        auto tot = hit->TOT(ud, j);
	auto ltime = hit->GetTimeLeading(ud, j);
	auto ttime = hit->GetTimeTrailing(ud, j);
	event.rc_tot[plane][seg][ud][j]   = tot;
	event.rc_ltime[plane][seg][ud][j] = ltime;
	event.rc_ttime[plane][seg][ud][j] = ttime;
        // HF1(RCHid+plane*1000+5+ud, tot);
        HF2(RCHid+plane*1000+13+ud, seg, tot);
	HF2(RCHid+plane*1000+56+ud, seg, ttime);
	HF2(RCHid+plane*1000+58+ud, seg, ltime);
	HF1(RCHid+plane*1000+60+ud, ltime);
        // HF1(RCHid+plane*1000+seg+300+ud*100, tot);
      }
    }
    auto de_high = hit->DeltaEHighGain();
    auto de_low = hit->DeltaELowGain();
    event.rc_de_high[plane][seg] = de_high;
    event.rc_de_low[plane][seg] = de_low;
    HF1(RCHid+plane*1000+24, de_high);
    HF1(RCHid+plane*1000+25, de_low);
    HF2(RCHid+plane*1000+34, seg, de_high);
    HF2(RCHid+plane*1000+35, seg, de_low);
  }
  for(Int_t plane=0; plane<NumOfPlaneRC; ++plane){
    //HF1(RCHid+plane*1000+1, event.rc_nhits[plane]);
  }

  //RC cluster
#if TIME_CUT
  hodoAna.TimeCut("RC", MinTimeRC, MaxTimeRC);
#endif
#if DE_CUT
  hodoAna.DeCut("RC", MinDeRC, MaxDeRC);
#endif
  {
    int nclrc = hodoAna.GetNClusters("RC");
    event.rc_ncl = nclrc;
    double desum = 0;
    if (nclrc > MaxCluster) nclrc = MaxCluster;
    for(Int_t i=0; i<nclrc; ++i){
      const auto& cl = hodoAna.GetCluster("RC", i);
      if(!cl) continue;
      Int_t    plane  = cl->PlaneId();
      Double_t clsize = cl->ClusterSize();
      Double_t time   = cl->MeanTime();
      Double_t tot    = cl->TOT();
      Double_t pos    = cl->MeanPosition();
      Double_t seg    = cl->MeanSeg();
      Double_t de     = cl->DeltaE();
      desum = desum + de;
      event.rc_clsize[i]  = clsize;
      event.rc_clseg[i]   = seg;
      event.rc_cltot[i]   = tot;
      event.rc_clde[i]    = de;
      event.rc_cltime[i]  = time;
      event.rc_clplane[i] = plane;
      event.rc_clpos[i]   = pos;
      HF1(RCHid + 102, clsize);
      HF1(RCHid + 103, time);
      HF1(RCHid + 104, tot);
      HF2(RCHid + 105, time, tot);
      HF1(RCHid + 106, seg);
    }
    HF1(RCHid + 101, nclrc);
    HF1(RCHid + 107, desum);
    event.rc_desum = desum;
  }


  //rc_analysis
  // int multiplicity_pair[18] = { 0 };
  // for(int ud=0; ud<kUorD; ud++){
  //   for(int plane=0; plane<NumOfPlaneRC; plane++){
  //     std::vector<std::pair<int,int>> adc_seg_pair;
  //     int multiplicity = 0;
  //     for(int s=0; s<NumOfSegRCarr.at(plane); s++){
  // 	//--------------------------------------------------------
  // 	for(int depth=0; depth<MaxDepth; depth++){
  // 	  double ltime = event.rc_ltime[plane][s][ud][depth];
  // 	  double adc   = event.rc_adc_high[plane][s][ud];
  // 	  double tot   = event.rc_tot[plane][s][ud][depth];
  // 	  double mt    = event.rc_mt[plane][s][depth];
  // 	  double de_high = event.rc_de_high[plane][s];
  // 	  bool Timecut = ( MinTimeRC<ltime && ltime<MaxTimeRC);
  // 	  bool MeanTimecut = ( MinTimeRC<mt && mt<MaxTimeRC);
  // 	  bool Decut   = ( 0.2<de_high );
  // 	  if( Timecut )  HF2(RCHid+plane*1000+62+ud, adc, tot);
  // 	  if (MeanTimecut ) HF2(RCHid+plane*1000+64+ud, de_high, mt);
  // 	  if( MeanTimecut && Decut){
  // 	    multiplicity++;
  // 	    HF1(RCHid+plane*1000+2, s);
  // 	  }
  // 	}
  // 	//--------------------------------------------------------
  // 	if(std::isfinite(event.rc_tdc[plane][s][ud][0])){
  // 	  double adc = event.rc_adc_high[plane][s][ud];
  // 	  if(1000<adc) adc_seg_pair.push_back( {adc, s} );
  // 	}
  //     }//for seg
  //     if(ud==0) multiplicity_pair[plane/2] += multiplicity;
  //     if( ((plane%2)==1) && (ud==0) )  HF1(RCHid+(plane/2)*1000+1, multiplicity_pair[plane/2]);
  //     if(!adc_seg_pair.empty()){
  // 	std::sort(adc_seg_pair.rbegin(), adc_seg_pair.rend());
  // 	int adcmax = adc_seg_pair.at(0).first;
  // 	int adcmax_seg = adc_seg_pair.at(0).second;
  // 	int tdc_adcmax_seg = event.rc_tdc[plane][adcmax_seg][ud][0];
  // 	HF2(RCHid+plane*1000+50+ud, adcmax_seg, tdc_adcmax_seg);
  // 	HF1(RCHid+plane*1000+52+ud, adcmax_seg);
  // 	HF2(RCHid+plane*1000+54+ud, adcmax_seg, adcmax);
  // 	HF1(RCHid+plane*1000+56+ud, adcmax);
  //     }
  //   }//for plane
  // }//for ud

  return true;
}

//_____________________________________________________________________________
Bool_t
ProcessingEnd()
{
  tree->Fill();
  return true;
}

//_____________________________________________________________________________
Bool_t
ConfMan::InitializeHistograms()
{
  const Int_t    NbinAdc = 4000;
  const Double_t MinAdc  =    0.;
  const Double_t MaxAdc  = 4000.;

  const Int_t    NbinTdc = 1000;
  const Double_t MinTdc  =    0.;
  const Double_t MaxTdc  = 1000.;

  const Int_t    NbinTot =  170;
  const Double_t MinTot  =  -10.;
  const Double_t MaxTot  =  160.;

  const Int_t    NbinTime = 60;
  const Double_t MinTime  = -30.;
  const Double_t MaxTime  =  30.;

  const Int_t    NbinDe = 1000;
  const Double_t MinDe  =  0.;
  const Double_t MaxDe  = 10.;

  HB1( 1, "Status",  20,   0., 20.);
  HB1(10, "Trigger HitPat", NumOfSegTrig, 0., Double_t(NumOfSegTrig));
  for(Int_t i=0; i<NumOfSegTrig; ++i){
    HB1(10+i+1, Form("Trigger Flag %d", i+1), 0x1000, 0, 0x1000);
  }

  //BFT
  HB1(BFTHid + 1, "BFT Nhits U",   NumOfSegBFT, 0., (Double_t)NumOfSegBFT);
  HB1(BFTHid + 2, "BFT Nhits D",   NumOfSegBFT, 0., (Double_t)NumOfSegBFT);
  HB1(BFTHid + 3, "BFT Nhits",     NumOfSegBFT, 0., (Double_t)NumOfSegBFT);
  HB1(BFTHid + 4, "BFT Hitpat U",  NumOfSegBFT, 0., (Double_t)NumOfSegBFT);
  HB1(BFTHid + 5, "BFT Hitpat D",  NumOfSegBFT, 0., (Double_t)NumOfSegBFT);
  HB1(BFTHid + 6, "BFT Tdc U",      NbinTdc, MinTdc, MaxTdc);
  HB1(BFTHid + 7, "BFT Tdc D",      NbinTdc, MinTdc, MaxTdc);
  HB1(BFTHid + 8, "BFT Tot U",      NbinTot, MinTot, MaxTot);
  HB1(BFTHid + 9, "BFT Tot D",      NbinTot, MinTot, MaxTot);
  HB2(BFTHid +10, "BFT Tdc U%Seg",
      NumOfSegBFT, 0., (Double_t)NumOfSegBFT, NbinTdc, MinTdc, MaxTdc);
  HB2(BFTHid +11, "BFT Tdc D%Seg",
      NumOfSegBFT, 0., (Double_t)NumOfSegBFT, NbinTdc, MinTdc, MaxTdc);
  HB2(BFTHid +12, "BFT Tot U%Seg",
      NumOfSegBFT, 0., (Double_t)NumOfSegBFT, NbinTot, MinTot, MaxTot);
  HB2(BFTHid +13, "BFT Tot D%Seg",
      NumOfSegBFT, 0., (Double_t)NumOfSegBFT, NbinTot, MinTot, MaxTot);
  for(Int_t i=0; i<NumOfSegBFT; i++){
    HB1(BFTHid +1000+i+1, Form("BFT Tdc U-%d", i+1), NbinTdc, MinTdc, MaxTdc);
    HB1(BFTHid +2000+i+1, Form("BFT Tdc D-%d", i+1), NbinTdc, MinTdc, MaxTdc);
    HB1(BFTHid +3000+i+1, Form("BFT Tot U-%d", i+1), NbinTot, MinTot, MaxTot);
    HB1(BFTHid +4000+i+1, Form("BFT Tot D-%d", i+1), NbinTot, MinTot, MaxTot);
    HB2(BFTHid +5000+i+1, Form("BFT Time%%Tot U-%d", i+1),
        NbinTot, MinTot, MaxTot, NbinTime, MinTime, MaxTime);
    HB2(BFTHid +6000+i+1, Form("BFT Time%%Tot D-%d", i+1),
        NbinTot, MinTot, MaxTot, NbinTime, MinTime, MaxTime);
    HB2(BFTHid +7000+i+1, Form("BFT CTime%%Tot U-%d", i+1),
        NbinTot, MinTot, MaxTot, NbinTime, MinTime, MaxTime);
    HB2(BFTHid +8000+i+1, Form("BFT CTime%%Tot D-%d", i+1),
        NbinTot, MinTot, MaxTot, NbinTime, MinTime, MaxTime);
  }
  HB1(BFTHid +21, "BFT Time U",     NbinTime, MinTime, MaxTime);
  HB1(BFTHid +22, "BFT Time D",     NbinTime, MinTime, MaxTime);
  HB2(BFTHid +23, "BFT Time%Tot U", NbinTot, MinTot, MaxTot, NbinTime, MinTime, MaxTime);
  HB2(BFTHid +24, "BFT Time%Tot D", NbinTot, MinTot, MaxTot, NbinTime, MinTime, MaxTime);
  HB1(BFTHid +31, "BFT CTime U",     NbinTime, MinTime, MaxTime);
  HB1(BFTHid +32, "BFT CTime D",     NbinTime, MinTime, MaxTime);
  HB2(BFTHid +33, "BFT CTime%Tot U", NbinTot, MinTot, MaxTot, NbinTime, MinTime, MaxTime);
  HB2(BFTHid +34, "BFT CTime%Tot D", NbinTot, MinTot, MaxTot, NbinTime, MinTime, MaxTime);

  HB1(BFTHid +101, "BFT NCluster", 100, 0, 100);
  HB1(BFTHid +102, "BFT Cluster Size", 5, 0, 5);
  HB1(BFTHid +103, "BFT CTime (Cluster)", 100., -20., 30.);
  HB1(BFTHid +104, "BFT Tot (Cluster)", NbinTot, MinTot, MaxTot);
  HB2(BFTHid +105, "BFT CTime%Tot (Cluster)",
      NbinTot, MinTot, MaxTot, NbinTime, MinTime, MaxTime);
  HB1(BFTHid +106, "BFT Cluster Position",
      NumOfSegBFT, -0.5*(Double_t)NumOfSegBFT, 0.5*(Double_t)NumOfSegBFT);

  //RC
  for(Int_t plane=0; plane<NumOfPlaneRC; ++plane){
    HB1(RCHid+plane*1000+1, Form("RC Nhits Plane#%d", plane), NumOfSegRC, 0., NumOfSegRC);
    HB1(RCHid+plane*1000+2, Form("RC Hitpat Plane#%d", plane), NumOfSegRC, 0., NumOfSegRC);
    HB1(RCHid+plane*1000+21, Form("RC MeanTime Plane#%d", plane), NbinTime, MinTime, MaxTime);
    HB1(RCHid+plane*1000+22, Form("RC CMeanTime Plane#%d", plane), NbinTime, MinTime, MaxTime);
    HB1(RCHid+plane*1000+23, Form("RC MeanTot Plane#%d", plane), NbinTot, MinTot, MaxTot);
    HB1(RCHid+plane*1000+24, Form("RC DeltaE HighGain Plane#%d", plane), NbinDe, MinDe, MaxDe);
    HB1(RCHid+plane*1000+25, Form("RC DeltaE LowGain Plane#%d", plane), NbinDe, MinDe, MaxDe);
    HB2(RCHid+plane*1000+31, Form("RC MeanTime Plane#%d", plane),
        NumOfSegRC, 0., NumOfSegRC, NbinTime, MinTime, MaxTime);
    HB2(RCHid+plane*1000+32, Form("RC CMeanTime Plane#%d", plane),
        NumOfSegRC, 0., NumOfSegRC, NbinTime, MinTime, MaxTime);
    HB2(RCHid+plane*1000+33, Form("RC MeanTot Plane#%d", plane),
        NumOfSegRC, 0., NumOfSegRC, NbinTot, MinTot, MaxTot);
    HB2(RCHid+plane*1000+34, Form("RC DeltaE HighGain Plane#%d", plane),
        NumOfSegRC, 0., NumOfSegRC, NbinDe, MinDe, MaxDe);
    HB2(RCHid+plane*1000+35, Form("RC DeltaE LowGain Plane#%d", plane),
        NumOfSegRC, 0., NumOfSegRC, NbinDe, MinDe, MaxDe);
    for(Int_t ud=0; ud<kUorD; ++ud){
      const Char_t* s = (ud == kU) ? "U" : "D";
      HB1(RCHid+plane*1000+3+ud, Form("RC Tdc %s Plane#%d", s, plane), NbinTdc, MinTdc, MaxTdc);
      // HB1(RCHid+plane*1000+5+ud, Form("RC Tot %s Plane#%d", s, plane), NbinTdc, MinTdc, MaxTdc);
      HB2(RCHid+plane*1000+5+ud, Form("RC DepthTdc %s Plane#%d", s, plane),
	  NumOfSegRC, 0, NumOfSegRC, 16, 0, 16);
      HB1(RCHid+plane*1000+7+ud, Form("RC AdcHigh %s Plane#%d", s, plane), NbinAdc, MinAdc, MaxAdc);
      HB1(RCHid+plane*1000+9+ud, Form("RC AdcLow %s Plane#%d", s, plane), NbinAdc, MinAdc, MaxAdc);
      HB2(RCHid+plane*1000+11+ud, Form("RC Tdc %s%%Seg Plane#%d", s, plane),
          NumOfSegRC, 0., NumOfSegRC, NbinTdc, MinTdc, MaxTdc);
      HB2(RCHid+plane*1000+100+13+ud, Form("RC Tra %s%%Seg Plane#%d", s, plane),
          NumOfSegRC, 0., NumOfSegRC, NbinTdc, MinTdc, MaxTdc);
      HB2(RCHid+plane*1000+100+15+ud, Form("RC Tot %s%%Seg Plane#%d", s, plane),
          NumOfSegRC, 0., NumOfSegRC, 200, 0, 200);

      HB2(RCHid+plane*1000+15+ud, Form("RC Adchigh %s%%Seg Plane#%d", s, plane),
          NumOfSegRC, 0., NumOfSegRC, NbinAdc, MinAdc, MaxAdc);
      HB2(RCHid+plane*1000+17+ud, Form("RC AdcLow %s%%Seg Plane#%d", s, plane),
          NumOfSegRC, 0., NumOfSegRC, NbinAdc, MinAdc, MaxAdc);
      HB2(RCHid+plane*1000+81+ud, Form("RC AdcHigh w/ TDC %s%%Seg Plane%d", s, plane),
	  NumOfSegRC, 0., NumOfSegRC, NbinAdc, MinAdc, MaxAdc);
      HB2(RCHid+plane*1000+83+ud, Form("RC AdcLow w/ TDC %s%%Seg Plane%d", s, plane),
	  NumOfSegRC, 0., NumOfSegRC, NbinAdc, MinAdc, MaxAdc);

      // HB2(RCHid+plane*1000+50+ud, Form("RC Tdc w/maxadc %s%%Seg Plane#%d", s, plane),
      //     NumOfSegRCarr.at(plane), 0., NumOfSegRCarr.at(plane), NbinTdc, MinTdc, MaxTdc);
      // HB1(RCHid+plane*1000+52+ud, Form("RC HitPat %s Plane#%d", s, plane), NumOfSegRCarr.at(plane), 0, NumOfSegRCarr.at(plane));
      // HB2(RCHid+plane*1000+54+ud, Form("RC AdcHigh w/maxadc %s%%Seg Plane#%d", s, plane),
      //     NumOfSegRCarr.at(plane), 0., NumOfSegRCarr.at(plane), NbinAdc, MinAdc, MaxAdc);
      HB2(RCHid+plane*1000+13+ud, Form("RC Time TOT %s%%Seg Plane#%d", s, plane),
	  NumOfSegRC, 0., NumOfSegRC, NbinTime, MinTime, MaxTime);
      HB2(RCHid+plane*1000+56+ud, Form("RC Time Trailing %s%%Seg Plane#%d", s, plane),
	  NumOfSegRC, 0., NumOfSegRC, NbinTime, MinTime, MaxTime);
      HB2(RCHid+plane*1000+58+ud, Form("RC Time Leading %s%%Seg Plane#%d", s, plane),
	  NumOfSegRC, 0., NumOfSegRC, NbinTime, MinTime, MaxTime);
      HB1(RCHid+plane*1000+60+ud, Form("RC Time %s Plane#%d", s, plane), NbinTime, MinTime, MaxTime);
      HB2(RCHid+plane*1000+62+ud, Form("RC adc:tot w/timecut %s%%Plane#%d", s, plane),
          NbinAdc, MinAdc, MaxAdc, 200, 0, 200);
      HB2(RCHid+plane*1000+64+ud, Form("RC de_high:time  %s%%Plane#%d", s, plane),
		  NbinDe, MinDe, MaxDe, NbinTime, MinTime, MaxTime);
    }
  }

  HB1(RCHid +101, "RC NCluster", 100, 0, 100);
  HB1(RCHid +102, "RC Cluster Size", 5, 0, 5);
  HB1(RCHid +103, "RC CTime (Cluster)", 100., -20., 30.);
  HB1(RCHid +104, "RC Tot (Cluster)", NbinTot, MinTot, MaxTot);
  HB2(RCHid +105, "RC CTime%Tot (Cluster)",
      NbinTot, MinTot, MaxTot, NbinTime, MinTime, MaxTime);
  HB1(RCHid +106, "RC Cluster Position",
      NumOfSegRC, -0.5*(Double_t)NumOfSegRC, 0.5*(Double_t)NumOfSegRC);
  HB1(RCHid +107, "RC DeltaE Sum(Cluster)", 800, 0, 80);


  //Tree
  HBTree("ea0c", "tree of Easiroc");
  //Trig
  tree->Branch("evnum",     &event.evnum,     "evnum/I");
  tree->Branch("trigpat",    event.trigpat,   Form("trigpat[%d]/I", NumOfSegTrig));
  tree->Branch("trigflag",   event.trigflag,  Form("trigflag[%d]/I", NumOfSegTrig));

  //BH1
  tree->Branch("bh1_nhits",     &event.bh1_nhits,        "bh1_nhits/I");
  tree->Branch("bh1seg",     event.bh1seg,        Form("bh1seg[%d]/I", NumOfSegBH1));
  tree->Branch("bh1_ncl",       &event.bh1_ncl,          "bh1_ncl/I");
  tree->Branch("bh1_clsize",     event.bh1_clsize,       "bh1_clsize[bh1_ncl]/I");
  tree->Branch("bh1_clseg",      event.bh1_clseg,        "bh1_clseg[bh1_ncl]/D");

  //BFT
#if FHitBranch
  tree->Branch("bft_nhits",     &event.bft_nhits,        "bft_nhits/I");
  tree->Branch("bft_unhits",    &event.bft_unhits,       "bft_unhits/I");
  tree->Branch("bft_dnhits",    &event.bft_dnhits,       "bft_dnhits/I");
  tree->Branch("bft_uhitpat",    event.bft_uhitpat,      "bft_uhitpat[bft_unhits]/I");
  tree->Branch("bft_dhitpat",    event.bft_dhitpat,      "bft_dhitpat[bft_dnhits]/I");
  tree->Branch("bft_utdc",       event.bft_utdc,         Form("bft_utdc[%d][%d]/D",
							      NumOfSegBFT, MaxDepth));
  tree->Branch("bft_dtdc",       event.bft_dtdc,         Form("bft_dtdc[%d][%d]/D",
							      NumOfSegBFT, MaxDepth));
  tree->Branch("bft_utrailing",  event.bft_utrailing,    Form("bft_utrailing[%d][%d]/D",
							      NumOfSegBFT, MaxDepth));
  tree->Branch("bft_dtrailing",  event.bft_dtrailing,    Form("bft_dtrailing[%d][%d]/D",
							      NumOfSegBFT, MaxDepth));
  tree->Branch("bft_utot",       event.bft_utot,         Form("bft_utot[%d][%d]/D",
                                                              NumOfSegBFT, MaxDepth));
  tree->Branch("bft_dtot",       event.bft_dtot,         Form("bft_dtot[%d][%d]/D",
                                                              NumOfSegBFT, MaxDepth));
  tree->Branch("bft_udepth",     event.bft_udepth,       Form("bft_udepth[%d]/I", NumOfSegBFT));
  tree->Branch("bft_ddepth",     event.bft_ddepth,       Form("bft_ddepth[%d]/I", NumOfSegBFT));
#endif
  tree->Branch("bft_ncl",       &event.bft_ncl,          "bft_ncl/I");
  tree->Branch("bft_clsize",     event.bft_clsize,       "bft_clsize[bft_ncl]/I");
  tree->Branch("bft_ctime",      event.bft_ctime,        "bft_ctime[bft_ncl]/D");
  tree->Branch("bft_ctot",       event.bft_ctot,         "bft_ctot[bft_ncl]/D");
  tree->Branch("bft_clpos",      event.bft_clpos,        "bft_clpos[bft_ncl]/D");
  tree->Branch("bft_clseg",      event.bft_clseg,        "bft_clseg[bft_ncl]/D");

  //RC
#if RawHitRCBranch
  tree->Branch("rc_nhits", event.rc_nhits, Form("rc_nhits[%d]/I", NumOfPlaneRC));
  tree->Branch("rc_hitpat", event.rc_hitpat,
               Form("rc_hitpat[%d][%d]/I", NumOfPlaneRC, NumOfSegRC));
  tree->Branch("rc_adc_high", event.rc_adc_high,
               Form("rc_adc_high[%d][%d][%d]/D", NumOfPlaneRC, NumOfSegRC, kUorD));
  tree->Branch("rc_adc_low", event.rc_adc_low,
               Form("rc_adc_low[%d][%d][%d]/D", NumOfPlaneRC, NumOfSegRC, kUorD));
  tree->Branch("rc_tdc", event.rc_tdc,
               Form("rc_tdc[%d][%d][%d][%d]/D",
                    NumOfPlaneRC, NumOfSegRC, kUorD, MaxDepth));
  tree->Branch("rc_mt", event.rc_mt,
               Form("rc_mt[%d][%d][%d]/D", NumOfPlaneRC, NumOfSegRC, MaxDepth));
  tree->Branch("rc_cmt", event.rc_cmt,
               Form("rc_cmt[%d][%d][%d]/D", NumOfPlaneRC, NumOfSegRC, MaxDepth));
  tree->Branch("rc_mtot", event.rc_mtot,
               Form("rc_mtot[%d][%d][%d]/D", NumOfPlaneRC, NumOfSegRC, MaxDepth));
  tree->Branch("rc_de_high", event.rc_de_high,
               Form("rc_de_high[%d][%d]/D", NumOfPlaneRC, NumOfSegRC));
  tree->Branch("rc_de_low", event.rc_de_low,
               Form("rc_de_low[%d][%d]/D", NumOfPlaneRC, NumOfSegRC));
  tree->Branch("rc_ltime", event.rc_ltime,
               Form("rc_ltime[%d][%d][%d][%d]/D",
                    NumOfPlaneRC, NumOfSegRC, kUorD, MaxDepth));
  tree->Branch("rc_ttime", event.rc_ttime,
               Form("rc_ttime[%d][%d][%d][%d]/D",
                    NumOfPlaneRC, NumOfSegRC, kUorD, MaxDepth));
  tree->Branch("rc_tot", event.rc_tot,
               Form("rc_tot[%d][%d][%d][%d]/D",
                    NumOfPlaneRC, NumOfSegRC, kUorD, MaxDepth));
#endif

#if ClusterHitRCBranch
  tree->Branch("rc_ncl",       &event.rc_ncl,          "rc_ncl/I");
  tree->Branch("rc_desum",     &event.rc_desum,        "rc_desum/D");
  tree->Branch("rc_clsize",     event.rc_clsize,       "rc_clsize[rc_ncl]/I");
  tree->Branch("rc_cltime",     event.rc_cltime,       "rc_ctime[rc_ncl]/D");
  tree->Branch("rc_cltot",      event.rc_cltot,        "rc_ctot[rc_ncl]/D");
  tree->Branch("rc_clseg",      event.rc_clseg,        "rc_clseg[rc_ncl]/D");
  tree->Branch("rc_clde",       event.rc_clde,         "rc_clde[rc_ncl]/D");
  tree->Branch("rc_clplane",    event.rc_clplane,      "rc_clplane[rc_ncl]/I");
  tree->Branch("rc_clpos",      event.rc_clpos,        "rc_clpos[rc_ncl]/D");
#endif

  // HPrint();
  return true;
}

//_____________________________________________________________________________
Bool_t
ConfMan::InitializeParameterFiles()
{
  return
    (InitializeParameter<DCGeomMan>("DCGEO")    &&
     InitializeParameter<HodoParamMan>("HDPRM") &&
     InitializeParameter<HodoPHCMan>("HDPHC")   &&
     InitializeParameter<UserParamMan>("USER"));
}

//_____________________________________________________________________________
Bool_t
ConfMan::FinalizeProcess()
{
  return true;
}
