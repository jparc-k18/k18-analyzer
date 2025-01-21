// -*- C++ -*-

#include "VEvent.hh"

#include <cmath>
#include <iostream>
#include <sstream>

#include <UnpackerManager.hh>

// #include "BH2Cluster.hh"
#include "BH2Hit.hh"
#include "FiberCluster.hh"
#include "FiberHit.hh"
#include "ConfMan.hh"
#include "DetectorID.hh"
#include "RMAnalyzer.hh"
#include "HodoHit.hh"
#include "HodoAnalyzer.hh"
#include "HodoCluster.hh"
#include "HodoParamMan.hh"
#include "HodoPHCMan.hh"
#include "HodoRawHit.hh"
#include "S2sLib.hh"
#include "RawData.hh"
#include "RootHelper.hh"
#include "UserParamMan.hh"
#include "DCGeomMan.hh"

// #define TimeCut    1 // in cluster analysis
#define FHitBranch 0 // make FiberHit branches (becomes heavy)
#define HodoHitPos 0

namespace
{
  using namespace root;
  const auto qnan = TMath::QuietNaN();
  auto &gUnpacker = hddaq::unpacker::GUnpacker::get_instance();
  auto &gRM = RMAnalyzer::GetInstance();
  auto &gUser = UserParamMan::GetInstance();
}

//_____________________________________________________________________________
struct Event
{
  Int_t evnum;
  Int_t spill;

  // Caen V792
  Double_t v792_adc[NumOfSegCaenV792];
  Double_t v792_tdc[NumOfSegCaenV792][MaxDepth];
  Double_t v792_adc_w_tdc[NumOfSegCaenV792];

  // VME-EASIROC
  Int_t vmeeasiroc_nhits[4];
  Int_t vmeeasiroc_hitpat[4];
  Double_t vmeeasiroc_leading [4][NumOfSegVMEEASIROC][MaxDepth];
  Double_t vmeeasiroc_trailing[4][NumOfSegVMEEASIROC][MaxDepth];
  Double_t vmeeasiroc_tot[4][NumOfSegVMEEASIROC][MaxDepth];
  Double_t vmeeasiroc_highgain[4][NumOfSegVMEEASIROC];
  Double_t vmeeasiroc_lowgain [4][NumOfSegVMEEASIROC];
  Double_t vmeeasiroc_hg_w_tdc[4][NumOfSegVMEEASIROC];
  Double_t vmeeasiroc_leading_onehit [4][NumOfSegVMEEASIROC];

  void clear();
};

//_____________________________________________________________________________
void Event::clear()
{
  evnum = 0;
  spill = 0;
  for (Int_t i = 0; i < NumOfSegCaenV792; ++i)
  {
    v792_adc[i] = qnan;
    v792_adc_w_tdc[i] = qnan;
    for(Int_t m=0; m<MaxDepth; m++){
      v792_tdc[i][m] = qnan;
    }
  }

  for (Int_t i = 0; i < 4; ++i)
  {
    vmeeasiroc_nhits[i] = 0;
    for (Int_t j = 0; j < NumOfSegVMEEASIROC; ++j)
    {
      for(Int_t m=0; m<MaxDepth; m++){
	vmeeasiroc_hitpat[i] = -1;
	vmeeasiroc_leading [i][j][m] = qnan;
	vmeeasiroc_trailing[i][j][m] = qnan;
	vmeeasiroc_tot[i][j][m] = qnan;
	vmeeasiroc_highgain[i][j] = qnan;
	vmeeasiroc_lowgain [i][j] = qnan;
	vmeeasiroc_hg_w_tdc[i][j] = qnan;
      }
      vmeeasiroc_leading_onehit[i][j] = qnan;
    }
  }
}

//_____________________________________________________________________________
namespace root
{
  Event event;
  TH1 *h[MaxHist];
  TTree *tree;
  enum eDetHid
  {
    CaenV792Hid = 10000,
    VMEEASIROCHid = 20000,
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
  // static const auto MinTdcBH1 = gUser.GetParameter("TdcBH1", 0);
  // static const auto MaxTdcBH1 = gUser.GetParameter("TdcBH1", 1);

  RawData rawData;
  // HodoAnalyzer hodoAna(rawData);

  event.evnum = gUnpacker.get_event_number();

  HF1(1, 0);

  { ///// CaenV792
    static const auto device_id = gUnpacker.get_device_id("HODO");
    static const auto adc_id = gUnpacker.get_data_id("HODO", "adc");
    static const auto tdc_id = gUnpacker.get_data_id("HODO", "tdc");
    static const Int_t n_seg = NumOfSegCaenV792;
    for (Int_t seg = 0; seg < n_seg; ++seg)
    {
      auto nhit = gUnpacker.get_entries(device_id, 0, seg, 0, adc_id);
      UInt_t adc = 0;
      if (nhit != 0)
      {
        adc = gUnpacker.get(device_id, 0, seg, 0, adc_id);
        HF1(CaenV792Hid + seg, adc);
        event.v792_adc[seg] = adc;
      }
      Bool_t hit_flag = false;
      for (Int_t i = 0, j = gUnpacker.get_entries(device_id, 0, seg, 0, tdc_id); i < j; ++i)
      {
        UInt_t tdc = gUnpacker.get(device_id, 0, seg, 0, tdc_id, i);
        if (tdc != 0)
        {
          HF1(CaenV792Hid + seg + 1000, tdc);
          event.v792_tdc[seg][i] = tdc;
          hit_flag = true;
        }
        if (hit_flag)
        {
          HF1(CaenV792Hid + seg + 2000, adc);
          event.v792_adc_w_tdc[seg] = adc;
        }
      }
    }
    // gUnpacker.dump_data_device(k_device);
  }

  { ///// VME-EASIROC
    static const auto device_id = gUnpacker.get_device_id("VMEEASIROC");
    static const auto leading_id = gUnpacker.get_data_id("VMEEASIROC", "leading");
    static const auto trailing_id = gUnpacker.get_data_id("VMEEASIROC", "trailing");
    static const auto highgain_id = gUnpacker.get_data_id("VMEEASIROC", "highgain");
    static const auto lowgain_id = gUnpacker.get_data_id("VMEEASIROC", "lowgain");

    // TDC gate range
    static const int tdc_min = gUser.GetParameter("TdcVMEEASIROC", 0);
    static const int tdc_max = gUser.GetParameter("TdcVMEEASIROC", 1);

    for (int l = 0; l < 4; ++l)
    {
      Int_t plane = 11 + l;
      for (int seg = 0; seg < NumOfSegVMEEASIROC; ++seg)
      {
        {
          int nhit_l = gUnpacker.get_entries(device_id, plane, seg, 0, leading_id);
          bool hit_flag_wt = false;
          if (nhit_l != 0)
          { // hit pattern
            HF1(VMEEASIROCHid + l * 100, seg);
            event.vmeeasiroc_nhits[l] = nhit_l;
            event.vmeeasiroc_hitpat[l]= seg;
            for (int i = 0; i < nhit_l; ++i)
            { // tdc
              int tdc = gUnpacker.get(device_id, plane, seg, 0, leading_id, i);
              HF2(VMEEASIROCHid + l * 100 + 1000, seg, tdc);
              HF1(VMEEASIROCHid + l * 100 + 1000 + seg + 1, tdc);
              event.vmeeasiroc_leading[l][seg][i] = tdc;
	      if(nhit_l==1) event.vmeeasiroc_leading_onehit[l][seg] = tdc;
              if (tdc_min < tdc && tdc < tdc_max)
              {
                hit_flag_wt = true;
              }
            }
            if (hit_flag_wt)
            { // high gain w/ tdc
              int nhit_hg = gUnpacker.get_entries(device_id, plane, seg, 0, highgain_id);
              if (nhit_hg != 0)
              {
                int adc_hg = gUnpacker.get(device_id, plane, seg, 0, highgain_id);
                HF1(VMEEASIROCHid + l * 100 + 5000 + seg + 1, adc_hg);
                event.vmeeasiroc_hg_w_tdc[l][seg] = adc_hg;
              }
              // low gain w/ tdc
              int nhit_lg = gUnpacker.get_entries(device_id, plane, seg, 0, lowgain_id);
              if (nhit_lg != 0)
              {
                int adc_lg = gUnpacker.get(device_id, plane, seg, 0, lowgain_id);
                HF1(VMEEASIROCHid + l * 100 + 6000 + seg + 1, adc_lg);
              }
            }
          }

          { // tot
            int nhit_l = gUnpacker.get_entries(device_id, plane, seg, 0, leading_id);
            int nhit_t = gUnpacker.get_entries(device_id, plane, seg, 0, trailing_id);
            Int_t hit_l_max = 0;
            Int_t hit_t_max = 0;
            if (nhit_l != 0)
              hit_l_max = gUnpacker.get(device_id, plane, seg, 0, leading_id, nhit_l - 1);
            if (nhit_t != 0)
              hit_t_max = gUnpacker.get(device_id, plane, seg, 0, trailing_id, nhit_t - 1);

            if (nhit_l == nhit_t && hit_l_max > hit_t_max)
            {
              for (int i = 0; i < nhit_l; ++i)
              {
                int tdc = gUnpacker.get(device_id, plane, seg, 0, leading_id, i);
                int tdc_t = gUnpacker.get(device_id, plane, seg, 0, trailing_id, i);
                //int tot = tdc_t - tdc; 
		int tot = tdc - tdc_t; // in case of tdc, leading - trailing
                HF2(VMEEASIROCHid + l * 100 + 2000, seg, tot);
                HF1(VMEEASIROCHid + l * 100 + 2000 + seg + 1, tot);
                event.vmeeasiroc_trailing[l][seg][i] = tdc_t;
		event.vmeeasiroc_tot[l][seg][i] = tot;
              }
            }
          }

          { // high gain
            int nhit_hg = gUnpacker.get_entries(device_id, plane, seg, 0, highgain_id);
            for (int i = 0; i < nhit_hg; ++i)
            {
              int adc_hg = gUnpacker.get(device_id, plane, seg, 0, highgain_id, i);
              HF2(VMEEASIROCHid + l * 100 + 3000, seg, adc_hg);
              HF1(VMEEASIROCHid + l * 100 + 3000 + seg + 1, adc_hg);
              event.vmeeasiroc_highgain[l][seg] = adc_hg;
            }
          }

          { // low gain
            int nhit_lg = gUnpacker.get_entries(device_id, plane, seg, 0, lowgain_id);
            for (int i = 0; i < nhit_lg; ++i)
            {
              int adc_lg = gUnpacker.get(device_id, plane, seg, 0, lowgain_id, i);
              HF2(VMEEASIROCHid + l * 100 + 4000, seg, adc_lg);
              HF1(VMEEASIROCHid + l * 100 + 4000 + seg + 1, adc_lg);
              event.vmeeasiroc_lowgain[l][seg] = adc_lg;
            }
          }

          { // multi hit
            int nhit_l = gUnpacker.get_entries(device_id, plane, seg, 0, leading_id);
            HF2(VMEEASIROCHid + l * 100 + 5000, seg, nhit_l);
          }
        }
      } // for seg
    } // for board
  } // VME-EASIROC
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
namespace
{
  const Int_t NbinAdc = 4096;
  const Double_t MinAdc = 0.;
  const Double_t MaxAdc = 4096.;

  const Int_t NbinTdc = 4096;
  const Double_t MinTdc = 0.;
  const Double_t MaxTdc = 4096.;

  const Int_t    NbinTot =  170;
  const Double_t MinTot  =  -10.;
  const Double_t MaxTot  =  160.;

  const Int_t NbinTdcHr = 1e6 / 10;
  const Double_t MinTdcHr = 0.;
  const Double_t MaxTdcHr = 1e6;
}

//_____________________________________________________________________________
Bool_t
ConfMan::InitializeHistograms()
{
  HB1(1, "Status", 20, 0., 20.);

  // !!! DO NOT MAKE ADC vs TDC TH2 !!!
  // CaenV792
  for (Int_t i = 0; i < NumOfSegTrig; ++i)
  {
    HB1(CaenV792Hid + i, Form("CaenV792 ADC %d", i), NbinAdc, MinAdc, MaxAdc);
    HB1(CaenV792Hid + i + 2000, Form("Caen V792 ADC w/ TDC %d", i), NbinAdc, MinAdc, MaxAdc);
    if (i == 21)
    {
      HB1(CaenV792Hid + i + 1000, Form("CaenV792 TDC %d", i), NbinTdcHr, MinTdcHr, MaxTdcHr);
    }
    else
    {
      HB1(CaenV792Hid + i + 1000, Form("CaenV792 TDC %d", i), NbinTdc, MinTdc, MaxTdc);
    }
  }

  // VME-EASIROC
  for (int l = 0; l < 4; ++l)
  {
    Int_t plane = 11 + l;
    HB1(VMEEASIROCHid + l * 100, Form("VMEEASIROC %d Hitpattern", plane), NumOfSegVMEEASIROC, 0, NumOfSegVMEEASIROC);
    HB2(VMEEASIROCHid + l * 100 + 1000, Form("VMEEASIROC TDC %d", plane), NumOfSegVMEEASIROC, 0, NumOfSegVMEEASIROC, NbinTdc, MinTdc, MaxTdc);
    HB2(VMEEASIROCHid + l * 100 + 2000, Form("VMEEASIROC TOT %d", plane), NumOfSegVMEEASIROC, 0, NumOfSegVMEEASIROC, NbinTot, MinTot, MaxTot);
    HB2(VMEEASIROCHid + l * 100 + 3000, Form("VMEEASIROC HG %d", plane), NumOfSegVMEEASIROC, 0, NumOfSegVMEEASIROC, NbinAdc, MinAdc, MaxAdc);
    HB2(VMEEASIROCHid + l * 100 + 4000, Form("VMEEASIROC LG %d", plane), NumOfSegVMEEASIROC, 0, NumOfSegVMEEASIROC, NbinAdc, MinAdc, MaxAdc);
    HB2(VMEEASIROCHid + l * 100 + 5000, Form("VMEEASIROC Multihit %d", plane), NumOfSegVMEEASIROC, 0, NumOfSegVMEEASIROC, 20, 0, 20);
    for (Int_t seg = 0; seg < NumOfSegVMEEASIROC; ++seg)
    {
      HB1(VMEEASIROCHid + l * 100 + 1000 + seg + 1, Form("VMEEASIROC TDC %d-%d", plane, seg), NbinTdc, MinTdc, MaxTdc);
      HB1(VMEEASIROCHid + l * 100 + 2000 + seg + 1, Form("VMEEASIROC TOT %d-%d", plane, seg), NbinTot, MinTot, MaxTot);
      HB1(VMEEASIROCHid + l * 100 + 3000 + seg + 1, Form("VMEEASIROC HG %d-%d", plane, seg), NbinAdc, MinAdc, MaxAdc);
      HB1(VMEEASIROCHid + l * 100 + 4000 + seg + 1, Form("VMEEASIROC LG %d-%d", plane, seg), NbinAdc, MinAdc, MaxAdc);
      HB1(VMEEASIROCHid + l * 100 + 5000 + seg + 1, Form("VMEEASIROC HG w/ TDC %d-%d", plane, seg), NbinAdc, MinAdc, MaxAdc);
      HB1(VMEEASIROCHid + l * 100 + 6000 + seg + 1, Form("VMEEASIROC LG w/ TDC %d-%d", plane, seg), NbinAdc, MinAdc, MaxAdc);
    }
  }

  // Tree
  HBTree("tree", "tree of Counter");
  tree->Branch("evnum", &event.evnum, "evnum/I");
  tree->Branch("spill", &event.spill, "spill/I");
  // CAEN V792
  tree->Branch("v792_adc", &event.v792_adc, Form("v792_adc[%d]/D", NumOfSegCaenV792));
  tree->Branch("v792_tdc", event.v792_tdc, Form("v792_tdc[%d][%d]/D", NumOfSegCaenV792, MaxDepth));
  tree->Branch("v792_adc_w_tdc", &event.v792_adc_w_tdc, Form("v792_adc_w_tdc[%d]/D", NumOfSegCaenV792));
  // VME-EASIROC 
  tree->Branch("vmeeasiroc_nhits", event.vmeeasiroc_nhits, Form("vmeeasiroc_nhits[%d]/I", 4));
  tree->Branch("vmeeasiroc_hitpat", event.vmeeasiroc_hitpat, Form("vmeeasiroc_hitpat[%d]/I", 4));
  tree->Branch("vmeeasiroc_leading", event.vmeeasiroc_leading,
               Form("vmeeasiroc_leading[%d][%d][%d]/D", 4, NumOfSegVMEEASIROC, MaxDepth));
  tree->Branch("vmeeasiroc_trailing", event.vmeeasiroc_trailing,
               Form("vmeeasiroc_trailing[%d][%d][%d]/D", 4, NumOfSegVMEEASIROC, MaxDepth));
  tree->Branch("vmeeasiroc_tot", event.vmeeasiroc_tot,
               Form("vmeeasiroc_tot[%d][%d][%d]/D", 4, NumOfSegVMEEASIROC, MaxDepth));
  tree->Branch("vmeeasiroc_highgain", event.vmeeasiroc_highgain,
               Form("vmeeasiroc_highgain[%d][%d]/D", 4, NumOfSegVMEEASIROC));
  tree->Branch("vmeeasiroc_lowgain", event.vmeeasiroc_lowgain,
	       Form("vmeeasiroc_lowgain[%d][%d]/D", 4, NumOfSegVMEEASIROC));
  tree->Branch("vmeeasiroc_hg_w_tdc", event.vmeeasiroc_hg_w_tdc,
	       Form("vmeeasiroc_hg_w_tdc[%d][%d]/D", 4, NumOfSegVMEEASIROC));
  tree->Branch("vmeeasiroc_leading_onehit", event.vmeeasiroc_leading_onehit,
               Form("vmeeasiroc_leading[%d][%d]/D", 4, NumOfSegVMEEASIROC));

  return true;
}

//_____________________________________________________________________________
Bool_t
ConfMan::InitializeParameterFiles()
{
  return (InitializeParameter<DCGeomMan>("DCGEO") &&
          InitializeParameter<HodoParamMan>("HDPRM") &&
          InitializeParameter<HodoPHCMan>("HDPHC") &&
          InitializeParameter<UserParamMan>("USER"));
}

//_____________________________________________________________________________
Bool_t
ConfMan::FinalizeProcess()
{
  return true;
}
