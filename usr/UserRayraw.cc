// -*- C++ -*-
// Author: Rintaro Kurata

#include "VEvent.hh"

#include <iostream>
#include <sstream>
#include <cmath>

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

namespace
{
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

  std::vector<std::vector<Double_t>> waveform; // [seg][nsample]
  std::vector<Double_t>              pedestal;
  std::vector<Double_t>              max_adc;
  std::vector<Double_t>              adc_all;
  std::vector<Double_t>              integral;
  std::vector<std::vector<Double_t>> leading;  // [seg][depth]
  std::vector<std::vector<Double_t>> trailing; // [seg][depth]

  void clear();
};

//_____________________________________________________________________________
void
Event::clear()
{
  evnum = 0;

  waveform.clear();
  pedestal.clear();
  max_adc.clear();
  adc_all.clear();
  integral.clear();
  leading.clear();
  trailing.clear();
}

//_____________________________________________________________________________
namespace root
{
  Event  event;
  TH1   *h[MaxHist];
  TTree *tree;
  enum eDetHid
    {
      RayrawHid  = 100000,
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
  static const auto MinRange = gUser.GetParameter("RangeRAYRAW", 0);
  static const auto MaxRange = gUser.GetParameter("RangeRAYRAW", 1);
  static const Double_t V_per_ch = 3.3/1024.; // [V]
  static const Double_t T_per_ch = 13.33; // [ns]

  RawData rawData;
  rawData.DecodeHits();
  HodoAnalyzer hodoAna(rawData);

  event.evnum = gUnpacker.get_event_number();

  HF1(1, 0);

    {
    // from Rawdata
    const auto& cont = rawData.GetHodoRawHC("RAYRAW");
    Int_t nh = cont.size();
    Double_t true_hit = 0;
    Double_t baseline[NumOfSegRayraw] = {qnan};
    //    std::vector<std::vector<Int_t>> fadc_single;

    event.waveform.resize(nh);
    event.leading.resize(nh);
    event.trailing.resize(nh);

    for( int i=0; i<nh; ++i ){

      HodoRawHit *hit = cont[i];

      Int_t seg     = hit->SegmentId();
      Int_t max_adc = 0;
      // Int_t min_l   = 2450;
      // Int_t max_l   = 2600;
      // Int_t min_t   = 2400;
      // Int_t max_t   = 2550;
      // Int_t leading_hit_in = 0;
      // Int_t trailing_hit_in = 0;
      // Int_t leading_hit_out = 0;
      // Int_t trailing_hit_out = 0;
      baseline[seg] = gUser.GetParameter("PedRAYRAW", seg);
      Double_t integral    = 0;
      Double_t V_fadc      = qnan;
      Double_t V_fadc_calc = qnan;
      Double_t T_fadc      = qnan;
      Double_t charge      = 0;

      Int_t nsample        = 0;

      Int_t hid_wf         = RayrawHid + (seg)*1000 + 0; // Raw Waveform
      Int_t hid_wf_vt      = RayrawHid + (seg)*1000 + 1; // Waveform in V-T
      Int_t hid_wf_vt_calc = RayrawHid + (seg)*1000 + 2; // Waveform after remove baseline in V-T
      Int_t hid_ped        = RayrawHid + (seg)*1000 + 3; // Pedestal (w/o signal)
      Int_t hid_adc        = RayrawHid + (seg)*1000 + 4; // Max ADC
      Int_t hid_adc_all    = RayrawHid + (seg)*1000 + 5; // All FADC
      Int_t hid_integral   = RayrawHid + (seg)*1000 + 6; // Integral
      Int_t hid_charge     = RayrawHid + (seg)*1000 + 7; // Integral in pC
      Int_t hid_tdc_l      = RayrawHid + (seg)*1000 + 8; // TDC Leading
      Int_t hid_tdc_t      = RayrawHid + (seg)*1000 + 9; // TDC Trailing
      // Int_t hid_tot     = RayrawHid + (seg)*1000 + 6; // TOT
      // Int_t hid_single  = RayrawHid + (seg)*1000 + 100 + event.evnum; // Waveform of single event

      // TDC block
      // Leading
      for (const auto& tdc_l : hit->GetArrayTdcLeading() ) {
	HF1(hid_tdc_l, tdc_l);
	event.leading[seg].push_back(tdc_l);

	// if(tdc_l != 0 ){
	//   if(min_l < tdc_l && tdc_l < max_l){
	//     leading_hit_in += 1;
	//   }else{
	//     leading_hit_out += 1;
	//   }
	// }
      }

      // Trailing
      for (const auto& tdc_t : hit->GetArrayTdcTrailing() ) {
	HF1(hid_tdc_t, tdc_t);
	event.trailing[seg].push_back(tdc_t);

	// if(tdc_t != 0 ){
	//   if(min_t < tdc_t && tdc_t < max_t){
	//     trailing_hit_in += 1;
	//   }else{
	//     trailing_hit_out += 1;
	//   }
	// }
      }

      // ADC block
      for(const auto& fadc : hit->GetArrayAdc()){
	V_fadc      = fadc*V_per_ch;
	V_fadc_calc = (fadc - baseline[seg])*V_per_ch;
	// T_fadc      = nsample*T_per_ch;

	// Raw Waveform
	HF2(hid_wf, nsample, fadc);
	event.waveform[seg].push_back(fadc);

	// Waveform in V-T
	// HF2(hid_wf_vt, T_fadc, V_fadc);
	HF2(hid_wf_vt, nsample, V_fadc);

	// Waveform after remove baseline in V-T
	// HF2(hid_wf_vt_calc, T_fadc, V_fadc_calc);
	HF2(hid_wf_vt_calc, nsample, V_fadc_calc);

	// All ADC
	  HF1(hid_adc_all, fadc);
	  event.adc_all.push_back(fadc);

	// Pedestal (w/o signal)
	  if(nsample < MinRange || MaxRange < nsample){
	    HF1(hid_ped, fadc);
	    event.pedestal.push_back(fadc);
	  }

	// Max ADC
	if(MinRange <= nsample && nsample <= MaxRange){
	  integral += fadc - baseline[seg];
	  charge   += V_fadc_calc / 50 * T_per_ch*1000; // [pC]
	  if (fadc > max_adc)
	    max_adc = fadc;
	}

	// // efficiency?
	// if(28 < j && j < 32){
	//   if(fadc > 522 && fadc > max_adc)
	//     true_hit += 1;
	// }
	nsample++;
      } // for ADC block

      // Charge
      HF1(hid_integral, integral);
      event.integral.push_back(integral);

      // Charge in pC
      HF1(hid_charge, charge);

      // Max ADC
      HF1(hid_adc, max_adc);
      event.max_adc.push_back(max_adc);

    } // for nh

    //    std::cout << "Efficiency = " << true_hit/nh << std::endl;

    return true;
  }
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
  const Int_t  NbinFADC_X     = 100;
  const Int_t  NbinFADC_Y     = 1024;
  const double MinFADC        = 0.;
  const double MaxFADC        = 1024;
  const double NbinTDC        = 10000;
  const double MinTDC         = 0;
  const double MaxTDC         = 10000;
  // const double MinTime        = -20. * 1024;
  // const double MaxTime        = 20. * 1024;
  // const double MinPulseHeight = -1.5;
  // const double MaxPulseHeight =  1.5;
  // const double MaxEnergy      = 200.;

  HB1( 1, "Status",  20,   0., 20.);

  char buf[100];

  for (Int_t seg=0; seg<NumOfSegRayraw; seg++) {

    sprintf(buf, "RAYRAW - Raw Waveform (ch%d)", seg);
    Int_t hid = RayrawHid + (seg)*1000 + 0;
    HB2( hid, buf, NbinFADC_X, 0, NbinFADC_X, NbinFADC_Y, MinFADC, MaxFADC );

    sprintf(buf, "RAYRAW - Waveform(V-T) (ch%d)", seg);
    hid = RayrawHid + (seg)*1000 + 1;
    HB2( hid, buf, NbinFADC_X, 0, NbinFADC_X, NbinFADC_Y, MinFADC, 3.3 );
    // HB2( hid, buf, NbinFADC_X, 0, 13.33*NbinFADC_X, NbinFADC_Y, MinFADC, 3.3 );

    sprintf(buf, "RAYRAW - Waveform after calc (V-T) (ch%d)", seg);
    hid = RayrawHid + (seg)*1000 + 2;
    HB2( hid, buf, NbinFADC_X, 0, NbinFADC_X, NbinFADC_Y, -1.65, 1.65 );
    // HB2( hid, buf, NbinFADC_X, 0, 13.33*NbinFADC_X, NbinFADC_Y, -1.65, 1.65 );

    sprintf(buf, "RAYRAW - Pedestal (ch%d)", seg);
    hid = RayrawHid + (seg)*1000 + 3;
    HB1( hid, buf, NbinFADC_Y, MinFADC, MaxFADC);

    sprintf(buf, "RAYRAW - Max ADC (ch%d)", seg);
    hid = RayrawHid + (seg)*1000 + 4;
    HB1( hid, buf, NbinFADC_Y, MinFADC, MaxFADC);

    sprintf(buf, "RAYRAW - ADC All (ch%d)", seg);
    hid = RayrawHid + (seg)*1000 + 5;
    HB1( hid, buf, NbinFADC_Y, MinFADC, MaxFADC);

    sprintf(buf, "RAYRAW - Integral (ch%d)", seg);
    hid = RayrawHid + (seg)*1000 + 6;
    // HB1( hid, buf, 2100*100, -100, 2000);
    HB1( hid, buf, 5000, -1000, 4000);

    sprintf(buf, "RAYRAW - Charge (ch%d)", seg);
    hid = RayrawHid + (seg)*1000 + 7;
    // HB1( hid, buf, 10000, -100, 10000);
    HB1( hid, buf, 5000, -1000, 4000);

    sprintf(buf, "RAYRAW - TDC Leading (ch%d)", seg);
    hid = RayrawHid + (seg)*1000 + 8;
    HB1( hid, buf, NbinTDC, MinTDC,  MaxTDC);

    sprintf(buf, "RAYRAW - TDC Trailing (ch%d)", seg);
    hid = RayrawHid + (seg)*1000 + 9;
    HB1( hid, buf, NbinTDC, MinTDC,  MaxTDC);

    // sprintf(buf, "RAYRAW - TOT (ch%d)", seg);
    // hid = RayrawHid + (seg)*1000 + 12;
    // HB1( hid, buf, NbinTDC, MinTDC,  MaxTDC);

  }

  //Tree
  HBTree( "tree", "tree" );
  tree->Branch("evnum",       &event.evnum,     "evnum/I");
  tree->Branch("waveform",    &event.waveform);
  tree->Branch("pedestal",    &event.pedestal);
  tree->Branch("max_adc",     &event.max_adc);
  tree->Branch("adc_all",     &event.adc_all);
  tree->Branch("integral",    &event.integral);
  tree->Branch("leading",     &event.leading);
  tree->Branch("trailing",    &event.trailing);

  HPrint();
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
