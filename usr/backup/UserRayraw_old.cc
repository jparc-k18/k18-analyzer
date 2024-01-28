/**
 *  file: UserBGO.cc
 *  date: 2017.04.10
 *  
 */

// modified for RAYRAW (copied from UserMppcASIC.cc)
// date:    2024.01.13
// author:  Rintaro Kurata


#include <iostream>
#include <sstream>
#include <cmath>

// #include "HodoWave1Hit.hh"
// #include "HodoWave2Hit.hh"
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
#include "RawData.hh"
#include "UnpackerManager.hh"
#include "VEvent.hh"
#include "UserParamMan.hh"

#define HodoCut 0 // with BH1/BH2
#define TimeCut 0 // in cluster analysis

namespace
{
using namespace root;
const std::string& classname("EventK11BH1");
const auto qnan = TMath::QuietNaN();
RMAnalyzer&         gRM   = RMAnalyzer::GetInstance();
const UserParamMan& gUser = UserParamMan::GetInstance();
}

//______________________________________________________________________________


#ifndef MaxHits2
#define MaxHits2 60
#endif


struct Event
{
  int evnum;

  void clear();
};

void
Event::clear()
{
}

//______________________________________________________________________________
namespace root
{
  Event  event;
  TH1   *h[MaxHist];
  TTree *tree;

  enum eDetHid {
    RayrawHid  = 10000
  };


}

//______________________________________________________________________________
bool
ProcessingBegin( void )
{
  event.clear();
  return true;
}

//______________________________________________________________________________
bool
ProcessingNormal( void )
{
  static const std::string funcname("["+classname+"::"+__func__+"]");

  RawData rawData;
  rawData.DecodeHits();
  HodoAnalyzer hodoAna(rawData);

  gRM.Decode();

  event.evnum = gRM.EventNumber();

  HF1(1, 0);

  {
    // from Rawdata
    const auto& cont = rawData.GetHodoRawHC("RAYRAW");
    Int_t nh = cont.size();

    for( int i=0; i<nh; ++i ){

      HodoRawHit *hit = cont[i];
      Int_t seg     = hit->SegmentId();
      Int_t max_adc = 0;
      Int_t min_l   = 2400;
      Int_t max_l   = 2600;
      Int_t min_t   = 2500;
      Int_t max_t   = 2600;
      Int_t leading_hit_in = 0;
      Int_t trailing_hit_in = 0;
      Int_t leading_hit_out = 0;
      Int_t trailing_hit_out = 0;
      Int_t charge = 0;

      Int_t hid_wf_raw  = RayrawHid + (seg+1)*100 + 0;  // Raw Waveform
      Int_t hid_wf      = RayrawHid + (seg+1)*100 + 1;  // Waveform
      Int_t hid_wf_wl   = RayrawHid + (seg+1)*100 + 2;  // Waveform w/ Leading
      Int_t hid_wf_wol  = RayrawHid + (seg+1)*100 + 3;  // Waveform w/o Leading
      Int_t hid_wf_wt   = RayrawHid + (seg+1)*100 + 4;  // Waveform w/ Trailing
      Int_t hid_wf_wot  = RayrawHid + (seg+1)*100 + 5;  // Waveform w/o Trailing
      Int_t hid_charge  = RayrawHid + (seg+1)*100 + 6;  // Charge
      Int_t hid_adc     = RayrawHid + (seg+1)*100 + 7;  // ADC
      Int_t hid_adc_wl  = RayrawHid + (seg+1)*100 + 8;  // ADC w/ Leading
      Int_t hid_adc_wol = RayrawHid + (seg+1)*100 + 9;  // ADC w/o Leading
      Int_t hid_tdc_l   = RayrawHid + (seg+1)*100 + 10; // TDC Leading
      Int_t hid_tdc_t   = RayrawHid + (seg+1)*100 + 11; // TDC Trailing
      Int_t hid_tot     = RayrawHid + (seg+1)*100 + 12; // TOT

      // TDC block
      // Leading
      for (Int_t j=0, k=hit->GetSizeTdcLeading(); j<k; j++ ) {
	auto  tdc_l = hit->GetTdcLeading(0, j);
	HF1(hid_tdc_l, tdc_l);
	//	std::cout << "i= " << i << ", seg= " << seg << ", nhtdc_l= " << nhtdc_l << ", j= " << j << ", tdc_l= " << tdc_l << std::endl;

	if(tdc_l != 0 ){
	  if(min_l < tdc_l && tdc_l < max_l){
	    leading_hit_in += 1;
	  }else{ 
	    leading_hit_out += 1;
	  }
	}
      }

      // Trailing
      for (Int_t j=0, k=hit->GetSizeTdcTrailing(); j<k; j++ ) {
	auto tdc_t = hit->GetTdcTrailing(0, j);
	HF1(hid_tdc_t, tdc_t);
	//	std::cout << "i= " << i << ", seg= " << seg << ", nhtdc_t= " << nhtdc_t << ", j= " << j << ", tdc_t= " << tdc_t << std::endl;

	if(tdc_t != 0 ){
	  if(min_t < tdc_t && tdc_t < max_t){
	    trailing_hit_in += 1;
	  }else{
	    trailing_hit_out += 1;
	  }
	}
      }

      // ADC block
      for (Int_t j=0, k=hit->GetSizeAdcHigh(); j<k; j++ ) {
	// for(const auto& fadc : hit->GetArrayAdc())
	auto fadc = hit->GetAdcUp(j);

	// Waveform
	HF2(hid_wf_raw, j, fadc);
	//	std::cout << "i= " << i << ", seg= " << seg << ", nhfadc= " << nhfadc << ", j= " << j << ", fadc= " << fadc << std::endl;

	// Charge
	HF1(hid_charge, charge);

	// Max ADC
	if(45 < j && j < 70){
	  charge += fadc;
	  if (fadc > max_adc)
	    max_adc = fadc;
	}

	// Waveform w/ leading
	if( leading_hit_in > 0){
	  HF2(hid_wf_wl, j, fadc);
	}	

	// Waveform w/o leading
	if( leading_hit_out > 0){
	  HF2(hid_wf_wol, j, fadc);
	}	

	// Waveform w/ trailing
	if( trailing_hit_in > 0){
	  HF2(hid_wf_wt, j, fadc);
	}	

	// Waveform w/o trailing
	if( trailing_hit_in > 0){
	  HF2(hid_wf_wot, j, fadc);
	}	

      } // for nhfadc

      // Max ADC
      HF1(hid_adc, max_adc);

      // Max ADC w/ leading
      if( leading_hit_in > 0){
	HF1(hid_adc_wl, max_adc);
      }	

      // Max ADC w/o leading
      if( leading_hit_out > 0){
	HF1(hid_adc_wol, max_adc);
      }	

    } // for nh

// #if 0
//     // T2
//     hodoAna->DecodeRayrawHits(rawData);
//     {
//       Int_t nh = hodoAna->GetNHitsRayraw();
//       for(Int_t i=0; i<nh; ++i){
// 	HodoWaveHit *hit = hodoAna->GetHitsRayraw(i);
// 	if(!hit) continue;
// 	Int_t seg = hit->SegmentId();
// 	Int_t nhfadc1 = hit->GetNumOfWaveform();
// 	// TDC Leading
// 	Int_t nhit_tdc = hit->GetEntries();
// 	if (nhit_tdc>0) {
// 	  double t1 = hit->GetT();
// 	  Int_t hid = RayrawHid + (seg+1)*100+3;;
// 	  HF1(hid, t1);
// 	}
// 	HodoWave1Hit *hit = hodoAna->GetHitsRayraw(i);
// 	if(!hit) continue;
// 	Int_t seg = hit->SegmentId();
// 	Int_t nhfadc1 = hit->GetNumOfWaveform();
// 	// V-T
// 	for (Int_t j=0; j<nhfadc1; j++) {
// 	  double volt = hit->GetWaveformA(j);
// 	  double time = hit->GetWaveformT(j);
// 	  Int_t hid = RayrawHid + (seg+1)*100+2;;
// 	  HF2(hid, time, volt);
// 	}

// 	// TDC Leading
// 	Int_t nhit_tdc = hit->GetNumOfHit();
// 	if (nhit_tdc>0) {
// 	  double t1 = hit->GetT();
// 	  Int_t hid = RayrawHid + (seg+1)*100+3;;
// 	  HF1(hid, t1);
// 	}

// 	// Charge
// 	double charge1 = hit->GetQ();
// 	Int_t hid = RayrawHid + (seg+1)*100+4;      
// 	HF1(hid, charge1);

// 	// Charge w/ TDC
// 	if (nhit_tdc>0) {
// 	  Int_t hid = RayrawHid + (seg+1)*100+5;
// 	  HF1(hid, charge1);
// 	}

// 	// Charge of Baseline
// 	hid = RayrawHid + (seg+1)*100+6;
// 	double charge_ped = hit->GetQ_ped();
// 	HF1(hid, charge_ped);

// 	// dE
// 	double dE1 = hit->GetA();
// 	hid = RayrawHid + (seg+1)*100+7;
// 	HF1(hid, dE1);
//       }

//     }

//     // Int_t ref_sample_BH1[NumOfSegK11_BH1] = {
//     //   ref_sample1, ref_sample1, ref_sample1, ref_sample1, ref_sample1
//     // };

//     // Int_t ref_sample_T1[NumOfSegT1] = {
//     //   ref_sample1, ref_sample1, ref_sample2, ref_sample2, ref_sample2
//     // };

//     // Int_t ref_sample_T2[NumOfSegT2] = {
//     //   ref_sample2, ref_sample2, ref_sample2, ref_sample2, ref_sample2
//     // };

//     // Int_t ref_sample_V1743[NumOfSegV1743];
//     // for (Int_t i=0; i<NumOfSegV1743; i++) {
//     //   ref_sample_V1743[i]=0;
//     //  }

// #endif

    return true;
  }
}
//______________________________________________________________________________
bool
ProcessingEnd( void )
{
  tree->Fill();
  return true;
}

//______________________________________________________________________________
namespace
{
  const Int_t  NbinFADC_X     = 100;
  const Int_t  NbinFADC_Y     = 1024;
  const double MinFADC        = 0.;
  const double MaxFADC        = 1024;
  const double NbinTDC        = 10000;
  const double MinTDC         = 0;
  const double MaxTDC         = 10000;
  const double MinTime        = -20. * 1024;
  const double MaxTime        = 20. * 1024;
  const double MinPulseHeight = -1.5;
  const double MaxPulseHeight =  1.5;
  const double MaxEnergy      = 200.;
}

//______________________________________________________________________________
bool
ConfMan:: InitializeHistograms( void )
{
  HB1(  1, "Status",  20,   0., 20. );
  HB1( 10, "Trigger HitPat", NumOfSegTrig, 0., double(NumOfSegTrig) );
  for(Int_t i=0; i<NumOfSegTrig; ++i){
    HB1( 10+i+1, Form("Trigger Flag %d", i+1), 0x1000, 0, 0x1000 );
  }


  char buf[100];

  for (Int_t seg=0; seg<NumOfSegRayraw; seg++) {

    sprintf(buf, "RAYRAW - Raw Waveform (ch%d)", seg);
    Int_t hid = RayrawHid + (seg+1)*100 + 0;
    HB2( hid, buf, NbinFADC_X, 0, NbinFADC_X, NbinFADC_Y, MinFADC, MaxFADC );

    sprintf(buf, "RAYRAW - Waveform (ch%d)", seg);
    hid = RayrawHid + (seg+1)*100 + 1;
    HB2( hid, buf, NbinFADC_X, 0, NbinFADC_X, NbinFADC_Y, MinFADC, MaxFADC );

    sprintf(buf, "RAYRAW - Waveform w/ Leading (ch%d)", seg);
    hid = RayrawHid + (seg+1)*100 + 2;
    HB2( hid, buf, NbinFADC_X, 0, NbinFADC_X, NbinFADC_Y, MinFADC, MaxFADC );

    sprintf(buf, "RAYRAW - Waveform w/o Leading (ch%d)", seg);
    hid = RayrawHid + (seg+1)*100 + 3;
    HB2( hid, buf, NbinFADC_X, 0, NbinFADC_X, NbinFADC_Y, MinFADC, MaxFADC );

    sprintf(buf, "RAYRAW - Waveform w/ Trailing (ch%d)", seg);
    hid = RayrawHid + (seg+1)*100 + 4;
    HB2( hid, buf, NbinFADC_X, 0, NbinFADC_X, NbinFADC_Y, MinFADC, MaxFADC );

    sprintf(buf, "RAYRAW - Waveform w/o Trailing (ch%d)", seg);
    hid = RayrawHid + (seg+1)*100 + 5;
    HB2( hid, buf, NbinFADC_X, 0, NbinFADC_X, NbinFADC_Y, MinFADC, MaxFADC );

    sprintf(buf, "RAYRAW - Charge (ch%d)", seg);
    hid = RayrawHid + (seg+1)*100 + 6;
    HB1( hid, buf, NbinFADC_Y, MinFADC, MaxFADC);

    sprintf(buf, "RAYRAW - ADC (ch%d)", seg);
    hid = RayrawHid + (seg+1)*100 + 7;
    HB1( hid, buf, NbinFADC_Y, MinFADC, MaxFADC);

    sprintf(buf, "RAYRAW - ADC w/ Leading(ch%d)", seg);
    hid = RayrawHid + (seg+1)*100 + 8;
    HB1( hid, buf, NbinFADC_Y, MinFADC, MaxFADC);

    sprintf(buf, "RAYRAW - ADC w/o Leading(ch%d)", seg);
    hid = RayrawHid + (seg+1)*100 + 9;
    HB1( hid, buf, NbinFADC_Y, MinFADC, MaxFADC);

    sprintf(buf, "RAYRAW - TDC Leading (ch%d)", seg);
    hid = RayrawHid + (seg+1)*100 + 10;
    HB1( hid, buf, NbinTDC, MinTDC,  MaxTDC);

    sprintf(buf, "RAYRAW - TDC Trailing (ch%d)", seg);
    hid = RayrawHid + (seg+1)*100 + 11;
    HB1( hid, buf, NbinTDC, MinTDC,  MaxTDC);

    sprintf(buf, "RAYRAW - TOT (ch%d)", seg);
    hid = RayrawHid + (seg+1)*100 + 12;
    HB1( hid, buf, NbinTDC, MinTDC,  MaxTDC);

    // hid = RayrawHid + (seg+1)*100+2;
    // sprintf(buf, "RAYRAW - V-T (ch%d)", seg);
    // HB2( hid, buf, 2048, MinTime, MaxTime, 2000, MinPulseHeight, MaxPulseHeight );

    // hid = RayrawHid + (seg+1)*100+4;
    // sprintf(buf, "RAYRAW - Charge (ch%d)", seg);
    // HB1( hid, buf, 4000, -500, 9500 );

    // hid = RayrawHid + (seg+1)*100+5;
    // sprintf(buf, "RAYRAW - Charge w/ TDC (ch%d)", seg);
    // HB1( hid, buf, 4000, -500, 9500 );

    // hid = RayrawHid + (seg+1)*100+6;
    // sprintf(buf, "RAYRAW - BaseLine Charge (ch%d)", seg);
    // HB1( hid, buf, 4000, -500, 9500 );

    // hid = RayrawHid + (seg+1)*100+7;
    // sprintf(buf, "RAYRAW - dE (ch%d)", seg);
    // HB1( hid, buf, 200, 0, 5 );

  }

  //Tree
  HBTree( "tree", "tree" );
#if 0
  //Trig
  tree->Branch("evnum",     &event.evnum,     "evnum/I");

  tree->Branch("BH1Tdc_L",   event.BH1Tdc_L,   Form("BH1Tdc_L[%d]/I",  NumOfSegK11_BH1));
  tree->Branch("BH1Qdc_L",   event.BH1Qdc_L,   Form("BH1Qdc_L[%d]/D",  NumOfSegK11_BH1));
  tree->Branch("BH1Tdc_R",   event.BH1Tdc_R,   Form("BH1Tdc_R[%d]/I",  NumOfSegK11_BH1));
  tree->Branch("BH1Qdc_R",   event.BH1Qdc_R,   Form("BH1Qdc_R[%d]/D",  NumOfSegK11_BH1));
  tree->Branch("BH1Time_L",  event.BH1Time_L,  Form("BH1Time_L[%d]/D", NumOfSegK11_BH1));
  tree->Branch("BH1dE_L",    event.BH1dE_L,    Form("BH1dE_L[%d]/D",   NumOfSegK11_BH1));
  tree->Branch("BH1Time_R",  event.BH1Time_R,  Form("BH1Time_R[%d]/D", NumOfSegK11_BH1));
  tree->Branch("BH1dE_R",    event.BH1dE_R,    Form("BH1dE_R[%d]/D",   NumOfSegK11_BH1));
  tree->Branch("BH1MTime" ,  event.BH1MTime,   Form("BH1MTime[%d]/D",  NumOfSegK11_BH1));

  tree->Branch("T1Tdc_L",   event.T1Tdc_L,   Form("T1Tdc_L[%d]/I",  NumOfSegT1));
  tree->Branch("T1Qdc_L",   event.T1Qdc_L,   Form("T1Qdc_L[%d]/D",  NumOfSegT1));
  tree->Branch("T1Tdc_R",   event.T1Tdc_R,   Form("T1Tdc_R[%d]/I",  NumOfSegT1));
  tree->Branch("T1Qdc_R",   event.T1Qdc_R,   Form("T1Qdc_R[%d]/D",  NumOfSegT1));
  tree->Branch("T1Time_L",  event.T1Time_L,  Form("T1Time_L[%d]/D", NumOfSegT1));
  tree->Branch("T1dE_L",    event.T1dE_L,    Form("T1dE_L[%d]/D",   NumOfSegT1));
  tree->Branch("T1Time_R",  event.T1Time_R,  Form("T1Time_R[%d]/D", NumOfSegT1));
  tree->Branch("T1dE_R",    event.T1dE_R,    Form("T1dE_R[%d]/D",   NumOfSegT1));
  tree->Branch("T1MTime",   event.T1MTime,   Form("T1MTime[%d]/D",  NumOfSegT1));

  tree->Branch("T2Tdc",   event.T2Tdc,   Form("T2Tdc[%d]/I",  NumOfSegT2));
  tree->Branch("T2Qdc",   event.T2Qdc,   Form("T2Qdc[%d]/D",  NumOfSegT2));
  tree->Branch("T2Time",  event.T2Time,  Form("T2Time[%d]/D", NumOfSegT2));
  tree->Branch("T2dE",    event.T2dE,    Form("T2dE[%d]/D",   NumOfSegT2));

  tree->Branch("V1743Tdc",   event.V1743Tdc,   Form("V1743Tdc[%d]/I",  NumOfSegV1743));
  tree->Branch("V1743Qdc",   event.V1743Qdc,   Form("V1743Qdc[%d]/D",  NumOfSegV1743));
  tree->Branch("V1743Time",  event.V1743Time,  Form("V1743Time[%d]/D", NumOfSegV1743));
  tree->Branch("V1743dE",    event.V1743dE,    Form("V1743dE[%d]/D",   NumOfSegV1743));

  tree->Branch("RCTdc_L",   event.RCTdc_L,   Form("RCTdc_L[%d][16]/I",  NumOfSegE63RC));
  tree->Branch("RCAdc_L",   event.RCAdc_L,   Form("RCAdc_L[%d]/I",      NumOfSegE63RC));
  tree->Branch("RCTdc_R",   event.RCTdc_R,   Form("RCTdc_R[%d][16]/I",  NumOfSegE63RC));
  tree->Branch("RCAdc_R",   event.RCAdc_R,   Form("RCAdc_R[%d]/I",      NumOfSegE63RC));
#endif

  HPrint();
  return true;
}

//______________________________________________________________________________
bool
ConfMan::InitializeParameterFiles( void )
{
  return 
    ( //InitializeParameter<DCGeomMan>("DCGEO")    &&
      InitializeParameter<HodoParamMan>("HDPRM") &&
      InitializeParameter<HodoPHCMan>("HDPHC")   &&
      //InitializeParameter<BGOTemplateManager>("BGOTEMP") &&
      //InitializeParameter<BGOCalibMan>("BGOCALIB") &&
      InitializeParameter<UserParamMan>("USER")  );

}

//______________________________________________________________________________
bool
ConfMan::FinalizeProcess( void )
{
  return true;
}
