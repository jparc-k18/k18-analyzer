// -*- C++ -*-

#include <iostream>
#include <iterator>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include <TGFileBrowser.h>
#include <TH1.h>
#include <TH2.h>
#include <TMath.h>
#include <TStyle.h>

#include <DAQNode.hh>
#include <filesystem_util.hh>
#include <Unpacker.hh>
#include <UnpackerManager.hh>

#include "Controller.hh"
#include "user_analyzer.hh"

#include "ConfMan.hh"
#include "DCDriftParamMan.hh"
#include "DCGeomMan.hh"
#include "DCTdcCalibMan.hh"
#include "DetectorID.hh"
#include "GuiPs.hh"
#include "HistMaker.hh"
#include "HodoParamMan.hh"
#include "HodoPHCMan.hh"
#include "MacroBuilder.hh"
#include "MatrixParamMan.hh"
#include "MsTParamMan.hh"
#include "ProcInfo.hh"
#include "PsMaker.hh"
#include "SsdAnalyzer.hh"
#include "UserParamMan.hh"

#define DEBUG      0
#define FLAG_DAQ   1
#define TIME_STAMP 0

namespace
{
using hddaq::unpacker::GUnpacker;
using hddaq::unpacker::DAQNode;
const auto& gUnpacker = GUnpacker::get_instance();
      auto& gHist     = HistMaker::getInstance();
const auto& gMatrix   = MatrixParamMan::GetInstance();
const auto& gUser     = UserParamMan::GetInstance();
std::vector<TH1*> hptr_array;
Bool_t flag_event_cut = false;
Int_t event_cut_factor = 1; // for fast semi-online analysis
}

namespace analyzer
{
//____________________________________________________________________________
int
process_begin(const std::vector<std::string>& argv)
{
  ConfMan& gConfMan = ConfMan::GetInstance();
  gConfMan.Initialize(argv);
  //  GUnpacker::get_instance().set_decode_mode(false);
  gConfMan.InitializeParameter<UserParamMan>("USER");
  if( !gConfMan.IsGood() ) return -1;
  // unpacker and all the parameter managers are initialized at this stage

  if( argv.size()==4 ){
    Int_t factor = std::strtod( argv[3].c_str(), NULL );
    if( factor!=0 ) event_cut_factor = std::abs( factor );
    flag_event_cut = true;
    std::cout << "#D Event cut flag on : factor="
	      << event_cut_factor << std::endl;
  }

  // Make tabs
  hddaq::gui::Controller& gCon = hddaq::gui::Controller::getInstance();
  TGFileBrowser *tab_hist  = gCon.makeFileBrowser("Hist");
  TGFileBrowser *tab_macro = gCon.makeFileBrowser("Macro");

  // Add macros to the Macro tab
  tab_macro->Add(macro::Get("clear_all_canvas"));
  tab_macro->Add(macro::Get("clear_canvas"));
  //tab_macro->Add(macro::Get("split22"));
  //tab_macro->Add(macro::Get("dispBGOTdc"));
  //tab_macro->Add(macro::Get("BGOhitpat"));
  //tab_macro->Add(macro::Get("split32"));
  //tab_macro->Add(macro::Get("split33"));
  tab_macro->Add(macro::Get("dispGeAdc"));
  tab_macro->Add(macro::Get("dispSplitGeAdC"));
  tab_macro->Add(macro::Get("dispGeAdc_wTrigFlag"));
  tab_macro->Add(macro::Get("dispGeTdc"));
  tab_macro->Add(macro::Get("dispGeHitpat"));
  tab_macro->Add(macro::Get("dispGeTdcDepth"));
  tab_macro->Add(macro::Get("dispGe2dhist"));
  tab_macro->Add(macro::Get("dispBGOTdc"));
  tab_macro->Add(macro::Get("dispBGOHitpat"));
  tab_macro->Add(macro::Get("dispHBXTrigFlag"));
  tab_macro->Add(macro::Get("dispGeBGO_2D"));
  tab_macro->Add(macro::Get("dispHBXScaler"));
  //tab_macro->Add(macro::Get("dispEventDisplay"));
  
  //tab_macro->Add(macro::Get("dispGe2dhist"));
  //tab_macro->Add(macro::Get("dispGeAdc_60Co"));
  //tab_macro->Add(macro::Get("dispGeAdc_LSO"));
  //tab_macro->Add(macro::Get("dispGeAdc_LSO_off"));
  //tab_macro->Add(macro::Get("dispBGOTDC"));
  //tab_macro->Add(macro::Get("BGO hitpat"));


  // Add histograms to the Hist tab
  tab_hist->Add(gHist.createGe());
  tab_hist->Add(gHist.createBGO());
  tab_hist->Add(gHist.createTriggerFlag());

  // Set histogram pointers to the vector sequentially.
  // This vector contains both TH1 and TH2.
  // Then you need to do down cast when you use TH2.
  if(0 != gHist.setHistPtr(hptr_array)){ return -1; }

  gStyle->SetOptStat(1110);
  gStyle->SetTitleW(.400);
  gStyle->SetTitleH(.100);
  // gStyle->SetStatW(.420);
  // gStyle->SetStatH(.350);
  gStyle->SetStatW(.320);
  gStyle->SetStatH(.250);

  return 0;
}

//____________________________________________________________________________
int
process_end()
{
  hptr_array.clear();
  return 0;
}

//____________________________________________________________________________
int
process_event()
{
  //std::cout<<"debag1"<<std::endl;
  // TriggerFlag ---------------------------------------------------
  std::bitset<NumOfSegTFlag> trigger_flag;
  {
    static const Int_t k_device = gUnpacker.get_device_id("HBXTFlag");
    static const Int_t k_tdc    = gUnpacker.get_data_id("HBXTFlag", "tdc");
    //static const Int_t tdc_id   = gHist.getSequentialID( kTriggerFlag, 0, kTDC );
    //static const Int_t hit_id   = gHist.getSequentialID( kTriggerFlag, 0, kHitPat );
    static const Int_t tdc_id   = gHist.getSequentialID( kGe, 0, kFlagTDC );
    static const Int_t hit_id   = gHist.getSequentialID( kGe, 0, kFlagHitPat );
    //std::cout<<123<<std::endl;
    for( Int_t seg=0; seg<NumOfSegTFlag; ++seg ){
      for( Int_t m=0, n=gUnpacker.get_entries( k_device, 0, seg, 0, k_tdc );
           m<n; ++m ){
	auto tdc = gUnpacker.get( k_device, 0, seg, 0, k_tdc, m );
	if( tdc>0 ){
	  trigger_flag.set( seg );
	  hptr_array[tdc_id+seg]->Fill( tdc );
	}
      }
      if( trigger_flag[seg] ) hptr_array[hit_id]->Fill( seg );
    }
    if( !( trigger_flag[trigger::kSpillEnd] |
	   trigger_flag[trigger::kLevel1OR] ) |
        !( trigger_flag[trigger::kL1SpillOn] |
	   trigger_flag[trigger::kL1SpillOff] |
           trigger_flag[trigger::kSpillEnd] ) ){
      //hddaq::cerr << "#W Trigger flag is missing!!! "
      // << trigger_flag << std::endl;
    }
#if 0
    gUnpacker.dump_data_device(k_device);
#endif
    //std::cout<<"debag2"<<std::endl;
  }

#if DEBUG
  std::cout << __FILE__ << " " << __LINE__ << std::endl;
#endif

  //------------------------------------------------------------------
  // Hyperball-X'
  //------------------------------------------------------------------

  // Ge --------------------------------------------------------------
  {
#if 0
    static const Int_t tfa_min = gUser.GetParameter("TFA_TDC", 0);
    static const Int_t tfa_max = gUser.GetParameter("TFA_TDC", 1);
    static const Int_t crm_min = gUser.GetParameter("CRM_TDC", 0);
    static const Int_t crm_max = gUser.GetParameter("CRM_TDC", 1);
#endif
    static const Int_t tfa_min = 3500;
    static const Int_t tfa_max = 5000;
    static const Int_t crm_min = 3500;
    static const Int_t crm_max = 5000;
    // data type
    static const Int_t k_device = gUnpacker.get_device_id("Ge");
    static const Int_t k_adc    = gUnpacker.get_data_id("Ge","adc");
    static const Int_t k_tfa    = gUnpacker.get_data_id("Ge","tfa_leading");
    static const Int_t k_crm    = gUnpacker.get_data_id("Ge","crm_leading");
    static const Int_t k_rst    = gUnpacker.get_data_id("Ge","reset_time");

    static const Int_t ge_adc_id    = gHist.getSequentialID(kGe, 0, kADC);
    static const Int_t ge_adc_wt_id = gHist.getSequentialID(kGe, 0, kADCwTDC);
    static const Int_t ge_adc_wc_id = gHist.getSequentialID(kGe, 0, kADCwTDC, NumOfSegGe+1);
    static const Int_t ge_adc_lso_on_id = gHist.getSequentialID(kGe, 0, kADCwTDC, NumOfSegGe*2+1);
    static const Int_t ge_adc_lso_off_id = gHist.getSequentialID(kGe, 0, kADCwTDC, NumOfSegGe*3+1);
    static const Int_t ge_adc_gecoin_on_id = gHist.getSequentialID(kGe, 0, kADCwTDC, NumOfSegGe*4+1);
    static const Int_t ge_adc_gecoin_off_id = gHist.getSequentialID(kGe, 0, kADCwTDC, NumOfSegGe*5+1);
    static const Int_t ge_tdc_depth_id = gHist.getSequentialID(kGe, 0, kADCwTDC, NumOfSegGe*6+1);
    static const Int_t ge_adc_wL1trig_id = gHist.getSequentialID(kGe, 0, kADCwTDC, NumOfSegGe*7+1);
    static const Int_t ge_adc_wLSOtrig_id = gHist.getSequentialID(kGe, 0, kADCwTDC, NumOfSegGe*8+1);
    static const Int_t ge_scaler_id = gHist.getSequentialID(kGe, 0, kADCwTDC, NumOfSegGe*9+1);
    
    static const Int_t ge_tfa_id    = gHist.getSequentialID(kGe, 0, kTFA);
    static const Int_t ge_tfa_depth_id    = gHist.getSequentialID(kGe, 0, kTFA, NumOfSegGe+1);
    static const Int_t ge_crm_id    = gHist.getSequentialID(kGe, 0, kCRM);
    static const Int_t ge_rst_id    = gHist.getSequentialID(kGe, 0, kRST);

    // hitpat hist id 
    static const Int_t ge_hitpat_id = gHist.getSequentialID(kGe, 0, kHitPat);//671ADC
    static const Int_t ge_hitpat_id2 = gHist.getSequentialID(kGe, 0, kHitPat, 1+1);//973ADC
    //static const Int_t ge_hitpat_id3 = gHist.getSequentialID(kGe, 0, kHitPat, 1+2);//Ge vs BGO HitPat
    // 2d hist id
    static const Int_t ge_adc2d_id = gHist.getSequentialID(kGe, 0, kADC2D);
    static const Int_t ge_tfa2d_id = gHist.getSequentialID(kGe, 0, kTFA2D);
    static const Int_t ge_crm2d_id = gHist.getSequentialID(kGe, 0, kCRM2D);
    static const Int_t ge_rst2d_id = gHist.getSequentialID(kGe, 0, kRST2D);

    static const Int_t ge_tfa_adc_id = gHist.getSequentialID(kGe, 0, kTFA_ADC);
    static const Int_t ge_crm_adc_id = gHist.getSequentialID(kGe, 0, kCRM_ADC);
    static const Int_t ge_rst_adc_id = gHist.getSequentialID(kGe, 0, kRST_ADC);
    static const Int_t ge_tfa_crm_id = gHist.getSequentialID(kGe, 0, kTFA_CRM);

    //event slip flag//
    Int_t flag_eventslip = 0;
    ////////////////////////////
    
    for(Int_t seg = 0; seg<NumOfSegGe; ++seg){
      // ADC
      Int_t nhit_adc = gUnpacker.get_entries(k_device, 0, seg, 0, k_adc);
      Int_t adc = -9999;
      //4-1seg dead -> no event slip
      if(seg==12)flag_eventslip = 0;
      ///////////////////////////
      if(nhit_adc != 0){
	adc = gUnpacker.get(k_device, 0, seg, 0, k_adc);
	hptr_array[ge_adc_id + seg]->Fill(adc);
	if(trigger_flag[1]>0){
	  //hptr_array[ge_adc_wL1trig_id + seg]->Fill(adc);
	}
	if(trigger_flag[2]>0){
	  //hptr_array[ge_adc_wLSOtrig_id + seg]->Fill(adc);
	}
	hptr_array[ge_adc2d_id]->Fill(seg, adc);

	if(115 < adc && adc < 7500){
	  if(0<=seg&&seg<=15)hptr_array[ge_hitpat_id]->Fill(seg);
	  if(16<=seg&&seg<=31)hptr_array[ge_hitpat_id2]->Fill(seg-NumOfSegGe/2);
	}
      }

      // TFA
      Int_t nhit_tfa = gUnpacker.get_entries(k_device, 0, seg, 0, k_tfa);
      Int_t tfa_first = -9999;
      Int_t tfaflag = 0;
      Int_t flag_973 = 0;
      if(nhit_tfa != 0){
	tfa_first = gUnpacker.get(k_device, 0, seg, 0, k_tfa, 0);
	if(tfa_min<tfa_first&&tfa_first<tfa_max){
	  flag_973 = 1;
	}
	if(adc >= 0) hptr_array[ge_tfa_adc_id + seg]->Fill(tfa_first, adc);
	for(Int_t m = 0; m<nhit_tfa; ++m){
	  Int_t tfa = gUnpacker.get(k_device, 0, seg, 0, k_tfa, m);
	  //detect event slip//
	  if(tfa==-9999&&adc>0&&seg!=12){
	    flag_eventslip = 1;
	  }
	  else{
	    flag_eventslip = 0;
	  }
	  //std::cout<<"seg "<<seg<<" tfa "<<tfa<<" adc "<<adc<<" event slip flag = "<<flag_eventslip<<std::endl;
	  if(flag_eventslip == 1){
	    //std::cout<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<std::endl;
	    //std::cout<<"!!!!!!!!!!!!!!!!Event Slip Occur!!!!!!!!!!!!!!!"<<std::endl;
	    //std::cout<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<std::endl;
	    //getchar();
	  }
	  if(flag_eventslip == 0){
	    //std::cout<<"Event Slip does not occur"<<std::endl;
	  }
	  //fill tfa 
	  hptr_array[ge_tfa_id + seg]->Fill(tfa);
	  //fill depth
	  hptr_array[ge_tdc_depth_id + seg]->Fill(m);
	  hptr_array[ge_tfa2d_id]->Fill(seg, tfa);
	  if( tfa_min < tfa && tfa < tfa_max ) tfaflag = 1;
	}
	//ADCwTFA
	if(flag_973 == 1){
	  hptr_array[ge_adc_wt_id + seg]->Fill(adc);
	  //if(trigger_flag[trigger::kTrigEPS]){
	  //  if (trigger_flag[trigger::kL1SpillOn]) hptr_array[ge_adc_lso_on_id + seg]->Fill(adc);
	  //  if (trigger_flag[trigger::kL1SpillOff]) hptr_array[ge_adc_lso_off_id + seg]->Fill(adc);
	  //}
	  //if(trigger_flag[trigger::kTrigFPS]){
	  //  if (trigger_flag[trigger::kL1SpillOn]) hptr_array[ge_adc_gecoin_on_id + seg]->Fill(adc);
	  //  if (trigger_flag[trigger::kL1SpillOff]) hptr_array[ge_adc_gecoin_off_id + seg]->Fill(adc);
	  //}
	}
      }//nhit tfa
      
      // CRM
      Int_t nhit_crm = gUnpacker.get_entries(k_device, 0, seg, 0, k_crm);
      Int_t crmflag = 0;
      if(nhit_crm != 0){
	Int_t crm_first = gUnpacker.get(k_device, 0, seg, 0, k_crm, 0);
	if(adc >= 0) hptr_array[ge_crm_adc_id + seg]->Fill(crm_first, adc);
	if(tfa_first > 0) hptr_array[ge_tfa_crm_id + seg]->Fill(tfa_first, crm_first);
	for(Int_t m = 0; m<nhit_crm; ++m){
	  Int_t crm = gUnpacker.get(k_device, 0, seg, 0, k_crm, m);
	  hptr_array[ge_crm_id + seg]->Fill(crm);
	  hptr_array[ge_crm2d_id]->Fill(seg, crm);
	  if( crm_min < crm && crm < crm_max ) crmflag = 1;
	}
	//ADCwCRM
	if(crmflag == 1){
	  hptr_array[ge_adc_wc_id + seg]->Fill(adc);
	  if(trigger_flag[trigger::kTrigEPS]){
	    if (trigger_flag[trigger::kL1SpillOn]) hptr_array[ge_adc_lso_on_id + seg]->Fill(adc);
	    if (trigger_flag[trigger::kL1SpillOff]) hptr_array[ge_adc_lso_off_id + seg]->Fill(adc);
	  }
	  if(trigger_flag[trigger::kTrigFPS]){
	    if (trigger_flag[trigger::kL1SpillOn]) hptr_array[ge_adc_gecoin_on_id + seg]->Fill(adc);
	    if (trigger_flag[trigger::kL1SpillOff]) hptr_array[ge_adc_gecoin_off_id + seg]->Fill(adc);
	  }
	}
      }

      // RST
      Int_t nhit_rst = gUnpacker.get_entries(k_device, 0, seg, 0, k_rst);
      if(nhit_rst != 0){
	  Int_t rst = gUnpacker.get(k_device, 0, seg, 0, k_rst);
	  hptr_array[ge_rst_id + seg]->Fill(rst);
	  hptr_array[ge_rst2d_id]->Fill(seg, rst);
	  if(adc >= 0) hptr_array[ge_rst_adc_id + seg]->Fill(rst, adc);
      }
    }
#if 0
    // Debug, dump data relating this detector
    gUnpacker.dump_data_device(k_device);
#endif
    #if 0
    if(flag_eventslip==0){
      std::cout<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<std::endl;
      std::cout<<"!!!!!!!!!!!!Event slip ccurs!!!!!!!!!!!!!"<<std::endl;
      std::cout<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<std::endl;
      getchar();
    }
    if(flag_eventslip==1)std::cout<<"Event slip does not occur"<<std::endl;
    #endif
  }//Ge


  // BGO --------------------------------------------------------------
  {
    // data typep
    static const Int_t k_device = gUnpacker.get_device_id("BGO");
    static const Int_t k_tdc    = gUnpacker.get_data_id("BGO","leading");
    static const Int_t k_device_ge = gUnpacker.get_device_id("Ge");
    static const Int_t k_adc    = gUnpacker.get_data_id("Ge","adc");
    // sequential id
    // sequential hist
    static const Int_t bgo_tdc_id    = gHist.getSequentialID(kBGO, 0, kTDC);
    static const Int_t bgo_tdc2d_id  = gHist.getSequentialID(kBGO, 0, kTDC2D);
    static const Int_t bgo_hit_id    = gHist.getSequentialID(kBGO, 0, kHitPat);
    static const Int_t ge_hitpat_id3 = gHist.getSequentialID(kGe, 0, kHitPat, 1+1+1);//Ge vs BGO HitPat
    
    for(Int_t seg = 0; seg<NumOfSegBGO; ++seg){

      TH2* htdc2d = dynamic_cast<TH2*>(hptr_array[bgo_tdc2d_id]);

      // TDC
      Int_t nhit_tdc = gUnpacker.get_entries(k_device, 0, seg, 0, k_tdc);
      for(Int_t m = 0; m<nhit_tdc; ++m){
	Int_t tdc = gUnpacker.get(k_device, 0, seg, 0, k_tdc, m);
	Int_t tdc_first = gUnpacker.get(k_device, 0, seg, 0, k_tdc, 0);
	//hptr_array[bgo_tdc_id + seg]->Fill(tdc_first);
	hptr_array[bgo_tdc_id + seg]->Fill(tdc_first);
	htdc2d->Fill(seg, tdc);
      }
      
      // HitPat
      if(nhit_tdc != 0){
	//BGO HitPat
	hptr_array[bgo_hit_id]->Fill(seg);
	//BGO and Ge HitPat
	for(int seg_ge = 0; seg_ge<NumOfSegGe/2; seg_ge++){
	  #if 1
	  Int_t nhit_adc = gUnpacker.get_entries(k_device, 0, seg_ge, 0, k_adc);
	  if(nhit_adc>0){
	    //hptr_array[ge_hitpat_id3]->Fill(seg,seg_ge);
	    //std::cout<<"hit seg bgo ge "<<std::endl;
	    hptr_array[ge_hitpat_id3]->Fill(seg_ge,seg);
	  }
	  #endif
	}
      }
      
    }

#if 0
    // Debug, dump data relating this detector
    gUnpacker.dump_data_device(k_device);
#endif
  }

#if DEBUG
  std::cout << __FILE__ << " " << __LINE__ << std::endl;
#endif


  return 0;
} //process_event()

}//analyzer
