/**
 *  file: UserScaler.cc
 *  date: 2017.04.10
 *
 */

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

#include <filesystem_util.hh>
#include <lexical_cast.hh>

#include "ConfMan.hh"
#include "DetectorID.hh"
#include "FuncName.hh"
#include "RootHelper.hh"
#include "ScalerAnalyzer.hh"
#include "Unpacker.hh"
#include "UnpackerManager.hh"
#include "VEvent.hh"

#define USE_COMMA   0
#define SPILL_RESET 0
#define MAKE_LOG    1

namespace
{
  using namespace root;
  using namespace hddaq::unpacker;
  const TString& class_name("EventScaler");
  ScalerAnalyzer&  gScaler   = ScalerAnalyzer::GetInstance();
  UnpackerManager& gUnpacker = GUnpacker::get_instance();
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
class EventScaler : public VEvent
{
public:
          EventScaler( void );
         ~EventScaler( void );
  Bool_t  ProcessingBegin( void );
  Bool_t  ProcessingEnd( void );
  Bool_t  ProcessingNormal( void );
  Bool_t  InitializeHistograms( void );
};

//______________________________________________________________________________
struct Event
{
  Int_t evnum;
};

//______________________________________________________________________________
namespace root
{
  Event  event;
  TH1   *h[MaxHist];
  TTree *tree;
}

//______________________________________________________________________________
EventScaler::EventScaler( void )
  : VEvent()
{
}

//______________________________________________________________________________
EventScaler::~EventScaler( void )
{
}

//______________________________________________________________________________
Bool_t
EventScaler::ProcessingBegin( void )
{
  return true;
}

//______________________________________________________________________________
Bool_t
EventScaler::ProcessingNormal( void )
{
  event.evnum++;
  gScaler.Decode();

  //if( event.evnum%400==0 )
  //  gScaler.Print();
  
  if( event.evnum%100000==0 )
    gScaler.Print();
  
#if SPILL_RESET
  if( gScaler.SpillIncrement() )
    gScaler.Clear();
#endif

  return true;
}

//______________________________________________________________________________
Bool_t
EventScaler::ProcessingEnd( void )
{
  return true;
}

//______________________________________________________________________________
VEvent*
ConfMan::EventAllocator( void )
{
  return new EventScaler;
}

//______________________________________________________________________________
Bool_t
ConfMan::InitializeHistograms( void )
{
  HBTree( "scaler", "tree of Scaler" );
  event.evnum = 0;

  //////////////////// Set Channels
  // ScalerAnalylzer::Set( Int_t column,
  //                       Int_t raw,
  //                       ScalerInfo( name, module, channel ) );
  // scaler information is defined from here.
  // please do not use a white space character.
  {
    Int_t c = ScalerAnalyzer::kLeft;
    Int_t r = 0;
    gScaler.Set( c, r++, ScalerInfo( "BH1",        0, 16 ) );
    gScaler.Set( c, r++, ScalerInfo( "BH1-SUM",   -1, -1 ) );
    gScaler.Set( c, r++, ScalerInfo( "BH1-01",     1,  0 ) );
    gScaler.Set( c, r++, ScalerInfo( "BH1-02",     1,  1 ) );
    gScaler.Set( c, r++, ScalerInfo( "BH1-03",     1,  2 ) );
    gScaler.Set( c, r++, ScalerInfo( "BH1-04",     1,  3 ) );
    gScaler.Set( c, r++, ScalerInfo( "BH1-05",     1,  4 ) );
    gScaler.Set( c, r++, ScalerInfo( "BH1-06",     1,  5 ) );
    gScaler.Set( c, r++, ScalerInfo( "BH1-07",     1,  6 ) );
    gScaler.Set( c, r++, ScalerInfo( "BH1-08",     1,  7 ) );
    gScaler.Set( c, r++, ScalerInfo( "BH1-09",     1,  8 ) );
    gScaler.Set( c, r++, ScalerInfo( "BH1-10",     1,  9 ) );
    gScaler.Set( c, r++, ScalerInfo( "BH1-11",     1, 10 ) );
    gScaler.Set( c, r++, ScalerInfo( "BH2",        0, 17 ) );
    gScaler.Set( c, r++, ScalerInfo( "BH2-SUM",   -1, -1 ) );
    gScaler.Set( c, r++, ScalerInfo( "BH2-01",     0, 64 ) );
    gScaler.Set( c, r++, ScalerInfo( "BH2-02",     0, 65 ) );
    gScaler.Set( c, r++, ScalerInfo( "BH2-03",     0, 66 ) );
    gScaler.Set( c, r++, ScalerInfo( "BH2-04",     0, 67 ) );
    gScaler.Set( c, r++, ScalerInfo( "BH2-05",     0, 68 ) );
    gScaler.Set( c, r++, ScalerInfo( "BAC",        0, 18 ) );
    gScaler.Set( c, r++, ScalerInfo( "PVAC",       0, 20 ) );
    gScaler.Set( c, r++, ScalerInfo( "FAC",        0, 21 ) );
    gScaler.Set( c, r++, ScalerInfo( "SCH",        0, 11 ) );
    gScaler.Set( c, r++, ScalerInfo( "TOF",        0, 22 ) );
    gScaler.Set( c, r++, ScalerInfo( "LAC",        0, 23 ) );
    gScaler.Set( c, r++, ScalerInfo( "WC",         0, 24 ) );
    for( Int_t i=0; i<16; ++i ){
      gScaler.Set( c, r++, ScalerInfo( Form( "TFA-%02d", i+1 ), 2, i ) );
    }
  }

  {
    Int_t c = ScalerAnalyzer::kCenter;
    Int_t r = 0;
    gScaler.Set( c, r++, ScalerInfo( "10M-Clock",    0,  0 ) );
    gScaler.Set( c, r++, ScalerInfo( "BH1-1/100-PS", 1, 11 ) );
    gScaler.Set( c, r++, ScalerInfo( "BH1-1/1e5-PS", 1, 12 ) );
    gScaler.Set( c, r++, ScalerInfo( "TOF-24",       0, 29 ) );
    gScaler.Set( c, r++, ScalerInfo( "Mtx2D-1",      0, 32 ) );
    gScaler.Set( c, r++, ScalerInfo( "Mtx2D-2",      0, 33 ) );
    gScaler.Set( c, r++, ScalerInfo( "Mtx3D",        0, 34 ) );
    gScaler.Set( c, r++, ScalerInfo( "Other1",       0, 25 ) );
    gScaler.Set( c, r++, ScalerInfo( "Other2",       0, 26 ) );
    gScaler.Set( c, r++, ScalerInfo( "Other3",       0, 27 ) );
    gScaler.Set( c, r++, ScalerInfo( "Other4",       0, 28 ) );
    gScaler.Set( c, r++, ScalerInfo( "LSO1",         2, 48 ) );
    gScaler.Set( c, r++, ScalerInfo( "LSO2",         2, 49 ) );
    gScaler.Set( c, r++, ScalerInfo( "LSO1xGe13",    2, 50 ) );
    gScaler.Set( c, r++, ScalerInfo( "LSO2xGe24",    2, 51 ) );
    gScaler.Set( c, r++, ScalerInfo( "GeCoin1",      2, 52 ) );
    gScaler.Set( c, r++, ScalerInfo( "GeCoin2",      2, 53 ) );
    gScaler.Set( c, r++, ScalerInfo( "LSO1High",     2, 54 ) );
    gScaler.Set( c, r++, ScalerInfo( "LSO2High",     2, 55 ) );
    gScaler.Set( c, r++, ScalerInfo( "GeHigh1",      2, 56 ) );
    gScaler.Set( c, r++, ScalerInfo( "GeHigh2",      2, 57 ) );
    gScaler.Set( c, r++, ScalerInfo( "BEAM-A",       0, 35 ) );
    gScaler.Set( c, r++, ScalerInfo( "BEAM-B",       0, 36 ) );
    gScaler.Set( c, r++, ScalerInfo( "BEAM-C",       0, 37 ) );
    gScaler.Set( c, r++, ScalerInfo( "BEAM-D",       0, 38 ) );
    gScaler.Set( c, r++, ScalerInfo( "BEAM-E",       0, 39 ) );
    gScaler.Set( c, r++, ScalerInfo( "BEAM-F",       0, 40 ) );
    for( Int_t i=0; i<16; ++i ){
      gScaler.Set( c, r++, ScalerInfo( Form( "CRM-%02d", i+1 ), 2, i+16 ) );
    }
  }

  {
    Int_t c = ScalerAnalyzer::kRight;
    Int_t r = 0;
    gScaler.Set( c, r++, ScalerInfo( "Spill",        -1, -1 ) );
    gScaler.Set( c, r++, ScalerInfo( "TM",            0,  9 ) );
    gScaler.Set( c, r++, ScalerInfo( "SY",            0, 10 ) );
    gScaler.Set( c, r++, ScalerInfo( "Real-Time",     0,  1 ) );
    gScaler.Set( c, r++, ScalerInfo( "Live-Time",     0,  2 ) );
    gScaler.Set( c, r++, ScalerInfo( "L1-Req",        0,  3 ) );
    gScaler.Set( c, r++, ScalerInfo( "L1-Acc",        0,  4 ) );
    // gScaler.Set( c, r++, ScalerInfo( "MstClr",        0,  5 ) );
    gScaler.Set( c, r++, ScalerInfo( "Clear",         0,  6 ) );
    gScaler.Set( c, r++, ScalerInfo( "L2-Req",        0,  7 ) );
    gScaler.Set( c, r++, ScalerInfo( "L2-Acc",        0,  8 ) );
    gScaler.Set( c, r++, ScalerInfo( "TRIG-A",        0, 41 ) );
    gScaler.Set( c, r++, ScalerInfo( "TRIG-B",        0, 42 ) );
    gScaler.Set( c, r++, ScalerInfo( "TRIG-C",        0, 43 ) );
    gScaler.Set( c, r++, ScalerInfo( "TRIG-D",        0, 44 ) );
    gScaler.Set( c, r++, ScalerInfo( "TRIG-E",        0, 45 ) );
    gScaler.Set( c, r++, ScalerInfo( "TRIG-F",        0, 46 ) );
    gScaler.Set( c, r++, ScalerInfo( "TRIG-A-PS",     0, 48 ) );
    gScaler.Set( c, r++, ScalerInfo( "TRIG-B-PS",     0, 49 ) );
    gScaler.Set( c, r++, ScalerInfo( "TRIG-C-PS",     0, 50 ) );
    gScaler.Set( c, r++, ScalerInfo( "TRIG-D-PS",     0, 51 ) );
    gScaler.Set( c, r++, ScalerInfo( "TRIG-E-PS",     0, 52 ) );
    gScaler.Set( c, r++, ScalerInfo( "TRIG-F-PS",     0, 53 ) );
    gScaler.Set( c, r++, ScalerInfo( "TRIG-PSOR-A",   0, 54 ) );
    gScaler.Set( c, r++, ScalerInfo( "TRIG-PSOR-B",   0, 55 ) );
    gScaler.Set( c, r++, ScalerInfo( "Clock-PS",      0, 56 ) );
    gScaler.Set( c, r++, ScalerInfo( "Reserve2-PS",   0, 57 ) );
    gScaler.Set( c, r++, ScalerInfo( "Level1-PS",     0, 58 ) );
    for( Int_t i=0; i<16; ++i ){
      gScaler.Set( c, r++, ScalerInfo( Form( "Reset-%02d", i+1 ), 2, i+32 ) );
    }
  }

#if USE_COMMA
  gScaler.SetFlag( ScalerAnalyzer::kSeparateComma );
#endif

  gScaler.PrintFlags();

  return true;
}

//______________________________________________________________________________
Bool_t
ConfMan::InitializeParameterFiles( void )
{
  return true;
}

//______________________________________________________________________________
Bool_t
ConfMan::FinalizeProcess( void )
{
  if( event.evnum==0 ) return true;

  const Int_t run_number = gUnpacker.get_root()->get_run_number();
  gScaler.Print();

#if MAKE_LOG
  const TString& bin_dir( hddaq::dirname(hddaq::selfpath()) );
  const TString& data_dir( hddaq::dirname(gUnpacker.get_istream()) );

  std::stringstream run_number_ss; run_number_ss << run_number;
  const TString& recorder_log( data_dir+"/recorder.log" );
  std::ifstream ifs( recorder_log );
  if( !ifs.is_open() ){
    std::cerr << "#E " << FUNC_NAME << " "
	      << "cannot open recorder.log : "
	      << recorder_log << std::endl;
    return false;
  }

  const TString& scaler_dir( bin_dir+"/../scaler" );
  const TString& scaler_txt = Form( "%s/scaler_%05d.txt",
				    scaler_dir.Data(), run_number );

  std::ofstream ofs( scaler_txt );
  if( !ofs.is_open() ){
    std::cerr << "#E " << FUNC_NAME << " "
	      << "cannot open scaler.txt : "
	      << scaler_txt << std::endl;
    return false;
  }

  Int_t recorder_event_number = 0;
  Bool_t found_run_number = false;
  std::string line;
  while( ifs.good() && std::getline(ifs,line) ){
    if( line.empty() ) continue;
    std::istringstream input_line( line );
    std::istream_iterator<std::string> line_begin( input_line );
    std::istream_iterator<std::string> line_end;
    std::vector<std::string> log_column( line_begin, line_end );
    if( log_column.at(0) != "RUN" ) continue;
    if( log_column.at(1) != run_number_ss.str() ) continue;
    recorder_event_number = hddaq::a2i( log_column.at(15) );
    ofs << line << std::endl;
    found_run_number = true;
  }

  if( !found_run_number ){
    std::cerr << "#E " << FUNC_NAME << " "
	      << "not found run# " << run_number
	      << " in " << recorder_log << std::endl;
    return false;
  }

  ofs << std::endl;
  ofs << std::left  << std::setw(15) << "" << "\t"
      << std::right << std::setw(15) << "Integral" << std::endl;
  ofs << std::left  << std::setw(15) << "Event"    << "\t"
      << std::right << std::setw(15) << event.evnum << std::endl;

  if( recorder_event_number != event.evnum ){
    std::cerr << "#W " << FUNC_NAME << " "
	      << "event number mismatch" << std::endl
	      << "   recorder : " << recorder_event_number << std::endl
	      << "   decode   : " << event.evnum << std::endl;
  }

  {
    std::vector<Int_t> order = {
      ScalerAnalyzer::kRight,
      ScalerAnalyzer::kLeft,
      ScalerAnalyzer::kCenter
    };
    for( auto&& c : order ){
      for(Int_t i=0; i<ScalerAnalyzer::MaxRow; i++){
	TString name = gScaler.GetScalerName( c, i );
	if( name=="n/a" ) continue;
	ofs << std::left  << std::setw(15) << name << "\t"
	    << std::right << std::setw(15) << gScaler.Get(name) << std::endl;
      }
      ofs << std::endl;
    }
  }

  Double_t reallive = gScaler.Fraction("Live-Time", "Real-Time");
  Double_t daqeff   = gScaler.Fraction("L1-Acc", "L1-Req");
  Double_t l2eff    = gScaler.Fraction("L2-Acc", "L1-Acc");
  Double_t pitm     = gScaler.Fraction("pi-Beam", "TM");
  Double_t ptm      = gScaler.Fraction("p-Beam", "TM");
  Double_t bh2kbeam = gScaler.Fraction("(BH2,K)", "Beam");
  Double_t pbgobeam = gScaler.Fraction("(p,BGO)", "Beam");
  Double_t pirate   = gScaler.Fraction("pi-Beam", "Spill");
  Double_t prate    = gScaler.Fraction("p-Beam", "Spill");
  Double_t l1rate   = gScaler.Fraction("L1-Req", "Spill");
  Double_t l2rate   = gScaler.Fraction("L2-Acc", "Spill");

  ofs << std::fixed << std::setprecision(6)
      << std::left  << std::setw(18) << "Live/Real"     << "\t"
      << std::right << std::setw(12)<<  reallive        << std::endl
      << std::left  << std::setw(18) << "DAQ-Eff"       << "\t"
      << std::right << std::setw(12) << daqeff          << std::endl
      << std::left  << std::setw(18) << "L2-Eff"        << "\t"
      << std::right << std::setw(12) << l2eff           << std::endl
      << std::left  << std::setw(18) << "Duty-Factor"   << "\t"
      << std::right << std::setw(12) << gScaler.Duty()  << std::endl
      << std::left  << std::setw(18) << "pi-Beam/TM"    << "\t"
      << std::right << std::setw(12) << pitm            << std::endl
      << std::left  << std::setw(18) << "p-Beam/TM"     << "\t"
      << std::right << std::setw(12) << ptm             << std::endl
      << std::left  << std::setw(18) << "(BH2,K)/Beam"  << "\t"
      << std::right << std::setw(12) << bh2kbeam        << std::endl
      << std::left  << std::setw(18) << "(p,BGO)/Beam"  << "\t"
      << std::right << std::setw(12) << pbgobeam        << std::endl
      << std::setprecision(0)
      << std::left  << std::setw(18) << "pi-Beam/Spill" << "\t"
      << std::right << std::setw(12) << pirate          << std::endl
      << std::left  << std::setw(18) << "p-Beam/Spill"  << "\t"
      << std::right << std::setw(12) << prate           << std::endl
      << std::left  << std::setw(18) << "L1-Req/Spill"  << "\t"
      << std::right << std::setw(12) << l1rate          << std::endl
      << std::left  << std::setw(18) << "L2-Acc/Spill"  << "\t"
      << std::right << std::setw(12) << l2rate          << std::endl;

#endif

  return true;
}
