// -*- C++ -*-

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
  std::vector<ScalerAnalyzer> gScaler(2);
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
  TString ClassName( void ) { return "EventScaler"; }
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

  gScaler[0].Decode();
  gScaler[1].Decode();

  // if( gUnpacker.get_event_number()%500==0 )
  //   gScaler[0].Print();

#if SPILL_RESET
  if( gScaler[0].SpillIncrement() )
    gScaler[0].Clear();
  if( gScaler[1].SpillIncrement() )
    gScaler[1].Clear();
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
    gScaler[0].Set( c, r++, ScalerInfo( "K-in",       0, 39 ) );
    gScaler[0].Set( c, r++, ScalerInfo( "Pi-in",      0, 40 ) );
    gScaler[0].Set( c, r++, ScalerInfo( "BH1",        0, 16 ) );
    gScaler[0].Set( c, r++, ScalerInfo( "BH1-SUM",   -1, -1 ) );
    gScaler[0].Set( c, r++, ScalerInfo( "BH1-01",     1,  0 ) );
    gScaler[0].Set( c, r++, ScalerInfo( "BH1-02",     1,  1 ) );
    gScaler[0].Set( c, r++, ScalerInfo( "BH1-03",     1,  2 ) );
    gScaler[0].Set( c, r++, ScalerInfo( "BH1-04",     1,  3 ) );
    gScaler[0].Set( c, r++, ScalerInfo( "BH1-05",     1,  4 ) );
    gScaler[0].Set( c, r++, ScalerInfo( "BH1-06",     1,  5 ) );
    gScaler[0].Set( c, r++, ScalerInfo( "BH1-07",     1,  6 ) );
    gScaler[0].Set( c, r++, ScalerInfo( "BH1-08",     1,  7 ) );
    gScaler[0].Set( c, r++, ScalerInfo( "BH1-09",     1,  8 ) );
    gScaler[0].Set( c, r++, ScalerInfo( "BH1-10",     1,  9 ) );
    gScaler[0].Set( c, r++, ScalerInfo( "BH1-11",     1, 10 ) );
    gScaler[0].Set( c, r++, ScalerInfo( "BH2",        0, 17 ) );
    gScaler[0].Set( c, r++, ScalerInfo( "BH2-SUM",   -1, -1 ) );
    gScaler[0].Set( c, r++, ScalerInfo( "BH2-01",     0, 64 ) );
    gScaler[0].Set( c, r++, ScalerInfo( "BH2-02",     0, 65 ) );
    gScaler[0].Set( c, r++, ScalerInfo( "BH2-03",     0, 66 ) );
    gScaler[0].Set( c, r++, ScalerInfo( "BH2-04",     0, 67 ) );
    gScaler[0].Set( c, r++, ScalerInfo( "BH2-05",     0, 68 ) );
    gScaler[0].Set( c, r++, ScalerInfo( "BAC",        0, 18 ) );
    gScaler[0].Set( c, r++, ScalerInfo( "PVAC",       0, 20 ) );
    gScaler[0].Set( c, r++, ScalerInfo( "FAC",        0, 21 ) );
    gScaler[0].Set( c, r++, ScalerInfo( "SCH",        0, 11 ) );
    gScaler[0].Set( c, r++, ScalerInfo( "TOF",        0, 22 ) );
    gScaler[0].Set( c, r++, ScalerInfo( "LAC",        0, 23 ) );
    gScaler[0].Set( c, r++, ScalerInfo( "WC",         0, 24 ) );
    gScaler[0].Set( c, r++, ScalerInfo( "BH1-1/100-PS", 1, 11 ) );
    gScaler[0].Set( c, r++, ScalerInfo( "BH1-1/1e5-PS", 1, 12 ) );
    gScaler[0].Set( c, r++, ScalerInfo( "TOF-24",       0, 29 ) );
    gScaler[0].Set( c, r++, ScalerInfo( "Mtx2D-1",      0, 32 ) );
    gScaler[0].Set( c, r++, ScalerInfo( "Mtx2D-2",      0, 33 ) );
    gScaler[0].Set( c, r++, ScalerInfo( "Mtx3D",        0, 34 ) );
    gScaler[0].Set( c, r++, ScalerInfo( "Other1",       0, 25 ) );
    gScaler[0].Set( c, r++, ScalerInfo( "Other2",       0, 26 ) );
    gScaler[0].Set( c, r++, ScalerInfo( "Other3",       0, 27 ) );
    gScaler[0].Set( c, r++, ScalerInfo( "Other4",       0, 28 ) );
  }

  {
    Int_t c = ScalerAnalyzer::kCenter;
    Int_t r = 0;
    for( Int_t i=0; i<16; ++i ){
      gScaler[0].Set( c, r++, ScalerInfo( Form( "TFA-%02d", i+1 ), 2, i ) );
    }
    for( Int_t i=0; i<16; ++i ){
      gScaler[0].Set( c, r++, ScalerInfo( Form( "CRM-%02d", i+1 ), 2, i+16 ) );
    }
    for( Int_t i=0; i<16; ++i ){
      gScaler[0].Set( c, r++, ScalerInfo( Form( "Reset-%02d", i+1 ), 2, i+32 ) );
    }
    gScaler[0].Set( c, r++, ScalerInfo( "LSO1",         2, 48 ) );
    gScaler[0].Set( c, r++, ScalerInfo( "LSO2",         2, 49 ) );
    gScaler[0].Set( c, r++, ScalerInfo( "LSO1xGe13",    2, 50 ) );
    gScaler[0].Set( c, r++, ScalerInfo( "LSO2xGe24",    2, 51 ) );
    gScaler[0].Set( c, r++, ScalerInfo( "GeCoin1",      2, 52 ) );
    gScaler[0].Set( c, r++, ScalerInfo( "GeCoin2",      2, 53 ) );
    gScaler[0].Set( c, r++, ScalerInfo( "LSO1High",     2, 54 ) );
    gScaler[0].Set( c, r++, ScalerInfo( "LSO2High",     2, 55 ) );
    gScaler[0].Set( c, r++, ScalerInfo( "GeHigh1",      2, 56 ) );
    gScaler[0].Set( c, r++, ScalerInfo( "GeHigh2",      2, 57 ) );
  }

  {
    Int_t c = ScalerAnalyzer::kRight;
    Int_t r = 0;
    gScaler[0].Set( c, r++, ScalerInfo( "Spill",        -1, -1 ) );
    gScaler[0].Set( c, r++, ScalerInfo( "10M-Clock",    0,  0 ) );
    gScaler[0].Set( c, r++, ScalerInfo( "TM",            0,  9 ) );
    gScaler[0].Set( c, r++, ScalerInfo( "SY",            0, 10 ) );
    gScaler[0].Set( c, r++, ScalerInfo( "Real-Time",     0,  1 ) );
    gScaler[0].Set( c, r++, ScalerInfo( "Live-Time",     0,  2 ) );
    gScaler[0].Set( c, r++, ScalerInfo( "L1-Req",        0,  3 ) );
    gScaler[0].Set( c, r++, ScalerInfo( "L1-Acc",        0,  4 ) );
    // gScaler[0].Set( c, r++, ScalerInfo( "MstClr",        0,  5 ) );
    gScaler[0].Set( c, r++, ScalerInfo( "Clear",         0,  6 ) );
    gScaler[0].Set( c, r++, ScalerInfo( "L2-Req",        0,  7 ) );
    gScaler[0].Set( c, r++, ScalerInfo( "L2-Acc",        0,  8 ) );
    gScaler[0].Set( c, r++, ScalerInfo( "BEAM-A",       0, 35 ) );
    gScaler[0].Set( c, r++, ScalerInfo( "BEAM-B",       0, 36 ) );
    gScaler[0].Set( c, r++, ScalerInfo( "BEAM-C",       0, 37 ) );
    gScaler[0].Set( c, r++, ScalerInfo( "BEAM-D",       0, 38 ) );
    gScaler[0].Set( c, r++, ScalerInfo( "BEAM-E",       0, 39 ) );
    gScaler[0].Set( c, r++, ScalerInfo( "BEAM-F",       0, 40 ) );
    gScaler[0].Set( c, r++, ScalerInfo( "TRIG-A",        0, 41 ) );
    gScaler[0].Set( c, r++, ScalerInfo( "TRIG-B",        0, 42 ) );
    gScaler[0].Set( c, r++, ScalerInfo( "TRIG-C",        0, 43 ) );
    gScaler[0].Set( c, r++, ScalerInfo( "TRIG-D",        0, 44 ) );
    gScaler[0].Set( c, r++, ScalerInfo( "TRIG-E",        0, 45 ) );
    gScaler[0].Set( c, r++, ScalerInfo( "TRIG-F",        0, 46 ) );
    gScaler[0].Set( c, r++, ScalerInfo( "TRIG-A-PS",     0, 48 ) );
    gScaler[0].Set( c, r++, ScalerInfo( "TRIG-B-PS",     0, 49 ) );
    gScaler[0].Set( c, r++, ScalerInfo( "TRIG-C-PS",     0, 50 ) );
    gScaler[0].Set( c, r++, ScalerInfo( "TRIG-D-PS",     0, 51 ) );
    gScaler[0].Set( c, r++, ScalerInfo( "TRIG-E-PS",     0, 52 ) );
    gScaler[0].Set( c, r++, ScalerInfo( "TRIG-F-PS",     0, 53 ) );
    gScaler[0].Set( c, r++, ScalerInfo( "TRIG-PSOR-A",   0, 54 ) );
    gScaler[0].Set( c, r++, ScalerInfo( "TRIG-PSOR-B",   0, 55 ) );
    gScaler[0].Set( c, r++, ScalerInfo( "Clock-PS",      0, 56 ) );
    gScaler[0].Set( c, r++, ScalerInfo( "Reserve2-PS",   0, 57 ) );
    gScaler[0].Set( c, r++, ScalerInfo( "Level1-PS",     0, 58 ) );
  }

  for( Int_t i=0; i<ScalerAnalyzer::MaxColumn; ++i ){
    for( Int_t j=0; j<ScalerAnalyzer::MaxRow; ++j ){
      gScaler[1].Set( i, j, gScaler[0].GetScalerInfo( i, j ) );
    }
  }


#if USE_COMMA
  gScaler[0].SetFlag( ScalerAnalyzer::kSeparateComma );
  gScaler[1].SetFlag( ScalerAnalyzer::kSeparateComma );
#endif

  gScaler[0].SetFlag( ScalerAnalyzer::kSpillOn );
  gScaler[1].SetFlag( ScalerAnalyzer::kSpillOff );

  gScaler[0].PrintFlags();
  gScaler[1].PrintFlags();

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
  gScaler[0].Print();

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

  Int_t recorder_event_number = 0;
  Bool_t found_run_number = false;
  std::string line;
  std::string recorder_line;
  while( ifs.good() && std::getline(ifs,line) ){
    if( line.empty() ) continue;
    std::istringstream input_line( line );
    std::istream_iterator<std::string> line_begin( input_line );
    std::istream_iterator<std::string> line_end;
    std::vector<std::string> log_column( line_begin, line_end );
    if( log_column.at(0) != "RUN" ) continue;
    if( log_column.at(1) != run_number_ss.str() ) continue;
    recorder_event_number = hddaq::a2i( log_column.at(15) );
    recorder_line = line;
    found_run_number = true;
  }

  if( !found_run_number ){
    std::cerr << "#E " << FUNC_NAME << " "
	      << "not found run# " << run_number
	      << " in " << recorder_log << std::endl;
    return false;
  }

  const TString& scaler_dir( bin_dir+"/../scaler" );

  for( const auto& scaler : gScaler ){
    TString scaler_txt;
    if( scaler.GetFlag( ScalerAnalyzer::kSpillOn ) ){
      scaler_txt = Form( "%s/spillon_%05d.txt",
			 scaler_dir.Data(), run_number );
    } else if( scaler.GetFlag( ScalerAnalyzer::kSpillOff ) ){
      scaler_txt = Form( "%s/spilloff_%05d.txt",
			 scaler_dir.Data(), run_number );
    } else {
      return false;
    }
    std::ofstream ofs( scaler_txt );
    if( !ofs.is_open() ){
      std::cerr << "#E " << FUNC_NAME << " "
		<< "cannot open scaler.txt : "
		<< scaler_txt << std::endl;
      return false;
    }

    ofs << recorder_line << std::endl;

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
	  TString name = scaler.GetScalerName( c, i );
	  if( name=="n/a" ) continue;
	  ofs << std::left  << std::setw(15) << name << "\t"
	      << std::right << std::setw(15) << scaler.Get(name) << std::endl;
	}
	ofs << std::endl;
      }
    }

    Double_t reallive = scaler.Fraction("Live-Time", "Real-Time");
    Double_t daqeff   = scaler.Fraction("L1-Acc", "L1-Req");
    Double_t l2eff    = scaler.Fraction("L2-Acc", "L1-Acc");
    Double_t ktm      = scaler.Fraction("K-in", "TM");
    Double_t pitm     = scaler.Fraction("Pi-in", "TM");
    Double_t l1reqbh2 = scaler.Fraction("L1-Req", "BH2");
    Double_t krate    = scaler.Fraction("K-in", "Spill");
    Double_t pirate   = scaler.Fraction("Pi-in", "Spill");
    Double_t l1rate   = scaler.Fraction("L1-Req", "Spill");
    Double_t l2rate   = scaler.Fraction("L2-Acc", "Spill");

    ofs << std::fixed << std::setprecision(6)
	<< std::left  << std::setw(18) << "Live/Real"     << "\t"
	<< std::right << std::setw(12)<<  reallive        << std::endl
	<< std::left  << std::setw(18) << "DAQ-Eff"       << "\t"
	<< std::right << std::setw(12) << daqeff          << std::endl
	<< std::left  << std::setw(18) << "L2-Eff"        << "\t"
	<< std::right << std::setw(12) << l2eff           << std::endl
	<< std::left  << std::setw(18) << "Duty-Factor"   << "\t"
	<< std::right << std::setw(12) << scaler.Duty()  << std::endl
	<< std::left  << std::setw(18) << "K-in/TM"       << "\t"
	<< std::right << std::setw(12) << ktm             << std::endl
	<< std::left  << std::setw(18) << "Pi-in/TM"      << "\t"
	<< std::right << std::setw(12) << pitm            << std::endl
	<< std::left  << std::setw(18) << "L1-Req/BH2"  << "\t"
	<< std::right << std::setw(12) << l1reqbh2        << std::endl
	<< std::setprecision(0)
	<< std::left  << std::setw(18) << "K-in/Spill" << "\t"
	<< std::right << std::setw(12) << krate          << std::endl
	<< std::left  << std::setw(18) << "Pi-in/Spill"  << "\t"
	<< std::right << std::setw(12) << pirate           << std::endl
	<< std::left  << std::setw(18) << "L1-Req/Spill"  << "\t"
	<< std::right << std::setw(12) << l1rate          << std::endl
	<< std::left  << std::setw(18) << "L2-Acc/Spill"  << "\t"
	<< std::right << std::setw(12) << l2rate          << std::endl;
  }
#endif

  return true;
}
