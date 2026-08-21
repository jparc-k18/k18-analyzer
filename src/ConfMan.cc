/*
  ConfMan.cc

  2013/4/14
*/

#include "ConfMan.hh"

#include <algorithm>
#include <iostream>
#include <fstream>
#include <iterator>
#include <sstream>
#include <vector>
#include <boost/lexical_cast.hpp>

#include "UnpackerManager.hh"
#include "HodoParamMan.hh"
#include "HodoPHCMan.hh"
#include "DCGeomMan.hh"
#include "DCTdcCalibMan.hh"
#include "DCDriftParamMan.hh"
#include "K18TransMatrix.hh"

#include "EvDisp.hh"

#include <sys/stat.h>

#include <TDirectory.h>
#include <TFile.h>
#include <TMacro.h>
#include <TNamed.h>
#include <TSystem.h>

ConfMan * ConfMan::confManager_=0;

ConfMan::ConfMan()
  : m_ConfFileName(),
    m_UnpackerFileName(),
    m_DigitFileName(),
    m_CmapFileName(),
    DCGeomManager_(0),
    DCTdcCalibManager_(0),
    DCDriftParamManager_(0),
    K18Matrix_(),
    K18Momentum_(1.05),
    SKSFieldNMR_(2.2), SKSFieldCalc_(2.2),
    evDisp_(0), FlagEvDisp_(false),
    FlagMHTDC_(false)
{
}

ConfMan::ConfMan(const std::string& confFile)
  : m_ConfFileName(confFile),
    m_UnpackerFileName(),
    HodoParamManager_(0),
    HodoPHCManager_(0),
    DCGeomManager_(0),
    DCTdcCalibManager_(0),
    DCDriftParamManager_(0),
    K18Matrix_(0),
    K18Momentum_(1.05),
    SKSFieldNMR_(2.2), SKSFieldCalc_(2.2),
    ScalerAnalyzer_(0),
    evDisp_(0), FlagEvDisp_(false),
    FlagMHTDC_(false)
{
  if (confManager_)
    std::exit(-1);
    
  confManager_ = this;
}

ConfMan::~ConfMan()
{
}

// ConfMan&
// ConfMan::GetInstance()
// {
//   static ConfMan gConfMan;
//   return gConfMan;
// }

bool
ConfMan::Initialize()
{
  bool status = true;

  double val;
  std::ifstream f(m_ConfFileName.c_str());
  //m_ConfFileName = confFile;
  if (f.fail()){
    std::cerr << "#E " 
	      << std::endl;
    return false;
  }
  
  while (f.good()){
    std::string line;
    std::getline(f, line);
    if (line.empty())
      continue;
    std::istringstream input_line(line);
    std::istream_iterator<std::string> line_begin(input_line);
    std::istream_iterator<std::string> line_end;
    std::vector<std::string> param_list(line_begin, line_end);
    if (param_list.size()!=2)
      continue;
    std::cout << param_list[0] << " " << param_list[1] << std::endl;
    if (param_list[0][0] != '#')
      m_paramLines[param_list[0]] = param_list[1];
    //Conf
    if (param_list[0] == "UNPACKER:")
      m_UnpackerFileName = param_list[1];
    if (param_list[0] == "DIGIT:")
      m_DigitFileName = param_list[1];
    if (param_list[0] == "CMAP:")
      m_CmapFileName = param_list[1];
    //Hodo
    if (param_list[0] == "HDPRM:")
      HodoParamFileName_ = param_list[1];	
    if (param_list[0] == "HDPHC:")
      HodoPHCFileName_ = param_list[1];	
    //DC
    if (param_list[0] == "DCGEOM:")
      DCGeomFileName_ = param_list[1];
    if (param_list[0] == "TDCCALIB:")
      DCTdcCalibFileName_ = param_list[1];
    if (param_list[0] == "DRFTPM:")
      DCDriftParamFileName_ = param_list[1];
    if (param_list[0] == "K18TM:")
      K18MatrixFileName_ = param_list[1];
    if (param_list[0] == "PK18:")
      K18Momentum_ = boost::lexical_cast<double>(param_list[1]);
    if (param_list[0] == "FLDMAP:")
      FieldMapFileName_ = param_list[1];
    if (param_list[0] == "FLDNMR:")
      SKSFieldNMR_ = boost::lexical_cast<double>(param_list[1]);
    if (param_list[0] == "FLDCALC:")
      SKSFieldCalc_ = boost::lexical_cast<double>(param_list[1]);

    // h_Utility
    if(param_list[0] == "FILEMAN:"){
      FileManagerName_ = param_list[1];
    }

    // Event display 
    if( param_list[0] == "EVDISP:" )
      if(param_list[1]=="1") FlagEvDisp_=true;
      else         FlagEvDisp_=false;
      
    //MHTDC event suppression flag
    if( param_list[0] == "MHTDC:" )
      if(param_list[1]=="1") FlagMHTDC_=true;
      else         FlagMHTDC_=false;
  }
  
  status = InitializeParameterFiles();
  InitializeHistograms();
  hddaq::unpacker::UnpackerManager& gUnpacker
    = hddaq::unpacker::GUnpacker::get_instance();

  gUnpacker.set_config_file( m_UnpackerFileName, m_DigitFileName, m_CmapFileName);

  return status;
}

bool ConfMan::InitializeEvDisp( void )
{
  static const std::string funcname = "[ConfMan::InitializeEvDisp]";
  evDisp_ = & EvDisp::GetInstance();
  evDisp_->Initialize();
  return true;
}

namespace
{
  // conf-referenced files above this size (e.g. the multi-ten-MB
  // SKS field map, FLDMAP:) are not small calibration tables and
  // must not be embedded: each unit-split sub-job would write its
  // own full copy, and hadd cannot merge/dedupe TMacro objects, so
  // N units means N copies (seen in practice: 70 units x 72 MB
  // FLDMAP = ~5 GB of pure waste in a single merged output file).
  const long MAX_EMBED_BYTES = 1 * 1024 * 1024; // 1 MB

  bool isEmbeddableFile( const std::string& path )
  {
    struct stat sb;
    return ( !path.empty() && stat( path.c_str(), &sb ) == 0
	     && S_ISREG( sb.st_mode ) && sb.st_size <= MAX_EMBED_BYTES );
  }
}

void
ConfMan::WriteParameters()
{
  if( !gFile ) return;
  // a unit-split run spawns one process per event-range chunk, each
  // with its own ConfMan/WriteParameters() call; only the first
  // chunk (skip==0) should write metadata, otherwise hadd ends up
  // with one copy per chunk (non-mergeable TNamed/TMacro objects
  // just accumulate as extra cycles -- e72's ConfMan has the same
  // guard for the same reason).
  hddaq::unpacker::UnpackerManager& gUnpacker
    = hddaq::unpacker::GUnpacker::get_instance();
  if( gUnpacker.get_skip() != 0 ) return;
  TDirectory* prevDir = gDirectory;
  gFile->cd();
  TNamed git( "git", ( "\n" + gSystem->GetFromPipe( "git log -1" ) ).Data() );
  git.Write();
  TNamed conffile( "conffile", m_ConfFileName.c_str() );
  conffile.Write();
  gFile->mkdir( "param" );
  gFile->cd( "param" );
  for( std::map<std::string, std::string>::const_iterator
	 it = m_paramLines.begin(); it != m_paramLines.end(); ++it ){
    std::string key = it->first;
    if( !key.empty() && key[key.size()-1] == ':' )
      key.erase( key.size()-1 );
    if( isEmbeddableFile( it->second ) ){
      TMacro macro;
      macro.SetName( key.c_str() );
      macro.SetTitle( it->second.c_str() );
      macro.ReadFile( it->second.c_str() );
      macro.Write();
    } else {
      TNamed named( key.c_str(), it->second.c_str() );
      named.Write();
    }
  }
  gFile->cd();
  if( prevDir ) prevDir->cd();
}
