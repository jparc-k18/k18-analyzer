/*
  DstMetaUtil.hh

  Header-only helpers shared by the Dst-stage merger programs
  (DstKKAna.cc/DstKPiAna.cc, DstSksHodoscope.cc) to embed
  reproducibility metadata into their output ROOT files: this
  program's own conf-file contents (mirrors
  ConfMan::WriteParameters() used by the Stage-1 usr/ analyzers),
  plus a COPY of the same metadata already embedded in each
  upstream Stage-1 input file, so the final Dst output is fully
  self-documenting without having to chase the upstream files.
*/

#ifndef DstMetaUtil_hh
#define DstMetaUtil_hh 1

#include <sys/stat.h>

#include <fstream>
#include <sstream>
#include <string>

#include <TDirectory.h>
#include <TFile.h>
#include <TIterator.h>
#include <TKey.h>
#include <TMacro.h>
#include <TNamed.h>
#include <TObject.h>
#include <TSystem.h>

namespace DstMeta
{
  // conf-referenced files above this size (e.g. the multi-ten-MB
  // SKS field map, FLDMAP:) are not small calibration tables and
  // must not be embedded -- seen in practice: unit-split Stage-1
  // jobs each embedding a full 72 MB FLDMAP copy, and hadd unable
  // to merge/dedupe TMacro objects, ballooned a single merged
  // output file to ~5 GB of pure waste.
  const long MAX_EMBED_BYTES = 1 * 1024 * 1024; // 1 MB

  inline bool
  IsReadableFile( const std::string& path )
  {
    struct stat sb;
    return ( !path.empty() && stat( path.c_str(), &sb ) == 0
	     && S_ISREG( sb.st_mode ) && sb.st_size <= MAX_EMBED_BYTES );
  }

  // Parse this program's own conf file (every 2-token
  // "KEY VALUE" line) and write it into a "param" TDirectory +
  // top-level "git"/"conffile" TNamed in the CURRENT directory
  // (caller must fOut->cd() first).
  inline void
  WriteOwnParam( const std::string& conf )
  {
    TNamed git( "git", ( "\n" + gSystem->GetFromPipe( "git log -1" ) ).Data() );
    git.Write();
    TNamed conffile( "conffile", conf.c_str() );
    conffile.Write();

    TDirectory* here = gDirectory;
    here->mkdir( "param" );
    here->cd( "param" );
    std::ifstream ifs( conf.c_str() );
    std::string line;
    while( std::getline( ifs, line ) ){
      std::istringstream iss( line );
      std::string k, v;
      if( !(iss >> k >> v) ) continue;
      if( k[0] == '#' ) continue; // skip commented-out lines
      std::string key = k;
      if( !key.empty() && key[key.size()-1] == ':' ) key.erase( key.size()-1 );
      if( IsReadableFile( v ) ){
	TMacro macro;
	macro.SetName( key.c_str() );
	macro.SetTitle( v.c_str() );
	macro.ReadFile( v.c_str() );
	macro.Write();
      } else {
	TNamed named( key.c_str(), v.c_str() );
	named.Write();
      }
    }
    here->cd();
  }

  // Copy the "param"/"git"/"conffile" objects already embedded in
  // an upstream Stage-1 (usr/, via ConfMan::WriteParameters()) or
  // Dst-stage (via WriteOwnParam() above) input file into the
  // CURRENT directory, tagged with a suffix so multiple upstream
  // inputs don't collide (e.g. "k18tracking", "skstracking").
  inline void
  CopyUpstreamParam( TFile* src, const std::string& tag )
  {
    if( !src ) return;
    TNamed* git = dynamic_cast<TNamed*>( src->Get( "git" ) );
    if( git ){
      TNamed copy( ( "git_" + tag ).c_str(), git->GetTitle() );
      copy.Write();
    }
    TNamed* conffile = dynamic_cast<TNamed*>( src->Get( "conffile" ) );
    if( conffile ){
      TNamed copy( ( "conffile_" + tag ).c_str(), conffile->GetTitle() );
      copy.Write();
    }
    TDirectory* srcParam = dynamic_cast<TDirectory*>( src->Get( "param" ) );
    if( !srcParam ) return;

    TDirectory* here = gDirectory;
    const std::string dirName = "param_" + tag;
    here->mkdir( dirName.c_str() );
    here->cd( dirName.c_str() );
    TIter next( srcParam->GetListOfKeys() );
    TKey* key;
    while( ( key = static_cast<TKey*>( next() ) ) ){
      TObject* obj = key->ReadObj();
      if( obj ){ obj->Write(); delete obj; }
    }
    here->cd();
  }
}

#endif
