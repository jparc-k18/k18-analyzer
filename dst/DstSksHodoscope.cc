/*
  DstSksHodoscope.cc

  Dst-stage combination of SksTracking and SksHodoscope: the TOF
  times and dE are RECOMPUTED from the raw TDC/ADC stored in the
  SksHodoscope tree with the CURRENT HodoParam/HodoPHC files, so
  the TOF calibration can be iterated without re-running the
  Runge-Kutta tracking (cf. e07 dst/DstKuramaHodoscope).

  The two trees are matched by the unpacker event number (evnum);
  events skipped by either analyzer are dropped. The scattered-
  particle flight time is
    stof = cmtTof - time0 - L_bt/(beta_beam c)
  with L_bt from the surveyed DCGEOM positions (K18Target id 130,
  BH2 id 202) and beta_beam from PK18 (kaon).

  Usage:
    bin/DstSksHodoscope [conf] [SksTracking.root]
                        [SksHodoscope.root] [out.root]
*/

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

#include <TFile.h>
#include <TROOT.h>
#include <TTree.h>

#include "DCGeomMan.hh"
#include "HodoParamMan.hh"
#include "HodoPHCMan.hh"

#include "DstMetaUtil.hh"

namespace
{
  const int    NumOfSegTOF = 32;
  const int    NumOfMisc   = 20; // trigger flag array size
  const int    MaxHits     = 20;
  const int    CidTOF      = 5;
  const int    IdK18Target = 130;
  const int    IdBH2       = 202;
  const double KaonMass    = 0.493677;
  const double LightVel    = 299.792458; // mm/ns
}

//______________________________________________________________________
struct Src
{
  // SksTracking
  int    evnum;
  int    ntSks;
  double chisqrSks[MaxHits];
  double p[MaxHits];
  double path[MaxHits];
  int    ntOut;
  double y0Out[MaxHits];
  double v0Out[MaxHits];
  // trigger flag, read from SksTracking's copy (see DstKKAna.cc
  // for why: same raw-DAQ quantity duplicated in each Stage-1 tree)
  int    trigpat[NumOfMisc], trigflag[NumOfMisc];
  // SksHodoscope
  int    hevnum;
  double time0;
  int    tofua[NumOfSegTOF], tofda[NumOfSegTOF];
  int    tofut[NumOfSegTOF], tofdt[NumOfSegTOF];
};

struct Dst
{
  int    evnum;
  int    ntSks;
  double chisqrSks[MaxHits];
  double p[MaxHits];
  double path[MaxHits];
  int    ntOut;
  double yTof[MaxHits];
  double time0;
  int    nhTof;
  int    TofSeg[NumOfSegTOF];
  double tuTof[NumOfSegTOF], tdTof[NumOfSegTOF];
  double ctuTof[NumOfSegTOF], ctdTof[NumOfSegTOF];
  double cmtTof[NumOfSegTOF];
  double deuTof[NumOfSegTOF], dedTof[NumOfSegTOF];
  double deTof[NumOfSegTOF];
  double stof[NumOfSegTOF];
  double beta[NumOfSegTOF];
  double m2[NumOfSegTOF];
  int    trigpat[NumOfMisc], trigflag[NumOfMisc];
};

//______________________________________________________________________
static std::string
confValue( const std::string& conf, const std::string& key )
{
  std::ifstream ifs( conf.c_str() );
  std::string line;
  while( std::getline( ifs, line ) ){
    std::istringstream iss( line );
    std::string k, v;
    if( iss >> k >> v && k == key ) return v;
  }
  return "";
}

//______________________________________________________________________
int
main( int argc, char* argv[] )
{
  if( argc != 5 ){
    std::cout << "Usage: " << argv[0]
	      << " [conf] [SksTracking.root] [SksHodoscope.root]"
	      << " [out.root]" << std::endl;
    return EXIT_FAILURE;
  }
  const std::string conf = argv[1];

  const std::string prmFile = confValue( conf, "HDPRM:" );
  const std::string phcFile = confValue( conf, "HDPHC:" );
  const std::string geoFile = confValue( conf, "DCGEOM:" );
  const std::string pk18Str = confValue( conf, "PK18:" );
  if( prmFile.empty() || geoFile.empty() || pk18Str.empty() ){
    std::cerr << "#E missing HDPRM/DCGEOM/PK18 in " << conf
	      << std::endl;
    return EXIT_FAILURE;
  }
  const double pK18 = std::atof( pk18Str.c_str() );

  HodoParamMan hodoMan( prmFile );
  hodoMan.Initialize();
  HodoPHCMan phcMan( phcFile );
  phcMan.Initialize();
  DCGeomMan& geomMan = DCGeomMan::GetInstance();
  geomMan.Initialize( geoFile );
  const double lbt =
    geomMan.GetLocalZ( IdK18Target ) - geomMan.GetLocalZ( IdBH2 );
  const double bBeam = pK18/std::sqrt( pK18*pK18
				       + KaonMass*KaonMass );
  const double tBeam = lbt/(bBeam*LightVel);
  std::cout << "#D L(BH2->target) = " << lbt << " mm, "
	    << "beam leg = " << tBeam << " ns (pK18 = "
	    << pK18 << ")" << std::endl;

  TFile* fSks  = TFile::Open( argv[2] );
  TFile* fHodo = TFile::Open( argv[3] );
  if( !fSks || !fHodo || fSks->IsZombie() || fHodo->IsZombie() ){
    std::cerr << "#E cannot open the input files" << std::endl;
    return EXIT_FAILURE;
  }
  TTree* tSks  = dynamic_cast<TTree*>( fSks->Get( "tree" ) );
  TTree* tHodo = dynamic_cast<TTree*>( fHodo->Get( "tree" ) );

  Src src;
  tSks->SetBranchAddress( "evnum",     &src.evnum );
  tSks->SetBranchAddress( "ntSks",     &src.ntSks );
  tSks->SetBranchAddress( "chisqrSks",  src.chisqrSks );
  tSks->SetBranchAddress( "p",          src.p );
  tSks->SetBranchAddress( "path",       src.path );
  tSks->SetBranchAddress( "ntOut",     &src.ntOut );
  tSks->SetBranchAddress( "y0Out",      src.y0Out );
  tSks->SetBranchAddress( "v0Out",      src.v0Out );
  tSks->SetBranchAddress( "trigpat",    src.trigpat );
  tSks->SetBranchAddress( "trigflag",   src.trigflag );
  tHodo->SetBranchAddress( "evnum", &src.hevnum );
  tHodo->SetBranchAddress( "time0", &src.time0 );
  tHodo->SetBranchAddress( "tofua",  src.tofua );
  tHodo->SetBranchAddress( "tofda",  src.tofda );
  tHodo->SetBranchAddress( "tofut",  src.tofut );
  tHodo->SetBranchAddress( "tofdt",  src.tofdt );
  tHodo->BuildIndex( "evnum" );

  TFile* fOut = TFile::Open( argv[4], "recreate" );
  fOut->cd();
  DstMeta::WriteOwnParam( conf );
  DstMeta::CopyUpstreamParam( fSks,  "skstracking" );
  DstMeta::CopyUpstreamParam( fHodo, "skshodoscope" );

  TTree* tOut = new TTree( "tree", "DstSksHodoscope" );
  Dst dst;
  tOut->Branch( "evnum",     &dst.evnum,     "evnum/I" );
  tOut->Branch( "ntSks",     &dst.ntSks,     "ntSks/I" );
  tOut->Branch( "chisqrSks",  dst.chisqrSks, "chisqrSks[ntSks]/D" );
  tOut->Branch( "p",          dst.p,         "p[ntSks]/D" );
  tOut->Branch( "path",       dst.path,      "path[ntSks]/D" );
  tOut->Branch( "ntOut",     &dst.ntOut,     "ntOut/I" );
  tOut->Branch( "yTof",       dst.yTof,      "yTof[ntOut]/D" );
  tOut->Branch( "time0",     &dst.time0,     "time0/D" );
  tOut->Branch( "nhTof",     &dst.nhTof,     "nhTof/I" );
  tOut->Branch( "TofSeg",     dst.TofSeg,    "TofSeg[nhTof]/I" );
  tOut->Branch( "tuTof",      dst.tuTof,     "tuTof[nhTof]/D" );
  tOut->Branch( "tdTof",      dst.tdTof,     "tdTof[nhTof]/D" );
  tOut->Branch( "ctuTof",     dst.ctuTof,    "ctuTof[nhTof]/D" );
  tOut->Branch( "ctdTof",     dst.ctdTof,    "ctdTof[nhTof]/D" );
  tOut->Branch( "cmtTof",     dst.cmtTof,    "cmtTof[nhTof]/D" );
  tOut->Branch( "deuTof",     dst.deuTof,    "deuTof[nhTof]/D" );
  tOut->Branch( "dedTof",     dst.dedTof,    "dedTof[nhTof]/D" );
  tOut->Branch( "deTof",      dst.deTof,     "deTof[nhTof]/D" );
  tOut->Branch( "stof",       dst.stof,      "stof[nhTof]/D" );
  tOut->Branch( "beta",       dst.beta,      "beta[nhTof]/D" );
  tOut->Branch( "m2",         dst.m2,        "m2[nhTof]/D" );
  tOut->Branch( "trigpat",    dst.trigpat,   "trigpat[20]/I" );
  tOut->Branch( "trigflag",   dst.trigflag,  "trigflag[20]/I" );

  const Long64_t nev = tSks->GetEntries();
  Long64_t nout = 0, nmiss = 0;
  for( Long64_t iev = 0; iev < nev; ++iev ){
    tSks->GetEntry( iev );
    if( src.ntSks < 1 ) continue;
    if( tHodo->GetEntryWithIndex( src.evnum ) <= 0 ){
      ++nmiss;
      continue;
    }
    if( src.evnum != src.hevnum ) { ++nmiss; continue; }

    dst.evnum = src.evnum;
    for( int i=0; i<NumOfMisc; ++i ){
      dst.trigpat[i]  = src.trigpat[i];
      dst.trigflag[i] = src.trigflag[i];
    }
    dst.ntSks = src.ntSks > MaxHits ? MaxHits : src.ntSks;
    for( int it = 0; it < dst.ntSks; ++it ){
      dst.chisqrSks[it] = src.chisqrSks[it];
      dst.p[it]         = src.p[it];
      dst.path[it]      = src.path[it];
    }
    // the TOF plane is the origin of the SdcOut local frame
    // (DCGEOM id 51, l = 0), so y(TOF) = y0Out
    dst.ntOut = src.ntOut > MaxHits ? MaxHits : src.ntOut;
    for( int it = 0; it < dst.ntOut; ++it )
      dst.yTof[it] = src.y0Out[it];
    dst.time0 = src.time0;

    // recompute the TOF observables from the raw values with the
    // CURRENT parameter files
    dst.nhTof = 0;
    const double p0    = src.p[0];
    const double path0 = src.path[0];
    for( int seg = 0; seg < NumOfSegTOF; ++seg ){
      if( src.tofut[seg] <= 0 || src.tofdt[seg] <= 0 ) continue;
      double tu, td, deu = 0., ded = 0.;
      if( !hodoMan.GetTime( CidTOF, 0, seg, 0, src.tofut[seg], tu ) )
	continue;
      if( !hodoMan.GetTime( CidTOF, 0, seg, 1, src.tofdt[seg], td ) )
	continue;
      if( src.tofua[seg] >= 0 )
	hodoMan.GetDe( CidTOF, 0, seg, 0, src.tofua[seg], deu );
      if( src.tofda[seg] >= 0 )
	hodoMan.GetDe( CidTOF, 0, seg, 1, src.tofda[seg], ded );
      double ctu = tu, ctd = td;
      phcMan.doCorrection( CidTOF, 0, seg, 0, tu, deu, ctu );
      phcMan.doCorrection( CidTOF, 0, seg, 1, td, ded, ctd );
      const double cmt = 0.5*(ctu + ctd);
      const double de  = std::sqrt( std::max( deu, 0. )
				    *std::max( ded, 0. ) );
      const double st  = cmt - dst.time0 - tBeam;
      const int n = dst.nhTof;
      dst.TofSeg[n] = seg;
      dst.tuTof[n] = tu;   dst.tdTof[n] = td;
      dst.ctuTof[n] = ctu; dst.ctdTof[n] = ctd;
      dst.cmtTof[n] = cmt;
      dst.deuTof[n] = deu; dst.dedTof[n] = ded;
      dst.deTof[n] = de;
      dst.stof[n] = st;
      double beta = ( st > 0. ) ? path0/(st*LightVel) : -1.;
      dst.beta[n] = beta;
      dst.m2[n] = ( beta > 0. )
	? p0*p0*(1. - beta*beta)/(beta*beta) : -9999.;
      ++dst.nhTof;
    }
    tOut->Fill();
    ++nout;
  }
  std::cout << "#D processed " << nev << " tracking entries, "
	    << nout << " written, " << nmiss
	    << " without hodo match" << std::endl;

  fOut->Write();
  fOut->Close();
  return EXIT_SUCCESS;
}
