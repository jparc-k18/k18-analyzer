/*
  DstKKAna.cc / DstKPiAna.cc

  Dst-stage (K-, K+) / (K-, pi+) pair maker: K18Tracking x
  SksTracking x SksHodoscope matched by evnum. The TOF m2 is
  recomputed from the raw TDC/ADC with the current HodoParam/
  HodoPHC (as in DstSksHodoscope) and the reaction vertex is
  computed from the beam and scattered tracks at the target
  (closeDist, vertex x/y/z, lab angle). No physics cuts are
  applied to the output TREE: the selection (chisqr, m2,
  closeDist, vtz, theta) is left visible/variable at the
  plotting stage.

  In ADDITION, this program applies the momentum calibration
  (h(p)/k(u)/PO(p) from param/MOMCALIB/e05/mom_calib.json,
  hardcoded below -- keep in sync by hand whenever mom_calib.py
  is rerun) and fills the final Xi-/Sigma- missing-mass
  histograms (basic + vertex selection) directly in this
  compiled loop, writing them into the same output file. Do the
  looping here, not in a PyROOT `for e in tree:` loop over
  millions of entries: plotting scripts should just hadd/Draw()
  these histograms. These same histograms also require the
  trigger flag matching each hypothesis' own production reaction
  ((K,K) for Xi-, (K,pi) for Sigma-/Lambda -- see DetectorID.hh's
  TrigFlag enum), since this rejects most of the residual pi-
  beam contamination the m2-window cut alone lets through.

  The scattered-mass hypothesis only affects the plotting stage;
  this program stores the kinematic quantities for both. The two
  binary names exist for bookkeeping and produce identical trees.

  The target material (ch2/carbon, for the eloss calculation)
  is read from the conf file's "TARGET:" key (defaults to ch2
  if absent), the same way HDPRM:/DCGEO:/PK18: are read --
  NOT a command-line argument, since runmanager invokes Dst
  programs with a fixed argument count built from "dstin" and
  cannot pass an extra one.

  Usage:
    bin/DstKKAna [conf] [K18Tracking.root] [SksTracking.root]
                 [SksHodoscope.root] [out.root]
*/

#include <sys/stat.h>

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include <TDirectory.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TMacro.h>
#include <TNamed.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TTree.h>

#include "ThreeVector.hh"
#include "DCGeomMan.hh"
#include "DetectorID.hh"
#include "HodoParamMan.hh"
#include "HodoPHCMan.hh"
#include "Kinematics.hh"

#include "DstMetaUtil.hh"

namespace
{
  // NumOfSegTOF, NumOfMisc (trigflag array size), IdK18Target
  // all come from DetectorID.hh (same values: 32, 20, 130)
  const int    MaxHits     = 20;
  const int    CidTOF      = 5;
  const int    IdBH2       = 202;
  const double KaonMass    = 0.493677;
  const double PionMass    = 0.139570;
  const double ProtonMass  = 0.938272;
  const double LightVel    = 299.792458; // mm/ns

  // ---- energy loss (Bethe, no density correction) ----
  const double ME      = 0.510999e-3;
  const double K_BETHE = 0.307075e-3; // GeV cm^2/mol, per Z/A
  struct Mat { double za, exc; };
  const Mat MAT_POLYSTYRENE = { 0.53768, 68.7e-9 };
  const Mat MAT_CH2         = { 0.57034, 57.4e-9 };
  const Mat MAT_CARBON      = { 0.49955, 78.0e-9 };
  const Mat MAT_AIR         = { 0.49919, 85.7e-9 };
  const Mat MAT_AEROGEL     = { 0.50093, 139.2e-9 };
  const Mat MAT_MYLAR       = { 0.52037, 78.7e-9 };
  const double HALF_CH2    = 4.6;  // half of the effective
  const double HALF_CARBON = 4.4;  // target thickness [g/cm2]

  // ---- momentum calibration ----
  // NOT CALIBRATED in this repository: the coefficients below
  // are zeroed, i.e. h(p) = k(u) = po(pb) = 0, so the momenta
  // pass through unchanged. Fit them for your own data set and
  // fill them in here (see param/MOMCALIB/ for the file format).
  const double H_COEF[4] = { 0., 0., 0., 0. };
  const double K_COEF[2] = { 0., 0. };
  const double PO0 = 0.;
  const double PO1 = 0.;

  const double MAX_CHISQR    = 5.;
  const double MAX_CLOSEDIST = 8.;
  const double MAX_VTZ       = 120.;
  const double MIN_THETA     = 2.;
}

//______________________________________________________________________
double dedx( double mass, double p, const Mat& mat )
{
  double e = std::sqrt( p*p + mass*mass );
  double beta2 = (p/e)*(p/e);
  double gamma = e/mass;
  double bg2 = beta2*gamma*gamma;
  double wmax = 2.*ME*bg2
    /(1. + 2.*gamma*ME/mass + (ME/mass)*(ME/mass));
  double arg = 2.*ME*bg2*wmax/(mat.exc*mat.exc);
  return K_BETHE*mat.za/beta2*(0.5*std::log(arg) - beta2);
}

//______________________________________________________________________
double dpLoss( double mass, double p, double thickness,
	       const Mat& mat, int nstep = 10 )
{
  double dx = thickness/nstep;
  double e = std::sqrt( p*p + mass*mass );
  for( int i = 0; i < nstep; ++i ){
    double pp = std::sqrt( std::max( e*e - mass*mass, 1e-12 ) );
    e -= dedx( mass, pp, mat )*dx;
  }
  double pout = std::sqrt( std::max( e*e - mass*mass, 0. ) );
  return p - pout;
}

//______________________________________________________________________
// momentum loss from BFT to the target plane, beam-line budget
double beamLineDp( double mass, double p )
{
  double total = 0., pp = p;
  struct { Mat mat; double t; } budget[4] = {
    { MAT_POLYSTYRENE, 0.4 }, { MAT_AEROGEL, 0.6 },
    { MAT_MYLAR, 0.1 }, { MAT_AIR, 0.5 } };
  for( int i = 0; i < 4; ++i ){
    double d = dpLoss( mass, pp, budget[i].t, budget[i].mat );
    total += d;
    pp -= d;
  }
  return total;
}

//______________________________________________________________________
double hPoly( double p )
{
  double r = 0., pk = 1.;
  for( int i = 0; i < 4; ++i ){ r += H_COEF[i]*pk; pk *= p; }
  return r;
}

//______________________________________________________________________
double kU( double u ) { return K_COEF[0]*u + K_COEF[1]*u*u; }

//______________________________________________________________________
double po( double pb ) { return PO0 + PO1*(pb - 1.8); }

//______________________________________________________________________
struct Src
{
  // K18Tracking
  int    evnum;
  int    ntK18;
  double p_3rd[MaxHits];
  double xtgtb[MaxHits], ytgtb[MaxHits];
  double utgtb[MaxHits], vtgtb[MaxHits];
  // beam TOF (BH1 cluster mean time - BH2 time0), for the beam
  // Kaon-peak calibration-quality check (hBtof below); NOT part
  // of the Dst pair tree, since it's a beam-only quantity
  int    nhBH1;
  double btof[MaxHits];
  // trigger flag: read from K18Tracking's copy, since the trigger
  // pattern is the same raw-DAQ quantity duplicated in all three
  // Stage-1 trees, and K18Tracking's fill is bug-free (unlike the
  // historical off-by-one in SksTracking/SksHodoscope, fixed
  // 2026-07-18 but not yet reflected in already-processed files)
  int    trigpat[NumOfMisc], trigflag[NumOfMisc];
  // SksTracking
  int    sevnum;
  int    ntSks;
  double chisqrSks[MaxHits];
  double p[MaxHits];
  double path[MaxHits];
  double xtgts[MaxHits], ytgts[MaxHits];
  double utgts[MaxHits], vtgts[MaxHits];
  // SksHodoscope
  int    hevnum;
  double time0;
  int    tofua[NumOfSegTOF], tofda[NumOfSegTOF];
  int    tofut[NumOfSegTOF], tofdt[NumOfSegTOF];
};

struct Dst
{
  int    runnum, evnum;
  // beam
  double pb;
  double xb, yb, ub, vb;
  // scattered
  int    ntSks;
  double ps, chisqrSks, path;
  double xs, ys, us, vs;
  // pair
  double theta, cost;
  double closedist, vtx, vty, vtz;
  // TOF pid
  int    nhTof;
  int    TofSeg[NumOfSegTOF];
  double stof[NumOfSegTOF], m2[NumOfSegTOF];
  double deTof[NumOfSegTOF];
  // trigger flag (see Src comment)
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
  if( argc != 6 ){
    std::cout << "Usage: " << argv[0]
	      << " [conf] [K18Tracking.root] [SksTracking.root]"
	      << " [SksHodoscope.root] [out.root]" << std::endl;
    return EXIT_FAILURE;
  }
  const std::string conf = argv[1];
  std::string target = confValue( conf, "TARGET:" );
  if( target.empty() ) target = "ch2";
  Mat targetMat = MAT_CH2;
  double half = HALF_CH2;
  double mRecoilXi = ProtonMass, mRecoilSg = ProtonMass;
  if( target == "carbon" ){
    targetMat = MAT_CARBON;
    half = HALF_CARBON;
    // K- + n -> K+ + Xi- and K- + p -> pi+ + Sigma- each keep
    // the mass number at A=12 (strangeness-exchange on a bound
    // nucleon, not knockout), so the recoiling system's "free
    // hyperon, undisturbed core" mass is M(12C) with that one
    // nucleon's mass traded for the hyperon's -- handled at the
    // plotting stage; here the missing-mass calculation itself
    // just needs the recoil mass under the "target is 12C"
    // hypothesis, i.e. M(12C) itself (12.174862 GeV... see below)
    mRecoilXi = 11.174862;
    mRecoilSg = 11.174862;
  } else if( target != "ch2" ){
    std::cerr << "#E unknown target: " << target << std::endl;
    return EXIT_FAILURE;
  }

  const std::string prmFile = confValue( conf, "HDPRM:" );
  const std::string phcFile = confValue( conf, "HDPHC:" );
  std::string geoFile = confValue( conf, "DCGEO:" );
  // legacy E05 spelling, cf. ConfMan's alias table
  if( geoFile.empty() ) geoFile = confValue( conf, "DCGEOM:" );
  const std::string pk18Str = confValue( conf, "PK18:" );
  if( prmFile.empty() || geoFile.empty() || pk18Str.empty() ){
    std::cerr << "#E missing HDPRM/DCGEO/PK18 in " << conf
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

  TFile* fK18  = TFile::Open( argv[2] );
  TFile* fSks  = TFile::Open( argv[3] );
  TFile* fHodo = TFile::Open( argv[4] );
  if( !fK18 || !fSks || !fHodo
      || fK18->IsZombie() || fSks->IsZombie()
      || fHodo->IsZombie() ){
    std::cerr << "#E cannot open the input files" << std::endl;
    return EXIT_FAILURE;
  }
  TTree* tK18  = dynamic_cast<TTree*>( fK18->Get( "tree" ) );
  TTree* tSks  = dynamic_cast<TTree*>( fSks->Get( "tree" ) );
  TTree* tHodo = dynamic_cast<TTree*>( fHodo->Get( "tree" ) );

  Src src;
  tK18->SetBranchAddress( "evnum", &src.evnum );
  tK18->SetBranchAddress( "ntK18", &src.ntK18 );
  tK18->SetBranchAddress( "p_3rd",  src.p_3rd );
  tK18->SetBranchAddress( "xtgt",   src.xtgtb );
  tK18->SetBranchAddress( "ytgt",   src.ytgtb );
  tK18->SetBranchAddress( "utgt",   src.utgtb );
  tK18->SetBranchAddress( "vtgt",   src.vtgtb );
  tK18->SetBranchAddress( "trigpat",  src.trigpat );
  tK18->SetBranchAddress( "trigflag", src.trigflag );
  tK18->SetBranchAddress( "nhBH1", &src.nhBH1 );
  tK18->SetBranchAddress( "btof",   src.btof );
  tSks->SetBranchAddress( "evnum", &src.sevnum );
  tSks->SetBranchAddress( "ntSks", &src.ntSks );
  tSks->SetBranchAddress( "chisqrSks", src.chisqrSks );
  tSks->SetBranchAddress( "p",      src.p );
  tSks->SetBranchAddress( "path",   src.path );
  tSks->SetBranchAddress( "xtgt",   src.xtgts );
  tSks->SetBranchAddress( "ytgt",   src.ytgts );
  tSks->SetBranchAddress( "utgt",   src.utgts );
  tSks->SetBranchAddress( "vtgt",   src.vtgts );
  tHodo->SetBranchAddress( "evnum", &src.hevnum );
  tHodo->SetBranchAddress( "time0", &src.time0 );
  tHodo->SetBranchAddress( "tofua",  src.tofua );
  tHodo->SetBranchAddress( "tofda",  src.tofda );
  tHodo->SetBranchAddress( "tofut",  src.tofut );
  tHodo->SetBranchAddress( "tofdt",  src.tofdt );
  tSks->BuildIndex( "evnum" );
  tHodo->BuildIndex( "evnum" );

  int runnum = 0;
  {
    std::string s( argv[5] );
    std::size_t i = s.find( "run" );
    if( i != std::string::npos )
      runnum = std::atoi( s.substr( i + 3, 5 ).c_str() );
  }

  TFile* fOut = TFile::Open( argv[5], "recreate" );
  fOut->cd();
  DstMeta::WriteOwnParam( conf );
  DstMeta::CopyUpstreamParam( fK18,  "k18tracking" );
  DstMeta::CopyUpstreamParam( fSks,  "skstracking" );
  DstMeta::CopyUpstreamParam( fHodo, "skshodoscope" );

  TTree* tOut = new TTree( "tree", "DstKKAna" );
  Dst dst;
  tOut->Branch( "runnum",    &dst.runnum,    "runnum/I" );
  tOut->Branch( "evnum",     &dst.evnum,     "evnum/I" );
  tOut->Branch( "pb",        &dst.pb,        "pb/D" );
  tOut->Branch( "xb",        &dst.xb,        "xb/D" );
  tOut->Branch( "yb",        &dst.yb,        "yb/D" );
  tOut->Branch( "ub",        &dst.ub,        "ub/D" );
  tOut->Branch( "vb",        &dst.vb,        "vb/D" );
  tOut->Branch( "ntSks",     &dst.ntSks,     "ntSks/I" );
  tOut->Branch( "ps",        &dst.ps,        "ps/D" );
  tOut->Branch( "chisqrSks", &dst.chisqrSks, "chisqrSks/D" );
  tOut->Branch( "path",      &dst.path,      "path/D" );
  tOut->Branch( "xs",        &dst.xs,        "xs/D" );
  tOut->Branch( "ys",        &dst.ys,        "ys/D" );
  tOut->Branch( "us",        &dst.us,        "us/D" );
  tOut->Branch( "vs",        &dst.vs,        "vs/D" );
  tOut->Branch( "theta",     &dst.theta,     "theta/D" );
  tOut->Branch( "cost",      &dst.cost,      "cost/D" );
  tOut->Branch( "closedist", &dst.closedist, "closedist/D" );
  tOut->Branch( "vtx",       &dst.vtx,       "vtx/D" );
  tOut->Branch( "vty",       &dst.vty,       "vty/D" );
  tOut->Branch( "vtz",       &dst.vtz,       "vtz/D" );
  tOut->Branch( "nhTof",     &dst.nhTof,     "nhTof/I" );
  tOut->Branch( "TofSeg",     dst.TofSeg,    "TofSeg[nhTof]/I" );
  tOut->Branch( "stof",       dst.stof,      "stof[nhTof]/D" );
  tOut->Branch( "m2",         dst.m2,        "m2[nhTof]/D" );
  tOut->Branch( "deTof",      dst.deTof,     "deTof[nhTof]/D" );
  tOut->Branch( "trigpat",   dst.trigpat,   "trigpat[20]/I" );
  tOut->Branch( "trigflag",  dst.trigflag,  "trigflag[20]/I" );

  // final calibrated missing-mass histograms (basic + vertex
  // selection), so plotting scripts only have to Draw()
  const bool isCarbon = ( target == "carbon" );
  // Carbon range covers -B_Y in [-80, +540] MeV around the
  // Xi-+11B / Lambda+11Be thresholds (see pymacro/e05/
  // mm_carbon.py -- the Lambda channel's daughter core is 11Be,
  // not 11B: the Sigma- N -> Lambda N in-medium conversion also
  // turns a second proton into a neutron), wide enough for both
  // the bound-state zoom (-40/+10 MeV) and the quasi-free view
  // (-40/+500 MeV). The bounds below are the thresholds
  // (M(11B)+M_Xi- = 11.574258, M(11Be)+M_Lambda = 11.380248 GeV,
  // same constants as mm_carbon.py) offset by +-0.080/+0.540 GeV
  // -- exact multiples of the 2 MeV bin width -- so that after
  // mm_carbon.py subtracts the threshold, bin edges land on round
  // MeV values (0 falls on an edge, not mid-bin). Do not replace
  // with round M_X numbers like 11.50/12.10: that reintroduces a
  // fractional offset since the thresholds themselves are not
  // multiples of 0.002.
  const double xiLo = isCarbon ? 11.494258 : 1.24;
  const double xiHi = isCarbon ? 12.114258 : 1.42;
  const double sgLo = isCarbon ? 11.300248 : 1.10;
  const double sgHi = isCarbon ? 11.920248 : 1.31;
  const int xiNb = int( (xiHi - xiLo)/0.002 + 0.5 );
  const int sgNb = int( (sgHi - sgLo)/0.002 + 0.5 );
  TH1D* hXiBasic  = new TH1D( "xi_basic",  "", xiNb, xiLo, xiHi );
  TH1D* hXiVertex = new TH1D( "xi_vertex", "", xiNb, xiLo, xiHi );
  TH1D* hSgBasic  = new TH1D( "sg_basic",  "", sgNb, sgLo, sgHi );
  TH1D* hSgVertex = new TH1D( "sg_vertex", "", sgNb, sgLo, sgHi );
  // vertex-selected mass vs lab scattering angle, for 2-degree
  // angle-sliced spectra at the plotting stage (theta covers the
  // full SKS acceptance, 0-32 deg; the MIN_THETA>2 deg cut still
  // applies via passVertex before filling)
  const int ThNb = 16;
  const double ThLo = 0., ThHi = 32.;
  TH2D* hXiVertexTh = new TH2D( "xi_vertex_theta", "",
				 xiNb, xiLo, xiHi, ThNb, ThLo, ThHi );
  TH2D* hSgVertexTh = new TH2D( "sg_vertex_theta", "",
				 sgNb, sgLo, sgHi, ThNb, ThLo, ThHi );
  // beam Kaon-peak calibration check: BH1-BH2 TOF (ns), 20 ps
  // bins over +-5 ns (wide enough to see any pi- shoulder)
  TH1D* hBtof = new TH1D( "btof", "", 500, -5., 5. );
  const double m2refXi = KaonMass*KaonMass;
  const double m2refSg = PionMass*PionMass;

  const Long64_t nev = tK18->GetEntries();
  Long64_t nout = 0;
  for( Long64_t iev = 0; iev < nev; ++iev ){
    tK18->GetEntry( iev );
    if( src.ntK18 != 1 ) continue;
    if( src.nhBH1 >= 1 ){
      double best = src.btof[0];
      for( int i=1; i<src.nhBH1 && i<MaxHits; ++i )
	if( std::fabs( src.btof[i] ) < std::fabs( best ) ) best = src.btof[i];
      hBtof->Fill( best );
    }
    if( tSks->GetEntryWithIndex( src.evnum ) <= 0 ) continue;
    if( src.evnum != src.sevnum || src.ntSks < 1 ) continue;
    if( tHodo->GetEntryWithIndex( src.evnum ) <= 0 ) continue;
    if( src.evnum != src.hevnum ) continue;

    dst.runnum = runnum;
    dst.evnum  = src.evnum;
    for( int i=0; i<NumOfMisc; ++i ){
      dst.trigpat[i]  = src.trigpat[i];
      dst.trigflag[i] = src.trigflag[i];
    }
    dst.pb = src.p_3rd[0];
    dst.xb = src.xtgtb[0]; dst.yb = src.ytgtb[0];
    dst.ub = src.utgtb[0]; dst.vb = src.vtgtb[0];
    dst.ntSks = src.ntSks;
    dst.ps = src.p[0];
    dst.chisqrSks = src.chisqrSks[0];
    dst.path = src.path[0];
    dst.xs = src.xtgts[0]; dst.ys = src.ytgts[0];
    dst.us = src.utgts[0]; dst.vs = src.vtgts[0];

    // vertex of the two straight tracks in the target frame
    ThreeVector xk( dst.xb, dst.yb, 0. );
    ThreeVector pk( dst.ub, dst.vb, 1. );
    ThreeVector xs( dst.xs, dst.ys, 0. );
    ThreeVector ps( dst.us, dst.vs, 1. );
    ThreeVector vert = VertexPoint( xk, xs, pk, ps );
    dst.closedist = closeDist( xk, xs, pk, ps );
    dst.vtx = vert.x(); dst.vty = vert.y(); dst.vtz = vert.z();
    const double cost =
      (dst.ub*dst.us + dst.vb*dst.vs + 1.)
      / std::sqrt( (dst.ub*dst.ub + dst.vb*dst.vb + 1.)
		   *(dst.us*dst.us + dst.vs*dst.vs + 1.) );
    dst.cost = cost;
    dst.theta = std::acos( std::min( cost, 1. ) )*180./M_PI;

    // TOF m2 recomputed from the raw values (cf DstSksHodoscope)
    dst.nhTof = 0;
    for( int seg = 0; seg < NumOfSegTOF; ++seg ){
      if( src.tofut[seg] <= 0 || src.tofdt[seg] <= 0 ) continue;
      double tu, td, deu = 0., ded = 0.;
      if( !hodoMan.GetTime( CidTOF, 0, seg, 0, src.tofut[seg],
			    tu ) ) continue;
      if( !hodoMan.GetTime( CidTOF, 0, seg, 1, src.tofdt[seg],
			    td ) ) continue;
      if( src.tofua[seg] >= 0 )
	hodoMan.GetDe( CidTOF, 0, seg, 0, src.tofua[seg], deu );
      if( src.tofda[seg] >= 0 )
	hodoMan.GetDe( CidTOF, 0, seg, 1, src.tofda[seg], ded );
      double ctu = tu, ctd = td;
      phcMan.doCorrection( CidTOF, 0, seg, 0, tu, deu, ctu );
      phcMan.doCorrection( CidTOF, 0, seg, 1, td, ded, ctd );
      const double cmt = 0.5*(ctu + ctd);
      const double st  = cmt - src.time0 - tBeam;
      const int n = dst.nhTof;
      dst.TofSeg[n] = seg;
      dst.stof[n] = st;
      dst.deTof[n] = std::sqrt( std::max( deu, 0. )
				*std::max( ded, 0. ) );
      double beta = ( st > 0. ) ? dst.path/(st*LightVel) : -1.;
      dst.m2[n] = ( beta > 0. )
	? dst.ps*dst.ps*(1. - beta*beta)/(beta*beta) : -9999.;
      ++dst.nhTof;
    }
    tOut->Fill();
    ++nout;

    // ---- final calibrated missing mass (Xi- and Sigma-) ----
    if( dst.chisqrSks <= MAX_CHISQR && dst.nhTof >= 1 ){
      const bool passVertex =
	( dst.closedist < MAX_CLOSEDIST
	  && std::fabs( dst.vtz ) < MAX_VTZ
	  && dst.theta > MIN_THETA );
      const double pbCal = dst.pb + po( dst.pb );
      const double pb = pbCal - beamLineDp( KaonMass, pbCal )
	- dpLoss( KaonMass, pbCal, half, targetMat );
      const double eb = std::sqrt( pb*pb + KaonMass*KaonMass );

      // Xi- (via the (K,K) trigger) and Sigma-/Lambda (via the
      // (K,pi) trigger) each get the trigflag matching their own
      // production reaction, decoded in DetectorID.hh's TrigFlag
      // enum (see analysis-note.md 2026-07-22); this cuts the
      // residual pi- beam contamination on top of the m2-window
      // selection below.
      struct { double m2ref, mass, m2lo, m2hi, mRecoil;
	       TH1D *hb, *hv; TH2D *hv2; int trig; } hyp[2] = {
	{ m2refXi, KaonMass, 0.20, 0.32, mRecoilXi,
	  hXiBasic, hXiVertex, hXiVertexTh, TrigKK },
	{ m2refSg, PionMass, -0.06, 0.12, mRecoilSg,
	  hSgBasic, hSgVertex, hSgVertexTh, TrigKPi },
      };
      for( int ih = 0; ih < 2; ++ih ){
	if( !src.trigflag[hyp[ih].trig] ) continue;
	int best = -1;
	double bestd = 1e9;
	for( int j = 0; j < dst.nhTof; ++j ){
	  double d = std::fabs( dst.m2[j] - hyp[ih].m2ref );
	  if( d < bestd ){ bestd = d; best = j; }
	}
	if( best < 0 ) continue;
	const double m2 = dst.m2[best];
	if( !(hyp[ih].m2lo < m2 && m2 < hyp[ih].m2hi) ) continue;

	const double psCal = dst.ps + hPoly( dst.ps )
	  + kU( dst.us );
	const double ps = psCal
	  + dpLoss( hyp[ih].mass, psCal, half, targetMat );
	const double es = std::sqrt( ps*ps
				     + hyp[ih].mass*hyp[ih].mass );
	const double mm2 = (eb + hyp[ih].mRecoil - es)
	    *(eb + hyp[ih].mRecoil - es)
	  - (pb*pb + ps*ps - 2.*pb*ps*dst.cost);
	if( mm2 <= 0. ) continue;
	const double mm = std::sqrt( mm2 );
	hyp[ih].hb->Fill( mm );
	if( passVertex ){
	  hyp[ih].hv->Fill( mm );
	  hyp[ih].hv2->Fill( mm, dst.theta );
	}
      }
    }
  }
  std::cout << "#D " << argv[5] << ": " << nout << " pairs from "
	    << nev << " beam entries" << std::endl;

  fOut->cd();
  hXiBasic->Write();
  hXiVertex->Write();
  hXiVertexTh->Write();
  hSgBasic->Write();
  hSgVertex->Write();
  hSgVertexTh->Write();
  hBtof->Write();
  fOut->Write();
  fOut->Close();
  return EXIT_SUCCESS;
}
