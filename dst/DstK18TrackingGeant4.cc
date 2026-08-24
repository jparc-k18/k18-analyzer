// -*- C++ -*-

#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <algorithm>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include <TParticle.h>
#include <TNamed.h>
#include <TRandom.h>
#include <TString.h>
#include <TTreeReaderArray.h>
#include <TTreeReaderValue.h>
#include <UnpackerManager.hh>

#include "CatchSignal.hh"
#include "ConfMan.hh"
#include "DCAnalyzer.hh"
#include "DCGeomMan.hh"
#include "DCLocalTrack.hh"
#include "DetectorID.hh"
#include "DstHelper.hh"
#include "K18BeamlineRK.hh"
#include "K18TrackD2U.hh"
#include "K18TransMatrix.hh"
#include "RootHelper.hh"
#include "UserParamMan.hh"

namespace
{
using namespace root;
using namespace dst;
const std::string& class_name("DstK18TrackingGeant4");
const auto qnan = TMath::QuietNaN();
auto&       gConf = ConfMan::GetInstance();
const auto& gGeom = DCGeomMan::GetInstance();
const auto& gUnpacker = hddaq::unpacker::GUnpacker::get_instance();

constexpr Int_t kMaxBftCandidates = 16;
constexpr Int_t kBftHitMax = 64;

Bool_t
HasConfValue(const char* key)
{
  return !gConf.Get<TString>(key).IsNull();
}

Double_t
ConfDoubleOr(const char* key, Double_t fallback)
{
  return HasConfValue(key) ? gConf.Get<Double_t>(key) : fallback;
}

Bool_t
ConfFlagOr(const char* key, Bool_t fallback)
{
  return HasConfValue(key) ? std::abs(gConf.Get<Double_t>(key)) > 0.5 : fallback;
}

UInt_t
K18HitSmearSeed(Int_t evnum, std::uint64_t stream_tag)
{
  const auto base = static_cast<std::uint64_t>(
    std::llround(ConfDoubleOr("K18HitSmearSeed", 20260716.)));
  std::uint64_t z = base + static_cast<std::uint64_t>(evnum) + stream_tag +
                    0x9e3779b97f4a7c15ULL;
  z = (z ^ (z >> 30))*0xbf58476d1ce4e5b9ULL;
  z = (z ^ (z >> 27))*0x94d049bb133111ebULL;
  z ^= z >> 31;
  const UInt_t seed = static_cast<UInt_t>(z ^ (z >> 32));
  return seed == 0 ? 1U : seed;
}

void
SeedK18BftSmearing(Int_t evnum)
{
  static const Bool_t enabled =
    ConfDoubleOr("K18BFTPositionSmearSigma", 0.) > 0.;
  if(enabled) gRandom->SetSeed(K18HitSmearSeed(evnum, 0x424654ULL));
}

void
SeedK18DcSmearing(Int_t evnum)
{
  static const Bool_t enabled =
    ConfDoubleOr("G4DCSmearResolutionScale", 0.) > 0.;
  if(enabled) gRandom->SetSeed(K18HitSmearSeed(evnum, 0x4243ULL));
}

Double_t
CorrectK18NativeMomentum(const K18RKTrack& track)
{
  static const Bool_t use_corr = ConfFlagOr("K18NativeMomCorrEnable", false);
  if(!use_corr)
    return track.p;

  // Optional closure correction for the Geant4 native RK diagnostic.
  // Units: p [GeV/c], x/y/xbft [mm], u/v dimensionless.
  static const Double_t c0   = ConfDoubleOr("K18NativeMomCorrC0", 0.);
  static const Double_t c1   = ConfDoubleOr("K18NativeMomCorrC1", 1.);
  static const Double_t c2   = ConfDoubleOr("K18NativeMomCorrC2", 0.);
  static const Double_t cx   = ConfDoubleOr("K18NativeMomCorrX", 0.);
  static const Double_t cu   = ConfDoubleOr("K18NativeMomCorrU", 0.);
  static const Double_t cy   = ConfDoubleOr("K18NativeMomCorrY", 0.);
  static const Double_t cv   = ConfDoubleOr("K18NativeMomCorrV", 0.);
  static const Double_t cxx  = ConfDoubleOr("K18NativeMomCorrXX", 0.);
  static const Double_t cuu  = ConfDoubleOr("K18NativeMomCorrUU", 0.);
  static const Double_t cyy  = ConfDoubleOr("K18NativeMomCorrYY", 0.);
  static const Double_t cvv  = ConfDoubleOr("K18NativeMomCorrVV", 0.);
  static const Double_t cxu  = ConfDoubleOr("K18NativeMomCorrXU", 0.);
  static const Double_t cyv  = ConfDoubleOr("K18NativeMomCorrYV", 0.);
  static const Double_t cbft = ConfDoubleOr("K18NativeMomCorrXBFT", 0.);

  const Double_t p = track.p;
  const Double_t x = track.xvo;
  const Double_t u = track.uvo;
  const Double_t y = track.yvo;
  const Double_t v = track.vvo;
  return c0 + c1*p + c2*p*p
    + cx*x + cu*u + cy*y + cv*v
    + cxx*x*x + cuu*u*u + cyy*y*y + cvv*v*v
    + cxu*x*u + cyv*y*v + cbft*track.xbft;
}
}

namespace dst
{
enum kArgc
{
  kProcess, kConfFile,
  kK18Geant4, kOutFile, nArgc
};
std::vector<TString> ArgName =
{ "[Process]", "[ConfFile]", "[K18Geant4]", "[OutFile]" };
std::vector<TString> TreeName =
{ "", "", "g4s2s", "" };
std::vector<TFile*> TFileCont;
std::vector<TTree*> TTreeCont;
std::vector<TTreeReader*> TTreeReaderCont;
}

//_____________________________________________________________________
struct Event
{
  Int_t evnum;
  Double_t x0;
  Double_t y0;
  Double_t z0;
  Double_t xp0;
  Double_t yp0;
  Double_t phi0;
  Double_t theta0;
  Double_t p0;
  Double_t pB;
  Double_t xTgtTruth;
  Double_t yTgtTruth;
  Double_t zTgtTruth;
  Double_t uTgtTruth;
  Double_t vTgtTruth;
  Double_t pTgtTruth;
  Double_t pTruthK18;

  Int_t bh1_nhit;
  Int_t bh2_nhit;
  Bool_t beamTrigDefined;
  Bool_t beamTrig;

  Int_t bft_nhit;
  Int_t bft_ncl;
  Double_t bft_xmean;
  Double_t bft_xrms;
  Double_t bft_xcand[kMaxBftCandidates];

  Int_t nlBcOut;
  Int_t nhBcOutLayer[NumOfLayersBcOut];
  Int_t ntBcOut;
  Int_t nhBcOut[MaxHits];
  Double_t chisqrBcOut[MaxHits];
  Double_t x0BcOut[MaxHits];
  Double_t y0BcOut[MaxHits];
  Double_t u0BcOut[MaxHits];
  Double_t v0BcOut[MaxHits];

  Int_t ntK18;
  Int_t nhK18[MaxHits];
  Double_t chisqrK18[MaxHits];
  Double_t p_2nd[MaxHits];
  Double_t p_3rd[MaxHits];
  Double_t delta_2nd[MaxHits];
  Double_t delta_3rd[MaxHits];
  Double_t dpK18[MaxHits];
  Double_t relDpK18[MaxHits];

  Int_t ntK18RK;
  Int_t nhK18RK[MaxHits];
  Double_t chisqrK18RK[MaxHits];
  Double_t p_rk_raw[MaxHits];
  Double_t p_rk[MaxHits];
  Double_t p_rk_tgt_truthloss[MaxHits];
  Double_t delta_rk_raw[MaxHits];
  Double_t delta_rk[MaxHits];
  Double_t dpK18RK_raw[MaxHits];
  Double_t dpK18RK[MaxHits];
  Double_t dpK18RK_tgt_truthloss[MaxHits];
  Double_t relDpK18RK_raw[MaxHits];
  Double_t relDpK18RK[MaxHits];
  Double_t relDpK18RK_tgt_truthloss[MaxHits];
  Double_t bftResRK[MaxHits];
  Double_t xbftRK[MaxHits];
  Double_t xbftCalcRK[MaxHits];
  Double_t xoutRK[MaxHits];
  Double_t youtRK[MaxHits];
  Double_t uoutRK[MaxHits];
  Double_t voutRK[MaxHits];
  Double_t xtgtK18RK[MaxHits];
  Double_t ytgtK18RK[MaxHits];
  Double_t utgtK18RK[MaxHits];
  Double_t vtgtK18RK[MaxHits];

  Double_t xin[MaxHits];
  Double_t yin[MaxHits];
  Double_t uin[MaxHits];
  Double_t vin[MaxHits];
  Double_t xout[MaxHits];
  Double_t yout[MaxHits];
  Double_t uout[MaxHits];
  Double_t vout[MaxHits];

  Double_t xtgtK18[MaxHits];
  Double_t ytgtK18[MaxHits];
  Double_t utgtK18[MaxHits];
  Double_t vtgtK18[MaxHits];
  Double_t xbftK18[MaxHits];
  Double_t ybftK18[MaxHits];
  Double_t ubftK18[MaxHits];
  Double_t vbftK18[MaxHits];
  Double_t theta[MaxHits];
  Double_t phi[MaxHits];

  void clear();
};

//_____________________________________________________________________
void
Event::clear()
{
  evnum = 0;
  x0 = qnan;
  y0 = qnan;
  z0 = qnan;
  xp0 = qnan;
  yp0 = qnan;
  phi0 = qnan;
  theta0 = qnan;
  p0 = qnan;
  pB = qnan;
  xTgtTruth = qnan;
  yTgtTruth = qnan;
  zTgtTruth = qnan;
  uTgtTruth = qnan;
  vTgtTruth = qnan;
  pTgtTruth = qnan;
  pTruthK18 = qnan;

  bh1_nhit = 0;
  bh2_nhit = 0;
  beamTrigDefined = kFALSE;
  beamTrig = kFALSE;

  bft_nhit = 0;
  bft_ncl = 0;
  bft_xmean = qnan;
  bft_xrms = qnan;
  for(Int_t i=0; i<kMaxBftCandidates; ++i)
    bft_xcand[i] = qnan;

  nlBcOut = 0;
  ntBcOut = 0;
  ntK18 = 0;
  ntK18RK = 0;
  for(Int_t i=0; i<NumOfLayersBcOut; ++i)
    nhBcOutLayer[i] = 0;

  for(Int_t i=0; i<MaxHits; ++i){
    nhBcOut[i] = 0;
    chisqrBcOut[i] = qnan;
    x0BcOut[i] = qnan;
    y0BcOut[i] = qnan;
    u0BcOut[i] = qnan;
    v0BcOut[i] = qnan;

    nhK18[i] = 0;
    chisqrK18[i] = qnan;
    p_2nd[i] = qnan;
    p_3rd[i] = qnan;
    delta_2nd[i] = qnan;
    delta_3rd[i] = qnan;
    dpK18[i] = qnan;
    relDpK18[i] = qnan;

    nhK18RK[i] = 0;
    chisqrK18RK[i] = qnan;
    p_rk_raw[i] = qnan;
    p_rk[i] = qnan;
    p_rk_tgt_truthloss[i] = qnan;
    delta_rk_raw[i] = qnan;
    delta_rk[i] = qnan;
    dpK18RK_raw[i] = qnan;
    dpK18RK[i] = qnan;
    dpK18RK_tgt_truthloss[i] = qnan;
    relDpK18RK_raw[i] = qnan;
    relDpK18RK[i] = qnan;
    relDpK18RK_tgt_truthloss[i] = qnan;
    bftResRK[i] = qnan;
    xbftRK[i] = qnan;
    xbftCalcRK[i] = qnan;
    xoutRK[i] = qnan;
    youtRK[i] = qnan;
    uoutRK[i] = qnan;
    voutRK[i] = qnan;
    xtgtK18RK[i] = qnan;
    ytgtK18RK[i] = qnan;
    utgtK18RK[i] = qnan;
    vtgtK18RK[i] = qnan;

    xin[i] = qnan;
    yin[i] = qnan;
    uin[i] = qnan;
    vin[i] = qnan;
    xout[i] = qnan;
    yout[i] = qnan;
    uout[i] = qnan;
    vout[i] = qnan;

    xtgtK18[i] = qnan;
    ytgtK18[i] = qnan;
    utgtK18[i] = qnan;
    vtgtK18[i] = qnan;
    xbftK18[i] = qnan;
    ybftK18[i] = qnan;
    ubftK18[i] = qnan;
    vbftK18[i] = qnan;
    theta[i] = qnan;
    phi[i] = qnan;
  }
}

//_____________________________________________________________________
struct Src
{
  TTreeReaderValue<Int_t>* evnum;
  TTreeReaderValue<Double_t>* x0;
  TTreeReaderValue<Double_t>* y0;
  TTreeReaderValue<Double_t>* z0;
  TTreeReaderValue<Double_t>* xp0;
  TTreeReaderValue<Double_t>* yp0;
  TTreeReaderValue<Double_t>* phi0;
  TTreeReaderValue<Double_t>* theta0;
  TTreeReaderValue<Double_t>* p0;
  TTreeReaderValue<Double_t>* pB;
  TTreeReaderValue<Double_t>* xTgtTruth;
  TTreeReaderValue<Double_t>* yTgtTruth;
  TTreeReaderValue<Double_t>* zTgtTruth;
  TTreeReaderValue<Double_t>* uTgtTruth;
  TTreeReaderValue<Double_t>* vTgtTruth;
  TTreeReaderValue<Double_t>* pTgtTruth;
  TTreeReaderArray<TParticle>* K18BH1;
  TTreeReaderArray<TParticle>* K18BH2;
  TTreeReaderArray<TParticle>* K18BFT;
  TTreeReaderArray<TParticle>* K18BC;
};

namespace root
{
Event event;
Src src;
TH1 *h[MaxHist];
TTree *tree;
}

namespace
{
//_____________________________________________________________________
Double_t
TruthMomentumGeV(Double_t p0)
{
  if(!std::isfinite(p0)) return qnan;
  return (std::abs(p0) > 10.) ? p0/1000. : p0;
}

//_____________________________________________________________________
Double_t
ApplyTargetTruthLoss(Double_t p_reco, Double_t p0, Double_t p_tgt_truth)
{
  const Double_t p0_gev = TruthMomentumGeV(p0);
  if(!std::isfinite(p_reco) || !std::isfinite(p0_gev) ||
     !std::isfinite(p_tgt_truth))
    return qnan;
  return p_reco + (p_tgt_truth - p0_gev);
}

//_____________________________________________________________________
Bool_t
HasInputBranch(const char* name)
{
  return TTreeCont[kK18Geant4] && TTreeCont[kK18Geant4]->GetBranch(name);
}

//_____________________________________________________________________
Int_t
CountPrimaryBeamHits(const TTreeReaderArray<TParticle>* hits)
{
  if(!hits) return 0;
  const Int_t accepted_abs_pdg =
    static_cast<Int_t>(ConfDoubleOr("K18BFTAcceptedAbsPdg", 321.));
  Int_t count = 0;
  for(const auto& particle: *hits){
    if(particle.GetFirstMother() != 0) continue;
    if(accepted_abs_pdg > 0 &&
       std::abs(particle.GetPdgCode()) != accepted_abs_pdg) continue;
    ++count;
  }
  return count;
}

//_____________________________________________________________________
Bool_t
ExtractBftXCandidates(const TTreeReaderArray<TParticle>& bft,
                      std::vector<Double_t>& xCand)
{
  const Int_t id_bft_x  = gGeom.DetectorId("BFT-X");
  const Int_t id_bft_xp = gGeom.DetectorId("BFT-XP");
  const Int_t accepted_abs_pdg =
    static_cast<Int_t>(ConfDoubleOr("K18BFTAcceptedAbsPdg", 321.));
  static const Double_t cluster_gap =
    std::max(0., ConfDoubleOr("K18BFTClusterGap", 5.));
  static const Double_t smear_sigma =
    std::max(0., ConfDoubleOr("K18BFTPositionSmearSigma", 0.));
  std::vector<Double_t> xs;
  xs.reserve(kBftHitMax);

  for(const auto& particle: bft){
    if(particle.GetFirstMother()!=0) continue;
    const Int_t copy = particle.GetSecondMother();
    if(copy != id_bft_x && copy != id_bft_xp) continue;
    if(accepted_abs_pdg > 0 &&
       std::abs(particle.GetPdgCode()) != accepted_abs_pdg) continue;
    if(xs.size() >= kBftHitMax) break;
    xs.push_back(particle.Vx());
  }

  event.bft_nhit = xs.size();
  if(xs.empty()) return false;

  std::sort(xs.begin(), xs.end());

  Double_t sum = 0.;
  for(const auto x: xs) sum += x;
  event.bft_xmean = sum/Double_t(xs.size());

  Double_t sum2 = 0.;
  for(const auto x: xs) sum2 += (x - event.bft_xmean)*(x - event.bft_xmean);
  event.bft_xrms = TMath::Sqrt(sum2/Double_t(xs.size()));

  Double_t cluster_sum = 0.;
  Int_t cluster_n = 0;
  Double_t prev_x = qnan;
  auto flush_cluster = [&]() {
    if(cluster_n <= 0 || event.bft_ncl >= kMaxBftCandidates)
      return;
    Double_t x = cluster_sum/Double_t(cluster_n);
    if(smear_sigma > 0.) x = gRandom->Gaus(x, smear_sigma);
    event.bft_xcand[event.bft_ncl++] = x;
    xCand.push_back(x);
  };

  for(const auto x: xs){
    if(cluster_n > 0 && std::isfinite(prev_x) &&
       cluster_gap > 0. && std::abs(x - prev_x) > cluster_gap){
      flush_cluster();
      cluster_sum = 0.;
      cluster_n = 0;
    }
    cluster_sum += x;
    ++cluster_n;
    prev_x = x;
  }
  flush_cluster();

  if(xCand.empty()){
    event.bft_ncl = 1;
    event.bft_xcand[0] = event.bft_xmean;
    xCand.push_back(event.bft_xmean);
  }
  return true;
}

//_____________________________________________________________________
void
FillBcOutTracks(const DCAnalyzer& DCAna)
{
  Int_t ntBcOut = DCAna.GetNtracksBcOut();
  if(MaxHits < ntBcOut){
    std::cout << "#W too many ntBcOut " << ntBcOut << "/" << MaxHits << std::endl;
    ntBcOut = MaxHits;
  }
  event.ntBcOut = ntBcOut;
  HF1(21, Double_t(ntBcOut));

  for(Int_t it=0; it<ntBcOut; ++it){
    const auto& track = DCAna.GetTrackBcOut(it);
    if(!track) continue;
    const Int_t nh = track->GetNHit();
    event.nhBcOut[it] = nh;
    event.chisqrBcOut[it] = track->GetChiSquare();
    event.x0BcOut[it] = track->GetX(0.);
    event.y0BcOut[it] = track->GetY(0.);
    event.u0BcOut[it] = track->GetU0();
    event.v0BcOut[it] = track->GetV0();
    HF1(22, Double_t(nh));
    HF1(23, event.chisqrBcOut[it]);
  }
}

//_____________________________________________________________________
void
FillK18Tracks(const DCAnalyzer& DCAna)
{
  Int_t ntK18 = DCAna.GetNTracksK18D2U();
  if(MaxHits < ntK18){
    std::cout << "#W too many ntK18 " << ntK18 << "/" << MaxHits << std::endl;
    ntK18 = MaxHits;
  }
  event.ntK18 = ntK18;
  HF1(30, Double_t(ntK18));

  for(Int_t i=0; i<ntK18; ++i){
    const auto& track = DCAna.GetK18TrackD2U(i);
    if(!track) continue;
    const auto& ltrack = track->TrackOut();
    if(!ltrack) continue;

    event.nhK18[i] = ltrack->GetNHit();
    event.chisqrK18[i] = ltrack->GetChiSquare();
    event.p_2nd[i] = track->P();
    event.p_3rd[i] = track->P3rd();
    event.delta_2nd[i] = track->Delta();
    event.delta_3rd[i] = track->Delta3rd();
    event.dpK18[i] = event.p_3rd[i] - event.pTruthK18;
    event.relDpK18[i] = event.dpK18[i]/event.pTruthK18;

    event.xin[i] = track->Xin();
    event.yin[i] = track->Yin();
    event.uin[i] = track->Uin();
    event.vin[i] = track->Vin();
    event.xout[i] = track->Xout();
    event.yout[i] = track->Yout();
    event.uout[i] = track->Uout();
    event.vout[i] = track->Vout();

    event.xtgtK18[i] = track->Xtgt();
    event.ytgtK18[i] = track->Ytgt();
    event.utgtK18[i] = track->Utgt();
    event.vtgtK18[i] = track->Vtgt();
    event.xbftK18[i] = track->Xbft();
    event.ybftK18[i] = track->Ybft();
    event.ubftK18[i] = track->Ubft();
    event.vbftK18[i] = track->Vbft();
    event.theta[i] = ltrack->GetTheta();
    event.phi[i] = ltrack->GetPhi();

    HF1(31, event.p_3rd[i]);
    HF1(32, event.dpK18[i]);
    HF1(33, event.relDpK18[i]);
    HF2(34, event.pTruthK18, event.dpK18[i]);
  }
}

//_____________________________________________________________________
void
FillK18RKTracks(const DCAnalyzer& DCAna, const std::vector<Double_t>& xCand)
{
  static const K18BeamlineRK rk;
  static const Bool_t use_full_fit =
    !gConf.Get<TString>("K18NativeUseFullFit").IsNull() &&
    std::abs(gConf.Get<Double_t>("K18NativeUseFullFit")) > 0.5;
  const Double_t pk18 = std::abs(gConf.Get<Double_t>("PK18"));
  Int_t n = 0;

  for(const auto xbft: xCand){
    for(Int_t it=0, nt=DCAna.GetNtracksBcOut(); it<nt; ++it){
      if(n >= MaxHits) return;
      const auto& local = DCAna.GetTrackBcOut(it);
      if(!local) continue;
      K18RKTrack track;
      const Bool_t ok = use_full_fit ?
        rk.FitTrack(*local, xbft, track) :
        rk.FitMomentum(local->GetX(0.), local->GetY(0.),
                       local->GetU0(), local->GetV0(), xbft, track);
      if(!ok)
        continue;

      event.nhK18RK[n] = local->GetNHit();
      event.chisqrK18RK[n] = use_full_fit ? track.chisqr : local->GetChiSquare();
      event.p_rk_raw[n] = track.p;
      event.delta_rk_raw[n] = (pk18 > 0.) ? track.p/pk18 - 1. : qnan;
      event.dpK18RK_raw[n] = track.p - event.pTruthK18;
      event.relDpK18RK_raw[n] =
        (event.pTruthK18 != 0.) ? event.dpK18RK_raw[n]/event.pTruthK18 : qnan;
      event.p_rk[n] = CorrectK18NativeMomentum(track);
      event.delta_rk[n] = (pk18 > 0.) ? event.p_rk[n]/pk18 - 1. : qnan;
      event.dpK18RK[n] = event.p_rk[n] - event.pTruthK18;
      event.relDpK18RK[n] = (event.pTruthK18 != 0.) ? event.dpK18RK[n]/event.pTruthK18 : qnan;
      event.p_rk_tgt_truthloss[n] =
        ApplyTargetTruthLoss(event.p_rk[n], event.p0, event.pTruthK18);
      event.dpK18RK_tgt_truthloss[n] =
        event.p_rk_tgt_truthloss[n] - event.pTruthK18;
      event.relDpK18RK_tgt_truthloss[n] =
        (event.pTruthK18 != 0.) ?
        event.dpK18RK_tgt_truthloss[n]/event.pTruthK18 : qnan;
      event.bftResRK[n] = track.bft_residual;
      event.xbftRK[n] = track.xbft;
      event.xbftCalcRK[n] = track.xbft_calc;
      event.xoutRK[n] = track.xvo;
      event.youtRK[n] = track.yvo;
      event.uoutRK[n] = track.uvo;
      event.voutRK[n] = track.vvo;
      event.xtgtK18RK[n] = track.xtgt;
      event.ytgtK18RK[n] = track.ytgt;
      event.utgtK18RK[n] = track.utgt;
      event.vtgtK18RK[n] = track.vtgt;

      HF1(41, event.p_rk[n]);
      HF1(42, event.dpK18RK[n]);
      HF1(43, event.relDpK18RK[n]);
      HF1(44, event.bftResRK[n]);
      HF2(45, event.pTruthK18, event.dpK18RK[n]);
      ++n;
    }
  }

  event.ntK18RK = n;
  HF1(40, Double_t(n));
}
}

//_____________________________________________________________________
int
main(int argc, char **argv)
{
  std::vector<std::string> arg(argv, argv+argc);
  if(!CheckArg(arg))
    return EXIT_FAILURE;
  if(!DstOpen(arg))
    return EXIT_FAILURE;
  if(!gConf.Initialize(arg[kConfFile]))
    return EXIT_FAILURE;
  if(!gConf.InitializeUnpacker())
    return EXIT_FAILURE;

  Int_t skip = gUnpacker.get_skip();
  if(skip < 0) skip = 0;
  Int_t max_loop = gUnpacker.get_max_loop();
  Int_t nevent = GetEntries(TTreeCont);
  if(max_loop > 0) nevent = skip + max_loop;
  CatchSignal::Set();

  Int_t ievent = skip;
  for(; ievent<nevent && !CatchSignal::Stop(); ++ievent){
    InitializeEvent();
    if(DstRead(ievent)) tree->Fill();
  }
  std::cout << "#D Event Number: " << std::setw(6) << ievent << std::endl;

  DstClose();
  return EXIT_SUCCESS;
}

//_____________________________________________________________________
Bool_t
dst::InitializeEvent()
{
  event.clear();
  return true;
}

//_____________________________________________________________________
Bool_t
dst::DstOpen(std::vector<std::string> arg)
{
  Int_t open_file = 0;
  Int_t open_tree = 0;
  for(std::size_t i=0; i<nArgc; ++i){
    if(i==kProcess || i==kConfFile || i==kOutFile) continue;
    open_file += OpenFile(TFileCont[i], arg[i]);
    open_tree += OpenTree(TFileCont[i], TTreeCont[i], TreeName[i]);
  }

  if(open_file!=open_tree || open_file!=nArgc-3)
    return false;
  if(!CheckEntries(TTreeCont))
    return false;

  TFileCont[kOutFile] = new TFile(arg[kOutFile].c_str(), "recreate");
  return true;
}

//_____________________________________________________________________
Bool_t
dst::DstRead(Int_t ievent)
{
  if(ievent%10000==0){
    std::cout << "#D Event Number: " << std::setw(6) << ievent << std::endl;
  }

  GetEntry(ievent);

  HF1(1, 0.);
  event.evnum = **src.evnum;
  event.x0 = **src.x0;
  event.y0 = **src.y0;
  event.z0 = **src.z0;
  event.xp0 = **src.xp0;
  event.yp0 = **src.yp0;
  event.phi0 = **src.phi0;
  event.theta0 = **src.theta0;
  event.p0 = **src.p0;
  event.pB = **src.pB;
  if(src.xTgtTruth) event.xTgtTruth = **src.xTgtTruth;
  if(src.yTgtTruth) event.yTgtTruth = **src.yTgtTruth;
  if(src.zTgtTruth) event.zTgtTruth = **src.zTgtTruth;
  if(src.uTgtTruth) event.uTgtTruth = **src.uTgtTruth;
  if(src.vTgtTruth) event.vTgtTruth = **src.vTgtTruth;
  if(src.pTgtTruth) event.pTgtTruth = **src.pTgtTruth;
  event.pTruthK18 = TruthMomentumGeV(
    std::isfinite(event.pTgtTruth) ? event.pTgtTruth : event.p0);
  event.bh1_nhit = CountPrimaryBeamHits(src.K18BH1);
  event.bh2_nhit = CountPrimaryBeamHits(src.K18BH2);
  event.beamTrigDefined = src.K18BH1 && src.K18BH2;
  event.beamTrig = event.bh1_nhit > 0 && event.bh2_nhit > 0;

  std::vector<Double_t> xCand;
  // The merged runmanager ROOT keeps the chunk-local Geant4 evnum, so evnum
  // repeats every chunk.  Use the global input entry for independent,
  // split-invariant detector smearing across the merged production.
  SeedK18BftSmearing(ievent);
  ExtractBftXCandidates(*src.K18BFT, xCand);
  HF1(10, Double_t(event.bft_nhit));
  if(event.bft_ncl > 0) HF1(11, event.bft_xcand[0]);

  DCAnalyzer DCAna;
  SeedK18DcSmearing(ievent);
  DCAna.DecodeBcOutHitsGeant4(*src.K18BC);

  for(Int_t l=0; l<NumOfLayersBcOut; ++l){
    const Int_t nh = DCAna.GetBcOutHC(l).size();
    event.nhBcOutLayer[l] = nh;
    if(nh > 0) ++event.nlBcOut;
  }
  HF1(20, Double_t(event.nlBcOut));

  HF1(1, 10.);
  DCAna.TrackSearchBcOut();
  DCAna.ChiSqrCutBcOut(20.);
  FillBcOutTracks(DCAna);

  if(event.ntBcOut > 0 && !xCand.empty()){
    HF1(1, 20.);
    DCAna.TrackSearchK18D2U(xCand);
    FillK18Tracks(DCAna);
    FillK18RKTracks(DCAna, xCand);
  }

  HF1(1, 29.);
  return true;
}

//_____________________________________________________________________
Bool_t
dst::DstClose()
{
  TNamed("k18_hit_smearing",
         Form("G4DCSmearResolutionScale=%.12g; K18BFTPositionSmearSigma_mm=%.12g; "
              "K18HitSmearSeed=%.0f; per-event detector-specific deterministic streams; "
              "BC dl~Gaus(abs(s-wire),scale*DCGEO.Res) after truth wire assignment and before DL gate; "
              "BFT xcand=true-hit cluster mean+Gaus(0,sigma), bft_xmean/xrms remain unsmeared",
              ConfDoubleOr("G4DCSmearResolutionScale", 0.),
              ConfDoubleOr("K18BFTPositionSmearSigma", 0.),
              ConfDoubleOr("K18HitSmearSeed", 20260716.))).Write();
  TFileCont[kOutFile]->Write();
  std::cout << "#D Close : " << TFileCont[kOutFile]->GetName() << std::endl;
  TFileCont[kOutFile]->Close();

  const std::size_t n = TFileCont.size();
  for(std::size_t i=0; i<n; ++i){
    if(TTreeCont[i]) delete TTreeCont[i];
    if(TFileCont[i]) delete TFileCont[i];
  }
  return true;
}

//_____________________________________________________________________
Bool_t
ConfMan::InitializeHistograms()
{
  HB1(1,  "status", 40, 0., 40.);
  HB1(10, "BFT nhit", 20, 0., 20.);
  HB1(11, "BFT x candidate", 400, -200., 200.);
  HB1(20, "BcOut hit layers", NumOfLayersBcOut+1, 0., NumOfLayersBcOut+1.);
  HB1(21, "BcOut tracks", 20, 0., 20.);
  HB1(22, "BcOut track nhit", 20, 0., 20.);
  HB1(23, "BcOut track chisqr", 200, 0., 100.);
  HB1(30, "K18 D2U tracks", 20, 0., 20.);
  HB1(31, "K18 p3rd", 400, 1.0, 1.8);
  HB1(32, "K18 p3rd-ptruth", 400, -0.2, 0.2);
  HB1(33, "K18 relative dp", 400, -0.15, 0.15);
  HB2(34, "K18 dp%truth momentum", 200, 1.0, 1.8, 400, -0.2, 0.2);
  HB1(40, "K18 RK tracks", 20, 0., 20.);
  HB1(41, "K18 RK p", 400, 1.0, 1.8);
  HB1(42, "K18 RK p-ptruth", 400, -0.2, 0.2);
  HB1(43, "K18 RK relative dp", 400, -0.15, 0.15);
  HB1(44, "K18 RK BFT residual", 400, -20., 20.);
  HB2(45, "K18 RK dp%truth momentum", 200, 1.0, 1.8, 400, -0.2, 0.2);

  HBTree("k18track", "tree of K18 D2U tracking using Geant4 K18BFT/K18BC input");
  tree->Branch("evnum",  &event.evnum,  "evnum/I");
  tree->Branch("x0",     &event.x0,     "x0/D");
  tree->Branch("y0",     &event.y0,     "y0/D");
  tree->Branch("z0",     &event.z0,     "z0/D");
  tree->Branch("xp0",    &event.xp0,    "xp0/D");
  tree->Branch("yp0",    &event.yp0,    "yp0/D");
  tree->Branch("phi0",   &event.phi0,   "phi0/D");
  tree->Branch("theta0", &event.theta0, "theta0/D");
  tree->Branch("p0",     &event.p0,     "p0/D");
  tree->Branch("pB",     &event.pB,     "pB/D");
  tree->Branch("xTgtTruth", &event.xTgtTruth, "xTgtTruth/D");
  tree->Branch("yTgtTruth", &event.yTgtTruth, "yTgtTruth/D");
  tree->Branch("zTgtTruth", &event.zTgtTruth, "zTgtTruth/D");
  tree->Branch("uTgtTruth", &event.uTgtTruth, "uTgtTruth/D");
  tree->Branch("vTgtTruth", &event.vTgtTruth, "vTgtTruth/D");
  tree->Branch("pTgtTruth", &event.pTgtTruth, "pTgtTruth/D");
  tree->Branch("pTruthK18", &event.pTruthK18, "pTruthK18/D");

  tree->Branch("bh1_nhit", &event.bh1_nhit, "bh1_nhit/I");
  tree->Branch("bh2_nhit", &event.bh2_nhit, "bh2_nhit/I");
  tree->Branch("beamTrigDefined", &event.beamTrigDefined, "beamTrigDefined/O");
  tree->Branch("beamTrig", &event.beamTrig, "beamTrig/O");

  tree->Branch("bft_nhit", &event.bft_nhit, "bft_nhit/I");
  tree->Branch("bft_ncl", &event.bft_ncl, "bft_ncl/I");
  tree->Branch("bft_xmean", &event.bft_xmean, "bft_xmean/D");
  tree->Branch("bft_xrms", &event.bft_xrms, "bft_xrms/D");
  tree->Branch("bft_xcand", event.bft_xcand, "bft_xcand[bft_ncl]/D");

  tree->Branch("nlBcOut", &event.nlBcOut, "nlBcOut/I");
  tree->Branch("nhBcOutLayer", event.nhBcOutLayer,
               Form("nhBcOutLayer[%d]/I", NumOfLayersBcOut));
  tree->Branch("ntBcOut", &event.ntBcOut, "ntBcOut/I");
  tree->Branch("nhBcOut", event.nhBcOut, "nhBcOut[ntBcOut]/I");
  tree->Branch("chisqrBcOut", event.chisqrBcOut, "chisqrBcOut[ntBcOut]/D");
  tree->Branch("x0BcOut", event.x0BcOut, "x0BcOut[ntBcOut]/D");
  tree->Branch("y0BcOut", event.y0BcOut, "y0BcOut[ntBcOut]/D");
  tree->Branch("u0BcOut", event.u0BcOut, "u0BcOut[ntBcOut]/D");
  tree->Branch("v0BcOut", event.v0BcOut, "v0BcOut[ntBcOut]/D");

  tree->Branch("ntK18", &event.ntK18, "ntK18/I");
  tree->Branch("nhK18", event.nhK18, "nhK18[ntK18]/I");
  tree->Branch("chisqrK18", event.chisqrK18, "chisqrK18[ntK18]/D");
  tree->Branch("p_2nd", event.p_2nd, "p_2nd[ntK18]/D");
  tree->Branch("p_3rd", event.p_3rd, "p_3rd[ntK18]/D");
  tree->Branch("delta_2nd", event.delta_2nd, "delta_2nd[ntK18]/D");
  tree->Branch("delta_3rd", event.delta_3rd, "delta_3rd[ntK18]/D");
  tree->Branch("dpK18", event.dpK18, "dpK18[ntK18]/D");
  tree->Branch("relDpK18", event.relDpK18, "relDpK18[ntK18]/D");

  tree->Branch("ntK18RK", &event.ntK18RK, "ntK18RK/I");
  tree->Branch("nhK18RK", event.nhK18RK, "nhK18RK[ntK18RK]/I");
  tree->Branch("chisqrK18RK", event.chisqrK18RK, "chisqrK18RK[ntK18RK]/D");
  tree->Branch("p_rk_raw", event.p_rk_raw, "p_rk_raw[ntK18RK]/D");
  tree->Branch("p_rk", event.p_rk, "p_rk[ntK18RK]/D");
  tree->Branch("p_rk_tgt_truthloss", event.p_rk_tgt_truthloss,
               "p_rk_tgt_truthloss[ntK18RK]/D");
  tree->Branch("delta_rk_raw", event.delta_rk_raw, "delta_rk_raw[ntK18RK]/D");
  tree->Branch("delta_rk", event.delta_rk, "delta_rk[ntK18RK]/D");
  tree->Branch("dpK18RK_raw", event.dpK18RK_raw, "dpK18RK_raw[ntK18RK]/D");
  tree->Branch("dpK18RK", event.dpK18RK, "dpK18RK[ntK18RK]/D");
  tree->Branch("dpK18RK_tgt_truthloss", event.dpK18RK_tgt_truthloss,
               "dpK18RK_tgt_truthloss[ntK18RK]/D");
  tree->Branch("relDpK18RK_raw", event.relDpK18RK_raw, "relDpK18RK_raw[ntK18RK]/D");
  tree->Branch("relDpK18RK", event.relDpK18RK, "relDpK18RK[ntK18RK]/D");
  tree->Branch("relDpK18RK_tgt_truthloss", event.relDpK18RK_tgt_truthloss,
               "relDpK18RK_tgt_truthloss[ntK18RK]/D");
  tree->Branch("bftResRK", event.bftResRK, "bftResRK[ntK18RK]/D");
  tree->Branch("xbftRK", event.xbftRK, "xbftRK[ntK18RK]/D");
  tree->Branch("xbftCalcRK", event.xbftCalcRK, "xbftCalcRK[ntK18RK]/D");
  tree->Branch("xoutRK", event.xoutRK, "xoutRK[ntK18RK]/D");
  tree->Branch("youtRK", event.youtRK, "youtRK[ntK18RK]/D");
  tree->Branch("uoutRK", event.uoutRK, "uoutRK[ntK18RK]/D");
  tree->Branch("voutRK", event.voutRK, "voutRK[ntK18RK]/D");
  tree->Branch("xtgtK18RK", event.xtgtK18RK, "xtgtK18RK[ntK18RK]/D");
  tree->Branch("ytgtK18RK", event.ytgtK18RK, "ytgtK18RK[ntK18RK]/D");
  tree->Branch("utgtK18RK", event.utgtK18RK, "utgtK18RK[ntK18RK]/D");
  tree->Branch("vtgtK18RK", event.vtgtK18RK, "vtgtK18RK[ntK18RK]/D");

  tree->Branch("xin", event.xin, "xin[ntK18]/D");
  tree->Branch("yin", event.yin, "yin[ntK18]/D");
  tree->Branch("uin", event.uin, "uin[ntK18]/D");
  tree->Branch("vin", event.vin, "vin[ntK18]/D");
  tree->Branch("xout", event.xout, "xout[ntK18]/D");
  tree->Branch("yout", event.yout, "yout[ntK18]/D");
  tree->Branch("uout", event.uout, "uout[ntK18]/D");
  tree->Branch("vout", event.vout, "vout[ntK18]/D");

  tree->Branch("xtgtK18", event.xtgtK18, "xtgtK18[ntK18]/D");
  tree->Branch("ytgtK18", event.ytgtK18, "ytgtK18[ntK18]/D");
  tree->Branch("utgtK18", event.utgtK18, "utgtK18[ntK18]/D");
  tree->Branch("vtgtK18", event.vtgtK18, "vtgtK18[ntK18]/D");
  tree->Branch("xbftK18", event.xbftK18, "xbftK18[ntK18]/D");
  tree->Branch("ybftK18", event.ybftK18, "ybftK18[ntK18]/D");
  tree->Branch("ubftK18", event.ubftK18, "ubftK18[ntK18]/D");
  tree->Branch("vbftK18", event.vbftK18, "vbftK18[ntK18]/D");
  tree->Branch("theta", event.theta, "theta[ntK18]/D");
  tree->Branch("phi", event.phi, "phi[ntK18]/D");

  TTreeReaderCont[kK18Geant4] = new TTreeReader("g4s2s", TFileCont[kK18Geant4]);
  const auto& reader = TTreeReaderCont[kK18Geant4];
  src.evnum  = new TTreeReaderValue<Int_t>(*reader,     "evnum");
  src.x0     = new TTreeReaderValue<Double_t>(*reader,  "x0");
  src.y0     = new TTreeReaderValue<Double_t>(*reader,  "y0");
  src.z0     = new TTreeReaderValue<Double_t>(*reader,  "z0");
  src.xp0    = new TTreeReaderValue<Double_t>(*reader,  "xp0");
  src.yp0    = new TTreeReaderValue<Double_t>(*reader,  "yp0");
  src.phi0   = new TTreeReaderValue<Double_t>(*reader,  "phi0");
  src.theta0 = new TTreeReaderValue<Double_t>(*reader,  "theta0");
  src.p0     = new TTreeReaderValue<Double_t>(*reader,  "p0");
  src.pB     = new TTreeReaderValue<Double_t>(*reader,  "pB");
  src.xTgtTruth =
    HasInputBranch("xTgtTruth") ? new TTreeReaderValue<Double_t>(*reader, "xTgtTruth") : nullptr;
  src.yTgtTruth =
    HasInputBranch("yTgtTruth") ? new TTreeReaderValue<Double_t>(*reader, "yTgtTruth") : nullptr;
  src.zTgtTruth =
    HasInputBranch("zTgtTruth") ? new TTreeReaderValue<Double_t>(*reader, "zTgtTruth") : nullptr;
  src.uTgtTruth =
    HasInputBranch("uTgtTruth") ? new TTreeReaderValue<Double_t>(*reader, "uTgtTruth") : nullptr;
  src.vTgtTruth =
    HasInputBranch("vTgtTruth") ? new TTreeReaderValue<Double_t>(*reader, "vTgtTruth") : nullptr;
  src.pTgtTruth =
    HasInputBranch("pTgtTruth") ? new TTreeReaderValue<Double_t>(*reader, "pTgtTruth") : nullptr;
  src.K18BH1 =
    HasInputBranch("K18BH1") ? new TTreeReaderArray<TParticle>(*reader, "K18BH1") : nullptr;
  src.K18BH2 =
    HasInputBranch("K18BH2") ? new TTreeReaderArray<TParticle>(*reader, "K18BH2") : nullptr;
  src.K18BFT = new TTreeReaderArray<TParticle>(*reader, "K18BFT");
  src.K18BC  = new TTreeReaderArray<TParticle>(*reader, "K18BC");
  return true;
}

//_____________________________________________________________________
Bool_t
ConfMan::InitializeParameterFiles()
{
  return
    (InitializeParameter<DCGeomMan>("DCGEO")      &&
     InitializeParameter<K18TransMatrix>("K18TM") &&
     InitializeParameter<UserParamMan>("USER"));
}

//_____________________________________________________________________
Bool_t
ConfMan::FinalizeProcess()
{
  return true;
}
