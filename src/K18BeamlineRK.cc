// -*- C++ -*-

#include "K18BeamlineRK.hh"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

#include "ConfMan.hh"
#include "DCGeomMan.hh"
#include "DCLocalTrack.hh"
#include "DCLTrackHit.hh"
#include "DebugCounter.hh"
#include "K18FieldMap.hh"

namespace
{
const auto& gConf = ConfMan::GetInstance();
const auto& gGeom = DCGeomMan::GetInstance();

const Double_t kDeg2Rad = std::acos(-1.)/180.;

const Double_t kRhoD4       = 4000.;
const Double_t kBendAngleD4 = 64.000;
const Double_t kD4SpaceX    = 657.;
const Double_t kD4HalfGap   = 100.;
const Double_t kD4B0        = 11.5137;
const Double_t kAlphaD4     = 23.5;

const Double_t kQ10z  = 956.;
const Double_t kQ10B0 = -7.30734;
const Double_t kQ10a0 = 100.;

const Double_t kQ11z  = 660.;
const Double_t kQ11B0 = 5.60272;
const Double_t kQ11a0 = 100.;

const Double_t kQ12z  = 660.;
const Double_t kQ12B0 = -5.99934;
const Double_t kQ12a0 = 100.;

const Double_t kQ13z  = 660.;
const Double_t kQ13B0 = 5.00944;
const Double_t kQ13a0 = 100.;

const Double_t kDriftL0 = 272.; // VI to Q10 upstream EFB (magnet surface +102)
const Double_t kDriftL1 = 242.;
const Double_t kDriftL2 = 470.;
const Double_t kDriftL3 = 470.;
const Double_t kDriftL4 = 240.;
const Double_t kDriftL5 = 270.; // Q13 downstream EFB to VO (magnet surface +100)

inline Double_t d4_tan_half()
{
  return kRhoD4*std::tan(0.5*kBendAngleD4*kDeg2Rad);
}

Bool_t field_map_enabled(const TString& value)
{
  TString normalized(value);
  normalized.ToLower();
  return !normalized.IsNull() && normalized != "none" &&
    normalized != "skip" && normalized != "0";
}

}

//_____________________________________________________________________________
const TString&
K18BeamlineRK::ClassName()
{
  static TString s_name("K18BeamlineRK");
  return s_name;
}

//_____________________________________________________________________________
K18BeamlineRK::K18BeamlineRK()
  : m_field_map(nullptr),
    m_p_min(ConfDoubleOr("K18NativePMin", 0.6)),
    m_p_max(ConfDoubleOr("K18NativePMax", 1.8)),
    m_rk_step(std::max(0.1, ConfDoubleOr("K18NativeRKStep", 1.))),
    m_fit_max_iteration(std::max(
      8, static_cast<Int_t>(ConfDoubleOr("K18NativeFitMaxIteration", 60.)))),
    m_fit_momentum_tolerance(std::max(
      1.e-10, ConfDoubleOr("K18NativeFitMomentumTolerance", 1.e-7))),
    m_fit_residual_tolerance(std::max(0., ConfDoubleOr("K18NativeFitResidualTolerance", 1.e-4))),
    m_max_path(std::max(1., ConfDoubleOr("K18NativeMaxPath", 12000.))),
    m_bft_residual_max(ConfDoubleOr("K18NativeBFTResidualMax", 100.)),
    m_charge(ConfDoubleOr("K18RKCharge", -1.)),
    m_global_scale(ConfDoubleOr("K18GlobalScale", 1.0)),
    m_rk_drift_step(std::max(1., ConfDoubleOr("K18NativeRKDriftStep", 10.))),
    m_sect_cot_alpha(0.),
    m_sect_chord_x(0.),
    m_cos_theta(),
    m_sin_theta(),
    m_magnet()
{
  const TString field_map_file = gConf.Get<TString>("K18FLDMAP");
  if(field_map_enabled(field_map_file)){
    const Double_t value_nmr = ConfDoubleOr("K18FLDNMR", 1.);
    const Double_t value_calc = ConfDoubleOr("K18FLDCALC", 1.);
    m_field_map.reset(new K18FieldMap(field_map_file, value_nmr, value_calc));
    if(!m_field_map->Initialize())
      throw std::runtime_error("K18BeamlineRK failed to initialize K18FLDMAP");
  }

  m_magnet[0].type = kQuad;
  m_magnet[0].name = "K18Q10";
  m_magnet[0].scale = MagnetScale("K18Q10");
  m_magnet[0].theta = 0.;
  m_magnet[0].quad = {kQ10B0, kQ10a0, kQ10z,
    -(d4_tan_half() + kDriftL2 + kQ11z + kDriftL1 + 0.5*kQ10z)};

  m_magnet[1].type = kQuad;
  m_magnet[1].name = "K18Q11";
  m_magnet[1].scale = MagnetScale("K18Q11");
  m_magnet[1].theta = 0.;
  m_magnet[1].quad = {kQ11B0, kQ11a0, kQ11z,
    -(d4_tan_half() + kDriftL2 + 0.5*kQ11z)};

  m_magnet[2].type = kSect;
  m_magnet[2].name = "K18D4";
  m_magnet[2].scale = MagnetScale("K18D4");
  m_magnet[2].theta = 0.;
  m_magnet[2].sect = {kD4B0, kRhoD4, kD4SpaceX, kD4HalfGap,
    TVector3(-d4_tan_half(), kRhoD4, 0.), kBendAngleD4, kAlphaD4};

  m_magnet[3].type = kQuad;
  m_magnet[3].name = "K18Q12";
  m_magnet[3].scale = MagnetScale("K18Q12");
  m_magnet[3].theta = kBendAngleD4;
  m_magnet[3].quad = {kQ12B0, kQ12a0, kQ12z,
    d4_tan_half() + kDriftL3 + 0.5*kQ12z};

  m_magnet[4].type = kQuad;
  m_magnet[4].name = "K18Q13";
  m_magnet[4].scale = MagnetScale("K18Q13");
  m_magnet[4].theta = kBendAngleD4;
  m_magnet[4].quad = {kQ13B0, kQ13a0, kQ13z,
    d4_tan_half() + kDriftL3 + kQ12z + kDriftL4 + 0.5*kQ13z};

  m_sect_chord_x = kRhoD4*std::tan(0.5*kBendAngleD4*kDeg2Rad);
  m_sect_cot_alpha = 1./std::tan(kAlphaD4*kDeg2Rad);
  for(Int_t i=0; i<5; ++i){
    m_cos_theta[i] = std::cos(m_magnet[i].theta*kDeg2Rad);
    m_sin_theta[i] = std::sin(m_magnet[i].theta*kDeg2Rad);
  }

  debug::ObjectCounter::increase(ClassName());
}

//_____________________________________________________________________________
K18BeamlineRK::~K18BeamlineRK()
{
  debug::ObjectCounter::decrease(ClassName());
}

//_____________________________________________________________________________
Double_t
K18BeamlineRK::ConfDoubleOr(const TString& key, Double_t fallback) const
{
  return gConf.Get<TString>(key).IsNull() ? fallback : gConf.Get<Double_t>(key);
}

//_____________________________________________________________________________
Double_t
K18BeamlineRK::MagnetScale(const TString& name) const
{
  return m_global_scale*ConfDoubleOr(name + "Scale", 1.0);
}

//_____________________________________________________________________________
Double_t
K18BeamlineRK::BendAngle()
{
  return kBendAngleD4;
}

//_____________________________________________________________________________
TVector3
K18BeamlineRK::Axis(Double_t angle_deg)
{
  const Double_t a = angle_deg*kDeg2Rad;
  return TVector3(std::cos(a), std::sin(a), 0.);
}

//_____________________________________________________________________________
TVector3
K18BeamlineRK::XAxis(Double_t angle_deg)
{
  const Double_t a = angle_deg*kDeg2Rad;
  return TVector3(-std::sin(a), std::cos(a), 0.);
}

//_____________________________________________________________________________
TVector3
K18BeamlineRK::YAxis()
{
  return TVector3(0., 0., 1.);
}

//_____________________________________________________________________________
Double_t
K18BeamlineRK::VinL()
{
  return -(d4_tan_half() + kDriftL2 + kQ11z + kDriftL1 + kQ10z + kDriftL0);
}

//_____________________________________________________________________________
TVector3
K18BeamlineRK::Vin()
{
  return TVector3(VinL(), 0., 0.);
}

//_____________________________________________________________________________
Double_t
K18BeamlineRK::VoutL()
{
  return d4_tan_half() + kDriftL3 + kQ12z + kDriftL4 + kQ13z + kDriftL5;
}

//_____________________________________________________________________________
TVector3
K18BeamlineRK::Vout()
{
  return VoutL()*Axis(BendAngle());
}

//_____________________________________________________________________________
Double_t
K18BeamlineRK::TargetL()
{
  try {
    return gGeom.GetLocalZ("K18Target");
  } catch(...) {
    return gConf.Get<Double_t>("K18TargetL") != 0. ?
      gConf.Get<Double_t>("K18TargetL") : 1503.;
  }
}

//_____________________________________________________________________________
Double_t
K18BeamlineRK::BftL()
{
  try {
    return 0.5*(gGeom.GetLocalZ("BFT-X") + gGeom.GetLocalZ("BFT-XP"));
  } catch(...) {
    return gConf.Get<Double_t>("K18BFTL") != 0. ?
      gConf.Get<Double_t>("K18BFTL") : -33.5;
  }
}

//_____________________________________________________________________________
Double_t
K18BeamlineRK::UpstreamL(const TVector3& pos)
{
  return (pos - Vin()).Dot(Axis(0.));
}

//_____________________________________________________________________________
Double_t
K18BeamlineRK::UpstreamX(const TVector3& pos)
{
  return (pos - Vin()).Dot(XAxis(0.));
}

//_____________________________________________________________________________
Bool_t
K18BeamlineRK::InMagnet(const TVector3& pos, Int_t i) const
{
  const auto& mag = m_magnet[i];
  if(mag.type == kQuad){
    const Double_t ct = m_cos_theta[i];
    const Double_t st = m_sin_theta[i];
    const Double_t gx =  ct*pos.X() + st*pos.Y();
    const Double_t gy = -st*pos.X() + ct*pos.Y();
    const Double_t lx = gy;
    const Double_t ly = pos.Z();
    const Double_t lz = gx;
    return (std::abs(lz - mag.quad.z) <= 0.5*mag.quad.l
            && lx*lx + ly*ly <= kQuadFieldRadiusFactor*kQuadFieldRadiusFactor
                                *mag.quad.a0*mag.quad.a0);
  }

  const auto& sec = mag.sect;
  if(std::abs(pos.Z()) > sec.half_gap) return false;

  const Double_t lvecx = pos.X() - sec.center.X();
  const Double_t lvecy = pos.Y() - sec.center.Y();
  const Double_t r = std::sqrt(lvecx*lvecx + lvecy*lvecy);
  if(r <= 0.) return false;
  const Double_t ctheta = std::max(-1., std::min(1., lvecx/r));
  Double_t theta = std::acos(ctheta)/kDeg2Rad;
  if(lvecy < 0.) theta = 360. - theta;

  const Double_t entrance_line =
    -m_sect_cot_alpha*(pos.X() + m_sect_chord_x);

  if(r >= sec.rho - 0.5*sec.width
     && r <= sec.rho + 0.5*sec.width
     && theta >= 270.
     && theta <= 270. + sec.bend_angle){
    return pos.Y() >= entrance_line;
  }

  if(pos.Y() >= entrance_line
     && pos.X() <= -m_sect_chord_x
     && pos.Y() <= 0.5*sec.width)
    return true;

  return false;
}

//_____________________________________________________________________________
TVector3
K18BeamlineRK::CalcQuad(Int_t i, const TVector3& pos) const
{
  const auto& mag = m_magnet[i];
  const Double_t ct = m_cos_theta[i];
  const Double_t st = m_sin_theta[i];
  const Double_t lx = -st*pos.X() + ct*pos.Y();
  const Double_t ly = pos.Z();
  const Double_t g = mag.scale*mag.quad.b0/mag.quad.a0;
  const Double_t bgy = g*ly;
  const Double_t bgz = g*lx;
  return TVector3(st*bgy, ct*bgy, bgz);
}

//_____________________________________________________________________________
TVector3
K18BeamlineRK::CalcSect(Int_t i) const
{
  return TVector3(0., 0., m_magnet[i].scale*m_magnet[i].sect.b0);
}

//_____________________________________________________________________________
TVector3
K18BeamlineRK::Field(const TVector3& pos) const
{
  if(m_field_map){
    TVector3 field_tesla;
    m_field_map->GetFieldValue(pos, field_tesla);
    return field_tesla;
  }

  TVector3 b_kgauss;
  for(Int_t i=0; i<5; ++i){
    if(!InMagnet(pos, i)) continue;
    b_kgauss = (m_magnet[i].type == kQuad) ? CalcQuad(i, pos) : CalcSect(i);
    break;
  }
  return 0.1*b_kgauss; // kG -> T
}

//_____________________________________________________________________________
K18RKState
K18BeamlineRK::Derivative(const K18RKState& state) const
{
  K18RKState d;
  const Double_t p_mag = state.mom.Mag();
  if(p_mag <= 0.) return d;

  const TVector3 dir = state.mom.Unit();
  const TVector3 b_tesla = Field(state.pos);
  d.pos = dir;
  d.mom = m_charge*0.000299792458*dir.Cross(b_tesla);
  return d;
}

//_____________________________________________________________________________
void
K18BeamlineRK::StepRK4(K18RKState& state, Double_t ds) const
{
  const Double_t p_mag = state.mom.Mag();
  const auto k1 = Derivative(state);
  const K18RKState s2{state.pos + 0.5*ds*k1.pos,
                      state.mom + 0.5*ds*k1.mom};
  const auto k2 = Derivative(s2);
  const K18RKState s3{state.pos + 0.5*ds*k2.pos,
                      state.mom + 0.5*ds*k2.mom};
  const auto k3 = Derivative(s3);
  const K18RKState s4{state.pos + ds*k3.pos,
                      state.mom + ds*k3.mom};
  const auto k4 = Derivative(s4);

  state.pos += (ds/6.)*(k1.pos + 2.*k2.pos + 2.*k3.pos + k4.pos);
  state.mom += (ds/6.)*(k1.mom + 2.*k2.mom + 2.*k3.mom + k4.mom);
  if(state.mom.Mag2() > 0.)
    state.mom = state.mom.Unit()*p_mag;
}

//_____________________________________________________________________________
Bool_t
K18BeamlineRK::PropagateToBft(Double_t xvo, Double_t yvo,
                              Double_t uvo, Double_t vvo,
                              Double_t p,
                              K18RKState& state,
                              Double_t* path_length,
                              Int_t* n_steps) const
{
  const auto axis = Axis(BendAngle());
  const auto ex = XAxis(BendAngle());
  const auto ey = YAxis();
  state.pos = Vout() + xvo*ex + yvo*ey;
  state.mom = p*(axis + uvo*ex + vvo*ey).Unit();

  const Double_t target_bft_l = BftL();
  K18RKState prev = state;
  Double_t prev_delta = UpstreamL(prev.pos) - target_bft_l;
  if(prev_delta <= 0.) return true;

  const Double_t step_field = std::abs(m_rk_step);
  const Double_t step_drift = std::abs(m_rk_drift_step);
  const Int_t max_steps = static_cast<Int_t>(std::ceil(m_max_path/std::max(step_field, 1.e-6)));
  Double_t traveled = 0.;

  for(Int_t istep=0; istep<max_steps && traveled<m_max_path; ++istep){
    Double_t cur_step = (Field(prev.pos).Mag2() > 0.)
      ? step_field : step_drift;
    K18RKState cur = prev;
    StepRK4(cur, -cur_step);
    if(cur_step > step_field && Field(cur.pos).Mag2() > 0.){
      cur = prev;
      cur_step = step_field;
      StepRK4(cur, -cur_step);
    }
    traveled += cur_step;
    const Double_t cur_delta = UpstreamL(cur.pos) - target_bft_l;
    if(cur_delta <= 0.){
      const Double_t denom = prev_delta - cur_delta;
      const Double_t f = (denom != 0.) ? prev_delta/denom : 0.;
      state.pos = prev.pos + f*(cur.pos - prev.pos);
      state.mom = prev.mom + f*(cur.mom - prev.mom);
      if(state.mom.Mag2() > 0.) state.mom = state.mom.Unit()*p;
      if(path_length) *path_length = traveled - cur_step + f*cur_step;
      if(n_steps) *n_steps = istep + 1;
      return true;
    }
    prev = cur;
    prev_delta = cur_delta;
  }

  return false;
}

//_____________________________________________________________________________
Bool_t
K18BeamlineRK::FitMomentum(Double_t xvo, Double_t yvo,
                           Double_t uvo, Double_t vvo,
                           Double_t xbft,
                           K18RKTrack& track) const
{
  struct FitPoint {
    Bool_t ok = false;
    Double_t p = 0.;
    Double_t residual = std::numeric_limits<Double_t>::quiet_NaN();
    Double_t xcalc = std::numeric_limits<Double_t>::quiet_NaN();
    Double_t path = 0.;
    Int_t steps = 0;
    K18RKState state;

    Double_t Chi2() const { return residual*residual; }
  };

  const Double_t p_lo = std::min(m_p_min, m_p_max);
  const Double_t p_hi = std::max(m_p_min, m_p_max);
  if(!(p_lo > 0.) || !(p_hi > p_lo))
    return false;

  const Double_t p_seed_conf = std::abs(ConfDoubleOr("PK18", 0.5*(p_lo + p_hi)));
  const Double_t p_seed = std::clamp(p_seed_conf, p_lo, p_hi);
  const Int_t max_iteration = m_fit_max_iteration;
  const Double_t p_tolerance = m_fit_momentum_tolerance;
  const Double_t residual_tolerance = m_fit_residual_tolerance;

  auto evaluate = [&](Double_t p) -> FitPoint {
    FitPoint point;
    point.p = std::clamp(p, p_lo, p_hi);
    if(!std::isfinite(point.p))
      return point;
    if(!PropagateToBft(xvo, yvo, uvo, vvo, point.p,
                       point.state, &point.path, &point.steps))
      return point;
    point.xcalc = UpstreamX(point.state.pos);
    point.residual = point.xcalc - xbft;
    point.ok = std::isfinite(point.residual);
    return point;
  };

  auto better = [](const FitPoint& lhs, const FitPoint& rhs) -> FitPoint {
    if(!lhs.ok) return rhs;
    if(!rhs.ok) return lhs;
    return rhs.Chi2() < lhs.Chi2() ? rhs : lhs;
  };

  FitPoint best = evaluate(p_seed);
  FitPoint low = evaluate(p_lo);
  FitPoint high = evaluate(p_hi);
  best = better(best, low);
  best = better(best, high);
  if(!best.ok)
    return false;

  const auto has_crossing = [](const FitPoint& a, const FitPoint& b) -> Bool_t {
    return a.ok && b.ok &&
      (a.residual == 0. || b.residual == 0. || a.residual*b.residual < 0.);
  };

  Bool_t bracketed = false;
  FitPoint left;
  FitPoint right;
  if(has_crossing(low, best)){
    left = low;
    right = best;
    bracketed = true;
  }
  if(has_crossing(best, high)){
    const Double_t width = std::abs(high.p - best.p);
    if(!bracketed || width < std::abs(right.p - left.p)){
      left = best;
      right = high;
      bracketed = true;
    }
  }
  if(!bracketed && has_crossing(low, high)){
    left = low;
    right = high;
    bracketed = true;
  }

  Int_t fit_status = 0;

  if(bracketed){
    fit_status = 1;
    if(left.p > right.p)
      std::swap(left, right);
    for(Int_t iter=0; iter<max_iteration; ++iter){
      const Double_t width = right.p - left.p;
      if(width <= p_tolerance || std::abs(best.residual) <= residual_tolerance)
        break;

      Double_t p_next = 0.5*(left.p + right.p);
      const Double_t denom = right.residual - left.residual;
      if(std::abs(denom) > std::numeric_limits<Double_t>::epsilon()){
        const Double_t p_secant =
          left.p - left.residual*(right.p - left.p)/denom;
        if(std::isfinite(p_secant) && left.p < p_secant && p_secant < right.p)
          p_next = p_secant;
      }

      FitPoint next = evaluate(p_next);
      if(!next.ok)
        break;
      best = better(best, next);
      if(std::abs(next.residual) <= residual_tolerance)
        break;

      if(left.residual == 0.){
        right = left;
        best = better(best, left);
        break;
      }
      if(next.residual == 0.){
        left = right = next;
        best = better(best, next);
        break;
      }
      if(left.residual*next.residual < 0.)
        right = next;
      else
        left = next;
    }
  }
  else {
    fit_status = 2;
    const Double_t inv_phi = 0.6180339887498948482;
    Double_t a = p_lo;
    Double_t b = p_hi;
    Double_t c = b - inv_phi*(b - a);
    Double_t d = a + inv_phi*(b - a);
    FitPoint fc = evaluate(c);
    FitPoint fd = evaluate(d);
    best = better(best, fc);
    best = better(best, fd);

    for(Int_t iter=0; iter<max_iteration; ++iter){
      if(b - a <= p_tolerance || std::abs(best.residual) <= residual_tolerance)
        break;
      if(!fc.ok || !fd.ok)
        break;
      if(fc.Chi2() < fd.Chi2()){
        b = d;
        d = c;
        fd = fc;
        c = b - inv_phi*(b - a);
        fc = evaluate(c);
        best = better(best, fc);
      }
      else {
        a = c;
        c = d;
        fc = fd;
        d = a + inv_phi*(b - a);
        fd = evaluate(d);
        best = better(best, fd);
      }
    }
  }

  if(!best.ok || std::abs(best.residual) > m_bft_residual_max)
    return false;

  if(std::abs(best.p - p_lo) < p_tolerance || std::abs(best.p - p_hi) < p_tolerance)
    fit_status = 3;

  K18RKTrack fitted;
  fitted.status = true;
  fitted.p = best.p;
  fitted.xvo = xvo;
  fitted.yvo = yvo;
  fitted.uvo = uvo;
  fitted.vvo = vvo;
  const Double_t l = TargetL();
  fitted.xtgt = xvo + uvo*l;
  fitted.ytgt = yvo + vvo*l;
  fitted.utgt = uvo;
  fitted.vtgt = vvo;
  fitted.xbft = xbft;
  fitted.xbft_calc = best.xcalc;
  fitted.bft_residual = best.residual;
  fitted.chisqr = best.residual*best.residual;
  fitted.ndf = 1;
  fitted.n_iteration = 0;
  fitted.path_length = best.path;
  fitted.n_steps = best.steps;
  fitted.fit_status = fit_status;
  fitted.bft_state = best.state;
  track = fitted;
  return true;
}

//_____________________________________________________________________________
Bool_t
K18BeamlineRK::FitTrack(const DCLocalTrack& local, Double_t xbft,
                        K18RKTrack& track) const
{
  K18RKTrack seed;
  if(!FitMomentum(local.GetX(0.), local.GetY(0.),
                  local.GetU0(), local.GetV0(), xbft, seed))
    return false;

  const auto& hits = local.GetHitArray();
  if(hits.size() < 5)
    return false;

  struct FitEval {
    Bool_t ok = false;
    std::array<Double_t, 5> par{};
    std::vector<Double_t> res;
    Double_t chisqr = std::numeric_limits<Double_t>::infinity();
    Double_t raw_chisqr = std::numeric_limits<Double_t>::infinity();
    Int_t ndf = 0;
    Double_t xcalc = 0.;
    Double_t bft_residual = 0.;
    Double_t path = 0.;
    Int_t steps = 0;
    K18RKState state;
  };

  const Double_t p_lo = std::min(m_p_min, m_p_max);
  const Double_t p_hi = std::max(m_p_min, m_p_max);
  const Double_t bft_sigma =
    std::max(1.e-6, ConfDoubleOr("K18NativeBFTFitResolution", 0.2));
  const Int_t max_iteration =
    std::max(1, static_cast<Int_t>(ConfDoubleOr("K18NativeFullFitMaxIteration", 25.)));
  const Double_t lambda0 =
    std::max(1.e-12, ConfDoubleOr("K18NativeFullFitLambda", 1.e-3));

  auto evaluate = [&](const std::array<Double_t, 5>& par) -> FitEval {
    FitEval out;
    out.par = par;
    if(!std::isfinite(par[0]) || !std::isfinite(par[1]) ||
       !std::isfinite(par[2]) || !std::isfinite(par[3]) ||
       !std::isfinite(par[4]) || par[4] < p_lo || par[4] > p_hi)
      return out;

    out.res.reserve(hits.size() + 1);
    Double_t raw = 0.;
    for(const auto* hit: hits){
      if(!hit) continue;
      const Double_t a = hit->GetTiltAngle()*kDeg2Rad;
      const Double_t ct = std::cos(a);
      const Double_t st = std::sin(a);
      const Double_t z = hit->GetZ();
      const Double_t x = par[0] + par[2]*z;
      const Double_t y = par[1] + par[3]*z;
      const Double_t dsdz = par[2]*ct + par[3]*st;
      const Double_t coss = hit->IsHoneycomb() ? std::cos(std::atan(dsdz)) : 1.;
      const Double_t wp = hit->GetWirePosition();
      const Double_t ss = wp + (hit->GetLocalHitPos() - wp)/coss;
      const Double_t cal = x*ct + y*st;
      const Double_t residual = (ss - cal)*coss;
      const Double_t sigma = std::max(1.e-6, hit->GetResolution());
      out.res.push_back(residual/sigma);
      raw += (residual/sigma)*(residual/sigma);
    }

    if(!PropagateToBft(par[0], par[1], par[2], par[3], par[4],
                       out.state, &out.path, &out.steps))
      return out;
    out.xcalc = UpstreamX(out.state.pos);
    out.bft_residual = out.xcalc - xbft;
    out.res.push_back(out.bft_residual/bft_sigma);
    raw += (out.bft_residual/bft_sigma)*(out.bft_residual/bft_sigma);
    out.raw_chisqr = raw;
    out.ndf = static_cast<Int_t>(out.res.size()) - 5;
    out.chisqr = raw/std::max(1, out.ndf);
    out.ok = out.ndf > 0 && std::isfinite(out.chisqr);
    return out;
  };

  auto solve5 = [](std::array<std::array<Double_t, 5>, 5> a,
                   std::array<Double_t, 5> b,
                   std::array<Double_t, 5>& x) -> Bool_t {
    for(Int_t i=0; i<5; ++i){
      Int_t pivot = i;
      Double_t best = std::abs(a[i][i]);
      for(Int_t r=i+1; r<5; ++r){
        const Double_t v = std::abs(a[r][i]);
        if(v > best){ best = v; pivot = r; }
      }
      if(best < 1.e-14)
        return false;
      if(pivot != i){
        std::swap(a[pivot], a[i]);
        std::swap(b[pivot], b[i]);
      }
      const Double_t diag = a[i][i];
      for(Int_t c=i; c<5; ++c) a[i][c] /= diag;
      b[i] /= diag;
      for(Int_t r=0; r<5; ++r){
        if(r == i) continue;
        const Double_t f = a[r][i];
        if(f == 0.) continue;
        for(Int_t c=i; c<5; ++c) a[r][c] -= f*a[i][c];
        b[r] -= f*b[i];
      }
    }
    x = b;
    return true;
  };

  std::array<Double_t, 5> par = {
    seed.xvo, seed.yvo, seed.uvo, seed.vvo, seed.p
  };
  FitEval best = evaluate(par);
  if(!best.ok)
    return false;

  Double_t lambda = lambda0;
  Int_t accepted = 0;
  const std::array<Double_t, 5> step = {0.01, 0.01, 1.e-5, 1.e-5, 1.e-5};
  const std::array<Double_t, 5> max_delta = {2.0, 2.0, 0.002, 0.002, 0.01};

  for(Int_t iter=0; iter<max_iteration; ++iter){
    std::vector<std::array<Double_t, 5>> jac(best.res.size());
    Bool_t deriv_ok = true;
    for(Int_t ip=0; ip<5; ++ip){
      const Double_t h = step[ip]*std::max(1., std::abs(par[ip]));
      auto plus = par;
      auto minus = par;
      plus[ip] += h;
      minus[ip] -= h;
      plus[4] = std::clamp(plus[4], p_lo, p_hi);
      minus[4] = std::clamp(minus[4], p_lo, p_hi);
      const FitEval ep = evaluate(plus);
      const FitEval em = evaluate(minus);
      if(!ep.ok || !em.ok || ep.res.size() != best.res.size() ||
         em.res.size() != best.res.size() || plus[ip] == minus[ip]){
        deriv_ok = false;
        break;
      }
      const Double_t inv = 1./(plus[ip] - minus[ip]);
      for(std::size_t ir=0; ir<best.res.size(); ++ir)
        jac[ir][ip] = (ep.res[ir] - em.res[ir])*inv;
    }
    if(!deriv_ok)
      break;

    std::array<std::array<Double_t, 5>, 5> normal{};
    std::array<Double_t, 5> rhs{};
    for(std::size_t ir=0; ir<best.res.size(); ++ir){
      for(Int_t i=0; i<5; ++i){
        rhs[i] -= jac[ir][i]*best.res[ir];
        for(Int_t j=0; j<5; ++j)
          normal[i][j] += jac[ir][i]*jac[ir][j];
      }
    }
    for(Int_t i=0; i<5; ++i)
      normal[i][i] *= (1. + lambda);

    std::array<Double_t, 5> delta{};
    if(!solve5(normal, rhs, delta)){
      lambda *= 10.;
      continue;
    }

    Double_t scale = 1.;
    for(Int_t i=0; i<5; ++i){
      if(std::abs(delta[i]) > max_delta[i])
        scale = std::min(scale, max_delta[i]/std::abs(delta[i]));
    }
    auto trial_par = par;
    for(Int_t i=0; i<5; ++i)
      trial_par[i] += scale*delta[i];
    trial_par[4] = std::clamp(trial_par[4], p_lo, p_hi);

    FitEval trial = evaluate(trial_par);
    if(trial.ok && trial.raw_chisqr < best.raw_chisqr){
      const Double_t improvement = best.raw_chisqr - trial.raw_chisqr;
      par = trial_par;
      best = trial;
      ++accepted;
      lambda = std::max(1.e-12, 0.3*lambda);
      if(improvement < 1.e-6 || scale < 1.e-4)
        break;
    }
    else {
      lambda *= 10.;
    }
  }

  if(!best.ok || accepted < 1
     || std::abs(best.bft_residual) > m_bft_residual_max)
    return false;

  K18RKTrack fitted;
  fitted.status = true;
  fitted.p = best.par[4];
  fitted.xvo = best.par[0];
  fitted.yvo = best.par[1];
  fitted.uvo = best.par[2];
  fitted.vvo = best.par[3];
  const Double_t l = TargetL();
  fitted.xtgt = fitted.xvo + fitted.uvo*l;
  fitted.ytgt = fitted.yvo + fitted.vvo*l;
  fitted.utgt = fitted.uvo;
  fitted.vtgt = fitted.vvo;
  fitted.xbft = xbft;
  fitted.xbft_calc = best.xcalc;
  fitted.bft_residual = best.bft_residual;
  fitted.chisqr = best.chisqr;
  fitted.ndf = best.ndf;
  fitted.n_iteration = accepted;
  fitted.path_length = best.path;
  fitted.n_steps = best.steps;
  fitted.fit_status = 1;
  fitted.bft_state = best.state;
  track = fitted;
  return true;
}
