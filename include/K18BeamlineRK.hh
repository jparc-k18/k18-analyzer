// -*- C++ -*-

#ifndef K18_BEAMLINE_RK_HH
#define K18_BEAMLINE_RK_HH

#include <memory>

#include <TString.h>
#include <TVector3.h>

class DCLocalTrack;
class K18FieldMap;

struct K18RKState
{
  TVector3 pos;
  TVector3 mom; // GeV/c in K18 native internal coordinates
};

struct K18RKTrack
{
  Bool_t status = false;
  Double_t p = 0.;
  Double_t xvo = 0.;
  Double_t yvo = 0.;
  Double_t uvo = 0.;
  Double_t vvo = 0.;
  Double_t xtgt = 0.;
  Double_t ytgt = 0.;
  Double_t utgt = 0.;
  Double_t vtgt = 0.;
  Double_t xbft = 0.;
  Double_t xbft_calc = 0.;
  Double_t bft_residual = 0.;
  Double_t chisqr = 0.;
  Int_t ndf = 0;
  Int_t n_iteration = 0;
  Double_t path_length = 0.;
  Int_t n_steps = 0;
  Int_t fit_status = 0; // 0=none, 1=bracketed, 2=minimized, 3=boundary
  K18RKState bft_state;
};

// Fits K1.8 momentum by propagating a VO track backward through QQDQQ until
// its calculated BFT coordinate matches the measured cluster position.
class K18BeamlineRK
{
public:
  static const TString& ClassName();
  K18BeamlineRK();
  ~K18BeamlineRK();

  Bool_t FitMomentum(Double_t xvo, Double_t yvo,
                     Double_t uvo, Double_t vvo,
                     Double_t xbft,
                     K18RKTrack& track) const;
  Bool_t FitTrack(const DCLocalTrack& local, Double_t xbft,
                  K18RKTrack& track) const;
  Bool_t PropagateToBft(Double_t xvo, Double_t yvo,
                        Double_t uvo, Double_t vvo,
                        Double_t p,
                        K18RKState& state,
                        Double_t* path_length=nullptr,
                        Int_t* n_steps=nullptr) const;

  Double_t PMin() const { return m_p_min; }
  Double_t PMax() const { return m_p_max; }
  Double_t RKStep() const { return m_rk_step; }
  Double_t Charge() const { return m_charge; }
  Double_t GlobalScale() const { return m_global_scale; }
  Bool_t UsesFieldMap() const { return static_cast<Bool_t>(m_field_map); }

private:
  struct Quad {
    Double_t b0;
    Double_t a0;
    Double_t l;
    Double_t z;
  };
  struct Sect {
    Double_t b0;
    Double_t rho;
    Double_t width;
    Double_t half_gap;
    TVector3 center;
    Double_t bend_angle;
    Double_t alpha;
  };
  enum MagnetType { kQuad, kSect };
  struct Magnet {
    MagnetType type;
    TString name;
    Double_t scale;
    Double_t theta;
    Quad quad;
    Sect sect;
  };

  std::unique_ptr<K18FieldMap> m_field_map;
  Double_t m_p_min;
  Double_t m_p_max;
  Double_t m_rk_step;
  Int_t m_fit_max_iteration;
  Double_t m_fit_momentum_tolerance;
  Double_t m_fit_residual_tolerance;
  Double_t m_max_path;
  Double_t m_bft_residual_max;
  Double_t m_charge;
  Double_t m_global_scale;
  Double_t m_rk_drift_step;
  Double_t m_sect_cot_alpha;
  Double_t m_sect_chord_x;
  Double_t m_cos_theta[5];
  Double_t m_sin_theta[5];
  Magnet m_magnet[5];

  static constexpr Double_t kQuadFieldRadiusFactor = 2.;

  static Double_t BendAngle();
  static TVector3 Axis(Double_t angle_deg);
  static TVector3 XAxis(Double_t angle_deg);
  static TVector3 YAxis();
  static Double_t VinL();
  static TVector3 Vin();
  static Double_t VoutL();
  static TVector3 Vout();
  static Double_t TargetL();
  static Double_t BftL();
  static Double_t UpstreamL(const TVector3& pos);
  static Double_t UpstreamX(const TVector3& pos);

  Double_t ConfDoubleOr(const TString& key, Double_t fallback) const;
  Double_t MagnetScale(const TString& name) const;
  Bool_t InMagnet(const TVector3& pos, Int_t i) const;
  TVector3 Field(const TVector3& pos) const;
  TVector3 CalcQuad(Int_t i, const TVector3& pos) const;
  TVector3 CalcSect(Int_t i) const;
  K18RKState Derivative(const K18RKState& state) const;
  void StepRK4(K18RKState& state, Double_t ds) const;
};

#endif
