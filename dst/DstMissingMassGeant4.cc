// DstMissingMassGeant4.cc
// Evaluate reaction missing-mass resolution with the data-like vertex + BB flow.
//
// The primary quantity uses K18 D2U p_3rd and S-2S pS2s, extrapolates both
// target-plane states to one Kinematics::VertexPoint, applies mean Bethe--Bloch
// corrections to that same reconstructed vertex, and then builds missing mass.
// The one-arm rows use the identical reconstructed vertex and correction path;
// only the opposite reconstructed four-vector is replaced by G4 reaction truth.

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <initializer_list>
#include <iostream>
#include <memory>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <system_error>
#include <utility>
#include <vector>

#include "Kinematics.hh"

#include <TFile.h>
#include <TF1.h>
#include <TFitResult.h>
#include <TFitResultPtr.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLorentzVector.h>
#include <TNamed.h>
#include <TParticle.h>
#include <TROOT.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TTreeReaderValue.h>
#include <TVector3.h>

namespace {
constexpr double kFwhm = 2.*std::sqrt(2.*std::log(2.));
constexpr double kElectronMassMeV = 0.51099895;


void
WriteDstMetadata(TDirectory* dir, const char* inputContract,
                 const char* selection)
{
  if(!dir) return;
  dir->cd();
  const char* revision = std::getenv("K18_ANALYZER_REVISION");
  if(!revision || !*revision) revision = "unknown";
  TNamed("git_hash", revision).Write("", TObject::kOverwrite);
  TNamed("conf", inputContract).Write("", TObject::kOverwrite);
  TNamed("user", selection).Write("", TObject::kOverwrite);
  TNamed("user_digest", "missing").Write("", TObject::kOverwrite);
  TNamed("dst_git_hash", revision).Write("", TObject::kOverwrite);
  TNamed("dst_input_contract", inputContract).Write("", TObject::kOverwrite);
  TNamed("dst_selection", selection).Write("", TObject::kOverwrite);
  TNamed("dst_program", "DstMissingMassGeant4").Write(
    "", TObject::kOverwrite);
}

struct Options {
  std::string g4, k18, s2s, root, csv, label, targetShape = "slab";
  int beamPdg = 211;
  int scatPdg = 321;
  double beamMass = 0.13957039;
  double scatMass = 0.493677;
  double targetMass = 11.17486323534;
  double targetThickness = 16.666666667;
  double targetDensity = 1.8;
  double targetIEv = 81.0;
  double targetZOverA = 6.0/12.011;
  double targetSizeX = 100.;
  double targetSizeY = 100.;
  double targetCenterX = 0.;
  double targetCenterY = 0.;
  double targetPlaneZ = -5023.5;
  double targetZMargin = 50.;
  double targetRadius = 0.;
  double targetHalfLengthY = 0.;
  double cellRadius = 0.;
  double cellThickness = 0.;
  double cellDensity = 0.;
  double cellIEv = 0.;
  double cellZOverA = 0.;
  double vacuumShellRadius = 0.;
  double vacuumShellThickness = 0.;
  double vacuumShellDensity = 0.;
  double vacuumShellIEv = 0.;
  double vacuumShellZOverA = 0.;
  double closeDistMax = 100.;
  double bh2Thickness = 5.;
  double chi2Max = 20.;
  double thetaMin = 2.;
  double thetaMax = 14.;
};

std::filesystem::path
ComparablePath(const std::string& value)
{
  std::error_code error;
  auto path = std::filesystem::absolute(value, error);
  if(error) path = std::filesystem::path(value);
  auto resolved = std::filesystem::weakly_canonical(path, error);
  return error ? path.lexically_normal() : resolved;
}

void
ValidateDistinctPaths(const Options& o)
{
  const std::vector<std::pair<std::string,std::string>> stages{
    {"g4",o.g4}, {"k18",o.k18}, {"s2s",o.s2s},
    {"root",o.root}, {"csv",o.csv}
  };
  for(std::size_t i=0; i<stages.size(); ++i){
    const auto left=ComparablePath(stages[i].second);
    for(std::size_t j=i+1; j<stages.size(); ++j){
      const auto right=ComparablePath(stages[j].second);
      std::error_code error;
      const bool samePath=left==right;
      const bool sameFile=std::filesystem::exists(left,error) && !error &&
        std::filesystem::exists(right,error) && !error &&
        std::filesystem::equivalent(left,right,error) && !error;
      if(samePath || sameFile)
        throw std::runtime_error(
          "stage path collision: "+stages[i].first+"="+left.string()+
          " and "+stages[j].first+"="+right.string());
    }
  }
}

void
RequireBranches(TTree* tree, const char* treeName,
                std::initializer_list<const char*> names)
{
  std::vector<std::string> missing;
  for(const char* name:names)
    if(!tree->GetBranch(name)) missing.emplace_back(name);
  if(missing.empty()) return;
  std::ostringstream message;
  message << treeName << " missing required branch";
  if(missing.size()!=1) message << "es";
  message << ": ";
  for(std::size_t i=0; i<missing.size(); ++i){
    if(i) message << ", ";
    message << missing[i];
  }
  throw std::runtime_error(message.str());
}

Options Parse(int argc, char** argv)
{
  Options o;
  for(int i=1; i<argc; ++i){
    const std::string arg = argv[i];
    auto takeString = [&](const char* key, std::string& value){
      const std::string prefix = std::string("--") + key + "=";
      if(arg.rfind(prefix, 0) != 0) return false;
      value = arg.substr(prefix.size());
      return true;
    };
    auto takeDouble = [&](const char* key, double& value){
      const std::string prefix = std::string("--") + key + "=";
      if(arg.rfind(prefix, 0) != 0) return false;
      value = std::stod(arg.substr(prefix.size()));
      return true;
    };
    auto takeInt = [&](const char* key, int& value){
      const std::string prefix = std::string("--") + key + "=";
      if(arg.rfind(prefix, 0) != 0) return false;
      value = std::stoi(arg.substr(prefix.size()));
      return true;
    };
    if(takeString("g4", o.g4) || takeString("k18", o.k18) ||
       takeString("s2s", o.s2s) || takeString("root", o.root) ||
       takeString("csv", o.csv) || takeString("label", o.label) ||
       takeString("target-shape", o.targetShape) ||
       takeInt("beam-pdg", o.beamPdg) || takeInt("scat-pdg", o.scatPdg) ||
       takeDouble("beam-mass", o.beamMass) ||
       takeDouble("scat-mass", o.scatMass) ||
       takeDouble("target-mass", o.targetMass) ||
       takeDouble("target-thickness-mm", o.targetThickness) ||
       takeDouble("target-density", o.targetDensity) ||
       takeDouble("target-i-ev", o.targetIEv) ||
       takeDouble("target-z-over-a", o.targetZOverA) ||
       takeDouble("target-size-x-mm", o.targetSizeX) ||
       takeDouble("target-size-y-mm", o.targetSizeY) ||
       takeDouble("target-center-x-mm", o.targetCenterX) ||
       takeDouble("target-center-y-mm", o.targetCenterY) ||
       takeDouble("target-plane-z-mm", o.targetPlaneZ) ||
       takeDouble("target-z-margin-mm", o.targetZMargin) ||
       takeDouble("target-radius-mm", o.targetRadius) ||
       takeDouble("target-half-length-y-mm", o.targetHalfLengthY) ||
       takeDouble("cell-radius-mm", o.cellRadius) ||
       takeDouble("cell-thickness-mm", o.cellThickness) ||
       takeDouble("cell-density", o.cellDensity) ||
       takeDouble("cell-i-ev", o.cellIEv) ||
       takeDouble("cell-z-over-a", o.cellZOverA) ||
       takeDouble("vacuum-shell-radius-mm", o.vacuumShellRadius) ||
       takeDouble("vacuum-shell-thickness-mm", o.vacuumShellThickness) ||
       takeDouble("vacuum-shell-density", o.vacuumShellDensity) ||
       takeDouble("vacuum-shell-i-ev", o.vacuumShellIEv) ||
       takeDouble("vacuum-shell-z-over-a", o.vacuumShellZOverA) ||
       takeDouble("close-dist-max-mm", o.closeDistMax) ||
       takeDouble("bh2-thickness-mm", o.bh2Thickness) ||
       takeDouble("chi2-max", o.chi2Max) ||
       takeDouble("theta-min-deg", o.thetaMin) ||
       takeDouble("theta-max-deg", o.thetaMax)) continue;
    throw std::runtime_error("unknown option: " + arg);
  }
  if(o.g4.empty() || o.k18.empty() || o.s2s.empty() || o.root.empty() ||
     o.csv.empty() || o.label.empty())
    throw std::runtime_error("required: --g4 --k18 --s2s --root --csv --label");
  if(o.targetShape != "slab" && o.targetShape != "cylinder")
    throw std::runtime_error("--target-shape must be slab or cylinder");
  if(o.targetShape == "cylinder" &&
     (!(o.targetRadius > 0.) || !(o.targetHalfLengthY > 0.)))
    throw std::runtime_error(
      "cylinder target requires --target-radius-mm and "
      "--target-half-length-y-mm");
  ValidateDistinctPaths(o);
  for(const auto* path : {&o.root, &o.csv}){
    const auto parent = std::filesystem::path(*path).parent_path();
    if(!parent.empty()) std::filesystem::create_directories(parent);
  }
  return o;
}

bool FiniteVector(const TVector3& v)
{
  return std::isfinite(v.X()) && std::isfinite(v.Y()) && std::isfinite(v.Z());
}

TLorentzVector FourVector(double p, double u, double v, double mass)
{
  const double norm = std::sqrt(1. + u*u + v*v);
  return {p*u/norm, p*v/norm, p/norm, std::hypot(p, mass)};
}

TLorentzVector ParticleVector(const TParticle& p)
{
  return {p.Px()*1e-3, p.Py()*1e-3, p.Pz()*1e-3, p.Energy()*1e-3};
}

double MissingMass(const TLorentzVector& beam, const TLorentzVector& scat,
                   double targetMass)
{
  const TLorentzVector target(0., 0., 0., targetMass);
  const double m2 = (beam + target - scat).M2();
  return m2 > 0. ? std::sqrt(m2) : NAN;
}

double Residual(const TLorentzVector& beam, const TLorentzVector& scat,
                double targetMass, double truthMass)
{
  const double mass = MissingMass(beam, scat, targetMass);
  return std::isfinite(mass) ? 1e3*(mass-truthMass) : NAN;
}

const TParticle* SinglePrimary(const TTreeReaderArray<TParticle>& particles)
{
  return particles.GetSize()==1 ? &particles[0] : nullptr;
}

bool HasPrimary(const TTreeReaderArray<TParticle>& particles)
{
  for(const auto& p : particles)
    if(p.GetFirstMother()==0) return true;
  return false;
}

double BetheBloch(double pGeV, double massGeV, double rho, double iEv,
                  double zOverA)
{
  if(!(pGeV>0.) || !(massGeV>0.) || !(rho>0.) || !(iEv>0.) || !(zOverA>0.))
    return NAN;
  const double energy = std::hypot(pGeV, massGeV);
  const double beta2 = pGeV*pGeV/(energy*energy);
  const double gamma2 = 1./(1.-beta2);
  const double massRatio = kElectronMassMeV/(1e3*massGeV);
  const double tmax = 2.*kElectronMassMeV*beta2*gamma2/
    (1. + 2.*std::sqrt(gamma2)*massRatio + massRatio*massRatio);
  const double iMeV = iEv*1e-6;
  const double argument = 2.*kElectronMassMeV*beta2*gamma2*tmax/(iMeV*iMeV);
  if(!(argument>1.)) return NAN;
  const double bracket = .5*std::log(argument)-beta2;
  return 0.307075*rho*zOverA*bracket/beta2; // MeV/cm
}

// Forward propagation loses energy; inverse propagation restores it. Momenta
// are GeV/c and path length is mm.
double PropagateBB(double pGeV, double massGeV, double pathMm, double rho,
                   double iEv, double zOverA, bool inverse)
{
  if(!(pathMm>0.) || !(rho>0.)) return pGeV;
  double pMeV = 1e3*pGeV;
  const double massMeV = 1e3*massGeV;
  double energy = std::hypot(pMeV, massMeV);
  const int nstep = std::max(1, static_cast<int>(std::ceil(pathMm/.2)));
  const double dsCm = pathMm/(10.*nstep);
  for(int i=0; i<nstep; ++i){
    const double dedx = BetheBloch(1e-3*pMeV, massGeV, rho, iEv, zOverA);
    if(!std::isfinite(dedx)) return NAN;
    energy += (inverse ? 1. : -1.)*dedx*dsCm;
    if(energy<=massMeV) return NAN;
    pMeV = std::sqrt(energy*energy-massMeV*massMeV);
  }
  return 1e-3*pMeV;
}

double TargetSlabPath(double vertexZ, double thickness, double u, double v,
                      bool incoming)
{
  if(!(thickness>0.) || !std::isfinite(vertexZ)) return 0.;
  const double half=.5*thickness;
  const double dz = incoming ? std::clamp(vertexZ+half, 0., thickness)
                             : std::clamp(half-vertexZ, 0., thickness);
  return dz*std::sqrt(1.+u*u+v*v);
}

// Distance from an interior point to a cylindrical radial boundary.  The
// cylinder axis is y, matching the vertical E90 liquid-target cell.  Slopes
// are dx/dz and dy/dz in the S-2S target coordinate system.
double CylinderRadialDistance(double vertexX, double vertexZ,
                              double u, double v, double radius,
                              bool forward)
{
  if(!(radius > 0.) || !std::isfinite(vertexX) ||
     !std::isfinite(vertexZ)) return NAN;
  const double norm=std::sqrt(1.+u*u+v*v);
  const double sign=forward ? 1. : -1.;
  const double dx=sign*u/norm, dz=sign/norm;
  const double a=dx*dx+dz*dz;
  const double b=2.*(vertexX*dx+vertexZ*dz);
  const double c=vertexX*vertexX+vertexZ*vertexZ-radius*radius;
  const double disc=b*b-4.*a*c;
  if(!(a>0.) || disc<0.) return NAN;
  const double root=(-b+std::sqrt(std::max(0.,disc)))/(2.*a);
  return root >= 0. ? root : NAN;
}

double CylinderShellPath(double vertexX, double vertexZ,
                         double u, double v, double innerRadius,
                         double thickness, bool forward)
{
  if(!(thickness>0.)) return 0.;
  const double inner=CylinderRadialDistance(vertexX,vertexZ,u,v,
                                             innerRadius,forward);
  const double outer=CylinderRadialDistance(vertexX,vertexZ,u,v,
                                             innerRadius+thickness,forward);
  return std::isfinite(inner) && std::isfinite(outer) && outer>=inner
    ? outer-inner : NAN;
}

bool InsideTarget(const TVector3& vertex, const Options& o, double margin)
{
  if(!FiniteVector(vertex)) return false;
  if(o.targetShape == "cylinder"){
    const double dx=vertex.X()-o.targetCenterX;
    const double radius=o.targetRadius+margin;
    return radius>0. && dx*dx+vertex.Z()*vertex.Z()<radius*radius &&
      std::abs(vertex.Y()-o.targetCenterY)<o.targetHalfLengthY+margin;
  }
  return std::abs(vertex.X()-o.targetCenterX)<.5*o.targetSizeX &&
    std::abs(vertex.Y()-o.targetCenterY)<.5*o.targetSizeY &&
    std::abs(vertex.Z())<.5*o.targetThickness+margin;
}

bool InsideSelection(const TVector3& vertex, double closeDist, const Options& o)
{
  return std::isfinite(closeDist) && InsideTarget(vertex,o,o.targetZMargin) &&
    closeDist<o.closeDistMax;
}

double Quantile(std::vector<double> values, double probability)
{
  if(values.empty()) return NAN;
  std::sort(values.begin(), values.end());
  const double pos=probability*(values.size()-1);
  const auto lo=static_cast<std::size_t>(std::floor(pos));
  const auto hi=static_cast<std::size_t>(std::ceil(pos));
  if(lo==hi) return values[lo];
  return values[lo]+(pos-lo)*(values[hi]-values[lo]);
}

double Quantile68(const std::vector<double>& values)
{
  return .5*(Quantile(values,.841344746)-Quantile(values,.158655254));
}

double Mean(const std::vector<double>& values)
{
  if(values.empty()) return NAN;
  return std::accumulate(values.begin(),values.end(),0.)/values.size();
}

double Variance(const std::vector<double>& values)
{
  if(values.empty()) return NAN;
  const double mean=Mean(values);
  double sum=0.;
  for(const double v:values) sum+=(v-mean)*(v-mean);
  return sum/values.size();
}

double Covariance(const std::vector<double>& x,const std::vector<double>& y)
{
  if(x.empty() || x.size()!=y.size()) return NAN;
  const double mx=Mean(x), my=Mean(y);
  double sum=0.;
  for(std::size_t i=0;i<x.size();++i) sum+=(x[i]-mx)*(y[i]-my);
  return sum/x.size();
}

struct Summary {
  long long filled=0, regular=0, underflow=0, overflow=0;
  double mean=NAN, rms=NAN, q68=NAN, crossingFwhm=NAN;
  double fitMean=NAN, fitSigma=NAN, fitFwhm=NAN, fitFwhmError=NAN;
  double fitLow=NAN, fitHigh=NAN;
  double chi2=NAN;
  int ndf=0, fitStatus=-999, covariance=-1;
  std::string quality="insufficient_statistics";
};

Summary Summarize(TH1D& h, const std::vector<double>& values)
{
  // Core FWHM strategy (E90 fine-bin spectra):
  // 1) half-maximum crossing remains a diagnostic (crossing_fwhm)
  // 2) primary Gaussian core uses the central 68% window [q16,q84]
  //    so a single-bin Poisson spike cannot collapse fit_fwhm to ~bin width
  // 3) if [q16,q84] fails, fall back to a validated half-max window only when
  //    its width is compatible with q68 (not a spike)
  Summary out;
  out.filled=static_cast<long long>(values.size());
  out.regular=static_cast<long long>(h.Integral(1,h.GetNbinsX()));
  out.underflow=static_cast<long long>(h.GetBinContent(0));
  out.overflow=static_cast<long long>(h.GetBinContent(h.GetNbinsX()+1));
  out.mean=Mean(values); out.rms=std::sqrt(Variance(values));
  out.q68=Quantile68(values);
  if(out.regular<200 || h.GetMaximum()<=0. || !(out.q68>0.)) return out;

  const double q16=Quantile(values,.158655254);
  const double q84=Quantile(values,.841344746);
  const int peak=h.GetMaximumBin();
  const double half=.5*h.GetBinContent(peak);
  int left=peak-1, right=peak+1;
  while(left>=1 && h.GetBinContent(left)>=half) --left;
  while(right<=h.GetNbinsX() && h.GetBinContent(right)>=half) ++right;
  bool haveHalf=false;
  if(left>=1 && right<=h.GetNbinsX()){
    auto cross=[&](int b1, int b2){
      const double x1=h.GetBinCenter(b1), x2=h.GetBinCenter(b2);
      const double y1=h.GetBinContent(b1), y2=h.GetBinContent(b2);
      return y2==y1 ? .5*(x1+x2) : x1+(half-y1)*(x2-x1)/(y2-y1);
    };
    out.crossingFwhm=cross(right-1,right)-cross(left,left+1);
    haveHalf=std::isfinite(out.crossingFwhm) && out.crossingFwhm>0.;
  }

  auto runFit=[&](double low, double high, const char* tag)->bool{
    if(!(high>low)) return false;
    TF1 fit((std::string("fit_")+tag+"_"+h.GetName()).c_str(),"gaus",low,high);
    fit.SetParameters(h.GetMaximum(),.5*(low+high),
                      std::max(h.GetBinWidth(1),(high-low)/kFwhm));
    const TFitResultPtr result=h.Fit(&fit,"QRNSI0");
    out.fitStatus=int(result);
    if(result.Get()) out.covariance=result->CovMatrixStatus();
    out.fitLow=low; out.fitHigh=high;
    out.fitMean=fit.GetParameter(1); out.fitSigma=std::abs(fit.GetParameter(2));
    out.fitFwhm=kFwhm*out.fitSigma; out.fitFwhmError=kFwhm*fit.GetParError(2);
    out.chi2=fit.GetChisquare(); out.ndf=fit.GetNDF();
    auto* saved=static_cast<TF1*>(fit.Clone());
    saved->SetName((std::string("fit_")+h.GetName()).c_str());
    h.GetListOfFunctions()->Clear();
    h.GetListOfFunctions()->Add(saved);
    return out.fitStatus==0 && out.covariance>=2 && out.ndf>=1 &&
           out.fitFwhm>0. && out.fitFwhm/out.q68>=1. && out.fitFwhm/out.q68<=4.;
  };

  // Prefer quantile-window core (stable on fine bins / non-Gaussian cores).
  if(runFit(q16,q84,"q68")){
    out.quality = (out.chi2/out.ndf>3.) ? "shape_chi2_failed" : "ok";
    return out;
  }
  // Fall back to half-max only when it is not a narrow spike.
  if(haveHalf && out.crossingFwhm>=.5*out.q68 &&
     runFit(h.GetBinCenter(left),h.GetBinCenter(right),"half")){
    out.quality = (out.chi2/out.ndf>3.) ? "shape_chi2_failed" : "ok";
    return out;
  }
  if(out.fitStatus!=0) out.quality="fit_failed";
  else if(out.covariance<2 || out.ndf<1) out.quality="binning_limited";
  else if(!(out.q68>0.) || !(out.fitFwhm>0.) ||
          out.fitFwhm/out.q68<1. || out.fitFwhm/out.q68>4.)
    out.quality="shape_width_failed";
  else out.quality="shape_chi2_failed";
  return out;
}

struct Item {
  std::string name, title, unit;
  std::unique_ptr<TH1D> hist;
  std::vector<double> values;
  Summary stat;
};
}

int main(int argc, char** argv)
{
  try{
    const auto o=Parse(argc,argv); gROOT->SetBatch();
    TFile fg(o.g4.c_str(),"READ"), fk(o.k18.c_str(),"READ"), fs(o.s2s.c_str(),"READ");
    auto* tg=dynamic_cast<TTree*>(fg.Get("g4s2s"));
    auto* tk=dynamic_cast<TTree*>(fk.Get("k18track"));
    auto* ts=dynamic_cast<TTree*>(fs.Get("s2s"));
    if(!tg || !tk || !ts) throw std::runtime_error("missing g4s2s/k18track/s2s tree");
    RequireBranches(tg,"g4s2s",{
      "evnum","ReactionBeamVertex","ReactionScat"
    });
    if(o.bh2Thickness>0.) RequireBranches(tg,"g4s2s",{"K18BH2"});
    RequireBranches(tk,"k18track",{
      "evnum","ntK18","p_3rd","xtgtK18","ytgtK18","utgtK18",
      "vtgtK18","chisqrK18"
    });
    RequireBranches(ts,"s2s",{
      "evnum","ntS2s","pS2s","xtgtS2s","ytgtS2s","utgtS2s",
      "vtgtS2s","chisqrS2s"
    });
    const long long generated=tg->GetEntries();
    if(generated!=tk->GetEntries() || generated!=ts->GetEntries())
      throw std::runtime_error("entry count mismatch");

    TTreeReader rg(tg), rk(tk), rs(ts);
    TTreeReaderValue<Int_t> eg(rg,"evnum"), ek(rk,"evnum"), es(rs,"evnum");
    TTreeReaderValue<Int_t> nk(rk,"ntK18"), ns(rs,"ntS2s");
    TTreeReaderArray<TParticle> beamTruth(rg,"ReactionBeamVertex"),
      scatTruth(rg,"ReactionScat"), bh2Truth(rg,"K18BH2");
    TTreeReaderArray<Double_t> pk(rk,"p_3rd"), xk(rk,"xtgtK18"),
      yk(rk,"ytgtK18"), uk(rk,"utgtK18"), vk(rk,"vtgtK18"),
      chiK(rk,"chisqrK18");
    TTreeReaderArray<Double_t> ps(rs,"pS2s"), xs(rs,"xtgtS2s"),
      ys(rs,"ytgtS2s"), us(rs,"utgtS2s"), vs(rs,"vtgtS2s"),
      chiS(rs,"chisqrS2s");

    std::vector<Item> items;
    auto add=[&](const char* name,const char* title,const char* unit,
                 int bins,double low,double high){
      Item item; item.name=name; item.title=title; item.unit=unit;
      item.hist=std::make_unique<TH1D>(name,title,bins,low,high);
      item.hist->SetDirectory(nullptr); item.hist->Sumw2();
      items.push_back(std::move(item));
    };
    add("h_mm_raw_full","No BB, both reconstructed;M_{X}-M_{X,true} [MeV];events","MeV",80,-10.,10.);
    add("h_mm_bb_full_reco_vtx","BB to reconstructed vertex, both reconstructed;M_{X}-M_{X,true} [MeV];events","MeV",80,-10.,10.);
    add("h_mm_bb_full_truth_vtx","BB to truth vertex, both reconstructed;M_{X}-M_{X,true} [MeV];events","MeV",80,-10.,10.);
    add("h_mm_bb_k18_only_reco_vtx","BB to reconstructed vertex, K18 reconstructed only;M_{X}-M_{X,true} [MeV];events","MeV",80,-10.,10.);
    add("h_mm_bb_s2s_only_reco_vtx","BB to reconstructed vertex, S-2S reconstructed only;M_{X}-M_{X,true} [MeV];events","MeV",80,-10.,10.);
    add("h_mm_bb_linear_sum","K18-only + S-2S-only;sum [MeV];events","MeV",80,-10.,10.);
    add("h_mm_bb_interaction","Full - K18-only - S-2S-only;nonlinear closure [MeV];events","MeV",40,-1.,1.);
    add("h_mm_bb_vertex_path_delta","Reco-vertex BB - truth-vertex BB;#Delta M_{X} [MeV];events","MeV",80,-10.,10.);
    add("h_k18_dp_vertex_bb","K18 BB-corrected momentum at vertex;p_{K18}^{BB}-p_{truth,vtx} [MeV/c];events","MeV/c",400,-200.,200.);
    add("h_s2s_dp_vertex_bb","S-2S BB-corrected momentum at vertex;p_{S2S}^{BB}-p_{truth,vtx} [MeV/c];events","MeV/c",400,-200.,200.);
    add("h_s2s_dp_raw","Raw S-2S plane mismatch;p_{S2S}-p_{truth,vtx} [MeV/c];events","MeV/c",400,-200.,200.);
    add("h_vtx_z_residual","Reaction vertex z residual;z_{reco}-z_{truth} [mm];events","mm",200,-200.,200.);
    TH2D hContrib("h_mm_k18_vs_s2s_contribution",
      "One-arm missing-mass contributions;K18-only #DeltaM_{X} [MeV];S-2S-only #DeltaM_{X} [MeV]",
      400,-10.,10.,400,-10.,10.);
    TH2D hS2sRawZ("h_s2s_raw_dp_vs_truth_z",
      "Raw S-2S reference-plane mismatch;z_{truth} [mm];p_{S2S}-p_{truth,vtx} [MeV/c]",
      400,-20.,20.,400,-20.,20.);
    TH2D hS2sBBZ("h_s2s_bb_dp_vs_truth_z",
      "S-2S after BB to reconstructed vertex;z_{truth} [mm];p_{S2S}^{BB}-p_{truth,vtx} [MeV/c]",
      400,-20.,20.,400,-20.,20.);
    TH2D hMMFullZ("h_mm_bb_full_vs_truth_z",
      "Full BB-flow missing mass;z_{truth} [mm];M_{X}-M_{X,true} [MeV]",
      400,-20.,20.,400,-10.,10.);
    TH2D hMMFullX("h_mm_bb_full_vs_truth_x",
      "Full BB-flow missing mass;x_{truth} [mm];M_{X}-M_{X,true} [MeV]",
      400,-40.,40.,400,-10.,10.);
    TH2D hTruthXZ("h_truth_vertex_xz",
      "Generated reaction vertices inside target;x_{truth} [mm];z_{truth} [mm]",
      160,-40.,40.,160,-40.,40.);
    TH1D hMMTruthAbs("h_mm_truth_selected_abs",
      "Selected truth missing mass;M_{X,true} [MeV/c^{2}];events / 0.2 MeV",
      1250,2000.,2250.);
    TH1D hMMTruthPreselectionAbs("h_mm_truth_preselection_abs",
      "Truth missing mass before tracking selection;M_{X,true} [MeV/c^{2}];events / 0.2 MeV",
      1250,2000.,2250.);
    TH1D hMMRecoAbs("h_mm_reco_full_abs",
      "Selected reconstructed missing mass;M_{X,reco} [MeV/c^{2}];events / 0.2 MeV",
      1250,2000.,2250.);
    TH1D hTargetPathIn("h_target_path_in",
      "Incoming target-material path;path [mm];events",240,0.,60.);
    TH1D hTargetPathOut("h_target_path_out",
      "Outgoing target-material path;path [mm];events",240,0.,60.);
    TH1D hCellPathIn("h_cell_path_in",
      "Incoming target-cell path;path [mm];events",200,0.,2.);
    TH1D hCellPathOut("h_cell_path_out",
      "Outgoing target-cell path;path [mm];events",200,0.,2.);
    TH1D hVacuumShellPathIn("h_vacuum_shell_path_in",
      "Incoming target-vacuum-shell path;path [mm];events",1000,0.,10.);
    TH1D hVacuumShellPathOut("h_vacuum_shell_path_out",
      "Outgoing target-vacuum-shell path;path [mm];events",1000,0.,10.);
    for(auto* h:{&hContrib,&hS2sRawZ,&hS2sBBZ,&hMMFullZ,&hMMFullX,
                 &hTruthXZ}) h->SetDirectory(nullptr);
    for(auto* h:{&hMMTruthAbs,&hMMTruthPreselectionAbs,&hMMRecoAbs,
                 &hTargetPathIn,&hTargetPathOut,&hCellPathIn,&hCellPathOut,
                 &hVacuumShellPathIn,&hVacuumShellPathOut})
      h->SetDirectory(nullptr);
    auto fill=[&](std::size_t i,double value){
      if(!std::isfinite(value)) return;
      items[i].hist->Fill(value); items[i].values.push_back(value);
    };
    std::vector<double> truthZValues;

    TFile out(o.root.c_str(),"RECREATE");
    if(out.IsZombie()) throw std::runtime_error("cannot create output ROOT file");
    out.cd();

    int evnum=0; Long64_t sourceEntry=0;
    double xReco=0.,yReco=0.,zReco=0.,xTruth=0.,yTruth=0.,zTruth=0.;
    double closeDist=0., theta=0.;
    double mmTruthAbs=0., mmRecoAbs=0.;
    double mmRaw=0., mmFull=0., mmTruthVtx=0., mmK18=0., mmS2s=0.;
    double pK18BB=0., pS2sBB=0., pS2sBBTruthVtx=0.;
    double targetPathIn=0., targetPathOut=0., cellPathIn=0., cellPathOut=0.;
    double vacuumShellPathIn=0., vacuumShellPathOut=0.;
    TTree flat("bbflow","branch-limited reaction vertex + BB diagnostics");
    flat.Branch("evnum",&evnum,"evnum/I");
    flat.Branch("source_entry",&sourceEntry,"source_entry/L");
    flat.Branch("x_reco",&xReco,"x_reco/D");
    flat.Branch("y_reco",&yReco,"y_reco/D");
    flat.Branch("z_reco",&zReco,"z_reco/D");
    flat.Branch("x_truth",&xTruth,"x_truth/D");
    flat.Branch("y_truth",&yTruth,"y_truth/D");
    flat.Branch("z_truth",&zTruth,"z_truth/D");
    flat.Branch("close_dist",&closeDist,"close_dist/D");
    flat.Branch("theta",&theta,"theta/D");
    flat.Branch("mm_truth_abs",&mmTruthAbs,"mm_truth_abs/D");
    flat.Branch("mm_reco_abs",&mmRecoAbs,"mm_reco_abs/D");
    flat.Branch("mm_raw",&mmRaw,"mm_raw/D");
    flat.Branch("mm_bb_full",&mmFull,"mm_bb_full/D");
    flat.Branch("mm_bb_truth_vtx",&mmTruthVtx,"mm_bb_truth_vtx/D");
    flat.Branch("mm_bb_k18_only",&mmK18,"mm_bb_k18_only/D");
    flat.Branch("mm_bb_s2s_only",&mmS2s,"mm_bb_s2s_only/D");
    flat.Branch("p_k18_bb",&pK18BB,"p_k18_bb/D");
    flat.Branch("p_s2s_bb",&pS2sBB,"p_s2s_bb/D");
    flat.Branch("p_s2s_bb_truth_vtx",&pS2sBBTruthVtx,
                "p_s2s_bb_truth_vtx/D");
    flat.Branch("target_path_in",&targetPathIn,"target_path_in/D");
    flat.Branch("target_path_out",&targetPathOut,"target_path_out/D");
    flat.Branch("cell_path_in",&cellPathIn,"cell_path_in/D");
    flat.Branch("cell_path_out",&cellPathOut,"cell_path_out/D");
    flat.Branch("vacuum_shell_path_in",&vacuumShellPathIn,
                "vacuum_shell_path_in/D");
    flat.Branch("vacuum_shell_path_out",&vacuumShellPathOut,
                "vacuum_shell_path_out/D");

    long long evnumMismatch=0, primaryFail=0, truthTargetFail=0, pidFail=0,
      bh2Fail=0,
      multiplicityFail=0, trackFail=0, chi2Fail=0, thetaFail=0,
      vertexFail=0, selected=0;
    while(true){
      const bool ag=rg.Next(), ak=rk.Next(), as=rs.Next();
      if(ag!=ak || ag!=as) throw std::runtime_error("reader length mismatch");
      if(!ag) break;
      if(*eg!=*ek || *eg!=*es){ ++evnumMismatch; continue; }
      const auto* bt=SinglePrimary(beamTruth);
      const auto* st=SinglePrimary(scatTruth);
      if(!bt || !st){ ++primaryFail; continue; }
      const TVector3 truthVertex(st->Vx(),st->Vy(),st->Vz()-o.targetPlaneZ);
      if(o.targetShape == "cylinder" && !InsideTarget(truthVertex,o,0.)){
        ++truthTargetFail; continue;
      }
      if(bt->GetPdgCode()!=o.beamPdg || st->GetPdgCode()!=o.scatPdg){
        ++pidFail; continue;
      }
      if(o.bh2Thickness>0. && !HasPrimary(bh2Truth)){ ++bh2Fail; continue; }
      const auto bTruth=ParticleVector(*bt), sTruth=ParticleVector(*st);
      const double truthMass=MissingMass(bTruth,sTruth,o.targetMass);
      if(!std::isfinite(truthMass)){ ++primaryFail; continue; }
      mmTruthAbs=1e3*truthMass;
      hMMTruthPreselectionAbs.Fill(mmTruthAbs);
      if(*nk!=1 || *ns!=1){ ++multiplicityFail; continue; }
      if(pk.GetSize()<1 || xk.GetSize()<1 || yk.GetSize()<1 ||
         uk.GetSize()<1 || vk.GetSize()<1 || ps.GetSize()<1 || xs.GetSize()<1 ||
         ys.GetSize()<1 || us.GetSize()<1 || vs.GetSize()<1 || chiK.GetSize()<1 ||
         chiS.GetSize()<1 ||
         !std::isfinite(pk[0]) || !(pk[0]>0.) || !std::isfinite(ps[0]) ||
         !(ps[0]>0.)){
        ++trackFail; continue;
      }
      if(!std::isfinite(chiK[0]) || !std::isfinite(chiS[0]) ||
         chiK[0]>=o.chi2Max || chiS[0]>=o.chi2Max){
        ++chi2Fail; continue;
      }
      const TVector3 xK(xk[0],yk[0],0.), xS(xs[0],ys[0],0.);
      const auto bRecoRaw=FourVector(pk[0],uk[0],vk[0],o.beamMass);
      const auto sRecoRaw=FourVector(ps[0],us[0],vs[0],o.scatMass);
      const TVector3 pK=bRecoRaw.Vect(), pS=sRecoRaw.Vect();
      const double denom=pK.Mag()*pS.Mag();
      if(!(denom>0.)){ ++trackFail; continue; }
      const double cosTheta=std::clamp(pK.Dot(pS)/denom,-1.,1.);
      theta=std::acos(cosTheta)*180./std::acos(-1.);
      if(!std::isfinite(theta) || !(theta>o.thetaMin && theta<o.thetaMax)){
        ++thetaFail; continue;
      }
      closeDist=NAN;
      const TVector3 vertex=Kinematics::VertexPoint(xK,xS,pK,pS,closeDist);
      if(!InsideSelection(vertex,closeDist,o)){ ++vertexFail; continue; }
      xReco=vertex.X(); yReco=vertex.Y(); zReco=vertex.Z();
      xTruth=truthVertex.X(); yTruth=truthVertex.Y(); zTruth=truthVertex.Z();

      const double bh2Path=o.bh2Thickness*std::sqrt(1.+uk[0]*uk[0]+vk[0]*vk[0]);
      constexpr double bh2Rho=1.032, bh2IEv=64.7, bh2ZOverA=.54141;
      const double pAfterBH2=PropagateBB(pk[0],o.beamMass,bh2Path,
                                         bh2Rho,bh2IEv,bh2ZOverA,false);
      struct Paths { double targetIn=0., targetOut=0., cellIn=0., cellOut=0.;
        double shellIn=0., shellOut=0.; };
      auto pathsFor=[&](const TVector3& point){
        Paths p;
        if(o.targetShape == "cylinder"){
          const double x=point.X()-o.targetCenterX, z=point.Z();
          p.targetIn=CylinderRadialDistance(x,z,uk[0],vk[0],
                                             o.targetRadius,false);
          p.targetOut=CylinderRadialDistance(x,z,us[0],vs[0],
                                              o.targetRadius,true);
          p.cellIn=CylinderShellPath(x,z,uk[0],vk[0],o.cellRadius,
                                      o.cellThickness,false);
          p.cellOut=CylinderShellPath(x,z,us[0],vs[0],o.cellRadius,
                                       o.cellThickness,true);
          p.shellIn=CylinderShellPath(x,z,uk[0],vk[0],o.vacuumShellRadius,
                                       o.vacuumShellThickness,false);
          p.shellOut=CylinderShellPath(x,z,us[0],vs[0],o.vacuumShellRadius,
                                        o.vacuumShellThickness,true);
        }else{
          p.targetIn=TargetSlabPath(point.Z(),o.targetThickness,
                                    uk[0],vk[0],true);
          p.targetOut=TargetSlabPath(point.Z(),o.targetThickness,
                                     us[0],vs[0],false);
        }
        return p;
      };
      auto correct=[&](const Paths& path){
        double pin=pAfterBH2;
        pin=PropagateBB(pin,o.beamMass,path.shellIn,o.vacuumShellDensity,
                        o.vacuumShellIEv,o.vacuumShellZOverA,false);
        pin=PropagateBB(pin,o.beamMass,path.cellIn,o.cellDensity,
                        o.cellIEv,o.cellZOverA,false);
        pin=PropagateBB(pin,o.beamMass,path.targetIn,o.targetDensity,
                        o.targetIEv,o.targetZOverA,false);
        double pout=ps[0];
        pout=PropagateBB(pout,o.scatMass,path.shellOut,o.vacuumShellDensity,
                         o.vacuumShellIEv,o.vacuumShellZOverA,true);
        pout=PropagateBB(pout,o.scatMass,path.cellOut,o.cellDensity,
                         o.cellIEv,o.cellZOverA,true);
        pout=PropagateBB(pout,o.scatMass,path.targetOut,o.targetDensity,
                         o.targetIEv,o.targetZOverA,true);
        return std::make_pair(pin,pout);
      };
      const auto recoPaths=pathsFor(vertex), truthPaths=pathsFor(truthVertex);
      const auto ppReco=correct(recoPaths), ppTruth=correct(truthPaths);
      if(!std::isfinite(ppReco.first) || !std::isfinite(ppReco.second) ||
         !std::isfinite(ppTruth.first) || !std::isfinite(ppTruth.second)){
        ++trackFail; continue;
      }
      const auto bReco=FourVector(ppReco.first,uk[0],vk[0],o.beamMass);
      const auto sReco=FourVector(ppReco.second,us[0],vs[0],o.scatMass);
      const auto bRecoTruthVtx=FourVector(ppTruth.first,uk[0],vk[0],o.beamMass);
      const auto sRecoTruthVtx=FourVector(ppTruth.second,us[0],vs[0],o.scatMass);
      mmRaw=Residual(bRecoRaw,sRecoRaw,o.targetMass,truthMass);
      mmFull=Residual(bReco,sReco,o.targetMass,truthMass);
      mmRecoAbs=mmTruthAbs+mmFull;
      mmTruthVtx=Residual(bRecoTruthVtx,sRecoTruthVtx,o.targetMass,truthMass);
      mmK18=Residual(bReco,sTruth,o.targetMass,truthMass);
      mmS2s=Residual(bTruth,sReco,o.targetMass,truthMass);
      pK18BB=ppReco.first; pS2sBB=ppReco.second;
      pS2sBBTruthVtx=ppTruth.second; evnum=*eg;
      sourceEntry=rg.GetCurrentEntry();
      targetPathIn=recoPaths.targetIn; targetPathOut=recoPaths.targetOut;
      cellPathIn=recoPaths.cellIn; cellPathOut=recoPaths.cellOut;
      vacuumShellPathIn=recoPaths.shellIn;
      vacuumShellPathOut=recoPaths.shellOut;
      fill(0,mmRaw); fill(1,mmFull); fill(2,mmTruthVtx); fill(3,mmK18);
      fill(4,mmS2s); fill(5,mmK18+mmS2s); fill(6,mmFull-mmK18-mmS2s);
      fill(7,mmFull-mmTruthVtx); fill(8,1e3*(pK18BB-bTruth.P()));
      fill(9,1e3*(pS2sBB-sTruth.P())); fill(10,1e3*(ps[0]-sTruth.P()));
      fill(11,zReco-zTruth);
      truthZValues.push_back(zTruth);
      hContrib.Fill(mmK18,mmS2s);
      hS2sRawZ.Fill(zTruth,1e3*(ps[0]-sTruth.P()));
      hS2sBBZ.Fill(zTruth,1e3*(pS2sBB-sTruth.P()));
      hMMFullZ.Fill(zTruth,mmFull);
      hMMFullX.Fill(truthVertex.X(),mmFull);
      hTruthXZ.Fill(truthVertex.X(),zTruth);
      hMMTruthAbs.Fill(mmTruthAbs); hMMRecoAbs.Fill(mmRecoAbs);
      hTargetPathIn.Fill(targetPathIn); hTargetPathOut.Fill(targetPathOut);
      hCellPathIn.Fill(cellPathIn); hCellPathOut.Fill(cellPathOut);
      hVacuumShellPathIn.Fill(vacuumShellPathIn);
      hVacuumShellPathOut.Fill(vacuumShellPathOut);
      flat.Fill(); ++selected;
    }

    out.cd();
    WriteDstMetadata(&out,
      "entry-aligned g4s2s + k18track + s2s; D2U p_3rd is the primary K18 momentum",
      "truth beam/scatter species; truth reaction vertex inside the configured target; primary K18BH2 hit; ntK18==1; ntS2s==1; both track chi2<max; opening-angle window; reconstructed vertex inside the configured target plus margin; closeDist<max; detector m2 PID unavailable in tracking-only MC");
    std::ostringstream definition;
    definition << std::setprecision(12)
      << "Kinematics::VertexPoint shared by K18 and S-2S; incoming BB: BH2 then target entrance to vertex; outgoing inverse BB: target exit back to vertex; "
      << "beamMass=" << o.beamMass << " scatMass=" << o.scatMass
      << " targetMass=" << o.targetMass << " targetThicknessMm=" << o.targetThickness
      << " targetDensity=" << o.targetDensity << " targetIEv=" << o.targetIEv
      << " targetZOverA=" << o.targetZOverA << " targetSizeXYmm="
      << o.targetSizeX << "x" << o.targetSizeY << " targetPlaneZmm=" << o.targetPlaneZ
      << " targetShape=" << o.targetShape << " targetRadiusMm=" << o.targetRadius
      << " targetHalfLengthYmm=" << o.targetHalfLengthY
      << " cellRadiusThicknessMm=" << o.cellRadius << ":" << o.cellThickness
      << " cellDensityIEvZOverA=" << o.cellDensity << ":" << o.cellIEv
      << ":" << o.cellZOverA
      << " vacuumShellRadiusThicknessMm=" << o.vacuumShellRadius << ":"
      << o.vacuumShellThickness
      << " vacuumShellDensityIEvZOverA=" << o.vacuumShellDensity << ":"
      << o.vacuumShellIEv << ":" << o.vacuumShellZOverA
      << " targetZMarginMm=" << o.targetZMargin << " closeDistMaxMm=" << o.closeDistMax
      << " bh2ThicknessMm=" << o.bh2Thickness << " chi2Max=" << o.chi2Max
      << " thetaDeg=" << o.thetaMin << ":" << o.thetaMax
      << " truthPdg=" << o.beamPdg << ":" << o.scatPdg
      << " selectionName=GoodPiKEquivalent";
    TNamed("dst_definition",definition.str().c_str()).Write();
    TNamed("dst_counts",("generated="+std::to_string(generated)+
      " selected="+std::to_string(selected)+" evnum_mismatch="+
      std::to_string(evnumMismatch)+" primary_fail="+std::to_string(primaryFail)+
      " truth_target_fail="+std::to_string(truthTargetFail)+
      " pid_fail="+std::to_string(pidFail)+" bh2_fail="+std::to_string(bh2Fail)+
      " multiplicity_fail="+std::to_string(multiplicityFail)+" track_fail="+
      std::to_string(trackFail)+" chi2_fail="+std::to_string(chi2Fail)+
      " theta_fail="+std::to_string(thetaFail)+" vertex_fail="+
      std::to_string(vertexFail)).c_str()).Write();
    flat.Write();
    hContrib.Write(); hS2sRawZ.Write(); hS2sBBZ.Write(); hMMFullZ.Write();
    hMMFullX.Write(); hTruthXZ.Write(); hMMTruthAbs.Write();
    hMMTruthPreselectionAbs.Write(); hMMRecoAbs.Write();
    hTargetPathIn.Write(); hTargetPathOut.Write(); hCellPathIn.Write();
    hCellPathOut.Write(); hVacuumShellPathIn.Write();
    hVacuumShellPathOut.Write();

    for(auto& item:items){
      item.stat=Summarize(*item.hist,item.values); item.hist->Write();
    }
    const double varFull=Variance(items[1].values);
    const double varK18=Variance(items[3].values);
    const double varS2s=Variance(items[4].values);
    const double covK18S2s=Covariance(items[3].values,items[4].values);
    const double corrK18S2s=covK18S2s/std::sqrt(varK18*varS2s);
    const double shapleyK18=.5*(varK18+varFull-varS2s);
    const double shapleyS2s=.5*(varS2s+varFull-varK18);
    const double shapleyK18Fraction=shapleyK18/varFull;
    const double shapleyS2sFraction=shapleyS2s/varFull;
    const auto correlation=[&](const std::vector<double>& x,
                               const std::vector<double>& y){
      return Covariance(x,y)/std::sqrt(Variance(x)*Variance(y));
    };
    const double corrTruthZS2sRaw=correlation(truthZValues,items[10].values);
    const double corrTruthZS2sBB=correlation(truthZValues,items[9].values);
    const double corrTruthZMMRaw=correlation(truthZValues,items[0].values);
    const double corrTruthZMMBB=correlation(truthZValues,items[1].values);
    std::ostringstream budget;
    budget << std::setprecision(12)
      << "var_full_MeV2=" << varFull << " var_k18_only_MeV2=" << varK18
      << " var_s2s_only_MeV2=" << varS2s << " cov_k18_s2s_MeV2=" << covK18S2s
      << " corr_k18_s2s=" << corrK18S2s << " shapley_k18_MeV2=" << shapleyK18
      << " shapley_s2s_MeV2=" << shapleyS2s
      << " shapley_k18_fraction=" << shapleyK18Fraction
      << " shapley_s2s_fraction=" << shapleyS2sFraction
      << " corr_truth_z_s2s_raw=" << corrTruthZS2sRaw
      << " corr_truth_z_s2s_bb=" << corrTruthZS2sBB
      << " corr_truth_z_mm_raw=" << corrTruthZMMRaw
      << " corr_truth_z_mm_bb=" << corrTruthZMMBB;
    TNamed("dst_arm_variance_budget",budget.str().c_str()).Write();

    std::ofstream csv(o.csv);
    csv << "label,quantity,unit,generated,selected,evnum_mismatch,primary_fail,truth_target_fail,pid_fail,bh2_fail,multiplicity_fail,track_fail,chi2_fail,theta_fail,vertex_fail,filled,regular,underflow,overflow,in_range_fraction,mean,rms,q68,crossing_fwhm,fit_low,fit_high,fit_mean,fit_sigma,fit_fwhm,fit_fwhm_error,chi2,ndf,fit_status,covariance,quality,var_full,var_k18_only,var_s2s_only,cov_k18_s2s,corr_k18_s2s,shapley_k18_variance,shapley_s2s_variance,shapley_k18_fraction,shapley_s2s_fraction,corr_truth_z_s2s_raw,corr_truth_z_s2s_bb,corr_truth_z_mm_raw,corr_truth_z_mm_bb\n";
    csv << std::setprecision(12);
    for(const auto& item:items){
      const auto& s=item.stat;
      const double fraction=s.filled>0 ? double(s.regular)/s.filled : NAN;
      csv << o.label << ',' << item.name.substr(2) << ',' << item.unit << ','
          << generated << ',' << selected << ',' << evnumMismatch << ','
          << primaryFail << ',' << truthTargetFail << ',' << pidFail << ','
          << bh2Fail << ','
          << multiplicityFail << ',' << trackFail << ',' << chi2Fail << ','
          << thetaFail << ',' << vertexFail << ',' << s.filled << ',' << s.regular
          << ',' << s.underflow << ',' << s.overflow << ',' << fraction << ','
          << s.mean << ',' << s.rms << ',' << s.q68 << ',' << s.crossingFwhm
          << ',' << s.fitLow << ',' << s.fitHigh << ',' << s.fitMean << ','
          << s.fitSigma << ',' << s.fitFwhm << ','
          << s.fitFwhmError << ',' << s.chi2 << ',' << s.ndf << ','
          << s.fitStatus << ',' << s.covariance << ',' << s.quality << ','
          << varFull << ',' << varK18 << ',' << varS2s << ',' << covK18S2s
          << ',' << corrK18S2s << ',' << shapleyK18 << ',' << shapleyS2s
          << ',' << shapleyK18Fraction << ',' << shapleyS2sFraction
          << ',' << corrTruthZS2sRaw << ',' << corrTruthZS2sBB
          << ',' << corrTruthZMMRaw << ',' << corrTruthZMMBB << '\n';
    }
    out.Close();
    std::cout << "generated=" << generated << " selected=" << selected
              << " evnum_mismatch=" << evnumMismatch << '\n'
              << o.root << '\n' << o.csv << '\n';
    return evnumMismatch ? 2 : 0;
  }catch(const std::exception& e){
    std::cerr << "#E " << e.what() << '\n'; return 1;
  }
}
