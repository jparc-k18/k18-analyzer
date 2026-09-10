// Link against libK18Analyzer with the normal analyzer libraries.
// Usage: check_dcgeo_smearing DCGEO config (only G4DCSmearSignedPosition optional).
#include <cmath>
#include <iostream>
#include <TMath.h>
#include <TMemFile.h>
#include <TRandom3.h>
#include <TVector3.h>
#include "ConfMan.hh"
#include "DCGeomMan.hh"
#include "DCHit.hh"

Bool_t ConfMan::InitializeParameterFiles() { return true; }
Bool_t ConfMan::InitializeHistograms() { return true; }
Bool_t ConfMan::FinalizeProcess() { return true; }

int main(int argc, char** argv)
{
  if(argc != 3) return 2;
  TMemFile output("smearing-test", "RECREATE");
  auto& geom = DCGeomMan::GetInstance();
  if(!geom.Initialize(argv[1]) || !ConfMan::GetInstance().Initialize(argv[2])) return 2;
  const bool signed_position = ConfMan::Get<Double_t>("G4DCSmearSignedPosition") != 0.;
  TRandom3 actual(908231), reference(908231);
  TRandom* saved = gRandom;
  gRandom = &actual;
  int points = 0;
  for(int layer = 1; layer <= 124; ++layer){
    if(!(layer <= 10 || (31 <= layer && layer <= 42) || layer >= 113)) continue;
    const auto original_res = geom.GetResolution(layer);
    const auto wire = geom.CalcWireNumber(layer, 0.);
    const auto wpos = geom.CalcWirePosition(layer, wire);
    const auto angle = geom.GetTiltAngle(layer)*TMath::DegToRad();
    for(double res : {0.1, 0.2, 0.45}){
      geom.SetResolution(layer, res);
      for(double offset : {-0.7, 0., 0.7}){
        const double target = wpos + offset;
        TVector3 position(target*std::cos(angle), target*std::sin(angle), 0.);
        const double signed_dl = DCHit::CalcGeant4ReadoutPosition(layer, position)-wpos;
        const bool use_signed = layer < 113 && signed_position;
        const double draw = reference.Gaus(use_signed ? signed_dl : std::abs(signed_dl), res);
        const double expected_dl = use_signed ? std::abs((wpos+draw)-wpos) : draw;
        const double expected_s = wpos + (use_signed ? draw : std::copysign(draw, signed_dl));
        DCHit hit(0, layer, wire);
        hit.SetDCDataGeant4(position, 1.);
        if(!hit.CalcDCObservablesGeant4() || hit.GetEntries() != 1 ||
           std::abs(hit.GetDriftLength(0)-expected_dl) > 1.e-12 ||
           std::abs(hit.GetGeant4SmearedReadoutPosition()-expected_s) > 1.e-12){
          std::cerr << "FAIL layer=" << layer << " Res=" << res << '\n';
          gRandom = saved;
          return 1;
        }
        ++points;
      }
    }
    geom.SetResolution(layer, original_res);
  }
  gRandom = saved;
  std::cout << "PASS DCGEO.Res direct Gaussian samples=" << points
            << " signed_position=" << signed_position << '\n';
  return points == 306 ? 0 : 1;
}
