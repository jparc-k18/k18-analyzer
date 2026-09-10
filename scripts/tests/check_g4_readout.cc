// Link with libK18Analyzer and the normal analyzer libraries.
#include <cmath>
#include <iostream>
#include <TMath.h>
#include <TMemFile.h>
#include <TVector3.h>
#include "ConfMan.hh"
#include "DCGeomMan.hh"
#include "DCHit.hh"

Bool_t ConfMan::InitializeParameterFiles() { return true; }
Bool_t ConfMan::InitializeHistograms() { return true; }
Bool_t ConfMan::FinalizeProcess() { return true; }

int main(int argc, char** argv)
{
  if(argc == 3 && std::string(argv[1]) == "--config"){
    TMemFile output("readout-test", "RECREATE");
    return ConfMan::GetInstance().Initialize(argv[2]) ? 0 : 1;
  }
  if(argc != 2) return 2;
  auto& geom = DCGeomMan::GetInstance();
  if(!geom.Initialize(argv[1])) return 2;
  int count = 0;
  double max_difference = 0.;
  for(int layer = 1; layer <= 124; ++layer){
    if(!(layer <= 10 || (31 <= layer && layer <= 42) || layer >= 113)) continue;
    for(double x : {-40., 0., 40.}) for(double y : {-30., 0., 30.}){
      TVector3 chamber(x, y, 0.);
      TVector3 old_wire_local(chamber);
      old_wire_local.RotateZ(-geom.GetTiltAngle(layer)*TMath::DegToRad());
      const double measured = DCHit::CalcGeant4ReadoutPosition(layer, chamber);
      const double delta = std::abs(measured - old_wire_local.X());
      if(!std::isfinite(measured) || delta > 1.e-12){
        std::cerr << "FAIL layer=" << layer << " delta=" << delta << '\n';
        return 1;
      }
      max_difference = std::max(max_difference, delta);
      ++count;
    }
  }
  std::cout << "PASS readout points=" << count
            << " max_abs_difference_mm=" << max_difference << '\n';
  return count == 306 ? 0 : 1;
}
