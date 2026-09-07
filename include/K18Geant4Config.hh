// Configuration contract for DstK18TrackingGeant4. No ROOT dependency.
#ifndef K18_GEANT4_CONFIG_HH
#define K18_GEANT4_CONFIG_HH

#include <cmath>
#include <fstream>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>

namespace k18geant4
{
using Config = std::map<std::string, std::string>;

inline double Number(const Config& config, const std::string& key)
{
  const auto found = config.find(key);
  if(found == config.end()) throw std::runtime_error("missing required key " + key);
  try{
    std::size_t used = 0;
    const double value = std::stod(found->second, &used);
    if(used == found->second.size() && std::isfinite(value)) return value;
  }catch(const std::exception&){}
  throw std::runtime_error(key + " must be a finite number, got " + found->second);
}

inline Config ReadConfig(const std::string& path)
{
  std::ifstream input(path);
  if(!input) throw std::runtime_error("cannot read config " + path);
  Config config;
  std::string line;
  int lineno = 0;
  while(std::getline(input, line)){
    ++lineno;
    line = line.substr(0, line.find('#'));
    // Match ConfMan's supported spelling: KEY: value or KEY value.
    std::istringstream fields(line);
    std::string key, value, extra;
    if(!(fields >> key)) continue;
    if(key.back() == ':') key.pop_back();
    if(!(fields >> value) || fields >> extra || key.empty())
      throw std::runtime_error(path + ":" + std::to_string(lineno) +
                               ": expected KEY: value");
    if(value.size() >= 2 && value.front() == '"' && value.back() == '"')
      value = value.substr(1, value.size()-2);
    if(!config.emplace(key, value).second)
      throw std::runtime_error(path + ":" + std::to_string(lineno) +
                               ": duplicate key " + key);
  }
  return config;
}

inline void ValidateConfig(const Config& config)
{
  // Application-scoped retirement: S-2S still owns the G4DC controls.
  // No compatibility interpretation of old BC response switches is allowed.
  const std::set<std::string> files = {
    "UNPACK", "DIGIT", "CMAP", "DCGEO", "K18TM", "USER", "K18FLDMAP"
  };
  const std::set<std::string> numbers = {
    "K18TrackingConfigVersion", "PK18", "K18BFTAcceptedAbsPdg", "K18BFTClusterGap",
    "K18BFTPositionSmearSigma", "K18HitSmearSeed",
    "K18NativeUseFullFit", "K18NativeMomCorrEnable",
    "K18NativeMomCorrC0", "K18NativeMomCorrC1", "K18NativeMomCorrC2",
    "K18NativeMomCorrX", "K18NativeMomCorrU", "K18NativeMomCorrY",
    "K18NativeMomCorrV", "K18NativeMomCorrXX", "K18NativeMomCorrUU",
    "K18NativeMomCorrYY", "K18NativeMomCorrVV", "K18NativeMomCorrXU",
    "K18NativeMomCorrYV", "K18NativeMomCorrXBFT",
    "K18NativePMin", "K18NativePMax", "K18NativeRKStep",
    "K18NativeFitMaxIteration", "K18NativeFitMomentumTolerance",
    "K18NativeFitResidualTolerance", "K18NativeMaxPath",
    "K18NativeBFTResidualMax", "K18RKCharge", "K18GlobalScale",
    "K18NativeRKDriftStep", "K18FLDNMR", "K18FLDCALC",
    "K18Q10Scale", "K18Q11Scale", "K18D4Scale", "K18Q12Scale", "K18Q13Scale",
    "K18TargetL", "K18BFTL", "K18NativeBFTFitResolution",
    "K18NativeFullFitMaxIteration", "K18NativeFullFitLambda"
  };
  for(const auto& item: config){
    const auto& key = item.first;
    if(key.compare(0, 4, "G4DC") == 0)
      throw std::runtime_error("retired/unsupported K18 key " + key +
        "; BC uses plane-local x and DCGEO.Res for both smearing and fit. "
        "Use scripts/migrate_k18_tracking_config.py for an equivalent config; "
        "old diagnostic response models require their historical checkout.");
    if(files.count(key)){
      if(item.second.empty()) throw std::runtime_error("empty file key " + key);
    }else if(numbers.count(key)){
      Number(config, key);
    }else{
      throw std::runtime_error("unknown DstK18TrackingGeant4 key " + key);
    }
  }
  for(const auto key: {"UNPACK", "DIGIT", "CMAP", "DCGEO", "K18TM", "USER"})
    if(!config.count(key)) throw std::runtime_error(std::string("missing required key ") + key);
  // An unversioned legacy config with no BC scale meant ideal hits. Without
  // an explicit format boundary it would silently acquire noise in this build.
  if(Number(config, "K18TrackingConfigVersion") != 2.)
    throw std::runtime_error("K18TrackingConfigVersion must be 2; migrate legacy configs explicitly");
  if(Number(config, "PK18") <= 0.)
    throw std::runtime_error("PK18 must be positive (GeV/c)");
  const double pdg = Number(config, "K18BFTAcceptedAbsPdg");
  if(pdg <= 0. || pdg != std::floor(pdg) || pdg > std::numeric_limits<int>::max())
    throw std::runtime_error("K18BFTAcceptedAbsPdg must be a positive integer; shared by BFT and BH1/BH2");
  if(Number(config, "K18BFTPositionSmearSigma") < 0.)
    throw std::runtime_error("K18BFTPositionSmearSigma must be nonnegative (mm)");
  if(config.count("K18BFTClusterGap") && Number(config, "K18BFTClusterGap") < 0.)
    throw std::runtime_error("K18BFTClusterGap must be nonnegative (mm)");
  if(config.count("K18HitSmearSeed")){
    const double seed = Number(config, "K18HitSmearSeed");
    if(seed < 0. || seed != std::floor(seed) || seed > 4294967295.)
      throw std::runtime_error("K18HitSmearSeed must be an integer in [0, 4294967295]");
  }
  for(const auto key: {"K18NativeUseFullFit", "K18NativeMomCorrEnable"}){
    if(config.count(key) && Number(config, key) != 0. && Number(config, key) != 1.)
      throw std::runtime_error(std::string(key) + " must be 0 or 1");
  }
}

inline void ValidateConfigFile(const std::string& path)
{
  ValidateConfig(ReadConfig(path));
}
}
#endif
