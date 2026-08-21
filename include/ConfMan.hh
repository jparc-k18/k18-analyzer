/*
  ConfMan.hh

  2012/1/24
*/

#ifndef Confman_h
#define Confman_h 1

#include <map>
#include <string>

class VEvent;
class EvDisp;

class HodoParamMan;
class HodoPHCMan;

class ScalerAnalyzer_;
class ScalerAna;

class DCGeomMan;
class DCTdcCalibMan;
class DCDriftParamMan;
class K18TransMatrix;
class FieldMan;
class SimuData;

class ConfMan
{

private:
  //Conf
  std::string m_ConfFileName;
  std::string m_UnpackerFileName;
  std::string m_DigitFileName;
  std::string m_CmapFileName;
  static ConfMan *confManager_;

  //Hodo
  std::string   HodoParamFileName_;
  HodoParamMan *HodoParamManager_;
  std::string   HodoPHCFileName_;
  HodoPHCMan   *HodoPHCManager_;

  //DC
  std::string      DCGeomFileName_;
  DCGeomMan       *DCGeomManager_;
  std::string      DCTdcCalibFileName_;
  DCTdcCalibMan   *DCTdcCalibManager_;
  std::string      DCDriftParamFileName_;
  DCDriftParamMan *DCDriftParamManager_;
  std::string      K18MatrixFileName_;
  K18TransMatrix  *K18Matrix_;
  std::string      FieldMapFileName_;
  std::string      bh1FilterFileName_;
  std::string      bh2FilterFileName_;



  // h_Utility
  std::string      FileManagerName_;

  double K18Momentum_;
  double SKSFieldNMR_;
  double SKSFieldCalc_;

  //Scaler
  std::string ScalerDefinitionFileName_;
  ScalerAna  *ScalerAnalyzer_;

  EvDisp *evDisp_;
  bool FlagEvDisp_;

  //MHTDC event suppression
  bool FlagMHTDC_;

  // raw KEY -> VALUE lines from the conf file (every 2-token
  // line, not just the whitelisted keys above), used to embed a
  // full copy of the run's parameter files / conditions into the
  // output ROOT file for reproducibility (cf. e72 ConfMan)
  std::map<std::string, std::string> m_paramLines;

public:
  ~ConfMan();
  ConfMan(const std::string &confFile);
  VEvent*         EventAllocator();
  //  static ConfMan& GetInstance();
  bool            Initialize();
  bool InitializeEvDisp( void );
  static ConfMan *GetConfManager( void ) { return confManager_; }

  //Hodo
  HodoParamMan *GetHodoParamManager( void ) { return HodoParamManager_;}
  HodoPHCMan   *GetHodoPHCManager( void ) { return HodoPHCManager_;}

  //DC
  DCGeomMan       *GetDCGeomManager( void ) { return DCGeomManager_; }
  DCTdcCalibMan   *GetDCTdcCalibManager( void ) { return DCTdcCalibManager_; }
  DCDriftParamMan *GetDCDriftParamManager( void ) { return DCDriftParamManager_; } 
  K18TransMatrix  *GetK18Matrix( void ) { return K18Matrix_; }

  const std::string& GetFileManagerName( void ) const
  { return FileManagerName_; }
  double K18Momentum( void ) const { return K18Momentum_; }
  double SKSFieldNMR( void ) const { return SKSFieldNMR_; }
  double SKSFieldCalc( void ) const { return SKSFieldCalc_; }

  bool GetEvDispFlag( void ) const { return FlagEvDisp_; }
  EvDisp *GetEvDisp( void ) { return evDisp_; }

  bool GetMHTDCFlag( void ) const { return FlagMHTDC_; }

private:
  ConfMan();
  ConfMan(const ConfMan&);
  ConfMan& operator=(const ConfMan&);

public:
  bool InitializeHistograms();
  bool InitializeParameterFiles();

  // Write a "param" TDirectory (one TMacro per referenced param
  // file, one TNamed per scalar conf value) plus a top-level
  // "git" TNamed into the CURRENTLY OPEN output ROOT file
  // (gDirectory/gFile), so every run's parameter conditions are
  // self-documented inside its own output file.
  void WriteParameters();
};

#endif
