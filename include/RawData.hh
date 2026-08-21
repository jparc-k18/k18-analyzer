/*
  RawData.hh

  2012/1/24
*/

#ifndef RawData_h
#define RawData_h

#include "DetectorID.hh"
#include <vector>

#ifdef MemoryLeak
#include "DebugCounter.hh"
#endif

class HodoRawHit;
class DCRawHit;

typedef std::vector<HodoRawHit*> HodoRHitContainer;
typedef std::vector <DCRawHit*> DCRHitContainer;

class RawData
{

private:
  HodoRHitContainer GCRawHC;
  HodoRHitContainer BH1RawHC;
  HodoRHitContainer BH2RawHC;
  HodoRHitContainer BFTRawHC[NumOfPlaneBFT];
  //HodoRHitContainer SFTRawHC[NumOfPlaneSFT];

  HodoRHitContainer TOFRawHC;
  HodoRHitContainer ACRawHC;
  HodoRHitContainer AC1RawHC;
  HodoRHitContainer LCRawHC;

  // E13 counters
  HodoRHitContainer BMWRawHC;
  HodoRHitContainer BACRawHC;
  HodoRHitContainer SFVRawHC;
  HodoRHitContainer TOFMTRawHC;

  DCRHitContainer BcInRawHC[NumOfLayersBcIn+1];
  DCRHitContainer BcOutRawHC[NumOfLayersBcOut+1];
  DCRHitContainer SdcInRawHC[NumOfLayersSdcIn+1];
  DCRHitContainer SdcOutRawHC[NumOfLayersSdcOut+1];
  DCRHitContainer SsdRawHC[NumOfLayersSsd+1];

  HodoRHitContainer ScalerRawHC;
  HodoRHitContainer MiscRawHC;
  HodoRHitContainer MatrixRawHC;

#ifdef MemoryLeak
  static debug::Counter sm_counter;
#endif

public:
  RawData();
  ~RawData();

  void clearAll();
  bool DecodeHits();

private:
  RawData(const RawData&);
  RawData& operator=(const RawData&);

  bool AddHodoRawHit(HodoRHitContainer& cont,
 		     int DetId,
 		     int Plane,
 		     int Seg,
 		     int AorT,
 		     int UorD,
 		     int Data);

  bool AddDCRawHit( DCRHitContainer &cont, 
		    int Plane, 
		    int Wire, 
		    int Tdc,
		    int type=0); 

public:
  //const HodoRHitContainer& GetGCRawHC() const;
  const HodoRHitContainer& GetBH1RawHC() const;
  const HodoRHitContainer& GetBH2RawHC() const;
  const HodoRHitContainer & GetBFTRawHC( int plane ) const;
  //const HodoRHitContainer & GetSFTRawHC( int plane ) const;

  const HodoRHitContainer& GetTOFRawHC() const;
  //const HodoRHitContainer& GetACRawHC() const;
  //const HodoRHitContainer& GetAC1RawHC() const;
  const HodoRHitContainer& GetLCRawHC() const;
  const HodoRHitContainer& GetACRawHC() const;

  // E13 counters
  const HodoRHitContainer& GetBMWRawHC() const;
  const HodoRHitContainer& GetBACRawHC() const;
  const HodoRHitContainer& GetSFVRawHC() const;
  const HodoRHitContainer& GetTOFMTRawHC() const;

  const DCRHitContainer & GetBcInRawHC( int layer ) const;
  const DCRHitContainer & GetBcOutRawHC( int layer ) const;
  const DCRHitContainer & GetSsdRawHC( int layer ) const;
  const DCRHitContainer & GetSdcInRawHC( int layer ) const;
  const DCRHitContainer & GetSdcOutRawHC( int layer ) const;

  const HodoRHitContainer& GetScalerRawHC() const;
  const HodoRHitContainer& GetMiscRawHC() const;
  const HodoRHitContainer& GetMatrixRawHC() const;
};

#endif
