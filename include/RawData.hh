// -*- C++ -*-

#ifndef RAW_DATA_HH
#define RAW_DATA_HH

#include <TString.h>

#include <vector>

#include "DetectorID.hh"

class HodoRawHit;
class DCRawHit;

typedef std::vector<HodoRawHit*> HodoRHitContainer;
typedef std::vector<DCRawHit*>   DCRHitContainer;
typedef std::vector<Int_t>       FADCRHitContainer;

//_____________________________________________________________________________
class RawData
{
public:
  static TString ClassName( void );
  RawData( void );
  ~RawData( void );

private:
  RawData( const RawData& );
  RawData& operator=( const RawData& );

private:
  Bool_t            	  	 m_is_decoded;
  HodoRHitContainer              m_BH1RawHC;
  HodoRHitContainer              m_BH2RawHC;
  HodoRHitContainer              m_BACRawHC;
  HodoRHitContainer              m_PVACRawHC;
  HodoRHitContainer              m_FACRawHC;
  HodoRHitContainer              m_TOFRawHC;
  HodoRHitContainer              m_LACRawHC;
  HodoRHitContainer              m_WCRawHC;
  HodoRHitContainer              m_WCSUMRawHC;
  std::vector<HodoRHitContainer> m_BFTRawHC;
  HodoRHitContainer              m_SCHRawHC;
  std::vector<DCRHitContainer>   m_BcInRawHC;
  std::vector<DCRHitContainer>   m_BcOutRawHC;
  std::vector<DCRHitContainer>   m_SdcInRawHC;
  std::vector<DCRHitContainer>   m_SdcOutRawHC;
  HodoRHitContainer              m_ScalerRawHC;
  HodoRHitContainer              m_TrigRawHC;
  HodoRHitContainer              m_VmeCalibRawHC;

public:
  void                     ClearAll( void );
  Bool_t                   DecodeHits( void );
  Bool_t                   DecodeCalibHits( void );
  const HodoRHitContainer& GetBH1RawHC( void ) const;
  const HodoRHitContainer& GetBH2RawHC( void ) const;
  const HodoRHitContainer& GetBACRawHC( void ) const;
  const HodoRHitContainer& GetPVACRawHC( void ) const;
  const HodoRHitContainer& GetFACRawHC( void ) const;
  const HodoRHitContainer& GetTOFRawHC( void ) const;
  const HodoRHitContainer& GetLACRawHC( void ) const;
  const HodoRHitContainer& GetWCRawHC( void ) const;
  const HodoRHitContainer& GetWCSUMRawHC( void ) const;
  const HodoRHitContainer& GetBFTRawHC( Int_t plane ) const;
  const HodoRHitContainer& GetSCHRawHC( void ) const;
  const DCRHitContainer&   GetBcInRawHC( Int_t layer ) const;
  const DCRHitContainer&   GetBcOutRawHC( Int_t layer ) const;
  const DCRHitContainer&   GetSdcInRawHC( Int_t layer ) const;
  const DCRHitContainer&   GetSdcOutRawHC( Int_t layer ) const;
  const HodoRHitContainer& GetScalerRawHC( void ) const;
  const HodoRHitContainer& GetTrigRawHC( void ) const;
  const HodoRHitContainer& GetVmeCalibRawHC( void ) const;

private:
  enum EDCDataType   { kDcLeading, kDcTrailing, kDcOverflow, kDcNDataType };
  Bool_t AddHodoRawHit( HodoRHitContainer& cont,
			Int_t id, Int_t plane, Int_t seg,
			Int_t UorD, Int_t type, Int_t data );
  Bool_t AddDCRawHit( DCRHitContainer& cont,
		      Int_t plane, Int_t wire, Int_t data,
		      Int_t type=kDcLeading );
  void   DecodeHodo( Int_t id, Int_t plane, Int_t nseg, Int_t nch,
		     HodoRHitContainer& cont );
  void   DecodeHodo( Int_t id, Int_t nseg, Int_t nch,
		     HodoRHitContainer& cont );
};

//_____________________________________________________________________________
inline TString
RawData::ClassName( void )
{
  static TString s_name("RawData");
  return s_name;
}

//_____________________________________________________________________________
inline const HodoRHitContainer&
RawData::GetBH1RawHC( void ) const
{
  return m_BH1RawHC;
}

//_____________________________________________________________________________
inline const HodoRHitContainer&
RawData::GetBH2RawHC( void ) const
{
  return m_BH2RawHC;
}

//_____________________________________________________________________________
inline const HodoRHitContainer&
RawData::GetBACRawHC( void ) const
{
  return m_BACRawHC;
}

//_____________________________________________________________________________
inline const HodoRHitContainer&
RawData::GetPVACRawHC( void ) const
{
  return m_PVACRawHC;
}
//_____________________________________________________________________________
inline const HodoRHitContainer&
RawData::GetFACRawHC( void ) const
{
  return m_FACRawHC;
}

//_____________________________________________________________________________
inline const HodoRHitContainer&
RawData::GetTOFRawHC( void ) const
{
  return m_TOFRawHC;
}

//_____________________________________________________________________________
inline const HodoRHitContainer&
RawData::GetLACRawHC( void ) const
{
  return m_LACRawHC;
}

//_____________________________________________________________________________
inline const HodoRHitContainer&
RawData::GetWCRawHC( void ) const
{
  return m_WCRawHC;
}

//_____________________________________________________________________________
inline const HodoRHitContainer&
RawData::GetWCSUMRawHC( void ) const
{
  return m_WCSUMRawHC;
}


//_____________________________________________________________________________
inline const HodoRHitContainer&
RawData::GetBFTRawHC( Int_t plane ) const
{
  if( plane<0 || plane>NumOfPlaneBFT-1 ) plane=0;
  return m_BFTRawHC[plane];
}

//_____________________________________________________________________________
inline const HodoRHitContainer&
RawData::GetSCHRawHC( void ) const
{
  return m_SCHRawHC;
}

//_____________________________________________________________________________
inline const DCRHitContainer&
RawData::GetBcInRawHC( Int_t layer ) const
{
  if( layer<0 || layer>NumOfLayersBcIn ) layer = 0;
  return m_BcInRawHC[layer];
}

//_____________________________________________________________________________
inline const DCRHitContainer&
RawData::GetBcOutRawHC( Int_t layer ) const
{
  if( layer<0 || layer>NumOfLayersBcOut ) layer = 0;
  return m_BcOutRawHC[layer];
}

//_____________________________________________________________________________
inline const DCRHitContainer&
RawData::GetSdcInRawHC( Int_t layer ) const
{
  if( layer<0 || layer>NumOfLayersSdcIn ) layer = 0;
  return m_SdcInRawHC[layer];
}

//_____________________________________________________________________________
inline const DCRHitContainer&
RawData::GetSdcOutRawHC( Int_t layer ) const
{
  if( layer<0 || layer>NumOfLayersSdcOut ) layer = 0;
  return m_SdcOutRawHC[layer];
}

//_____________________________________________________________________________
inline const HodoRHitContainer&
RawData::GetScalerRawHC( void ) const
{
  return m_ScalerRawHC;
}

//_____________________________________________________________________________
inline const HodoRHitContainer&
RawData::GetTrigRawHC( void ) const
{
  return m_TrigRawHC;
}

//_____________________________________________________________________________
inline const HodoRHitContainer&
RawData::GetVmeCalibRawHC( void ) const
{
  return m_VmeCalibRawHC;
}

#endif
