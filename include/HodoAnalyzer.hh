// -*- C++ -*-

#ifndef HODO_ANALYZER_HH
#define HODO_ANALYZER_HH

#include <vector>

#include <TString.h>

#include "DetectorID.hh"
#include "RawData.hh"

class RawData;
class Hodo1Hit;
class Hodo2Hit;
class BH2Hit;
class FiberHit;
class FLHit;
class HodoCluster;
class BH2Cluster;
class FiberCluster;

typedef std::vector<Hodo1Hit*> Hodo1HitContainer;
typedef std::vector<Hodo2Hit*> Hodo2HitContainer;
typedef std::vector<BH2Hit*>   BH2HitContainer;
typedef std::vector<FiberHit*> FiberHitContainer;
typedef std::vector<FLHit*>    FLHitContainer;
typedef std::vector< std::vector<FiberHit*> >
MultiPlaneFiberHitContainer;

typedef std::vector<HodoCluster*>  HodoClusterContainer;
typedef std::vector<BH2Cluster*>   BH2ClusterContainer;
typedef std::vector<FiberCluster*> FiberClusterContainer;
typedef std::vector< std::vector<FiberCluster*> >
MultiPlaneFiberClusterContainer;

//_____________________________________________________________________________
class HodoAnalyzer
{
public:
  static TString       ClassName( void );
  static HodoAnalyzer& GetInstance( void );
  HodoAnalyzer( void );
  ~HodoAnalyzer( void );

private:
  HodoAnalyzer( const HodoAnalyzer& );
  HodoAnalyzer& operator =( const HodoAnalyzer& );

private:
  Hodo2HitContainer           m_BH1Cont;
  BH2HitContainer             m_BH2Cont;
  Hodo1HitContainer           m_BACCont;
  Hodo1HitContainer           m_PVACCont;
  Hodo1HitContainer           m_FACCont;
  Hodo2HitContainer           m_TOFCont;
  Hodo1HitContainer           m_LACCont;
  Hodo2HitContainer           m_WCCont;
  MultiPlaneFiberHitContainer m_BFTCont;
  FiberHitContainer           m_SCHCont;

  HodoClusterContainer            m_BH1ClCont;
  BH2ClusterContainer             m_BH2ClCont;
  HodoClusterContainer            m_BACClCont;
  HodoClusterContainer            m_PVACClCont;
  HodoClusterContainer            m_FACClCont;
  HodoClusterContainer            m_TOFClCont;
  HodoClusterContainer            m_LACClCont;
  HodoClusterContainer            m_WCClCont;
  FiberClusterContainer           m_BFTClCont;
  FiberClusterContainer           m_SCHClCont;

public:
  Bool_t DecodeRawHits( RawData* rawData );
  Bool_t DecodeBH1Hits( RawData* rawData );
  Bool_t DecodeBH2Hits( RawData* rawData );
  Bool_t DecodeBACHits( RawData* rawData );
  Bool_t DecodePVACHits( RawData* rawData );
  Bool_t DecodeFACHits( RawData* rawData );
  Bool_t DecodeTOFHits( RawData* rawData );
  Bool_t DecodeLACHits( RawData* rawData );
  Bool_t DecodeWCHits( RawData* rawData );
  Bool_t DecodeBFTHits( RawData* rawData );
  Bool_t DecodeSCHHits( RawData* rawData );
  Int_t  GetNHitsBH1( void ) const { return m_BH1Cont.size(); };
  Int_t  GetNHitsBH2( void ) const { return m_BH2Cont.size(); };
  Int_t  GetNHitsBAC( void ) const { return m_BACCont.size(); };
  Int_t  GetNHitsPVAC( void ) const { return m_PVACCont.size(); };
  Int_t  GetNHitsFAC( void ) const { return m_FACCont.size(); };
  Int_t  GetNHitsTOF( void ) const { return m_TOFCont.size(); };
  Int_t  GetNHitsLAC( void ) const { return m_LACCont.size(); };
  Int_t  GetNHitsWC( void ) const { return m_WCCont.size(); };
  Int_t  GetNHitsBFT( Int_t plane) const
  { return m_BFTCont.at( plane ).size(); };
  Int_t  GetNHitsSCH( void )  const { return m_SCHCont.size(); };

  inline Hodo2Hit* GetHitBH1( UInt_t i ) const;
  inline BH2Hit*   GetHitBH2( UInt_t i ) const;
  inline Hodo1Hit* GetHitBAC( UInt_t i ) const;
  inline Hodo1Hit* GetHitPVAC( UInt_t i ) const;
  inline Hodo1Hit* GetHitFAC( UInt_t i ) const;
  inline Hodo2Hit* GetHitTOF( UInt_t i ) const;
  inline Hodo1Hit* GetHitLAC( UInt_t i ) const;
  inline Hodo2Hit* GetHitWC( UInt_t i ) const;
  inline FiberHit* GetHitBFT( Int_t plane, UInt_t seg ) const;
  inline FiberHit* GetHitSCH( UInt_t seg ) const;

  Int_t GetNClustersBH1( void ) const { return m_BH1ClCont.size(); }
  Int_t GetNClustersBH2( void ) const { return m_BH2ClCont.size(); }
  Int_t GetNClustersBAC( void ) const { return m_BACClCont.size(); }
  Int_t GetNClustersPVAC( void ) const { return m_PVACClCont.size(); }
  Int_t GetNClustersFAC( void ) const { return m_FACClCont.size(); }
  Int_t GetNClustersTOF( void ) const { return m_TOFClCont.size(); }
  Int_t GetNClustersLAC( void ) const { return m_LACClCont.size(); }
  Int_t GetNClustersWC( void ) const { return m_WCClCont.size(); }
  Int_t GetNClustersBFT( void ) const { return m_BFTClCont.size(); };
  Int_t GetNClustersSCH( void ) const { return m_SCHClCont.size(); };
  inline HodoCluster*  GetClusterBH1( UInt_t i ) const;
  inline BH2Cluster*   GetClusterBH2( UInt_t i ) const;
  inline HodoCluster*  GetClusterBAC( UInt_t i ) const;
  inline HodoCluster*  GetClusterPVAC( UInt_t i ) const;
  inline HodoCluster*  GetClusterFAC( UInt_t i ) const;
  inline HodoCluster*  GetClusterTOF( UInt_t i ) const;
  inline HodoCluster*  GetClusterLAC( UInt_t i ) const;
  inline HodoCluster*  GetClusterWC( UInt_t i ) const;
  inline FiberCluster* GetClusterBFT( UInt_t i ) const;
  inline FiberCluster* GetClusterSCH( UInt_t i ) const;

  Bool_t ReCalcBH1Hits( Bool_t applyRecursively=false );
  Bool_t ReCalcBH2Hits( Bool_t applyRecursively=false );
  Bool_t ReCalcBACHits( Bool_t applyRecursively=false );
  Bool_t ReCalcPVACHits( Bool_t applyRecursively=false );
  Bool_t ReCalcFACHits( Bool_t applyRecursively=false );
  Bool_t ReCalcTOFHits( Bool_t applyRecursively=false );
  Bool_t ReCalcLACHits( Bool_t applyRecursively=false );
  Bool_t ReCalcWCHits( Bool_t applyRecursively=false );
  Bool_t ReCalcBH1Clusters( Bool_t applyRecursively=false );
  Bool_t ReCalcBH2Clusters( Bool_t applyRecursively=false );
  Bool_t ReCalcBACClusters( Bool_t applyRecursively=false );
  Bool_t ReCalcPVACClusters( Bool_t applyRecursively=false );
  Bool_t ReCalcFACClusters( Bool_t applyRecursively=false );
  Bool_t ReCalcTOFClusters( Bool_t applyRecursively=false );
  Bool_t ReCalcLACClusters( Bool_t applyRecursively=false );
  Bool_t ReCalcWCClusters( Bool_t applyRecursively=false );
  Bool_t ReCalcAll( void );

  void TimeCutBH1( Double_t tmin, Double_t tmax );
  void TimeCutBH2( Double_t tmin, Double_t tmax );
  void TimeCutTOF( Double_t tmin, Double_t tmax );
  void TimeCutBFT( Double_t tmin, Double_t tmax );
  void TimeCutSCH( Double_t tmin, Double_t tmax );

  void WidthCutBFT( Double_t min_width, Double_t max_width );
  void WidthCutSCH( Double_t min_width, Double_t max_width );
  BH2Cluster*  GetTime0BH2Cluster( void );
  HodoCluster* GetBtof0BH1Cluster( Double_t time0 );

private:
  void ClearBH1Hits( void );
  void ClearBH2Hits( void );
  void ClearBACHits( void );
  void ClearPVACHits( void );
  void ClearFACHits( void );
  void ClearTOFHits( void );
  void ClearLACHits( void );
  void ClearWCHits( void );
  void ClearBFTHits( void );
  void ClearSCHHits( void );

  template<typename TypeCluster>
  void TimeCut( std::vector<TypeCluster>& cont, Double_t tmin, Double_t tmax );
  template<typename TypeCluster>
  void WidthCut( std::vector<TypeCluster>& cont,
		 Double_t min_width, Double_t max_width, Bool_t adopt_nan );
  template<typename TypeCluster>
  void WidthCutR( std::vector<TypeCluster>& cont,
		  Double_t min_width, Double_t max_width, Bool_t adopt_nan );
  template<typename TypeCluster>
  void AdcCut( std::vector<TypeCluster>& cont, Double_t amin, Double_t amax );
  static Int_t MakeUpClusters( const Hodo1HitContainer& HitCont,
			       HodoClusterContainer& ClusterCont,
			       Double_t maxTimeDif );
  static Int_t MakeUpClusters( const Hodo2HitContainer& HitCont,
			       HodoClusterContainer& ClusterCont,
			       Double_t maxTimeDif );
  static Int_t MakeUpClusters( const BH2HitContainer& HitCont,
			       BH2ClusterContainer& ClusterCont,
			       Double_t maxTimeDif );
  static Int_t MakeUpClusters( const FiberHitContainer& cont,
			       FiberClusterContainer& ClusterCont,
			       Double_t maxTimeDif,
			       Int_t DifPairId );
  static Int_t MakeUpClusters( const FLHitContainer& cont,
			       FiberClusterContainer& ClusterCont,
			       Double_t maxTimeDif,
			       Int_t DifPairId );
  static Int_t MakeUpCoincidence( const FiberHitContainer& cont,
				  FLHitContainer& CoinCont,
				  Double_t maxTimeDif);
};

//_____________________________________________________________________________
inline TString
HodoAnalyzer::ClassName( void )
{
  static TString s_name("HodoAnalyzer");
  return s_name;
}

//_____________________________________________________________________________
inline HodoAnalyzer&
HodoAnalyzer::GetInstance( void )
{
  static HodoAnalyzer g_instance;
  return g_instance;
}

//_____________________________________________________________________________
inline HodoCluster*
HodoAnalyzer::GetClusterBH1( UInt_t i ) const
{
  if( i<m_BH1ClCont.size() )
    return m_BH1ClCont[i];
  else
    return nullptr;
}

//_____________________________________________________________________________
inline BH2Cluster*
HodoAnalyzer::GetClusterBH2( UInt_t i ) const
{
  if( i<m_BH2ClCont.size() )
    return m_BH2ClCont[i];
  else
    return nullptr;
}

//_____________________________________________________________________________
inline HodoCluster*
HodoAnalyzer::GetClusterBAC( UInt_t i ) const
{
  if( i<m_BACClCont.size() )
    return m_BACClCont[i];
  else
    return nullptr;
}

//_____________________________________________________________________________
inline HodoCluster*
HodoAnalyzer::GetClusterPVAC( UInt_t i ) const
{
  if( i<m_PVACClCont.size() )
    return m_PVACClCont[i];
  else
    return nullptr;
}

//_____________________________________________________________________________
inline HodoCluster*
HodoAnalyzer::GetClusterFAC( UInt_t i ) const
{
  if( i<m_FACClCont.size() )
    return m_FACClCont[i];
  else
    return nullptr;
}

//_____________________________________________________________________________
inline HodoCluster*
HodoAnalyzer::GetClusterTOF( UInt_t i ) const
{
  if( i<m_TOFClCont.size() )
    return m_TOFClCont[i];
  else
    return nullptr;
}

//_____________________________________________________________________________
inline HodoCluster*
HodoAnalyzer::GetClusterLAC( UInt_t i ) const
{
  if( i<m_LACClCont.size() )
    return m_LACClCont[i];
  else
    return nullptr;
}

//_____________________________________________________________________________
inline HodoCluster*
HodoAnalyzer::GetClusterWC( UInt_t i ) const
{
  if( i<m_WCClCont.size() )
    return m_WCClCont[i];
  else
    return nullptr;
}

//_____________________________________________________________________________
inline FiberCluster*
HodoAnalyzer::GetClusterBFT( UInt_t i ) const
{
  if( i<m_BFTClCont.size() )
    return m_BFTClCont[i];
  else
    return nullptr;
}


//_____________________________________________________________________________
inline FiberCluster*
HodoAnalyzer::GetClusterSCH( UInt_t i ) const
{
  if( i<m_SCHClCont.size() )
    return m_SCHClCont[i];
  else
    return nullptr;
}

//_____________________________________________________________________________
inline Hodo2Hit*
HodoAnalyzer::GetHitBH1( UInt_t i ) const
{
  if( i<m_BH1Cont.size() )
    return m_BH1Cont[i];
  else
    return nullptr;
}

//_____________________________________________________________________________
inline BH2Hit*
HodoAnalyzer::GetHitBH2( UInt_t i ) const
{
  if( i<m_BH2Cont.size() )
    return m_BH2Cont[i];
  else
    return nullptr;
}

//_____________________________________________________________________________
inline Hodo1Hit*
HodoAnalyzer::GetHitBAC( UInt_t i ) const
{
  if( i<m_BACCont.size() )
    return m_BACCont[i];
  else
    return nullptr;
}

//_____________________________________________________________________________
inline Hodo1Hit*
HodoAnalyzer::GetHitPVAC( UInt_t i ) const
{
  if( i<m_PVACCont.size() )
    return m_PVACCont[i];
  else
    return nullptr;
}

//_____________________________________________________________________________
inline Hodo1Hit*
HodoAnalyzer::GetHitFAC( UInt_t i ) const
{
  if( i<m_FACCont.size() )
    return m_FACCont[i];
  else
    return nullptr;
}

//_____________________________________________________________________________
inline Hodo2Hit*
HodoAnalyzer::GetHitTOF( UInt_t i ) const
{
  if( i<m_TOFCont.size() )
    return m_TOFCont[i];
  else
    return nullptr;
}

//_____________________________________________________________________________
inline Hodo1Hit*
HodoAnalyzer::GetHitLAC( UInt_t i ) const
{
  if( i<m_LACCont.size() )
    return m_LACCont[i];
  else
    return nullptr;
}

//_____________________________________________________________________________
inline Hodo2Hit*
HodoAnalyzer::GetHitWC( UInt_t i ) const
{
  if( i<m_WCCont.size() )
    return m_WCCont[i];
  else
    return nullptr;
}

//_____________________________________________________________________________
inline FiberHit*
HodoAnalyzer::GetHitBFT( Int_t plane, UInt_t seg ) const
{
  if( seg<m_BFTCont.at(plane).size() )
    return m_BFTCont.at(plane).at(seg);
  else
    return nullptr;
}


//_____________________________________________________________________________
inline FiberHit*
HodoAnalyzer::GetHitSCH( UInt_t seg ) const
{
  if( seg<m_SCHCont.size() )
    return m_SCHCont[seg];
  else
    return nullptr;
}

#endif
