/*
  HodoAnalyzer.hh

  2012/1/24
*/

#ifndef HodoAnalyzer_h
#define HodoAnalyzer_h

#include <cstddef>
#include <vector>

#include "DetectorID.hh"
#include "RawData.hh"

class RawData;
class Hodo2Hit;
class Hodo1Hit;
class BH2Hit;
class FiberHit;
class HodoCluster;
class BH2Cluster;
class FiberCluster;

typedef std::vector <Hodo2Hit*> Hodo2HitContainer;
typedef std::vector <Hodo1Hit*> Hodo1HitContainer;
typedef std::vector <BH2Hit*>   BH2HitContainer;
typedef std::vector <FiberHit*> FiberHitContainer;

typedef std::vector <HodoCluster*>  HodoClusterContainer;
typedef std::vector <BH2Cluster*>   BH2ClusterContainer;
typedef std::vector <FiberCluster*> FiberClusterContainer;

class HodoAnalyzer
{
public:
  HodoAnalyzer();
  ~HodoAnalyzer();

  static HodoAnalyzer& GetInstance();
private:
  HodoAnalyzer(const HodoAnalyzer &);
  HodoAnalyzer & operator =(const HodoAnalyzer &);

private:
  Hodo2HitContainer BH1Cont;
  BH2HitContainer   BH2Cont;
  FiberHitContainer BFTCont;
  //FiberHitContainer SFTCont[NumOfLayersSFT];
  //Hodo1HitContainer TgtCont;
  Hodo2HitContainer TOFCont;
  //Hodo1HitContainer ACCont;
  Hodo2HitContainer LCCont;

  // E13 counters
  Hodo1HitContainer BMWCont;
  Hodo1HitContainer BACCont;
  Hodo1HitContainer SFVCont;
  Hodo1HitContainer TOFMTCont;

  HodoClusterContainer  BH1ClCont;
  BH2ClusterContainer   BH2ClCont;
  FiberClusterContainer BFTClCont;
  //FiberClusterContainer SFTClCont[NumOfLayersSFT];
  HodoClusterContainer  TOFClCont;
  HodoClusterContainer  LCClCont;

public:
  bool DecodeRawHits(RawData* rawData);

  //bool DecodeGCHits(RawData* rawData);  
  bool DecodeBH1Hits(RawData* rawData);
  bool DecodeBH2Hits(RawData* rawData);
  bool DecodeBFTHits(RawData* rawData);
  //bool DecodeSFTHits(RawData* rawData);
  //bool DecodeTGTHits(RawData* rawData);
  bool DecodeTOFHits(RawData* rawData);
  //bool DecodeACHits(RawData* rawData);
  bool DecodeLCHits(RawData* rawData);

  // E13 counters
  bool DecodeBMWHits(RawData* rawData);
  bool DecodeBACHits(RawData* rawData);
  bool DecodeSFVHits(RawData* rawData);
  bool DecodeTOFMTHits(RawData* rawData);

  //int GetNHitsGC( void ) const { return GCCont.size(); };
  int GetNHitsBH1( void ) const { return BH1Cont.size(); };
  int GetNHitsBH2( void ) const { return BH2Cont.size(); };
  int GetNHitsBFT( void ) const { return BFTCont.size(); };
  //int GetNHitsSFT( int layer ) const { return SFTCont[layer].size(); };
  int GetNHitsTOF( void ) const { return TOFCont.size(); };
  //int GetNHitsAC( void ) const { return ACCont.size(); };
  int GetNHitsLC( void ) const { return LCCont.size(); };

  // E13 counters
  int GetNHitsBMW( void ) const { return BMWCont.size(); };
  int GetNHitsBAC( void ) const { return BACCont.size(); };
  int GetNHitsSFV( void ) const { return SFVCont.size(); };
  int GetNHitsTOFMT( void ) const { return TOFMTCont.size(); };
  //inline Hodo1Hit * GetHitGC( int i ) const;
  inline Hodo2Hit * GetHitBH1( int i ) const;
  inline BH2Hit   * GetHitBH2( int i ) const;
  inline FiberHit * GetHitBFT(int seg) const;
  //inline FiberHit * GetHitSFT(int lyaer, int seg) const;
  inline Hodo2Hit * GetHitTOF( int i ) const;
  //inline Hodo1Hit * GetHitAC( int i ) const;
  inline Hodo2Hit * GetHitLC( int i ) const;

  // E13 counters
  inline Hodo1Hit * GetHitBMW( int i ) const;
  inline Hodo1Hit * GetHitBAC( int i ) const;
  inline Hodo1Hit * GetHitSFV( int i ) const;
  inline Hodo1Hit * GetHitTOFMT( int i ) const;

  int GetNClustersBH1( void ) const { return BH1ClCont.size(); };
  int GetNClustersBH2( void ) const { return BH2ClCont.size(); };
  int GetNClustersBFT( void ) const { return BFTClCont.size(); };
  //int GetNClustersSFT( int layer ) const { return SFTClCont[layer].size(); };
  int GetNClustersTOF( void ) const { return TOFClCont.size(); }
  int GetNClustersLC( void )  const { return LCClCont.size(); }

  inline HodoCluster  * GetClusterBH1( int i ) const;
  inline BH2Cluster   * GetClusterBH2( int i ) const;
  inline FiberCluster * GetClusterBFT( int i ) const;
  //inline FiberCluster * GetClusterSFT( int layer, int i ) const;
  inline HodoCluster  * GetClusterTOF( int i ) const;
  inline HodoCluster  * GetClusterLC( int i )  const;

  //bool ReCalcGCHits( bool applyRecursively=false );
  bool ReCalcBH1Hits( bool applyRecursively=false );
  bool ReCalcBH2Hits( bool applyRecursively=false );
  bool ReCalcTOFHits( bool applyRecursively=false );
  //bool ReCalcACHits( bool applyRecursively=false );
  bool ReCalcLCHits( bool applyRecursively=false );

  // E13 counters
  bool ReCalcBMWHits( bool applyRecursively=false );
  bool ReCalcBACHits( bool applyRecursively=false );
  bool ReCalcSFVHits( bool applyRecursively=false );
  bool ReCalcTOFMTHits( bool applyRecursively=false );

  bool ReCalcBH1Clusters( bool applyRecursively=false );
  bool ReCalcBH2Clusters( bool applyRecursively=false );
  bool ReCalcTOFClusters( bool applyRecursively=false );
  bool ReCalcLCClusters( bool applyRecursively=false );

  bool ReCalcAll( void );
  
  void TimeCutBH1(double tmin, double tmax);
  void TimeCutBH2(double tmin, double tmax);
  void TimeCutBFT(double tmin, double tmax);
  //void TimeCutSFT(int layer, double tmin, double tmax);

private:
  static HodoAnalyzer *MySelf_;

  //void clearGCHits();
  void clearBH1Hits();
  void clearBH2Hits();
  void clearBFTHits();
  //void clearSFTHits();
  void clearTOFHits( void );
  void clearACHits( void );
  void clearLCHits( void );

  // E13 counters
  void clearBMWHits();
  void clearBACHits();
  void clearSFVHits();
  void clearTOFMTHits();

  // bool fl_decoded_SFT;
  // void DecodeSFT(RawData* rawData,
  // 		 int Begin, int End,
  // 		 int layer,
  // 		 char* Name);

  template<typename TypeCluster>
  void TimeCut(std::vector<TypeCluster>& cont, double tmin, double tmax);

  static int MakeUpClusters( const Hodo2HitContainer & HitCont,
			     HodoClusterContainer & ClusterCont,
			     double maxTimeDif );

  static int MakeUpClusters( const BH2HitContainer & HitCont,
			     BH2ClusterContainer & ClusterCont,
			     double maxTimeDif );

  static int MakeUpClusters( const FiberHitContainer & cont,
			     FiberClusterContainer& ClusterCont,
			     double maxTimeDif);
};

inline HodoCluster * HodoAnalyzer::GetClusterBH1( int i ) const
{
  if( i>=0 && i<BH1ClCont.size() )
    return BH1ClCont[i];
  else
    return 0;
}

inline BH2Cluster * HodoAnalyzer::GetClusterBH2( int i ) const
{
  if( i>=0 && i<BH2ClCont.size() )
    return BH2ClCont[i];
  else
    return 0;
}

inline FiberCluster * HodoAnalyzer::GetClusterBFT( int i ) const
{
  if( i>=0 && i<BFTClCont.size() )
    return BFTClCont[i];
  else
    return 0;
}

// inline FiberCluster * HodoAnalyzer::GetClusterSFT( int layer, int i ) const
// {
//   if( i>=0 && i<SFTClCont[layer].size() )
//     return SFTClCont[layer][i];
//   else
//     return 0;
// }

inline HodoCluster * HodoAnalyzer::GetClusterTOF( int i ) const
{
  if( i>=0 && i<TOFClCont.size() )
    return TOFClCont[i];
  else
    return 0;
}

inline HodoCluster * HodoAnalyzer::GetClusterLC( int i ) const
{
  if( i>=0 && i<LCClCont.size() )
    return LCClCont[i];
  else
    return 0;
}

inline Hodo2Hit * HodoAnalyzer::GetHitBH1( int i ) const
{
  if( i>=0 && i<BH1Cont.size() )
    return BH1Cont[i];
  else
    return 0;
}

inline BH2Hit * HodoAnalyzer::GetHitBH2( int i ) const
{
  if( i>=0 && i<BH2Cont.size() )
    return BH2Cont[i];
  else
    return 0;
}

inline FiberHit* HodoAnalyzer::GetHitBFT(int seg) const
{
  if( seg>=0 && seg<BFTCont.size() )
    return BFTCont[seg];
  else
    return NULL;
}

// inline FiberHit* HodoAnalyzer::GetHitSFT(int layer, int seg) const
// {
//   if( seg>=0 && seg<SFTCont[layer].size() )
//     return SFTCont[layer][seg];
//   else
//     return NULL;
// }

inline Hodo2Hit * HodoAnalyzer::GetHitTOF( int i ) const
{
  if( i>=0 && i<TOFCont.size() )
    return TOFCont[i];
  else
    return 0;
}

// inline Hodo1Hit * HodoAnalyzer::GetHitAC( int i ) const
// {
//   if( i>=0 && i<ACCont.size() )
//     return ACCont[i];
//   else
//     return 0;
// }

inline Hodo2Hit * HodoAnalyzer::GetHitLC( int i ) const
{
  if( i>=0 && i<LCCont.size() )
    return LCCont[i];
  else
    return 0;
}

// E13 counters
inline Hodo1Hit * HodoAnalyzer::GetHitBMW( int i ) const
{
  if( i>=0 && i<BMWCont.size() )
    return BMWCont[i];
  else
    return 0;
}
inline Hodo1Hit * HodoAnalyzer::GetHitBAC( int i ) const
{
  if( i>=0 && i<BACCont.size() )
    return BACCont[i];
  else
    return 0;
}
inline Hodo1Hit * HodoAnalyzer::GetHitSFV( int i ) const
{
  if( i>=0 && i<SFVCont.size() )
    return SFVCont[i];
  else
    return 0;
}
inline Hodo1Hit * HodoAnalyzer::GetHitTOFMT( int i ) const
{
  if( i>=0 && i<TOFMTCont.size() )
    return TOFMTCont[i];
  else
    return 0;
}



inline HodoAnalyzer& HodoAnalyzer::GetInstance(){
  return *MySelf_;
}

#endif
