/*
  K18Track_BFT.hh
*/

#ifndef K18Track_BFT_h
#define K18Track_BFT_h 1

#include "ThreeVector.hh"

#include <vector>
#include <functional>

class DCLocalTrack;
class TrackHit;
class DCAnalyzer;
class DCHit;

class K18Track_BFT
{
public: 
  K18Track_BFT( double LocalX, DCLocalTrack *tout, double P0 );
  ~K18Track_BFT(); 

private:
  K18Track_BFT();
  K18Track_BFT( const K18Track_BFT & );
  K18Track_BFT & operator = ( const K18Track_BFT & );

private:
  double LocalX_;
  DCLocalTrack *TrOut_;
  double P0_;
  double Xi_;
  double Xo_, Yo_, Uo_, Vo_;

  bool StatusD2U_;
  double P0D2U_;
  double DeltaD2U_;
  double DeltaD2U_3rd_;

  bool gfastatus_;

public:
  bool CalcMomentumD2U( void );
  bool StatusD2U( void ) const { return StatusD2U_; }

  DCLocalTrack *TrackOut( void ) { return TrOut_; }

  double Xin( void ) const { return Xi_; }
  double Xout( void ) const { return Xo_; }
  double Yout( void ) const { return Yo_; }
  double Uout( void ) const { return Uo_; }
  double Vout( void ) const { return Vo_; }

  double DeltaD2U( void ) const { return DeltaD2U_; }
  double DeltaD2U_3rd( void ) const { return DeltaD2U_3rd_; }
  double PD2U( void ) const { return P0_*(1.+DeltaD2U_); }
  double PD2U_3rd( void ) const { return P0_*(1.+DeltaD2U_3rd_); }

  ThreeVector BeamMomentumD2U( void ) const;
  double Xtgt( void ) const;
  double Ytgt( void ) const;
  double Utgt( void ) const;
  double Vtgt( void ) const;

  bool GoodForAnalysis( void ) const { return gfastatus_; }
  bool GoodForAnalysis( bool status )
  { bool ret=gfastatus_; gfastatus_=status; return ret; } 

private:
  void deleteHits( void );
  void addHits( void );
};


#endif
