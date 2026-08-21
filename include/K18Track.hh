/*
  K18Track.hh

  2012/1/24
*/

#ifndef K18Track_h
#define K18Track_h 1

#include "ThreeVector.hh"

#include <vector>
#include <functional>

class DCLocalTrack;
class TrackHit;
class DCAnalyzer;

class K18Track
{
public: 
  K18Track( DCLocalTrack *tin, DCLocalTrack *tout, double P0 );
  ~K18Track(); 

private:
  K18Track( const K18Track & );
  K18Track & operator = ( const K18Track & );

private:
  DCLocalTrack *TrIn_, *TrOut_;
  double P0_;
  bool Status_;
  double chisqr_;
  double Xi_, Yi_, Ui_, Vi_;
  double Xo_, Yo_, Uo_, Vo_;
  double Delta_; 
  bool gfastatus_;

  std::vector <TrackHit *> hitContIn, hitContOut;
  
public:
  bool doFit( void );
  bool Status( void ) const { return Status_; }

  DCLocalTrack *TrackIn( void )  { return TrIn_; }
  DCLocalTrack *TrackOut( void ) { return TrOut_; }

  double Xin( void ) const { return Xi_; }
  double Yin( void ) const { return Yi_; }
  double Uin( void ) const { return Ui_; }
  double Vin( void ) const { return Vi_; }
  double Xout( void ) const { return Xo_; }
  double Yout( void ) const { return Yo_; }
  double Uout( void ) const { return Uo_; }
  double Vout( void ) const { return Vo_; }
  double Delta( void ) const { return Delta_; }
  double P( void ) const { return P0_*(1.+Delta_); }
  double chisquare() const { return chisqr_; }

  ThreeVector BeamMomentum( void ) const;
  double Xtgt( void ) const;
  double Ytgt( void ) const;
  double Utgt( void ) const;
  double Vtgt( void ) const;

  int GetNHitsIn( void )  const { return hitContIn.size(); }
  int GetNHitsOut( void ) const { return hitContOut.size(); }
  int GetNHitsTotal( void ) const
  { return hitContIn.size()+hitContOut.size(); }

  TrackHit * GetK18HitIn( int i );
  TrackHit * GetK18HitOut( int i );
  TrackHit * GetK18HitTotal( int i );
  TrackHit * GetK18HitByPlaneId( int PlId );

  bool GoodForAnalysis( void ) const { return gfastatus_; }
  bool GoodForAnalysis( bool status )
  { bool ret=gfastatus_; gfastatus_=status; return ret; } 

  bool ReCalc( bool ApplyRecursively=false );

private:
  void deleteHits( void );
  void addHits( void );
};

struct K18TrackComp
  : public std::binary_function <K18Track *, K18Track *,bool>
{
  bool operator()( const K18Track * const p1,
		   const K18Track * const p2 ) const
  {
    return (p1->chisquare())<(p2->chisquare());
  }

};

#endif
