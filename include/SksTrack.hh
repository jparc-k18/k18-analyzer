/*
  SksTrack.hh

  2012/1/24
*/

#ifndef SksTrack_h
#define SksTrack_h 1

#include "RungeKuttaUtilities.hh"
#include "ThreeVector.hh"

#include <vector>
#include <iosfwd>
#include <functional>

class DCLocalTrack;
class TrackHit;
class DCAnalyzer;

class SksTrack
{
public:
  SksTrack( DCLocalTrack *In, DCLocalTrack *Out );
  ~SksTrack();

private:
  SksTrack( const SksTrack & );
  SksTrack & operator= ( const SksTrack & );

public:
  DCLocalTrack *GetLocalTrackIn( void ) { return In_;}
  DCLocalTrack *GetLocalTrackOut( void ) { return Out_; }
  bool doFit( void );
  bool Status( void ) const { return status_; }
  int Niteration( void ) const { return nIteration_; }
  void SetInitialMomentum( double Pini ) { iniP_=Pini; }

  const ThreeVector &
  PrimaryPosition( void ) const { return priPos_; }
  const ThreeVector &
  PrimaryMomentum( void ) const { return priMom_; }
  double PathLengthToTOF( void ) const { return pLenTtoT_; }
  double PathLengthTotal( void ) const { return pLenTot_; }
  double chisqr( void ) const { return chisqr_; }

  std::size_t GetNHits( void ) const { return hitArray_.size(); } 
  TrackHit * GetHit( std::size_t nth ) const;
  TrackHit * GetHitOfLayerNumber( int lnum ) const;
  bool GoodForAnalysis( void ) const { return gfastatus_; }
  bool GoodForAnalysis( bool status )
  { bool ret=gfastatus_; gfastatus_=status; return ret; } 

  bool GetTrajectoryLocalPosition( int layer, double & x, double & y ) const;

  bool ReCalc( bool applyRecursively=false );

private:
  DCLocalTrack *In_;
  DCLocalTrack *Out_;
  double iniP_;


  std::vector <TrackHit *> hitArray_;
  RKHitPointContainer hPntCont_;

  bool status_;
  int nIteration_;
  double chisqr_;

  ThreeVector priPos_, priMom_;
  double pLenTtoT_, pLenTot_;
  RKCordParameter CPval_;

  bool gfastatus_;

private:
  void fillHitArray( void );
  void clearHitArray( void );
  double calcChiSqr( const RKHitPointContainer &hpCont ) const;
  bool guessNextParameters( const RKHitPointContainer &hpCont,
			    RKCordParameter &Cord,
			    double &estDeltaChisqr,
			    double &lambdaCri, double dmp=0.0 ) const;

  void saveCalcPosition( const RKHitPointContainer &hpCont );
  void PrintCalcHits( const RKHitPointContainer &hpCont,
		      std::ostream &ost ) const;

  bool saveTrackParameters( const RKCordParameter &cp );  
  bool doFitAgain( RKCordParameter iniCord );
};

struct SksTrackComp
  : public std::binary_function <SksTrack *, SksTrack *, bool>
{
  bool operator()( const SksTrack * const p1,
		   const SksTrack * const p2 ) const
  {
    int n1=p1->GetNHits(), n2=p2->GetNHits();
    if( n1>n2+4 ) return true;
    else if( n2>n1+4 ) return false;
    else
      return (p1->chisqr())<(p2->chisqr());
  }
};

#endif
