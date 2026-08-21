/*
  Hodo1Hit.hh

  2012/1/24
*/

#ifndef Hodo1Hit_h
#define Hodo1Hit_h 1

#include "HodoRawHit.hh"

class RawData;

class Hodo1Hit
{
public:
  explicit Hodo1Hit( HodoRawHit *rhit );
  virtual ~Hodo1Hit();

private:
  Hodo1Hit( const Hodo1Hit & );
  Hodo1Hit & operator = ( const Hodo1Hit & );

protected:
  HodoRawHit *raw_;
  bool Status_;
  int multi_hit;
  std::vector<double> a_, t_, ct_;
  //  double a_, t_, ct_;

public:
  HodoRawHit * GetRawHit( void ) { return raw_; }
  bool calculate( void );

  bool status( void ) const { return Status_; } 

  int  GetNumOfHit(void) const {return multi_hit;};
  
  double GetA( void )  const { return a_.at(0); }
  double GetA( int n )  const { return a_.at(n); }
  double GetT( void )  const { return t_.at(0); }
  double GetT( int n )  const { return t_.at(n); }
  double GetCT( void ) const { return ct_.at(0); }
  double GetCT( int n )  const { return ct_.at(n); }

  double Time( void ) const { return GetT(); }
  double Time( int n ) const { return GetT(n); }
  double CTime( void ) const { return GetCT(); }
  double CTime( int n ) const { return GetCT(n); }
  double DeltaE( void ) const { return GetA(); }
  double DeltaE( int n ) const { return GetA(n); }

  int DetectorId( void ) const { return raw_->DetectorId(); }
  int PlaneId( void ) const { return raw_->PlaneId(); }
  int SegmentId( void ) const { return raw_->SegmentId(); }

  virtual bool ReCalc( bool applyRecursively=false ) 
  { return calculate(); }

};

#endif
