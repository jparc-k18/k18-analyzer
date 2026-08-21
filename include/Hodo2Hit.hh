/*
  Hodo2Hit.hh

  2012/1/24
*/

#ifndef HODO2_HIT_H 
#define HODO2_HIT_H
 
#include "HodoRawHit.hh"
#include "ThreeVector.hh"
#include <cstddef>

class RawData;

class Hodo2Hit
{
public:
  explicit Hodo2Hit( HodoRawHit *rhit, int index = 0);
  virtual ~Hodo2Hit();

private:
  Hodo2Hit( const Hodo2Hit & );
  Hodo2Hit & operator = ( const Hodo2Hit & );

protected:
  HodoRawHit *raw_;
  bool Status_;
  double a1_, a2_, t1_, t2_, ct1_, ct2_;
  int   index_;

public:
  HodoRawHit * GetRawHit( void ) { return raw_; }
  bool calculate( void );

  bool status( void ) const { return Status_; } 

  double GetAUp( void )    const { return a1_; }
  double GetALeft( void )  const { return a1_; }
  double GetADown( void )  const { return a2_; }
  double GetARight( void ) const { return a2_; }

  double GetTUp( void )    const { return t1_; }
  double GetTLeft( void )  const { return t1_; }
  double GetTDown( void )  const { return t2_; }
  double GetTRight( void ) const { return t2_; }

  double GetCTUp( void )    const { return ct1_; }
  double GetCTLeft( void )  const { return ct1_; }
  double GetCTDown( void )  const { return ct2_; }
  double GetCTRight( void ) const { return ct2_; }

  double MeanTime( void )  const { return 0.5*(t1_+t2_); }
  double CMeanTime( void ) const { return 0.5*(ct1_+ct2_); }
  double DeltaE( void ) const;

  int DetectorId( void ) const { return raw_->DetectorId(); }
  int PlaneId( void ) const { return raw_->PlaneId(); }
  int SegmentId( void ) const { return raw_->SegmentId(); }

  virtual bool ReCalc( bool applyRecursively=false )
  { return calculate(); }

};


#endif
