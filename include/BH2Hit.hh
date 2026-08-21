/*
  BH2Hit.hh

  2012/1/24
*/

#ifndef BH2HIT_H
#define BH2HIT_H

#include "Hodo1Hit.hh"

class BH2Hit : public Hodo1Hit
{
public:
  BH2Hit( HodoRawHit *rhit )
    : Hodo1Hit(rhit), Tofs_(0.)
  {};
  ~BH2Hit(){};
  
public:
  bool calculate( void );
  double Time0( void )  const { return t_.at(0)+Tofs_; };
  double CTime0( void ) const { return ct_.at(0)+Tofs_; };
  
  virtual bool ReCalc( bool applyRecursively=false )
  { return calculate(); };
  
private:
  double Tofs_;
};
#endif
