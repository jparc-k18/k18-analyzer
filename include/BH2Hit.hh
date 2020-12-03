/**
 *  file: BH2Hit.hh
 *  date: 2017.04.10
 *
 */

#ifndef BH2_HIT_HH
#define BH2_HIT_HH

#include "DebugCounter.hh"
#include "Hodo1Hit.hh"

//______________________________________________________________________________
class BH2Hit : public Hodo1Hit
{
public:
  BH2Hit( HodoRawHit *rhit, double max_time_diff=10. );
  ~BH2Hit( void );

private:
  double m_time_offset;

public:
  bool   Calculate( void );
  double UTime0( int n=0 )  const { return m_t.at(n)  +m_time_offset; }
  double UCTime0( int n=0 ) const { return m_t.at(n)  +m_time_offset; }
  double DTime0( int n=0 )  const { return m_ct.at(n) +m_time_offset; }
  double DCTime0( int n=0 ) const { return m_ct.at(n) +m_time_offset; }

  double MeanTime( int n = 0 )  const { return m_t.at(n) ; }
  double CMeanTime( int n = 0 ) const { return m_ct.at(n); }
  double Time0( int n=0 )   const { return m_t.at(n)  +m_time_offset; }
  double CTime0( int n=0 )  const { return m_ct.at(n) +m_time_offset; }

  virtual bool ReCalc( bool applyRecursively=false )
  { return Calculate(); };

};

#endif
