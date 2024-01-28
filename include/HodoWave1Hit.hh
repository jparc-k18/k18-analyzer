/**
 *  file: Hodo1Hit.hh
 *  date: 2017.04.10
 *
 */

#ifndef HODOWAVE1_HIT_HH
#define HODOWAVE1_HIT_HH

#include "HodoRawHit.hh"
#include "HodoHit.hh"

class RawData;

//______________________________________________________________________________
class HodoWave1Hit : public HodoHit
{
public:
  explicit HodoWave1Hit( HodoRawHit *rhit, int ref_sample=0 );
  virtual ~HodoWave1Hit();

private:
  HodoWave1Hit( const HodoWave1Hit& );
  HodoWave1Hit& operator =( const HodoWave1Hit& );

protected:
  HodoRawHit          *m_raw;
  bool                 m_is_calculated;
  int                  m_multi_hit_l;
  int                  m_multi_hit_t;
  std::vector<double>  m_a_wf;
  std::vector<double>  m_t_wf;
  std::vector<double>  m_q;
  std::vector<double>  m_q_pedestal;
  std::vector<double>  m_a;
  std::vector<double>  m_t;
  std::vector<double>  m_ct;
  std::vector<bool>    m_flag_join;
  int                  m_ref_sample;

public:
  HodoRawHit* GetRawHit( void ) { return m_raw; }
  //bool   Calculate( void );
  bool   Calculate( bool tdc_flag = true );
  bool   IsCalculated( void ) const { return m_is_calculated; }
  //  int    GetNumOfHit(int sel=0) const { return sel==0 ? m_multi_hit_l : m_multi_hit_t; };
  int    GetNumOfHit(int sel=0) const { return sel==0 ? m_multi_hit_t : m_multi_hit_l; };
  int    GetNumOfWaveform( void ) const { return m_a_wf.size() == m_t_wf.size() ? m_a_wf.size() : 0;};
  double GetWaveformA( int n=0 )   const { if (n<m_a_wf.size()) return m_a_wf.at(n); 
                                     else return -9999.;}
  double GetWaveformT( int n=0 )   const { if (n<m_t_wf.size()) return m_t_wf.at(n); 
                                     else return -9999.;}

  double GetQ( int n=0 )   const { if (n<m_q.size()) return m_q.at(n); 
                                   else return -9999.;}
  double GetQ_ped( int n=0)const { if (n<m_q_pedestal.size()) return m_q_pedestal.at(n);
                                   else return -9999.;}

  double GetA( int n=0 )   const { if (n<m_a.size()) return m_a.at(n); 
                                   else return -9999.;}
  double GetT( int n=0 )   const { return m_t.at(n); }
  double GetCT( int n=0 )  const { return m_ct.at(n); }

  double Time( int n=0 )   const { return GetT(n); }
  double CTime( int n=0 )  const { return GetCT(n); }
  virtual
  double DeltaE( int n=0 ) const { return GetA(n); }

  double GetAUp( int n=0 )     const { return m_a.at(n); }
  double GetALeft( int n=0 )   const { return m_a.at(n); }
  double GetADown( int n=0 )   const { return m_a.at(n); }
  double GetARight( int n=0 )  const { return m_a.at(n); }

  virtual
  double GetTUp( int n=0 )     const { return m_t.at(n); }
  double GetTLeft( int n=0 )   const { return m_t.at(n); }
  virtual
  double GetTDown( int n=0 )   const { return m_t.at(n); }
  double GetTRight( int n=0 )  const { return m_t.at(n); }

  double GetCTUp( int n=0 )    const { return m_ct.at(n); }
  double GetCTLeft( int n=0 )  const { return m_ct.at(n); }
  double GetCTDown( int n=0 )  const { return m_ct.at(n); }
  double GetCTRight( int n=0 ) const { return m_ct.at(n); }

  virtual
  double MeanTime( int n=0 )  const { return GetT(n); }
  virtual
  double CMeanTime( int n=0 ) const { return GetCT(n); }

  virtual
  int DetectorId( void ) const { return m_raw->DetectorId(); }
  virtual 
  int PlaneId( void )    const { return m_raw->PlaneId(); }
  virtual
  int SegmentId( void )  const { return m_raw->SegmentId(); }

  // For BGO
  void   ClearACont( void )  { m_a.clear(); }
  void   SetE( double energy )  { m_a.push_back(energy); }

  void   SetJoined( int m )           { m_flag_join.at(m) = true;         }
  bool   Joined( int m )        const { return m_flag_join.at(m);         }
  bool   JoinedAllMhit();

  virtual bool ReCalc( bool applyRecursively=false )
  { return Calculate(); }

};

#endif
