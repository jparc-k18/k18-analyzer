/*
  DCHit.hh

  2012/1/24
*/

#ifndef DCHit_h 
#define DCHit_h

#include <vector>
#include <deque>
#include <string>

//typedef std::vector <bool> BoolVec;
typedef std::deque <bool> BoolVec;
typedef std::vector <int> IntVec;
typedef std::vector <double> DoubleVec;

class DCLTrackHit;

#ifdef MemoryLeak
#include "DebugCounter.hh"
#endif

class DCHit
{
public:
  DCHit();
  DCHit( int layer );
  DCHit( int layer, double wire ); 

  ~DCHit();

private:
  DCHit( const DCHit & );
  DCHit & operator = ( const DCHit & );

private:
  int layer_;
  double wire_;
  IntVec tdc_;
  IntVec adc_;
  IntVec trailing_;
  DoubleVec dt_, dl_;
  DoubleVec trailingTime_;

  double wpos_;
  double angle_;
  BoolVec belongTrack_;
  BoolVec dlRange_;

  int clsize_;
  bool mwpcflag_;
  double mwire_;
  double mwpos_;

  mutable std::vector <DCLTrackHit *> Cont_;

#ifdef MemoryLeak
  static debug::Counter sm_counter;
#endif

public:
  bool CalcDCObservables( void );
  bool CalcMWPCObservables( void );
  bool CalcFiberObservables( void );
  bool CalcSsdObservables( void );
  //  bool CalcObservablesSimulation( double dlength);

  void SetLayer( int layer ) { layer_=layer; }
  void SetWire( double wire ) { wire_=wire; }  
  void SetTdcVal( int tdc );
  void SetAdcVal( int adc );
  void SetTdcTrailing(int tdc) {trailing_.push_back(tdc); }
  void SetDriftTime( double dt ) { dt_.push_back(dt); }
  void ClearDriftTime( void ) { dt_.clear(); }
  void SetDriftLength( double dl ) { dl_.push_back(dl); }
  void ClearDriftLength( void ) { dl_.clear();} 
  void SetTiltAngle( double angleDegree ) { angle_=angleDegree; }
  void SetTrailingTime( double t ) { trailingTime_.push_back(t); }

  void SetClusterSize( int size ) { clsize_=size; }
  void SetMWPCFlag( bool flag ) { mwpcflag_=flag; }
  void SetMeanWire( double mwire ) { mwire_=mwire; }
  void SetMeanWirePosition( double mwpos ) { mwpos_=mwpos; }
  void SetWirePosition( double wpos ) { wpos_ = wpos; }

  int GetLayer( void ) const { return layer_; }
  double GetWire( void )  const { 
    if( mwpcflag_ ) return mwire_; 
    else return int(wire_); 
  }

  int GetTdcSize( void ) const { return tdc_.size(); }
  int GetAdcSize( void ) const { return adc_.size(); }
  int GetDriftTimeSize( void ) const { return dt_.size(); }
  int GetDriftLengthSize( void ) const { return dl_.size(); }
  int GetTdcVal( int nh=0 ) const { return tdc_[nh]; }
  int GetAdcVal( int nh=0 ) const { return tdc_[nh]; }
  int GetTdcTrailing( int nh=0 ) const { return trailing_[nh]; }
  double GetDriftTime( int nh=0 ) const { return dt_[nh]; }
  double GetDriftLength( int nh=0 ) const { return dl_[nh]; }
  double GetTrailingTime( int nh=0 ) const { return trailingTime_[nh]; }

  double GetTiltAngle( void ) const { return angle_; }
  double GetWirePosition( void ) const { 
    if( mwpcflag_ ) return mwpos_; 
    else return wpos_; 
  }

  int GetClusterSize( void ) const { return clsize_; }
  double GetMeamWire( void ) const { return mwire_; }
  double GetMeamWirePosition( void ) const { return mwpos_; }

  void setFlags( int nh=0 ) { belongTrack_[nh]=true; }
  void clearFlags( int nh=0 ) { belongTrack_[nh]=false; }
  bool showFlags( int nh=0 ) const { return belongTrack_[nh]; }
  bool rangecheck( int nh=0 ) const { return dlRange_[nh]; }
  void setRangeCheckStatus( bool status, 
			    int nh ) 
  { dlRange_[nh] = (status); }

  void RegisterHits( DCLTrackHit *hit ) const
  { Cont_.push_back(hit); }

  bool ReCalcDC( bool applyRecursively=false ) 
  { return CalcDCObservables(); }
  bool ReCalcMWPC( bool applyRecursively=false ) 
  { return CalcMWPCObservables(); }

  void print(const std::string& arg="") const;

private:
  void clearRegisteredHits( void );
};


#endif
