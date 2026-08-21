/*
  DCRawHit.hh

  2012/1/24
*/

#ifndef DCRawHit_h
#define DCRawHit_h 1

#include <cstddef>
#include <vector>

#ifdef MemoryLeak
#include "DebugCounter.hh"
#endif

typedef std::vector <int> IntVec;

class DCRawHit
{
private:
  int PlaneId_, WireId_;
  IntVec Tdc_;
  IntVec trailing_;

#ifdef MemoryLeak
  static debug::Counter sm_counter;
#endif

public:
  DCRawHit( int pid, int wid )
    :PlaneId_(pid), 
     WireId_(wid),
     Tdc_(0),
     trailing_(0)
  {
#ifdef MemoryLeak
    ++sm_counter;
#endif
}
  ~DCRawHit() {
#ifdef MemoryLeak
    --sm_counter;
#endif
}

  int PlaneId( void ) const { return PlaneId_; }
  int WireId( void ) const { return WireId_; }

public:
  void SetTdc( int tdc ) { Tdc_.push_back(tdc); }
  void SetTrailing( int tdc ) { trailing_.push_back(tdc); }

  int  GetTdc( int nh ) const { return Tdc_[nh]; }
  int  GetTdcSize( void ) const { return Tdc_.size(); }
  int  GetTrailing(int nh ) const { return trailing_[nh]; }
  };

#endif
