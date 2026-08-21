/*
  SksParticle.hh

  2012/1/24
*/

#ifndef SksParticle_h
#define SksParticle_h 1

#include "ThreeVector.hh"

class SksTrack;
class HodoCluster;

class SksParticle
{
public:
  SksParticle( SksTrack *track, HodoCluster *Tof, HodoCluster *Lc )
    : Track_(track), Tof_(Tof), Lc_(Lc)
  {}
  ~SksParticle() {}
private:
  SksTrack *Track_;
  HodoCluster *Tof_, *Lc_;

public:
  SksTrack * GetTrack( void ) { return Track_; }
  HodoCluster * GetTof( void ) { return Tof_; }
  HodoCluster * GetLc( void ) { return Lc_; }

  ThreeVector Momentum( void ) const;
  ThreeVector Position( void ) const;
  double PathLengthToTOF( void ) const;
  double PathLengthTotal( void ) const;
  double MassSquare( double time0 ) const; 
};  

#endif
