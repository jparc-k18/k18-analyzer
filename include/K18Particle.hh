/*
  K18Particle.hh

  2012/1/24
*/

#ifndef K18Particle_h
#define K18Particle_h 1

#include "ThreeVector.hh"

class K18Track;
class HodoCluster;
class BH2Cluster;

class K18Particle
{
public:
  K18Particle( K18Track *track, HodoCluster *Bh1, BH2Cluster *Bh2 )
    : Track_(track), Bh1_(Bh1), Bh2_(Bh2)
  {}
  ~K18Particle() {}
private:
  K18Track *Track_;
  HodoCluster *Bh1_;
  BH2Cluster *Bh2_;

public:
  K18Track * GetTrack( void ) { return Track_; }
  HodoCluster *GetBH1( void ) { return Bh1_; }
  BH2Cluster *GetBH2( void ) { return Bh2_; }
  
  ThreeVector Momentum( void ) const;
  ThreeVector Position( void ) const;
  double TimeOfFlight( void ) const;
};



#endif
