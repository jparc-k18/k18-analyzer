/*
  K18Particle.cc

  2012/1/24
*/

#include "K18Particle.hh"
#include "K18Track.hh"
#include "HodoCluster.hh"
#include "BH2Cluster.hh"

ThreeVector K18Particle::Momentum( void ) const
{
  return Track_->BeamMomentum();
}

ThreeVector K18Particle::Position( void ) const
{
  return ThreeVector( Track_->Xtgt(), Track_->Ytgt(), 0.0 );
}

double K18Particle::TimeOfFlight( void ) const
{
  double t1=Bh1_->CMeanTime();
  double t2=Bh2_->CTime0();

  return t2-t1;
}
