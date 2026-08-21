/*
  SksParticle.cc

  2012/1/24
*/

#include "SksParticle.hh"
#include "SksTrack.hh"
#include "HodoCluster.hh"
#include "Kinematics.hh"

ThreeVector SksParticle::Momentum( void ) const
{
  return Track_->PrimaryMomentum();
}

ThreeVector SksParticle::Position( void ) const
{
  return Track_->PrimaryPosition();
}

double SksParticle::PathLengthToTOF( void ) const
{
  return Track_->PathLengthToTOF();
}

double SksParticle::PathLengthTotal( void ) const
{
  return Track_->PathLengthTotal();
}

double SksParticle::MassSquare( double time0 ) const
{
  double ttof=Tof_->CMeanTime();
  return ::MassSquare( Momentum().Mag(), PathLengthToTOF(), 
		       ttof-time0 );
}

