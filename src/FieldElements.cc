/*
  FieldElements.cc

  2012/1/24
*/

#include "FieldElements.hh"

FieldElements::FieldElements( const char *name,
                              const ThreeVector &pos,
                              double ta, double ra1, double ra2 )
  : record_( -1, name, pos, ta, ra1, ra2, 0., 0., 0., 0., 0. )
{}

ThreeVector FieldElements::
Local2GlobalPos( const ThreeVector &in ) const
{
  ThreeVector gPos = record_.Position();
  double x = gPos.x() 
    + record_.dxds()*in.x() 
    + record_.dxdt()*in.y() 
    + record_.dxdu()*in.z();
  double y = gPos.y() 
    + record_.dyds()*in.x() 
    + record_.dydt()*in.y() 
    + record_.dydu()*in.z();
  double z = gPos.z() 
    + record_.dzds()*in.x() 
    + record_.dzdt()*in.y() 
    + record_.dzdu()*in.z();

  return ThreeVector( x, y, z );
}

ThreeVector FieldElements::
Local2GlobalDir( const ThreeVector &in ) const
{
  double x 
    = record_.dxds()*in.x() 
    + record_.dxdt()*in.y() 
    + record_.dxdu()*in.z();
  double y 
    = record_.dyds()*in.x() 
    + record_.dydt()*in.y()
    + record_.dydu()*in.z();
  double z 
    = record_.dzds()*in.x() 
    + record_.dzdt()*in.y()
    + record_.dzdu()*in.z();

  return ThreeVector( x, y, z );
}

ThreeVector FieldElements::
Global2LocalPos( const ThreeVector &in ) const
{
  ThreeVector gPos = record_.Position();
  double x 
    = record_.dsdx()*(in.x()-gPos.x())
    + record_.dsdy()*(in.y()-gPos.y())
    + record_.dsdz()*(in.z()-gPos.z());
  double y 
    = record_.dtdx()*(in.x()-gPos.x())
    + record_.dtdy()*(in.y()-gPos.y())
    + record_.dtdz()*(in.z()-gPos.z());
  double z 
    = record_.dudx()*(in.x()-gPos.x())
    + record_.dudy()*(in.y()-gPos.y())
    + record_.dudz()*(in.z()-gPos.z());

  return ThreeVector( x, y, z );
}

ThreeVector FieldElements::
Global2LocalDir( const ThreeVector &in ) const
{
  double x 
    = record_.dsdx()*in.x()
    + record_.dsdy()*in.y()
    + record_.dsdz()*in.z();
  double y 
    = record_.dtdx()*in.x()
    + record_.dtdy()*in.y()
    + record_.dtdz()*in.z();
  double z 
    = record_.dudx()*in.x()
    + record_.dudy()*in.y()
    + record_.dudz()*in.z();

  return ThreeVector( x, y, z );
}
