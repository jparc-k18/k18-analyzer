/*
  FieldElements.hh

  2012/1/24
*/

#ifndef FieldElements_h 
#define FieldElements_h

#include "DCGeomRecord.hh"
#include "ThreeVector.hh"

enum FldElemReg { FERSurface=0, FERInside, FEROutside };

class FieldElements 
{
public:
  FieldElements( const char *name, const ThreeVector &pos,
                 double ta, double ra1, double ra2 );
  virtual ~FieldElements() {}

private:
  DCGeomRecord record_;

public:
  ThreeVector Local2GlobalPos( const ThreeVector &in ) const;
  ThreeVector Local2GlobalDir( const ThreeVector &in ) const;
  ThreeVector Global2LocalPos( const ThreeVector &in ) const;
  ThreeVector Global2LocalDir( const ThreeVector &in ) const;

  virtual ThreeVector GetField( const ThreeVector &gPos ) const = 0;
  virtual bool ExistField( const ThreeVector &gPos ) const = 0;

  virtual FldElemReg checkRegion( const ThreeVector &gPos, 
                                  double Tolerance ) const = 0;
};

#endif 
