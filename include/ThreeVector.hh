/*
  ThreeVector.hh

  2012/1/24
*/

#ifndef ThreeVector_h
#define ThreeVector_h 1

#include <iostream>

#include <TVector3.h>

typedef TVector3 ThreeVector;

// CLHEP Hep3Vector compatibility helpers (not provided by TVector3)
inline TVector3 operator/(const TVector3& v, double a)
{
  return TVector3(v.X()/a, v.Y()/a, v.Z()/a);
}

inline std::ostream& operator<<(std::ostream& ost, const TVector3& v)
{
  ost << " (" << v.x() << "," << v.y() << "," << v.z() << ")";
  return ost;
}

#endif
