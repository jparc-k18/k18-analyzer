/*
  DCTrackUtilities.cc

  2012/1/24
*/

#include "DCTrackUtilities.hh"
#include "DCLocalTrack.hh"
#include "DCGeomMan.hh"
#include "DCLTrackHit.hh"

#include <string>
#include <iostream>
#include <iomanip>
#include <cmath>

const double Deg2Rad = acos(-1.)/180.;
const double Rad2Deg = 180./acos(-1.);

const double CellSize[] = {
  0.0,
  5.0, 5.0, 5.0, 5.0, 10.0,
  3.0, 3.0, 3.0, 3.0, 3.0, 3.0
};

bool HasHitsAlongStraightLine( DCLocalTrack *TrSdcIn,
			       const ThreeVector & bPos,
			       const ThreeVector & bMom )
{
  static const std::string funcname = "[HasHitsAlongStraightLine]";

  double x0=bPos.x(), y0=bPos.y();
  double u0=bMom.x()/bMom.z(), v0=bMom.y()/bMom.z();

  const DCGeomMan & geomMan = DCGeomMan::GetInstance();
  
  int nh=0;

  int n=TrSdcIn->GetNHit();
  for( int i=0; i<n; ++i ){
    DCLTrackHit *thit=TrSdcIn->GetHit(i);
    if(!thit) continue;
    int layer=thit->GetLayer();
    double z=geomMan.GetLocalZ(layer);
    double ta=thit->GetTiltAngle()*Deg2Rad;
    double ct=cos(ta), st=sin(ta);
    double lpcal=(x0+u0*z)*ct+(y0+v0*z)*st;
    double wp=thit->GetWirePosition();
    double lp=thit->GetLocalHitPos();

#if 0
    std::cout << funcname << ": Layer=" << layer
	      << " " << lpcal << " " << wp 
	      << " " << lp << " --> " << wp-lpcal
	      << " " << lp-lpcal <<  std::endl;
#endif

    //    if( fabs(lpcal-wp)<CellSize[layer] ){
    if( fabs(lpcal-lp)<1.0 ){
      ++nh;
    }
  }
#if 0
  int it;
  std::cout << funcname << "#Hits along line = " << nh << std::endl;
  std::cout << "#? ";
  std::cin >> it;
#endif


  if( nh>0 ) return true;
  else       return false; 
}


