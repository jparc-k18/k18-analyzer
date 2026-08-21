/*
  DCTrackUtilities.hh

  2012/1/24
*/

#ifndef DCTrackUtilities_h

#define DCTrackUtilities_h

#include "ThreeVector.hh"

class DCLocalTrack;


bool HasHitsAlongStraightLine( DCLocalTrack * TrSdcIn, 
			       const ThreeVector & bPos,
			       const ThreeVector & bMom );

#endif
