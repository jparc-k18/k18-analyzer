/*
  TrackHit.cc

  2012/1/24
*/

#include "TrackHit.hh"
#include "DCGeomMan.hh"
#include "DCLocalTrack.hh"

#include <cstring>
#include <stdexcept>
#include <sstream>

TrackHit::TrackHit( DCLTrackHit *hit )
  : dchitp_(hit)
{

}

TrackHit::~TrackHit()
{
} 
