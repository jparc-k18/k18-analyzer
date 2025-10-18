// -*- C++ -*-

#ifndef UTILITY_HH
#define UTILITY_HH

#include <TTimeStamp.h>

#include "UnpackerManager.hh"

//_____________________________________________________________________________
namespace utility
{
UInt_t     EBDataSize( void );
UInt_t     UnixTime( void );
TTimeStamp TimeStamp( void );
}

#endif
