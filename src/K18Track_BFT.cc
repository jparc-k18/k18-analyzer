/*
  K18Track_BFT.cc
*/

#include "K18Track_BFT.hh"
#include "DCLocalTrack.hh"
#include "K18TransMatrix.hh"
#include "DCGeomMan.hh"
#include "DetectorID.hh"
#include "DCHit.hh"
#include "TrackHit.hh"
#include "TemplateLib.hh"
#include "Minuit.hh"
#include "ConfMan.hh"
#include "K18Parameters.hh"
#include "DCAnalyzer.hh"

#include <string>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <cstring>
#include <stdexcept>
#include <sstream>

const double LowBand[5] = 
  { MinK18InX, MinK18InY, MinK18InU, MinK18InV, MinK18Delta };

const double UpperBand[5] = 
  { MaxK18InX, MaxK18InY, MaxK18InU, MaxK18InV, MaxK18Delta };

const double LowBandOut[5] =
  { MinK18OutX, MinK18OutY, MinK18OutU, MinK18OutV, MinK18Delta };

const double UpperBandOut[5] =
  { MaxK18OutX, MaxK18OutY, MaxK18OutU, MaxK18OutV, MaxK18Delta };


K18Track_BFT::K18Track_BFT( double  LocalX, DCLocalTrack *tout, double P0 )
  : LocalX_(LocalX), TrOut_(tout), P0_(P0), 
    StatusD2U_(false), gfastatus_(true)
{}

K18Track_BFT::~K18Track_BFT()
{
}

bool K18Track_BFT::CalcMomentumD2U( void )
{
  static const std::string funcname = "[K18Track_BFT::CalcMomentumD2U]";
  StatusD2U_=false;

  K18TransMatrix *mat=ConfMan::GetConfManager()->GetK18Matrix();

  double xi, yi, ui, vi, xo, delta1, delta2;
  xi= -TrOut_->GetX0();
  yi=  TrOut_->GetY0();
  ui=  TrOut_->GetU0();
  vi= -TrOut_->GetV0();
  xo= -LocalX_;

  StatusD2U_ = mat->CalcDeltaD2U(xi, yi, ui, vi, xo, delta1, delta2);

#if 0
  std::cout << funcname << ": after calculation. "
            << " StatusD2U=" << StatusD2U_  << std::endl;
#endif

  if( StatusD2U_ ){
    //std::cout << "delta1 = " << delta1 << ", delta2 = " << delta2 << std::endl;
    DeltaD2U_     =delta1;
    DeltaD2U_3rd_ =delta2;
  }

  return StatusD2U_;
}

ThreeVector K18Track_BFT::BeamMomentumD2U( void ) const
{
  double u=TrOut_->GetU0(), v=TrOut_->GetV0();
  double pz=PD2U()/sqrt(1.+u*u+v*v);

  return ThreeVector( pz*u, pz*v, pz );
}

double K18Track_BFT::Xtgt( void ) const
{
  double z=DCGeomMan::GetInstance().GetLocalZ( IdK18Target );
  return TrOut_->GetU0()*z+TrOut_->GetX0();
}

double K18Track_BFT::Ytgt( void ) const
{
  double z=DCGeomMan::GetInstance().GetLocalZ( IdK18Target );
  return TrOut_->GetV0()*z+TrOut_->GetY0();
}

double K18Track_BFT::Utgt( void ) const
{
  return TrOut_->GetU0();
}

double K18Track_BFT::Vtgt( void ) const
{
  return TrOut_->GetV0();
}
