/*
  K18Track.cc

  2012/1/24
*/

#include "K18Track.hh"
#include "DCLocalTrack.hh"
#include "K18TransMatrix.hh"
#include "DCGeomMan.hh"
#include "DetectorID.hh"
#include "TrackHit.hh"
#include "TemplateLib.hh"
#include "K18TrackFCN.hh"
#include "Minuit.hh"
#include "ConfMan.hh"
#include "K18Parameters.hh"
#include "DCAnalyzer.hh"
#include "UnpackerManager.hh"

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

const double IniError[5] =
  { 0.4, 0.4, 1.E-5, 1.E-4, 1.E-3 };

const int MaxFCNCall = 300;
const double EPS = 1.;  


K18Track::K18Track( DCLocalTrack *tin, DCLocalTrack *tout, double P0 )
  : TrIn_(tin), TrOut_(tout), P0_(P0), 
    Status_(false), gfastatus_(true)
{}

K18Track::~K18Track()
{
  deleteHits();
}

bool K18Track::doFit( void )
{
  static const std::string funcname = "[K18Track::doFit]";
  Status_=false;
  deleteHits();
  addHits();

  K18TransMatrix *mat=ConfMan::GetConfManager()->GetK18Matrix();

  K18TrackFCN FCN( this, mat );

  double param[5], error[5];
  param[0]=TrIn_->GetX0(); 
  param[1]=TrIn_->GetY0();
  param[2]=TrIn_->GetU0();
  param[3]=TrIn_->GetV0();
  param[4]=0.0;
  error[0]=IniError[0];
  error[1]=IniError[1];
  error[2]=IniError[2];

  error[3]=IniError[3];
  error[4]=IniError[4];

  double LowBand_[5], UpperBand_[5];
  for( int i=0; i<5; ++i ){
    LowBand_[i]=LowBand[i]; 
    UpperBand_[i]=UpperBand[i];
    if( param[i]<LowBand_[i] )   LowBand_[i]=param[i];
    if( param[i]>UpperBand_[i] ) UpperBand_[i]=param[i];
  }

#if 0
  std::cout << funcname << ": before fitting." << std::endl;
#endif
  
  Status_ = Minuit(&FCN).Fit( 5, param, error, LowBand_, UpperBand_,
			      MaxFCNCall, EPS, chisqr_ );

#if 0
  std::cout << funcname << ": after fitting. " 
	    << " Status=" << Status_  << std::endl;
#endif

  if( Status_ ){
    Xi_=param[0]; Yi_=param[1]; Ui_=param[2]; Vi_=param[3];
    Delta_=param[4];
    mat->Transport( Xi_, Yi_, Ui_, Vi_, Delta_,
		    Xo_, Yo_, Uo_, Vo_ );
  }

  return Status_;
}

ThreeVector K18Track::BeamMomentum( void ) const
{
  double u=TrOut_->GetU0(), v=TrOut_->GetV0();
  double pz=P()/sqrt(1.+u*u+v*v);

  return ThreeVector( pz*u, pz*v, pz );
}

double K18Track::Xtgt( void ) const
{
  double z=DCGeomMan::GetInstance().GetLocalZ( IdK18Target );
  return TrOut_->GetU0()*z+TrOut_->GetX0();
}

double K18Track::Ytgt( void ) const
{
  double z=DCGeomMan::GetInstance().GetLocalZ( IdK18Target );
  return TrOut_->GetV0()*z+TrOut_->GetY0();
}

double K18Track::Utgt( void ) const
{
  return TrOut_->GetU0();
}

double K18Track::Vtgt( void ) const
{
  return TrOut_->GetV0();
}

void K18Track::deleteHits( void )
{
  for_each( hitContIn.begin(),  hitContIn.end(),  DeleteObject() );
  for_each( hitContOut.begin(), hitContOut.end(), DeleteObject() );
  hitContIn.clear(); hitContOut.clear();
}

void K18Track::addHits( void )
{
  static const std::string funcname="K18Track::addHits";

  int nhIn=TrIn_->GetNHit(), nhOut=TrOut_->GetNHit();

  for( int i=0; i<nhIn; ++i ){
    DCLTrackHit *lhit = TrIn_->GetHit(i);
    if(!lhit) continue;
    TrackHit *hit = new TrackHit(lhit);
    if(hit){
      hitContIn.push_back(hit);
    }
    else{
      std::cerr << funcname << ": new fail" << std::endl;
    }      
  }

  for( int i=0; i<nhOut; ++i ){
    DCLTrackHit *lhit = TrOut_->GetHit(i);
    if(!lhit) continue;
    TrackHit *hit = new TrackHit(lhit);
    if(hit){
      hitContOut.push_back(hit);
    }
    else{
      std::cerr << funcname << ": new fail" << std::endl;
    }      
  }
}

TrackHit * K18Track::GetK18HitIn( int i )
{
  if( i>=0 && i<hitContIn.size() )
    return hitContIn[i];
  else
    return 0;
}

TrackHit * K18Track::GetK18HitOut( int i )
{
  if( i>=0 && i<hitContOut.size() )
    return hitContOut[i];
  else
    return 0;
}

TrackHit * K18Track::GetK18HitTotal( int i )
{
  int nin=hitContIn.size();
  if( i<0 ) 
    return 0;
  else if( i<nin )
    return hitContIn[i];
  else if( i<nin+hitContOut.size() )
    return hitContOut[i-nin];
  else
    return 0;
} 


TrackHit *K18Track::GetK18HitByPlaneId( int PlId )
{
  int id = PlId-PlOffsBc;
  TrackHit *hit=0;
  if( id>=PlMinBcIn && id<=PlMaxBcIn ){
    for( int i=0; i<hitContIn.size(); ++i ){
      if( hitContIn[i]->GetLayer()==PlId ){
	hit=hitContIn[i]; break;
      }
    }
  }
  else if( id>=PlMinBcOut && id<=PlMaxBcOut ){
    for( int i=0; i<hitContOut.size(); ++i ){
      if( hitContOut[i]->GetLayer()==PlId ){
	hit=hitContOut[i]; break;
      }
    }
  }
  return hit;
}

bool K18Track::ReCalc( bool applyRecursively )
{
  static const std::string funcname = "[K18Track::ReCalc]"; 

  bool ret1=true, ret2=true, ret=false;
  if( applyRecursively ){
    ret1=TrIn_->ReCalc(applyRecursively);
    ret2=TrOut_->ReCalc(applyRecursively);
  }
  if( ret1 && ret2 ){
    ret=doFit();
  }
  if(!ret){
    std::cerr << funcname << ": ReCalculation fails" << std::endl;
  }

  return ret;
}
