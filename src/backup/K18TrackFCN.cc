/*
  K18TrackFCN.cc

  2012/1/24
*/

#include "K18TrackFCN.hh"
#include "K18Track.hh"
#include "K18TransMatrix.hh"
#include "TrackHit.hh"
#include "DCGeomMan.hh"

#include <string>
#include <cmath>

const double Deg2Rad = acos(-1.)/180.;
const double Rad2Deg = 180./acos(-1.);

K18TrackFCN::K18TrackFCN( K18Track *track, K18TransMatrix *trMatrix )
  : Tr_(track), trM_(trMatrix)
{}

K18TrackFCN::~K18TrackFCN()
{}

double K18TrackFCN::operator()( int np, double *g, double *u, int flag )
{
  static const std::string funcname = "[K18TrackFCN::operator()]";
  const DCGeomMan & geomMan=DCGeomMan::GetInstance();

  double xi=u[0], yi=u[1], ui=u[2], vi=u[3], delta=u[4];
  double xo, yo, uo, vo; 

  trM_->Transport( xi, yi, ui, vi, delta, xo, yo, uo, vo );

  int nin=Tr_->GetNHitsIn(), nout=Tr_->GetNHitsOut();

  double chi=0.0;
  int nh=0;

  for( int i=0; i<nin; ++i ){
    TrackHit *hit=Tr_->GetK18HitIn(i);
    if(!hit) continue;
    int lnum=hit->GetLayer();
    double z=geomMan.GetLocalZ( lnum );
    double dd=geomMan.GetResolution( lnum );
    double pos=hit->GetLocalHitPos();
    double aa=hit->GetHit()->GetTiltAngle()*Deg2Rad;
    double cpx=xi+ui*z, cpy=yi+vi*z;
    double c=cpx*cos(aa)+cpy*sin(aa);
    chi += ((pos-c)*(pos-c)/(dd*dd)); ++nh;
    if( flag==3 ) hit->SetCalLPos(c);
  }

  for( int i=0; i<nout; ++i ){
    TrackHit *hit=Tr_->GetK18HitOut(i);
    if(!hit) continue;
    int lnum=hit->GetLayer();
    double z=geomMan.GetLocalZ( lnum );
    double dd=geomMan.GetResolution( lnum );
    double pos=hit->GetLocalHitPos();
    double aa=hit->GetHit()->GetTiltAngle()*Deg2Rad;
    double cpx=xo+uo*z, cpy=yo+vo*z;
    double c=cpx*cos(aa)+cpy*sin(aa);
    chi += ((pos-c)*(pos-c)/(dd*dd)); ++nh;
    if( flag==3 ) hit->SetCalLPos(c);
  }

  if( nh>5 ) chi /= double(nh-5);

  return chi;
}
