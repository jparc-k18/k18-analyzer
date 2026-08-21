/*
  K18TrackFCN.hh

  2012/1/24
*/

#ifndef K18TrackFCN_h
#define K18TrackFCN_h 1

#include "Minuit.hh"

class K18Track;
class K18TransMatrix;

class K18TrackFCN : public MinuitFCN
{
  explicit K18TrackFCN( K18Track *track, K18TransMatrix *trMatrix );
  ~K18TrackFCN();

public:
  double operator()( int np, double *g, double *u, int flag );

private:  
  K18Track *Tr_;
  K18TransMatrix *trM_;

  friend class K18Track;
};


#endif
