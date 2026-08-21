/*
  DopplerCorr.hh
*/

#ifndef DopplerCorr_h 
#define DopplerCorr_h 1

#include "ThreeVector.hh"
#include <vector>

class DopplerCorr
{
private:
  ThreeVector RecoilMom;
  ThreeVector Vertex;
  
public:
  DopplerCorr();
  DopplerCorr( const ThreeVector recp, const ThreeVector vert )
    :RecoilMom(recp), Vertex(vert)
  {}
  ~DopplerCorr(){}

  void Initialize( ThreeVector recp, ThreeVector vert );
  double CalcDopplerFactor(const int slot);
  double CalcGtheta(const int slot);
  double CalcBeta (void);
  ThreeVector GetGePosition(const int slot);
};

#endif
