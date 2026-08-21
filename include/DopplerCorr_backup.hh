/*
  DopplerCorr.hh
*/

#ifndef DopplerCorr_h 
#define DopplerCorr_h 1

#include "ThreeVector.hh"
#include <vector>

typedef std::vector <ThreeVector> ThreeVectorContainer;

class DopplerCorr
{
private:
  ThreeVectorContainer RecoilMom;
  ThreeVectorContainer Vertex;
  
public:
  DopplerCorr();
  DopplerCorr( ThreeVector recp, ThreeVector vert )
    :RecoilMom(1, recp), Vertex(1, vert)
  {};
  ~DopplerCorr(){};

  void Initialize( void )
  {
    RecoilMom.clear();
    Vertex.clear();
    ThreeVector a(-100., -100., -100.);
    RecoilMom.push_back(a);
    Vertex.push_back(a);
  }

  void SetDopplerCorr( ThreeVector recp, ThreeVector vert)
  {
    ThreeVector a(-100., -100., -100.);
    if (a == RecoilMom.at(0) && a == Vertex.at(0))
      {
	RecoilMom.at(0) = recp;
	Vertex.at(0) = vert;
	//std::cout<<"Doppler Correction: First recoil momentum and vertex point set"<<std::endl;
      }
    else
      {
	RecoilMom.push_back(recp);
	Vertex.push_back(vert);
      }
  }

  int SizeDopplerCorr(void) const {
    ThreeVector a(-100., -100., -100.);
    if ( a == RecoilMom.at(0) && a == Vertex.at(0))
      {
	return 0;
      }
    else
      {
	int recpsize = RecoilMom.size();
	int vertsize = Vertex.size();
	if (recpsize == vertsize) return recpsize;
	else
	  {
	    std::cerr<<"Doppler Correction: Recoil momentum and vertex points size are not matched"<<std::endl;
	    return 0;
	  }
      }

}

  double CalcDopplerFactor(int slot, int id);
  double CalcGtheta(int slot, int id);
  double CalcBeta (int id);
  ThreeVector GetGePosition(int slot);
};

#endif
