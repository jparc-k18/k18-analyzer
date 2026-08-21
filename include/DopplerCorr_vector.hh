/*
  DopplerCorr.hh
*/

#ifndef DopplerCorr_h 
#define DopplerCorr_h

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
  DopplerCorr(ThreeVector recp, ThreeVector vert )
    :RecoilMom(1, recp), Vertex(1, vert)
  {}
  ~DopplerCorr(){}

  double CalcDopplerFactor(int slot, int id);
  ThreeVector GetGePosition(int slot);
  void initiallize(void)
  {
    RecoilMom.clear();
    Vertex.clear();
    ThreeVector a(-1.,-1.,-1.);
    RecoilMom.push_back(a);
    Vertex.push_back(a);
    /*
    ThreeVector a(1.,1.,1.);
    RecoilMom.at(0)= a;
    //   Vertex.at(0) = null;
    */
  }

  void SetDopplerCorr( ThreeVector recp, ThreeVector vert )
  {
    ThreeVector dummy(-1.,-1.,-1.);
    if (dummy == RecoilMom.at(0) && dummy == Vertex.at(0))
      {
	RecoilMom.at(0) = recp;
	Vertex.at(0) = vert;
	std::cout<<"doppler correction set at firsitposition"<<std::endl;
      }
    else
      {
	RecoilMom.push_back(recp);
	Vertex.push_back(vert);
      }
  }
  int SizeDopplerCorr(void) const {
    ThreeVector null(-1.,-1.,-1.);
    if (null == RecoilMom.at(0) && null == Vertex.at(0))
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
	    std::cerr<<"Recoil Size not match with Vertex Size"<<std::endl;
	    return 0;
	  }
      }
  }

};

#endif
