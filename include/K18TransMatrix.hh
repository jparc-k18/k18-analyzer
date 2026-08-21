/*
  K18TransMatrix.hh
*/

#ifndef K18TransMatrix_h
#define K18TransMatrix_h 1

#include <string>

class K18TransMatrix 
{
public:
  K18TransMatrix( const std::string &filename )
    : filename_(filename)
  {}
  ~K18TransMatrix() {} 

private:
  std::string filename_;

  enum NameX{X, A, T, XX, XA, XT, AA, AT, TT, YY,
	     YB, BB, XXX, XXA, XXT, XAA, XAT, XTT, XYY, XYB,
	     XBB, AAA, AAT, ATT, AYY, AYB, ABB, TTT, TYY, TYB,
	     TBB,
	     size_NameX};
  enum NameY{Y, B, YX, YA, YT, BX, BA, BT,
	     size_NameY};
  double Xpar[size_NameX], Ypar[size_NameY];
  double Upar[size_NameX], Vpar[size_NameY];

public:
  bool Initialize( void );

  bool Transport( double xin, double yin, double uin, double vin, double delta,
		  double & xout, double & yout, double & uout, double & vout ) const;
  bool CalcDeltaD2U( double xin, double yin, double uin, double vin, double xout, double & delta1, double & delta2 ) const;

};


#endif
