/*
  K18TransMatrix.cc
*/

#include "K18TransMatrix.hh"

#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <iterator>
#include <vector>
#include <cmath>
#include <cstdlib>

#include "lexical_cast.hh"

const double MMtoM = 1.E-3;
const double MtoMM = 1000.;

const double PI = 2.*std::asin(1);

#define DEBUG 0

using namespace hddaq;

bool K18TransMatrix::Initialize( void )
{
  static const std::string funcname = "[K18TansMatrix::Initialize]";
  std::vector<std::vector<double> > parameters;
  parameters.resize(4);
  std::ifstream file( filename_.c_str() );

  if(file.fail()){
    std::cerr << funcname << ": file open fail" << std::endl;
    exit(-1);
  }

  while(file.good())
    {
      std::string buf;
      std::getline(file,buf);

      std::istringstream is(buf);
      std::istream_iterator<std::string> issBegin(is);
      std::istream_iterator<std::string> issEnd;
      std::vector<std::string> param(issBegin,issEnd);
      if(param.empty())continue;
      if(param[0].empty())continue;
      if(param[0]=="X" && param[1]!="A")
	{
	  for(int i=1;i<param.size();++i)
	    {
	      double p = a2d(param[i]);
	      parameters[0].push_back(p);
	    }
	}
      if(param[0]=="A")
	{
	  for(int i=1;i<param.size();++i)
	    {
	      double p = a2d(param[i]);
	      parameters[1].push_back(p);
	    }
	}
      if(param[0]=="Y" && param[1]!="B")
	{
	  for(int i=1;i<param.size();++i)
	    {
	      double p = a2d(param[i]);
	      parameters[2].push_back(p);
	    }
	}
      if(param[0]=="B")
	{
	  for(int i=1;i<param.size();++i)
	    {
	      double p = a2d(param[i]);
	      parameters[3].push_back(p);
	    }
	}
    }
  file.close();

  for(unsigned int i=0;i<parameters[0].size();++i)
    {
      Xpar[i] = parameters[0][i];
    }
  for(unsigned int i=0;i<parameters[1].size();++i)
    {
      Upar[i] = parameters[1][i];
    }
  for(unsigned int i=0;i<parameters[2].size();++i)
    {
      Ypar[i] = parameters[2][i];
    }
  for(unsigned int i=0;i<parameters[3].size();++i)
    {
      Vpar[i] = parameters[3][i];
    }
#if DEBUG
  for(unsigned int i=0;i<parameters[0].size();++i)
    {
      std::cout<<"Xpar[i]"<<" ";
    }
  std::cout<<std::endl;
  for(unsigned int i=0;i<parameters[1].size();++i)
    {
      std::cout<<"Upar[i]"<<" ";
    }
  std::cout<<std::endl;
  for(unsigned int i=0;i<parameters[2].size();++i)
    {
      std::cout<<"Ypar[i]"<<" ";
    }
  std::cout<<std::endl;
  for(unsigned int i=0;i<parameters[3].size();++i)
    {
      std::cout<<"Vpar[i]"<<" ";
    }
  std::cout<<std::endl;
#endif
  std::cout << funcname << ": Initialization finished" << std::endl;
  return true;
}

bool K18TransMatrix::
Transport( double xin, double yin, double uin, double vin, double delta,
	   double & xout, double & yout, double & uout, double & vout ) const
{
  xin *= -MMtoM; yin *= -MMtoM; uin = -uin; vin = -vin;
  xout=yout=uout=vout=0.0;

  double bj1[31], bj2[8];
  bj1[ 0] = xin;     
  bj1[ 1] = uin; 
  bj1[ 2] = delta; 
  bj1[ 3] = xin*xin;
  bj1[ 4] = xin*uin; 
  bj1[ 5] = xin*delta;
  bj1[ 6] = uin*uin;
  bj1[ 7] = uin*delta;
  bj1[ 8] = delta*delta;
  bj1[ 9] = yin*yin;
  bj1[10] = yin*vin;
  bj1[11] = vin*vin;
  bj1[12] = xin*xin*xin;
  bj1[13] = xin*xin*uin;
  bj1[14] = xin*xin*delta;
  bj1[15] = xin*uin*uin;
  bj1[16] = xin*uin*delta;
  bj1[17] = xin*delta*delta;
  bj1[18] = xin*yin*yin;
  bj1[19] = xin*yin*vin;
  bj1[20] = xin*vin*vin;
  bj1[21] = uin*uin*uin;
  bj1[22] = uin*uin*delta;
  bj1[23] = uin*delta*delta;
  bj1[24] = uin*yin*yin;
  bj1[25] = uin*yin*vin;
  bj1[26] = uin*vin*vin;
  bj1[27] = delta*delta*delta;
  bj1[28] = delta*yin*yin;
  bj1[29] = delta*yin*vin;
  bj1[30] = delta*vin*vin;

  bj2[ 0] = yin;
  bj2[ 1] = vin;
  bj2[ 2] = yin*xin;
  bj2[ 3] = yin*uin;
  bj2[ 4] = yin*delta;
  bj2[ 5] = vin*xin;
  bj2[ 6] = vin*uin;
  bj2[ 7] = vin*delta;

  for( int i=0; i<31; ++i ){
    xout += Xpar[i]*bj1[i];
    uout += Upar[i]*bj1[i];
  }

  for( int i=0; i<8; ++i ){
    yout += Ypar[i]*bj2[i];
    vout += Vpar[i]*bj2[i];
  }

  xout *= -MtoMM; yout *= -MtoMM; uout = -uout; vout = -vout;
  return true;
}

bool K18TransMatrix::
CalcDeltaD2U( double xin, double yin, double uin, double vin, double xout, double & delta1, double & delta2 ) const
{
  xin *= MMtoM; yin *= MMtoM; uin = uin; vin = vin;  xout *= MMtoM;

  // use untill 2nd order
  double A_, B_, C_, D_;
  A_ = Xpar[8];
  B_ = Xpar[2]+Xpar[5]*xin+Xpar[7]*uin;
  C_ = Xpar[0]*xin+Xpar[1]*uin+Xpar[3]*xin*xin+Xpar[4]*xin*uin+Xpar[6]*uin*uin
    + Xpar[9]*yin*yin+Xpar[10]*yin*vin+Xpar[11]*vin*vin -xout;

  D_ = B_*B_-4.*A_*C_;
  if (D_<0) {
    std::cerr << "[K18TransMatrix::CalcDeltaD2U] : negative Hanbetu-shiki" << std::endl;
    return false;
  }

  double ans1, ans2;

  ans1 = (-B_+sqrt(D_))/(2.*A_);
  ans2 = (-B_-sqrt(D_))/(2.*A_);

  delta1 = ans1;
  
  /*
  Calc cubic equation
  
  a*x^3 + b*x^2 + c^x + d = 0;
  ==>
  x^3 + alpha*x^2 + beta*x + gamma = 0;
  { alpha = b/a;  beta = c/a;  gamma = d/a; }
  
  x = y - alpha/3;
  ==>
  y^3 + py + q = 0; -- (eq 1)
  { p = (-alpha^2)/3 + beta;  q = (2*alpha^3)/27 - (alpha*beta)/3 + gamma}

  general solution of the equation below are,
  X^3 -3*u*v*X -(u^3 + v^3) = 0;
  X = 
    u + v         -- (1);
    w*u + (w^2)*v -- (2);
    (w^2)*u + w*v -- (3);
     { w is a solution of X^3 = 1},  w = (-1 + sqrt(3)*i)/2;}

  then,
  p = -3uv;  q = -(u^3 + v^3);

  u^3 = -q/2 + sqrt(q^2/4 + p^3/27)
  v^3 = -q/2 - sqrt(q^2/4 + p^3/27)

  D = -((q/2)^2 + (p/3)^3);
  D > 0 => 3 real number;
  D < 0 => 1 real number, 2 complex number;
  D = 0 => 2 real equal root;

  if(D > 0) (in this case, p < 0)
   y = r*cos(t) (r > 0, 0 < t <pi);
   (r^3)*cos^3(t) +p*cos(t) +q = 0;
   To use 4*cos^3(t) -3*cos(t) = cos(3t), r should be,
   r = sqrt(-4*p/3);
   then, cos(3t) = 3q/(rp);
   since 0 < t < pi, 0 < (1/3.)*arccos(3q/(rp)) < pi/3;
   therefore, these are real number solutions of (eq 1)
   y = 
    r*cos(t), r*cos(t + 2pi/3), r*cos(t - 2pi/3);
    (t = (1/3.)*arccos(3q/(rp)), x = y -alpha/3, r = sqrt(-4p/3));

  if(D < 0)
   this is a unique real number solution of (eq 1)
   y = cbrt(-q/2 + sqrt(-D)) + cbrt(-q/2 - sqrt(-D)); (x = y - alpha/3);
  */
  {
    double alpha, beta, gamma;
    double x, a, y, b;
    x = xin; a = uin; y = yin; b = vin;
    alpha = Xpar[TT] + Xpar[XTT]*x + Xpar[ATT]*a;
    beta  = Xpar[T] + Xpar[XT]*x + Xpar[AT]*a + Xpar[XXT]*x*x + Xpar[XAT]*x*a
      + Xpar[AAT]*a*a + Xpar[TYY]*y*y + Xpar[TYB]*y*b + Xpar[TBB]*b*b;
    gamma = Xpar[X]*x + Xpar[A]*a + Xpar[XX]*x*x + Xpar[XA]*x*a + Xpar[AA]*a*a
      + Xpar[YY]*y*y + Xpar[YB]*y*b + Xpar[BB]*b*b + Xpar[XXX]*x*x*x
      + Xpar[XXA]*x*x*a + Xpar[XAA]*x*a*a + Xpar[XYY]*x*y*y + Xpar[XYB]*x*y*b
      + Xpar[XBB]*x*b*b + Xpar[AAA]*a*a*a + Xpar[AYY]*a*y*y + Xpar[AYB]*a*y*b
      + Xpar[ABB]*a*b*b - xout;
    
    double p, q;
    p = -(alpha*alpha)/3 + beta;
    q = (2*alpha*alpha*alpha)/27 - (alpha*beta)/3 + gamma;

    double hanbetu = -((q/2)*(q/2) + (p/3)*(p/3)*(p/3));
    if(hanbetu > 0){
      double r = sqrt(-4*p/3);
      double s1, s2, s3;
      double theta = (1/3.)*std::acos(3*q/(r*p));
      s1 = r*cos(theta) -alpha/3;
      s2 = r*cos(theta + 2*PI/3.) -alpha/3;
      s3 = r*cos(theta - 2*PI/3.) -alpha/3;

      double tmpans = s1;
      if(fabs(tmpans-ans1) > fabs(s2-ans1)){
	tmpans = s2;
      }
      if(fabs(tmpans-ans1) > fabs(s3-ans1)){
	tmpans = s3;
      }
      delta2 = tmpans;
      //      std::cout << "ans1 " << ans1 << std::endl;
      //      std::cout << delta2 << std::endl;
    }else{
      double s1;
      s1 = cbrt(-q/2 + sqrt(-hanbetu)) + cbrt(-q/2 - sqrt(-hanbetu))
	-alpha/3;
      delta2 = s1;
      //      std::cout << "ans1 " << ans1 << std::endl;
      //      std::cout << delta2 << std::endl;
    }
  }
  
  return true;
}

