/*
  Kinematics.hh

  2012/1/24
*/

#include "Kinematics.hh"
#include <cmath>
#include <iostream>
#include <cstdlib>
#include <stdio.h>

const double LightVel = 299.792458;
const double z_offset =  0.0;

double MassSquare( double p, double pathL, double flightTime )
{
  double beta=pathL/flightTime/LightVel;
  return p*p*(1.-beta*beta)/beta/beta;
}

ThreeVector VertexPoint( const ThreeVector & Xin, const ThreeVector & Xout,
                         const ThreeVector & Pin, const ThreeVector & Pout )
{
  double xi=Xin.x(), yi=Xin.y(), xo=Xout.x(), yo=Xout.y();
  double ui=Pin.x()/Pin.z(), vi=Pin.y()/Pin.z();
  double uo=Pout.x()/Pout.z(), vo=Pout.y()/Pout.z();

  double z=((xi-xo)*(uo-ui)+(yi-yo)*(vo-vi))/
    ((uo-ui)*(uo-ui)+(vo-vi)*(vo-vi));
  double x1=xi+ui*z, y1=yi+vi*z;
  double x2=xo+uo*z, y2=yo+vo*z;
  return ThreeVector( 0.5*(x1+x2), 0.5*(y1+y2), z+z_offset );

}
 
ThreeVector VertexPointByHonly( const ThreeVector & Xin, 
				const ThreeVector & Xout,
				const ThreeVector & Pin, 
				const ThreeVector & Pout )
{
  double xi=Xin.x(), yi=Xin.y(), xo=Xout.x(), yo=Xout.y();
  double ui=Pin.x()/Pin.z(), vi=Pin.y()/Pin.z();
  double uo=Pout.x()/Pout.z(), vo=Pout.y()/Pout.z();

  double z=(xi-xo)/(uo-ui);
  return ThreeVector( xi+ui*z, yi+vi*z, z+z_offset );
}

ThreeVector VertexPoint( const ThreeVector & Xin, const ThreeVector & Xout,
                         const ThreeVector & Pin, const ThreeVector & Pout,
			 double & dist )
{
  double xi=Xin.x(), yi=Xin.y(), xo=Xout.x(), yo=Xout.y();
  double ui=Pin.x()/Pin.z(), vi=Pin.y()/Pin.z();
  double uo=Pout.x()/Pout.z(), vo=Pout.y()/Pout.z();

  double z=((xi-xo)*(uo-ui)+(yi-yo)*(vo-vi))/
    ((uo-ui)*(uo-ui)+(vo-vi)*(vo-vi));
  double x1=xi+ui*z, y1=yi+vi*z;
  double x2=xo+uo*z, y2=yo+vo*z;
  dist=sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2));

  return ThreeVector( 0.5*(x1+x2), 0.5*(y1+y2), z+z_offset );
}
 
double closeDist( const ThreeVector & Xin, const ThreeVector & Xout,
		  const ThreeVector & Pin, const ThreeVector & Pout )
{
  double xi=Xin.x(), yi=Xin.y(), xo=Xout.x(), yo=Xout.y();
  double ui=Pin.x()/Pin.z(), vi=Pin.y()/Pin.z();
  double uo=Pout.x()/Pout.z(), vo=Pout.y()/Pout.z();

  double z=((xi-xo)*(uo-ui)+(yi-yo)*(vo-vi))/
    ((uo-ui)*(uo-ui)+(vo-vi)*(vo-vi));
  double x1=xi+ui*z, y1=yi+vi*z;
  double x2=xo+uo*z, y2=yo+vo*z;

  return sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2));
}

/*
 Correct Energy loss (Upstream part of Vertex)
*/
const double TARGETcenter  = 6.0;//Z position
const double TARGEThw      = 120.0/2.0;
const double TARGETradius  = 67.3/2.0;
const double TARGETcenterX = 0.0;
const double TARGETcenterY = 0.0;

ThreeVector CorrElossIn(const ThreeVector & Pin, const ThreeVector & Xin, 
			const ThreeVector & vtx, double mass)
{
  double length, dE;
  double mom_new, energy_new;

  ThreeVector CorPin = Pin;
  double mom=Pin.Mag();

  mom_new = mom;
  energy_new = sqrt(mass*mass+mom*mom);

  if (fabs(vtx.z()-TARGETcenter) <= TARGEThw)     /* in LH2 target */
    length = calcLengthBeam(Pin, Xin, vtx);
  else if (vtx.z() > TARGEThw+TARGETcenter)      /* Downstream of LH2 target */
    {
      double u=Pin.x()/Pin.z();
      double v=Pin.y()/Pin.z();
      double ztgtOut=TARGETcenter+TARGEThw;
      ThreeVector tgtOutPos(u*(ztgtOut)+Xin.x(), v*(ztgtOut)+Xin.y(), ztgtOut);
      length = calcLengthBeam(Pin, Xin, tgtOutPos);
    }
  else if (vtx.z() < (-1)*(TARGEThw+TARGETcenter)) /* Upstream of LH2 target */
    return CorPin;
  
  if (caldE(mom, mass, length, &mom_new, &energy_new)) {
    CorPin = mom_new/mom*Pin;
    //std::cout << "CorrElossIn:: mom = " << mom << ", mom_new = " << mom_new
    //<< std::endl;
    return CorPin;
  } else {
    return Pin;
  }
  
}

double calcLengthBeam(const ThreeVector & Pin, const ThreeVector & Xin, 
		      const ThreeVector & vtx)
{
  double u=Pin.x()/Pin.z();
  double v=Pin.y()/Pin.z();
  
  //std::cout << "calcLengthBeam:: vtx (x,y,z)=(" << vtx.x() << ", " 
  //	    << vtx.y() << ", " << vtx.z() << ")" << std::endl;
  
  ThreeVector point1, point2;
  /* vertex point check */
  if (IsInsideTarget(vtx))
    point1=vtx;
  else {
    double z1,z2;
    double z;
    if (calcCrossingPoint(u, v, vtx, &z1, &z2)) {
      if (fabs(vtx.z()-z1)<fabs(vtx.z()-z2))
	z=z1;
      else
	z=z2;
      
      point1 = ThreeVector(u*(z-vtx.z())+vtx.x(), v*(z-vtx.z())+vtx.y(), z);
    } else {
      return 0.0;
    }
  }
  //std::cout << "calcLengthBeam:: Point1 (x,y,z)=(" << point1.x() << ", " 
  //	    << point1.y() << ", " << point1.z() << ")" << std::endl;
  /* target entrance point check */
  double ztgtIn=TARGETcenter-TARGEThw;

  ThreeVector tgtInPos(u*(ztgtIn)+Xin.x(), v*(ztgtIn)+Xin.y(), ztgtIn);
  if (IsInsideTarget(tgtInPos))
    point2=tgtInPos;
  else {
    double z1,z2;
    double z;
    if (calcCrossingPoint(u, v, tgtInPos, &z1, &z2)) {
      if (fabs(tgtInPos.z()-z1)<fabs(tgtInPos.z()-z2))
	z=z1;
      else
	z=z2;
      
      point2 = ThreeVector(u*(z)+Xin.x(), v*(z)+Xin.y(), z);
    } else {
      return 0.0;
    }
  }

  //std::cout << "calcLengthBeam:: Point2 (x,y,z)=(" << point2.x() << ", " 
  //	    << point2.y() << ", " << point2.z() << ")" << std::endl;
  //std::cout << "calcLengthBeam:: length=" << (point1-point2).Mag() << std::endl;
  //
  return (point1-point2).Mag();
}

/*
  Correct Energy loss (Downstream part of Vertex)
*/
ThreeVector CorrElossOut(const ThreeVector & Pout, const ThreeVector & Xout, 
			 const ThreeVector & vtx, double mass)
{
  double FL,FH,FTMP;
  double Elow, Ehigh, Elast, EPS;
  double E, length;
  
  ThreeVector CorPout = Pout;
  double mom = Pout.Mag();

  if (fabs(vtx.z()-TARGETcenter) < TARGEThw)     /* in LH2 target */
    length =  calcLengthScat(Pout, Xout, vtx);
  else if (vtx.z() < (-1)*(TARGEThw-TARGETcenter)) /* Upstream of LH2 target */
    {
      double u=Pout.x()/Pout.z();
      double v=Pout.y()/Pout.z();
      double ztgtIn=TARGETcenter-TARGEThw;
      ThreeVector tgtInPos(u*(ztgtIn)+Xout.x(), v*(ztgtIn)+Xout.y(), ztgtIn);
      length = calcLengthScat(Pout, Xout, tgtInPos);
    }
  else if (vtx.z() > (TARGEThw+TARGETcenter)) {     /* Downstream of SciFi */
    return CorPout;
  }

  Elow  = mass;
  Ehigh = 10.;
  Elast = sqrt(mass*mass+mom*mom);
  EPS = 0.001;
  FL  = diffE(mass, Elow, length, Elast);
  FH  = diffE(mass, Ehigh, length, Elast);
  while (fabs((Ehigh-Elow)/Ehigh) > EPS) {
    E = Ehigh-(Ehigh-Elow)*FH/(FH-FL);
    //printf("-------E=%f, Elow=%f, Ehigh=%f, FL=%f, FH=%f----------\n",E, Elow, Ehigh, FL, FH);
     if (fabs(FTMP=diffE(mass, E, length, Elast)) < 0.0000001){
      Elow = E;
      Ehigh = E;
      FH = FTMP;
    } if ((FTMP=diffE(mass, E, length, Elast)) < 0) {
      Elow = E;
      FL = FTMP;
    } else if ((FTMP=diffE(mass, E, length, Elast)) > 0){
      Ehigh = E;
      FH = FTMP;
    }
  }

  double energy_new = E;
  double mom_new = sqrt(energy_new*energy_new-mass*mass);
  
  CorPout = mom_new/mom*Pout;
  //std::cout << "CorrElossOut:: mom = " << mom << ", mom_new = " << mom_new
  //	    << std::endl;

  return CorPout;
}

double calcLengthScat(const ThreeVector & Pout, const ThreeVector & Xout, 
		      const ThreeVector & vtx)
{
  double u=Pout.x()/Pout.z();
  double v=Pout.y()/Pout.z();

  //std::cout << "calcLengthScat:: vtx (x,y,z)=(" << vtx.x() << ", " 
  //	    << vtx.y() << ", " << vtx.z() << ")" << std::endl;

  ThreeVector point1, point2;
  /* vertex point check */
  if (IsInsideTarget(vtx))
    point1=vtx;
  else {
    double z1,z2;
    double z;
    if (calcCrossingPoint(u, v, vtx, &z1, &z2)) {
      if (fabs(vtx.z()-z1)<fabs(vtx.z()-z2))
	z=z1;
      else
	z=z2;
      
      point1 = ThreeVector(u*(z-vtx.z())+vtx.x(), v*(z-vtx.z())+vtx.y(), z);
    } else {
      return 0.0;
    }
  }
  //std::cout << "calcLengthScat:: Point1 (x,y,z)=(" << point1.x() << ", " 
  //	    << point1.y() << ", " << point1.z() << ")" << std::endl;

  /* target exit point check */
  double ztgtOut=TARGETcenter+TARGEThw;

  ThreeVector tgtOutPos(u*(ztgtOut)+Xout.x(), v*(ztgtOut)+Xout.y(), ztgtOut);
  if (IsInsideTarget(tgtOutPos))
    point2=tgtOutPos;
  else {
    double z1,z2;
    double z;
    if (calcCrossingPoint(u, v, tgtOutPos, &z1, &z2)) {
      if (fabs(tgtOutPos.z()-z1)<fabs(tgtOutPos.z()-z2))
	z=z1;
      else
	z=z2;
      
      point2 = ThreeVector(u*(z)+Xout.x(), v*(z)+Xout.y(), z);
    } else {
      return 0.0;
    }
  }

  //std::cout << "calcLengthScat:: Point2 (x,y,z)=(" << point2.x() << ", " 
  //	    << point2.y() << ", " << point2.z() << ")" << std::endl;
  //std::cout << "calcLengthScata:: length=" << (point1-point2).Mag() << std::endl;

  return (point1-point2).Mag();
    
}

ThreeVector CorrElossOutCheck(const ThreeVector & Pout, const ThreeVector & Xout, 
			      const ThreeVector & vtx, double mass)
{
  double length, dE;
  double energy_new, mom_new;

  double mom=Pout.Mag();

  mom_new = mom;
  energy_new = sqrt(mass*mass+mom*mom);

  ThreeVector CorPout = Pout;

  if (fabs(vtx.z()-TARGETcenter) <= TARGEThw)     /* in LH2 target */
    length =  calcLengthBeam(Pout, Xout, vtx);
  else if (vtx.z() < (-1)*(TARGEThw-TARGETcenter)) /* Upstream of LH2 target */
    {
      double u=Pout.x()/Pout.z();
      double v=Pout.y()/Pout.z();
      double ztgtIn=TARGETcenter-TARGEThw;
      ThreeVector tgtInPos(u*(ztgtIn)+Xout.x(), v*(ztgtIn)+Xout.y(), ztgtIn);
      length = calcLengthBeam(Pout, Xout, tgtInPos);
    }
  else if (vtx.z() > TARGEThw+TARGETcenter)      /* Downstream of LH2 target */
    return CorPout;

  if (caldE(mom, mass, length, &mom_new, &energy_new)) {
    CorPout = mom_new/mom*Pout;
    return CorPout;
  } else {
    return Pout;
  }
}

bool IsInsideTarget(const ThreeVector &point)
{
  if ((point.x()-TARGETcenterX)*(point.x()-TARGETcenterX)+
      (point.y()-TARGETcenterY)*(point.y()-TARGETcenterY) <=
      TARGETradius*TARGETradius)
    return true;
  else 
    return false;
}

bool  calcCrossingPoint(double u, double v, ThreeVector Point, double *z1, double *z2)
{
  double x0=TARGETcenterX, y0=TARGETcenterY;
  double a=Point.x()-u*Point.z()-x0;
  double b=Point.y()-v*Point.z()-y0;
  double r=TARGETradius;

  double c=(a*u+b*v)*(a*u+b*v)-(u*u+v*v)*(a*a+b*b-r*r);

  if (c<0) {
    std::cerr << "This track does not cross target" << std::endl;
    return false;
  } else {
    *z1 = (-(a*u+b*v)+sqrt(c))/(u*u+v*v);
    *z2 = (-(a*u+b*v)-sqrt(c))/(u*u+v*v);
    return true;
  }
}

double diffE(double mass, double E, double length, double Elast)
{
  double p;
  double mom_new, energy_new;

  p = sqrt(E*E-mass*mass);

  caldE(p, mass, length, &mom_new, &energy_new);
  return (energy_new-Elast);
  //return (sqrt(mass*mass+p*p)-caldE(particle,p,length)-Elast);
}

int caldE(double momentum, double mass, double distance, 
	  double *momentum_cor, double *energy_cor)
{
  double dE_dx; /*MeV/cm*/
  double eloss; /*MeV*/

  double beta;
  double thickness = distance/10.0; /*cm*/
  double m = mass*1000.0;  /*mass of incident particle(Mev)*/ 

  double E_0;
  double E;
  double p = momentum*1000.0;
  //double delta=0.01; /*cm*/
  double delta=0.1; /*cm*/
  int i;
  double length=0.0;
  double total_eloss=0.0;

  //E_0=mygamma(beta)*m;
  //p=mygamma(beta)*beta*m;

  E_0= sqrt(m*m + p*p);
  E=E_0;
  beta = mybeta(E,p);
  //printf("beta=%f, E=%f, p=%f\n",beta, E, p);
  if (beta<=0.0) {
    *momentum_cor = p/1000.0;
    *energy_cor = E/1000.0;
    return 1;
  }
  dE_dx=0.0;
  eloss=0.0;
  for(i=0;i<=thickness/delta;i++){
    dE_dx=calc_dE_dx(beta);
    eloss=dE_dx*delta;
    E=E-eloss;
    if(E<m){
      fprintf(stderr,"particle stops in material at %5.3fcm\n",length);
      *momentum_cor = p/1000.0;
      *energy_cor = E/1000.0;
      return 0;
      //break;
    }
    p=sqrt(pow(E,2.0)-pow(m,2.0));
    beta=mybeta(E,p);
    length=length+delta;
    total_eloss=total_eloss+eloss;
    /*
      printf("beta:%5.3f\n",beta);
      printf("dE_dx:%5.3f\teloss:%5.3f\n",dE_dx,eloss);
      printf("E:%5.3f(MeV)\tp:%5.3f(MeV/c)\n",E,p);
      printf("length:%5.3f(cm)\n",length);
      printf("total energy loss:%5.3f(MeV)\n",total_eloss);
    */
    //getchar();
  }
  *momentum_cor = p/1000.0;
  *energy_cor = E/1000.0;

  return 1;
}
 
double calc_dE_dx(double beta)
{
  double value;
  const double C=0.1535; /*MeVcm^2/g*/
  const double m_e=0.511;
  double logterm;
  double gamma_2;
  double W_max;
  double gamma;
  double X;
  double delta;

  double rho=0.0709;   /*g/cm^3 (C)*/
  double I=21.8;     /*eV*/
  double Z_A=0.99216;
  int z=1;
  double C0=-3.2632;
  double X0=0.4759;
  double X1=1.9215;
  double a=0.13483;
  double M=5.6249;

  gamma = mygamma(beta);
  X = log10(beta*gamma);
  if (X<=X0)
    delta=0.0;
  else if (X0<X && X<X1) 
    delta=4.6052*X+C0+a*pow((X1-X),M);
  else if (X>=X1) 
    delta=4.6052*X+C0;

  gamma_2=1/(1-pow(beta,2.0));
  //printf("sqrt(gamma_2):%f\n",sqrt(gamma_2));

  W_max=2.0*m_e*pow(beta,2.0)*gamma_2;

  logterm=log(2.0*m_e*gamma_2*pow(beta,2.0)*W_max*pow(10.0,12.0)/pow(I,2.0))-2.0*pow(beta,2.0)-delta/2.0;

  value=C*rho*Z_A*pow((double)z,2.0)*logterm/pow(beta,2.0);

  return value;

}

double mygamma(double beta)
{
  double value;
  
  value=1.0/sqrt(1.0-pow(beta,2.0));

  return value;

}

double mybeta(double energy,double mormentum)
{
  double value;

  value=mormentum/energy;

  return value;
}

