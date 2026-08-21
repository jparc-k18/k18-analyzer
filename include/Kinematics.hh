/*
  Kinematics.hh

  2012/1/24
*/

#ifndef Kinematics_h
#define Kinematics_h 1

#include "ThreeVector.hh"

double MassSquare( double p, double pathL, double flightTime );

ThreeVector VertexPoint( const ThreeVector & Xin, const ThreeVector & Xout,
			 const ThreeVector & Pin, const ThreeVector & Pout );

ThreeVector VertexPointByHonly( const ThreeVector & Xin, 
				const ThreeVector & Xout,
				const ThreeVector & Pin, 
				const ThreeVector & Pout );

ThreeVector VertexPoint( const ThreeVector & Xin, const ThreeVector & Xout,
			 const ThreeVector & Pin, const ThreeVector & Pout,
			 double & dist );

double closeDist( const ThreeVector & Xin, const ThreeVector & Xout,
		  const ThreeVector & Pin, const ThreeVector & Pout );

ThreeVector CorrElossIn(const ThreeVector & Pin, 
			const ThreeVector & Xin, 
			const ThreeVector & vtx, double mass);

double calcLengthBeam(const ThreeVector & Pin, 
		      const ThreeVector & Xin, 
		      const ThreeVector & vtx);

ThreeVector CorrElossOut(const ThreeVector & Pout, 
			 const ThreeVector & Xout, 
			 const ThreeVector & vtx, double mass);

double calcLengthScat(const ThreeVector & Pout, 
		      const ThreeVector & Xout, 
		      const ThreeVector & vtx);

ThreeVector CorrElossOutCheck(const ThreeVector & Pout, 
			      const ThreeVector & Xout, 
			      const ThreeVector & vtx, double mass);

bool IsInsideTarget(const ThreeVector &point);

bool  calcCrossingPoint(double u, double v, ThreeVector Point, double *z1, double *z2);

double diffE(double mass, double E, double length, double Elast);

int caldE(double momentum, double mass, double distance, double *momentum_cor, double *energy_cor);

double calc_dE_dx(double beta);
double mygamma(double beta);
double mybeta(double energy,double mormentum);

#endif
