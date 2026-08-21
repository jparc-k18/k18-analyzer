/*
  BH2Cluster.hh

  2012/1/24
*/

#ifndef BH2CLUSTER_H
#define BH2CLUSTER_H

class BH2Hit; 
class HodoAnalyzer;


#include <cstddef>

class BH2Cluster  
{
public:
  BH2Cluster( BH2Hit *hitA, BH2Hit *hitB=0, BH2Hit *hitC=0 );
  ~BH2Cluster(){};

private:
  BH2Cluster( const BH2Cluster & );
  BH2Cluster & operator = ( const BH2Cluster & );

private:
  BH2Hit *hitA_, *hitB_, *hitC_;
  int csize_;
  double MeanTime_, dE_;
  double MeanSeg_;
  double Tdiff_;
  double time0_;
  bool gfastatus_;
public:
  int ClusterSize( void ) const { return csize_; };
  double CMeanTime( void ) const { return MeanTime_; };
  double DeltaE( void ) const { return dE_; };
  double MeanSeg( void ) const { return MeanSeg_; };
  double CTime0( void ) const  { return time0_; };
  double TimeDif( void ) const { return Tdiff_; };

  BH2Hit * GetHit( int i ) const;
  bool GoodForAnalysis( void ) const { return gfastatus_; };
  bool GoodForAnalysis( bool status ) 
  { bool ret=gfastatus_; gfastatus_=status; return ret; } ;

  bool ReCalc( bool applyRecusively=false );

private:
  void calculate( void );

};
#endif
