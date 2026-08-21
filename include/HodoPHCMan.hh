/*
  HodoPHCMan.hh

  2012/1/24
*/

#ifndef HODO_PHC_MAN_H
#define HODO_PHC_MAN_H

#include <map>
#include <string>

class HodoPHCParam
{
public:
  HodoPHCParam( int type, int np, double *parlist );
  ~HodoPHCParam();
private:
  HodoPHCParam( const HodoPHCParam &right );
  HodoPHCParam & operator = ( const HodoPHCParam & );
private:
  int Type, NPar;
  double *ParList;
public:
  double DoPHC( double time, double de );
  double DoRPHC( double time, double de );
private:
  double type1Correction( double time, double de );
  double type2Correction( double time, double de ); // For fiber
  double type1RCorrection( double time, double de );
};

class HodoPHCMan
{
private:
  std::string PHCFileName;
public:
  explicit HodoPHCMan( const std::string & filename );
  ~HodoPHCMan();
private:
  HodoPHCMan();
  HodoPHCMan( const HodoPHCMan & );
  HodoPHCMan & operator = ( const HodoPHCMan & );

public:
  void SetFileName( const std::string & filename ) { PHCFileName=filename; } 

private:
  typedef std::map <int, HodoPHCParam *> PhcPContainer;
  typedef std::map <int, HodoPHCParam *>::const_iterator PhcPIterator;

  PhcPContainer Container;

public:
  bool Initialize( void );
  bool doCorrection( int cid, int plid, int seg, int ud, double time,
                     double de, double &ctime );
  bool doRCorrection( int cid, int plid, int seg, int ud, double time,
                      double de, double &ctime );
private:
  HodoPHCParam * GetMap( int cid, int plid, int seg, int ud ); 
  void clearMap( void );
};


#endif
