/*
  DCGeomMan.hh

  2012/1/24
*/

#ifndef DCGeomMan_h
#define DCGeomMan_h 1

#include "ThreeVector.hh"
#include <string>
#include <vector>
#include <map>

class DCGeomRecord;

class DCGeomMan 
{
private:
  DCGeomMan();
public:
  ~DCGeomMan();

public:
  void SetFileName( const char *filename ) { filename_=filename; }
  void SetFileName( const std::string &filename ) { filename_=filename; }

  bool Initialize( void );
  bool Initialize( const char *filename )
  { filename_=filename; return Initialize(); }
  bool Initialize( const std::string &filename )
  { filename_=filename; return Initialize(); }

  static DCGeomMan & GetInstance( void );
  double GetLocalZ( int lnum ) const;
  double GetResolution( int lnum ) const;

  // Do not use this method except for special cases
  void SetResolution( int lnum, double res ) const;
  //

  double GetTiltAngle( int lnum ) const;
  double GetRotAngle1( int lnum ) const;
  double GetRotAngle2( int lnum ) const;
  const ThreeVector & GetGlobalPosition( int lnum ) const;
  ThreeVector NormalVector( int lnum ) const;
  ThreeVector UnitVector( int lnum ) const;
  const DCGeomRecord *GetRecord( int lnum ) const;

  ThreeVector Local2GlobalPos( int lnum, const ThreeVector &in ) const;
  ThreeVector Global2LocalPos( int lnum, const ThreeVector &in ) const;
  ThreeVector Local2GlobalDir( int lnum, const ThreeVector &in ) const;
  ThreeVector Global2LocalDir( int lnum, const ThreeVector &in ) const;

  double calcWirePosition( int lnum, double wire ) const;
  int calcWireNumber( int lnum, double position ) const;

  std::vector<int> GetDetectorIDList( void ) const;
  int GetTofId( void ) const { return TOFid_; }
  int GetLcId( void ) const { return LCid_; }
  int GetDetectorId( const std::string &detName ) const;
  
private:
  static DCGeomMan *geomMan_;
  std::string filename_;
  mutable std::map <int, DCGeomRecord *> geomRecord_;
  int TOFid_, LCid_;

private:
  void clearElements( void );
};


#endif
