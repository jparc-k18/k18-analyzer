/*
  SksFieldMap.hh

  2012/1/24
*/

#ifndef SksFieldMap_h
#define SksFieldMap_h 1

#include <string>
#include <vector>

class SksFieldMap
{
public:
  SksFieldMap( const char *filename=0 );
  ~SksFieldMap();

private:
  SksFieldMap( const SksFieldMap & );
  SksFieldMap & operator = ( const SksFieldMap & );

private:
  std::string filename_;

public:
  bool Initialize( void );
  bool GetFieldValue( const double pointCM[3], double *BfieldTesla ) const;
private:
  struct FD {
    float x, y, z;
  };

  typedef std::vector < std::vector < std::vector < FD > > > FDContainer;
  FDContainer B;
  double X0, Y0, Z0, dX, dY, dZ;
  int Nx, Ny, Nz;

  void cleanupMap( void );
};

#endif
