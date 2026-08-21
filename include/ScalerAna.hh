/*
  ScalerAna.hh

  2012/1/24
*/

#ifndef ScalerAna_h
#define ScalerAna_h 1

#include <string>
#include "DaqForm.hh"

class ScalerAna
{
private:
  std::string ScalerDefinitionFileName_;
public:
  ScalerAna()
    : ScalerDefinitionFileName_()
  { ResetScaler(); }
  ScalerAna( const std::string & deffilename )
    : ScalerDefinitionFileName_(deffilename)
  { ResetScaler(); }
  ~ScalerAna() {}
private:
  ScalerAna( const ScalerAna & );
  ScalerAna & operator = ( const ScalerAna & );
  
private:
  int Counts_;
  std::string Name_[NumOfScalers];
  double Sum_[NumOfScalers];
  double Mean_[NumOfScalers];

public:
  bool Initialize( void );
  void ResetScaler( void );
  bool Processing( const unsigned int *buf );
  void PrintFile( const std::string & filename ) const;
  void PrintFile( int runnum ) const;
};

#endif
