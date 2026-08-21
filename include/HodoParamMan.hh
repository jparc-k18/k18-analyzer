/*
  HodoParamMan.hh

  2012/1/24
*/

#ifndef HODO_PARAM_MAN_H
#define HODO_PARAM_MAN_H

#include <string>
#include <map>

//Hodo TDC to Time
class HodoTParam
{
public:
  inline HodoTParam( double offset, double gain )
  : Offset(offset), Gain(gain)
  {}

  ~HodoTParam() {}
private:
  HodoTParam();
  HodoTParam( const HodoTParam & );
  HodoTParam & operator = ( const HodoTParam & );

private:
  double Offset, Gain;
public:
  double offset( void ) const { return Offset; }
  double gain( void )   const { return Gain; }
  double time( int tdc ) const
  { return ((double)tdc-Offset)*Gain; }
  inline int tdc( double time ) const
  { return (int)(time/Gain+Offset); }
};

//Hodo ADC to DeltaE
class HodoAParam
{
public:
  inline HodoAParam( double pedestal, double gain )
  : Pedestal(pedestal), Gain(gain)
  {}
  ~HodoAParam() {}
private:
  HodoAParam();
  HodoAParam( const HodoAParam & );
  HodoAParam & operator = ( const HodoAParam & );
private:
  double Pedestal, Gain;
public:
  double pedestal( void ) const { return Pedestal; }
  double gain( void )     const { return Gain; }
  double de( int adc ) const
  { return ((double)adc-Pedestal)/(Gain-Pedestal); }
  int adc( double de ) const
  { return (int)(Gain*de+Pedestal*(1.-de)); }
};

//HodoParam Main Class
class HodoParamMan
{
private:
  std::string ParamFileName;
public:
  explicit HodoParamMan( const std::string & filename );
  ~HodoParamMan();
private:
  HodoParamMan();
  HodoParamMan( const HodoParamMan & );
  HodoParamMan & operator = ( const HodoParamMan & );
public:
  void SetFileName( const std::string & filename ) { ParamFileName=filename; } 

  typedef std::map <int, HodoTParam *> TContainer;
  typedef std::map <int, HodoAParam *> AContainer;
  typedef std::map <int, HodoTParam *>::const_iterator TIterator;
  typedef std::map <int, HodoAParam *>::const_iterator AIterator;

private:
  TContainer TPContainer;
  AContainer APContainer;

public:
  bool Initialize( void );
  bool GetTime( int cid, int plid, int seg, int ud, int tdc, double &time );
  bool GetDe( int cid, int plid, int seg, int ud, int adc, double &de );
  bool GetTdc(int cid, int plid, int seg, int ud, double time, int &tdc ); 
  bool GetAdc(int cid, int plid, int seg, int ud, double de, int &adc ); 

private:
  HodoTParam *GetTmap( int cid, int plid, int seg, int ud );
  HodoAParam *GetAmap( int cid, int plid, int seg, int ud );

  void clearACont( void );
  void clearTCont( void );
};

#endif
