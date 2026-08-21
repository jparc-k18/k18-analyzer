/*
  HodoParamMan.cc

  2012/1/24
*/

#include "HodoParamMan.hh"

#include <iostream>
#include <iomanip>
#include <cstdio>
#include <fstream>
#include <cstdlib>

HodoParamMan::HodoParamMan( const std::string & filename )
  : ParamFileName(filename)
{}

HodoParamMan::~HodoParamMan()
{
  clearACont(); clearTCont();
}

void HodoParamMan::clearACont( void )
{
  for( AIterator i=APContainer.begin(); i!=APContainer.end(); ++i )
    delete i->second;
  APContainer.clear();
}

void HodoParamMan::clearTCont( void )
{
  for( TIterator i=TPContainer.begin(); i!=TPContainer.end(); ++i )
    delete i->second;
  TPContainer.clear();
}

const int SegMask  = 0x03FF;
const int CidMask  = 0x00FF;
const int PlidMask = 0x00FF;
const int UdMask   = 0x0003;

const int SegShift  =  0;
const int CidShift  = 11;
const int PlidShift = 19;
const int UdShift   = 27;

inline int KEY( int cid, int pl, int seg, int ud )
{
  return (((cid&CidMask)<<CidShift) | ((pl&PlidMask)<<PlidShift)
          | ((seg&SegMask)<<SegShift) | ((ud&UdMask)<<UdShift) );
}

//const int BufSize = 144;

bool HodoParamMan::Initialize( void )
{
  static const std::string funcname = "[HodoParamMan::Initialize()]";

  clearACont(); clearTCont();

  std::ifstream f(ParamFileName.c_str());
  //  char buf[BufSize];
  //  std::string buf;
  if(f.fail())
    {
       std::cerr << funcname << ": file open fail" << std::endl;
       exit(-1);
    }

  int cid, plid, seg, at, ud;
  double p0, p1;
  int invalid=0;
  while( f.good() )
    {
      std::string buf;
      std::getline(f,buf);
      ++invalid;
      if(buf.empty())
	continue;
      if(buf[0]!='#'){
	if( sscanf(buf.c_str(),"%d %d %d %d %d %lf %lf",&cid,&plid,
		   &seg,&at,&ud,&p0,&p1)==7 ){
        int key=KEY(cid,plid,seg,ud);
        if(at==0){ /* ADC */
          HodoAParam *pa=new HodoAParam(p0,p1);
          if(pa) APContainer[key]=pa;
          else{
            std::cerr << funcname << ": new fail (A)"
                      << " CId=" << std::setw(3) << cid
                      << " PlId=" << std::setw(2) << plid
                      << " Seg=" << std::setw(3) << seg
                      << " Ud=" << std::setw(1) << ud << std::endl;
          }
        }
        else if(at==1){ /* TDC */
          HodoTParam *ta=new HodoTParam(p0,p1);
	  if(ta){
	    TPContainer[key]=ta;
	  }
          else{ 
            std::cerr << funcname << ": new fail (T)"
                      << " CId=" << std::setw(3) << cid
                      << " PlId=" << std::setw(2) << plid
                      << " Seg=" << std::setw(3) << seg
                      << " Ud=" << std::setw(1) << ud << std::endl;
          }
        }
        else{
          std::cerr << funcname << ": Invalid Input" << std::endl;
          std::cerr << " ===> " << std::string(buf) <<" "
		    << invalid <<"a"<< std::endl;
        } /* if(at...) */
      }
      else {
        std::cerr << funcname << ": Invalid Input" << std::endl;
        std::cerr << " ===> " << std::string(buf) <<" "
		  << invalid <<"b"<< std::endl;
      } /* if( sscanf... ) */ 
    } /* if(str[0]...) */   
  } /* while( fgets...) */

  f.close();
  std::cout << funcname << ": Initialization finished" << std::endl;

  return true;
}

bool HodoParamMan::GetTime( int cid, int plid, int seg, int ud,
                            int tdc, double &time )
{
  HodoTParam *map=GetTmap(cid,plid,seg,ud);
  if(!map) return false;
  time=map->time(tdc);
  return true;
}

bool HodoParamMan::GetDe( int cid, int plid, int seg, int ud, 
                          int adc, double &de )
{
  HodoAParam *map=GetAmap(cid,plid,seg,ud);
  if(!map) return false;
  de=map->de(adc);
  return true;
}

HodoTParam * HodoParamMan::GetTmap( int cid, int plid, int seg, int ud )
{
  int key=KEY(cid,plid,seg,ud);
  HodoTParam *map=0;
  TIterator i=TPContainer.find(key);
  if( i!=TPContainer.end() ) map=i->second;
  return map;
}

HodoAParam * HodoParamMan::GetAmap( int cid, int plid, int seg, int ud )
{
  int key=KEY(cid,plid,seg,ud);
  HodoAParam *map=0;
  AIterator i=APContainer.find(key);
  if( i!=APContainer.end() ) map=i->second;
  return map;
}
