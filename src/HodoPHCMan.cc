/*
  HodoPHCMan.cc

  2012/1/24
*/

#include "HodoPHCMan.hh"

#include <iostream>
#include <iomanip>
#include <cstdio>
#include <cmath>
#include <limits>
#include <fstream>
#include <cstdlib>

const double EPS = std::numeric_limits<double>::epsilon();

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

const int BufSize = 144;

HodoPHCParam::HodoPHCParam( int type, int np, double *parlist )
  : Type(type),NPar(np)
{
  if(np>0){
    ParList = new double [NPar];
    for( int i=0; i<NPar; ++i ) ParList[i]=parlist[i];
  }
  else
    ParList=0;
}

HodoPHCParam::~HodoPHCParam()
{
  if(ParList)
    delete [] ParList;
}

double HodoPHCParam::DoPHC( double time, double de )
{
  static const std::string funcname = "[HodoPHCParam::DoPHC]";
  double retval=time;

  switch(Type){
  case 0:
    retval=time; break;
  case 1:
    retval=type1Correction(time,de); break;
  case 2:
    retval=type2Correction(time,de); break; // fiber
  default:
    std::cerr << funcname << ": No Correction Method. Type=" 
              << Type << std::endl;
  }
  return retval;
}

double HodoPHCParam::DoRPHC( double time, double de )
{
  static const std::string funcname = "[HodoPHCParam::DoRPHC]";
  double retval=time;

  switch(Type){
  case 0:
    retval=time; break;
  case 1:
    retval=type1RCorrection(time,de); break;
  default:
    std::cerr << funcname << ": No Correction Method. Type=" 
              << Type << std::endl;
  }
  return retval;
}

double HodoPHCParam::type1Correction( double time, double de )
{
  if(fabs(de-ParList[1])<EPS) de=ParList[1]+EPS;
  return time-ParList[0]/sqrt(fabs(de-ParList[1]))+ParList[2];
}

double HodoPHCParam::type2Correction( double time, double w )
{
  // Correction function for fiber is quadratic function
  return time-(ParList[0]*w*w + ParList[1]*w + ParList[2]);
}

double HodoPHCParam::type1RCorrection( double time, double de )
{
  if(fabs(de-ParList[1])<EPS) de=ParList[1]+EPS;
  return time+ParList[0]/sqrt(fabs(de-ParList[1]))-ParList[2];
}


HodoPHCMan::HodoPHCMan( const std::string & filename )
  : PHCFileName(filename)
{
}

HodoPHCMan::~HodoPHCMan()
{
  clearMap();
}

void HodoPHCMan::clearMap( void )
{
 for(PhcPIterator i=Container.begin(); i!=Container.end(); i++)
    delete i->second;
 Container.clear();
}

bool HodoPHCMan::Initialize( void )
{
  static const std::string funcname = "[HodoPHCMan::Initialize]";

  //  FILE *fp;
  std::ifstream f(PHCFileName.c_str());
  //  char buf[BufSize];

  if(f.fail()){
    std::cerr << funcname << ": file open fail :" <<
      PHCFileName << std::endl;
    exit(-1);
  }

  int cid, plid, seg, ud, type, np;
  double par[10];

  while(f.good() ){
    std::string buf;
    std::getline(f,buf);
    if(buf.empty())
      continue;
    if(buf[0]!='#'){
      if( sscanf(buf.c_str(),
                 "%d %d %d %d %d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
                 &cid, &plid, &seg, &ud, &type, &np, &par[0], &par[1], &par[2],
                 &par[3], &par[4], &par[5], &par[6], &par[7], &par[8], 
                 &par[9] ) >= 6 ){
        if(np>10) np=10;
        int key=KEY(cid,plid,seg,ud);
        HodoPHCParam *p=new HodoPHCParam(type,np,par);
        if(p) Container[key]=p;
        else{
          std::cerr << funcname << ": new fail." << std::endl;
          std::cerr << " PlId=" << std::setw(2) << plid
                    << " Seg=" << std::setw(3) << seg
                    << " Ud=" << std::setw(1) << ud << std::endl;
        }
      }
      else{
        std::cerr << funcname << ": Invalid format" << std::endl;
        std::cerr << " ===> " << buf << std::endl;
      } /* if( sscanf... ) */ 
    } /* if(buf[0]...) */ 
  } /* while( fgets... ) */

  //  fclose(fp);
  f.close();
  std::cout << funcname << ": Initialization finished." << std::endl;
  return true;
}

bool HodoPHCMan::doCorrection( int cid, int plid, int seg, int ud,
                               double time, double de, double & ctime )
{
  HodoPHCParam *map=GetMap(cid,plid,seg,ud);
  if(!map){ ctime=time; return false; }
  ctime=map->DoPHC(time,de);
  return true;
}

bool HodoPHCMan::doRCorrection( int cid, int plid, int seg, int ud,
                                double time, double de, double & ctime )
{
  HodoPHCParam *map=GetMap(cid,plid,seg,ud);
  if(!map){ ctime=time; return false; }
  ctime=map->DoRPHC(time,de);
  return true;
}

HodoPHCParam * HodoPHCMan::GetMap( int cid, int plid, int seg, int ud )
{
  int key=KEY(cid,plid,seg,ud);
  HodoPHCParam *map=0;
  PhcPIterator i=Container.find(key);
  if( i != Container.end() ) map=i->second;
  return map;
}
