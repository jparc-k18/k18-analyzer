/*
  DCGeomMan.cc

  2012/1/24
*/

#include "DCGeomMan.hh"
#include "DCGeomRecord.hh"

#include <string>
#include <stdexcept>
#include <cstdio>
#include <cstring>
#include <cstdlib>

const int MaxChar = 200;

DCGeomMan *DCGeomMan::geomMan_=0;

DCGeomMan::DCGeomMan()
  : TOFid_(51), LCid_(54)
{}

DCGeomMan::~DCGeomMan()
{}


DCGeomMan & DCGeomMan::GetInstance( void )
{
  if( !geomMan_ ){
    geomMan_ = new DCGeomMan();
  }
  return *geomMan_;
}

double DCGeomMan::GetLocalZ( int lnum ) const
{
  static const std::string funcname = "[DCGeomMan::GetLocalZ(int)]"; 
  DCGeomRecord *pGeo = geomRecord_[lnum];
  if( pGeo ) return pGeo->length_;
  else{
    std::cerr << funcname << ": No record. Layer#=" 
	      << lnum << std::endl;
    throw std::out_of_range(funcname+": No record" );
  }
}

double DCGeomMan::GetResolution( int lnum ) const
{
  static const std::string funcname = "[DCGeomMan::GetResolution(int)]"; 
  DCGeomRecord *pGeo = geomRecord_[lnum];
  if( pGeo ) return pGeo->resol_;
  else{
    std::cerr << funcname << ": No record. Layer#=" 
	      << lnum << std::endl;
    throw std::out_of_range(funcname+": No record" );
  }
}

double DCGeomMan::GetTiltAngle( int lnum ) const
{
  static const std::string funcname = "[DCGeomMan::GetTiltAngle(int)]";
  DCGeomRecord *pGeo = geomRecord_[lnum];
  if( pGeo ) return pGeo->tiltAngle_;
  else{
    std::cerr << funcname << ": No record. Layer#=" 
	      << lnum << std::endl;
    throw std::out_of_range(funcname+": No record" );
  }
}

double DCGeomMan::GetRotAngle1( int lnum ) const
{
  static const std::string funcname = "[DCGeomMan::GetRotAngle1(int)]";
  DCGeomRecord *pGeo = geomRecord_[lnum];
  if( pGeo ) return pGeo->rotAngle1_;
  else{
    std::cerr << funcname << ": No record. Layer#=" 
	      << lnum << std::endl;
    throw std::out_of_range(funcname+": No record" );
  }
}

double DCGeomMan::GetRotAngle2( int lnum ) const
{
  static const std::string funcname = "[DCGeomMan::GetRotAngle2(int)]";
  DCGeomRecord *pGeo = geomRecord_[lnum];
  if( pGeo ) return pGeo->rotAngle2_;
  else{
    std::cerr << funcname << ": No record. Layer#=" 
	      << lnum << std::endl;
    throw std::out_of_range(funcname+": No record" );
  }
}

const ThreeVector & DCGeomMan::GetGlobalPosition( int lnum ) const
{
  static const std::string funcname = "[DCGeomMan::GetGlobalPosition(int)]";
  DCGeomRecord *pGeo = geomRecord_[lnum];
  if( pGeo ) return pGeo->pos_;
  else{
    std::cerr << funcname << ": No record. Layer#=" 
	      << lnum << std::endl;
    throw std::out_of_range(funcname+": No record" );
  }
}

ThreeVector DCGeomMan::NormalVector( int lnum ) const
{
  static const std::string funcname = "[DCGeomMan::NormalVector(int)]";
  DCGeomRecord *pGeo = geomRecord_[lnum];
  if( pGeo ) return pGeo->NormalVector();
  else{
    std::cerr << funcname << ": No record. Layer#=" 
	      << lnum << std::endl;
    throw std::out_of_range(funcname+": No record" );
  }
}

ThreeVector DCGeomMan::UnitVector( int lnum ) const
{
  static const std::string funcname = "[DCGeomMan::UnitVector(int)]";
  DCGeomRecord *pGeo = geomRecord_[lnum];
  if( pGeo ) return pGeo->UnitVector();
  else{
    std::cerr << funcname << ": No record. Layer#=" 
	      << lnum << std::endl;
    throw std::out_of_range(funcname+": No record" );
  }
}

const DCGeomRecord * DCGeomMan::GetRecord( int lnum ) const
{
  static const std::string funcname = "[DCGeomMan::GetRecord(int)]";
  DCGeomRecord *pGeo = geomRecord_[lnum];
  if( pGeo ) return pGeo;
  else{
    std::cerr << funcname << ": No record. Layer#=" 
	      << lnum << std::endl;
    throw std::out_of_range(funcname+": No record" );
  }
}
 
  
double DCGeomMan::calcWirePosition( int lnum, double wire ) const
{
  static const std::string funcname = "[DCGeomMan::calcWirePosition()]";
  DCGeomRecord *pGeo = geomRecord_[lnum];
  if( pGeo ){
    return pGeo->WirePos(wire);
  }
  else{
    std::cerr << funcname << ": No record. Layer#=" 
	      << lnum << std::endl;
    throw std::out_of_range(funcname+": No record" );
  }
}

int DCGeomMan::calcWireNumber( int lnum, double pos ) const
{
  static const std::string funcname = "[DCGeomMan::calcWireNumber()]";
  DCGeomRecord *pGeo = geomRecord_[lnum];
  if( pGeo ){
    return pGeo->WireNumber(pos);
  }
  else{
    std::cerr << funcname << ": No record. Layer#=" 
	      << lnum << std::endl;
    throw std::out_of_range(funcname+": No record" );
  }
}
 

void DCGeomMan::clearElements( void )
{
  //  for_each( geomRecord_.begin(), geomRecord_.end(), DeleteObject() );
  std::map <int, DCGeomRecord *>::iterator itr;
  for( itr=geomRecord_.begin(); itr!=geomRecord_.end(); ++itr ){
    delete itr->second;
  }
  geomRecord_.clear();
  TOFid_=51;
}


bool DCGeomMan::Initialize( void )
{
  static const std::string funcname = "[DCGeomMan::Initialize]";
  char str[MaxChar];
  char cname[MaxChar];
  int id;
  double xs, ys, zs, ta, ra1, ra2, l, res, w0, dd, ofs;
  double X_offset = 0.0;
  double Y_offset = 0.0;
  double Z_offset = 0.0;
  //double A = 0.8/750.;

  FILE *fp;

  if( ( fp = fopen( filename_.c_str(), "r" ) ) == 0 ){
    throw std::invalid_argument(funcname+": file open fail");
  }

  clearElements();

  while( fgets( str, MaxChar, fp ) != 0 ){
    if( str[0]!='#' ){
      if( sscanf( str, "%d %s %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
		  &id, cname, &xs, &ys, &zs, &ta, &ra1, &ra2, &l, &res,
		  &w0, &dd, &ofs )
	  == 13 ){
	DCGeomRecord *pRec = 
	  new DCGeomRecord( id, cname, xs, ys, zs, 
			    ta, ra1, ra2, l, res,
			    w0, dd, ofs );
	DCGeomRecord *pOld = geomRecord_[id];
	geomRecord_[id] = pRec;

	if(strcmp(cname,"TOF")==0) {
	  TOFid_=id;
	}
	
	if( pOld ){
	  std::cerr << funcname << ": duplicated id number. "
		    << " following record is deleted." << std::endl;
	  std::cerr << "Id=" << pOld->id_ << " " << pOld->pos_
		    << " ) ... " << std::endl;
	  delete pOld;
	}
      }
      else {
	std::string strtemp=str;
	std::cerr << funcname << ": Invalid format " << strtemp << std::endl;
      }
    }
  }

  fclose(fp);

  std::cout << funcname << " Initialization finished." << std::endl;

  return true;
}


std::vector <int> DCGeomMan::GetDetectorIDList( void ) const
{
  std::vector<int> vlist;
  vlist.reserve(geomRecord_.size());
  std::map <int, DCGeomRecord *>::const_iterator 
    itr=geomRecord_.begin(), end=geomRecord_.end();

  for(; itr!=end; ++itr ){
    vlist.push_back( itr->first );
  }

  return vlist;
}

ThreeVector DCGeomMan::Local2GlobalPos( int lnum, 
					const ThreeVector &in ) const
{
  static const std::string funcname = 
    "[DCGeomMan::Local2GlobalPos(ThreeVecor &)]";

  DCGeomRecord *pGeo = geomRecord_[lnum];
  if( !pGeo ) 
    throw std::out_of_range(funcname+": No record" );

  double x = pGeo->dxds_*in.x() + pGeo->dxdt_*in.y()
    + pGeo->dxdu_*in.z() + pGeo->pos_.x();
  double y = pGeo->dyds_*in.x() + pGeo->dydt_*in.y()
    + pGeo->dydu_*in.z() + pGeo->pos_.y();
  double z = pGeo->dzds_*in.x() + pGeo->dzdt_*in.y()
    + pGeo->dzdu_*in.z() + pGeo->pos_.z();

  return ThreeVector( x, y, z );
}

ThreeVector DCGeomMan::Global2LocalPos( int lnum,
					const ThreeVector &in ) const
{
  static const std::string funcname = 
    "[DCGeomMan::Global2LocalPos(ThreeVecor &)]";

  DCGeomRecord *pGeo = geomRecord_[lnum];
  if( !pGeo ){
    std::cerr << funcname << ": No record. Layer#=" 
	      << lnum << std::endl;
    throw std::out_of_range(funcname+": No record" );
  }

  double x 
    = pGeo->dsdx_*(in.x()-pGeo->pos_.x())
    + pGeo->dsdy_*(in.y()-pGeo->pos_.y())
    + pGeo->dsdz_*(in.z()-pGeo->pos_.z());
  double y 
    = pGeo->dtdx_*(in.x()-pGeo->pos_.x())
    + pGeo->dtdy_*(in.y()-pGeo->pos_.y())
    + pGeo->dtdz_*(in.z()-pGeo->pos_.z());
  double z 
    = pGeo->dudx_*(in.x()-pGeo->pos_.x())
    + pGeo->dudy_*(in.y()-pGeo->pos_.y())
    + pGeo->dudz_*(in.z()-pGeo->pos_.z());

  return ThreeVector( x, y, z );
}

ThreeVector DCGeomMan::Local2GlobalDir( int lnum,
					const ThreeVector &in ) const
{
  static const std::string funcname = 
    "[DCGeomMan::Local2GlobalDir(ThreeVecor &)]";

  DCGeomRecord *pGeo = geomRecord_[lnum];
  if( !pGeo ){
    std::cerr << funcname << ": No record. Layer#=" 
	      << lnum << std::endl;
    throw std::out_of_range(funcname+": No record" );
  }

  double x = pGeo->dxds_*in.x() + pGeo->dxdt_*in.y()
    + pGeo->dxdu_*in.z();
  double y = pGeo->dyds_*in.x() + pGeo->dydt_*in.y()
    + pGeo->dydu_*in.z();
  double z = pGeo->dzds_*in.x() + pGeo->dzdt_*in.y()
    + pGeo->dzdu_*in.z();

  return ThreeVector( x, y, z );
}

ThreeVector DCGeomMan::Global2LocalDir( int lnum,
					const ThreeVector &in ) const
{
  static const std::string funcname = 
    "[DCGeomMan::Global2LocalDir(ThreeVecor &)]";

  DCGeomRecord *pGeo = geomRecord_[lnum];
  if( !pGeo ){
    std::cerr << funcname << ": No record. Layer#=" 
	      << lnum << std::endl;
    throw std::out_of_range(funcname+": No record" );
  }

  double x = pGeo->dsdx_*in.x() + pGeo->dsdy_*in.y()+ pGeo->dsdz_*in.z();
  double y = pGeo->dtdx_*in.x() + pGeo->dtdy_*in.y()+ pGeo->dtdz_*in.z();
  double z = pGeo->dudx_*in.x() + pGeo->dudy_*in.y()+ pGeo->dudz_*in.z();

  return ThreeVector( x, y, z );
}

void DCGeomMan::SetResolution( int lnum,  double res) const
{
  static const std::string funcname = "[DCGeomMan::SetResolution(int, double)]"; 
  DCGeomRecord *pGeo = geomRecord_[lnum];
  if( pGeo ) {
    pGeo->resol_=res;
    return;
  }
  else{
    std::cerr << funcname << ": No record. Layer#=" 
              << lnum << std::endl;
    throw std::out_of_range(funcname+": No record" );
  }
}

int DCGeomMan::GetDetectorId( const std::string &detName ) const
{
  const std::string funcName = "DCGeomMan::GetDetectorId";

  std::map <int, DCGeomRecord *>::const_iterator 
    itr=geomRecord_.begin(), end=geomRecord_.end();

  for(; itr!=end; ++itr ){
    if (itr->second->name_ == detName)
      return itr->second->id_;
  }

  std::cerr << funcName << " : No such detector " << detName << std::endl;
  exit(-1);
}
