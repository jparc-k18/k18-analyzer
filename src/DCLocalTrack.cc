/*
  DCLocalTrack.cc

  2012/1/24
*/

#include "DCLocalTrack.hh"
#include "DCLTrackHit.hh"
#include "DCGeomMan.hh"
#include "DCAnalyzer.hh"
#include "DetectorID.hh"

#include "MathTools.hh"

#include <string>
#include <vector>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <cstring>
#include <stdexcept>
#include <sstream>

const double Deg2Rad = acos(-1.)/180.;
const double Rad2Deg = 180./acos(-1.);

const int ReservedNumOfHits  = 16;
const int DCLocalMinNHits    = 6;
const int DCLocalMinNHitsVXU = 3;
DCLocalTrack::DCLocalTrack()
  : status_(false), x0_(0.0), y0_(0.0), u0_(0.0), v0_(0.0),a_(0.0),b_(0.0),
    gftstatus_(true)
{
  hitArray.reserve( ReservedNumOfHits );
}

DCLocalTrack::~DCLocalTrack()
{
}

DCLTrackHit *DCLocalTrack::GetHit( std::size_t nth ) const
{
  if( nth<hitArray.size() )
    return hitArray[nth];
  else
    return 0;
}

DCLTrackHit *DCLocalTrack::GetHitOfLayerNumber( int lnum ) const
{
  for( std::size_t i=0; i<hitArray.size(); ++i )
    if( hitArray[i]->GetLayer()==lnum )
      return hitArray[i];
  return 0;
}

bool DCLocalTrack::DoFit( void )
{
  const std::string funcname = "[DCLocalTrack::DoFit()]";
  const DCGeomMan & geomMan=DCGeomMan::GetInstance();
  //std::cout<<"###here000"<<std::endl;
  
  std::size_t n = hitArray.size();
  //std::cout<<"#####hitarray size:    "<< n <<std::endl;
  if(n < DCLocalMinNHits ) return status_ = false;
  //std::cout<<"###here111"<<std::endl;
  std::vector <double> z, w, s, ct, st;
  z.reserve(n); w.reserve(n); s.reserve(n);
  ct.reserve(n); st.reserve(n);
  
  for( std::size_t i=0; i<n; ++i ){
    DCLTrackHit *hitp = hitArray[i];
    if( hitp ){
      int lnum = hitp->GetLayer();
      double ww = geomMan.GetResolution( lnum );
      double zz = geomMan.GetLocalZ( lnum );
      double aa = hitp->GetTiltAngle()*Deg2Rad;

      z.push_back( zz ); w.push_back( 1./(ww*ww) ); 
      s.push_back( hitp->GetLocalHitPos() );
      ct.push_back( cos(aa) ); st.push_back( sin(aa) );

#if 0
      std::cout << std::setw(10) << "layer = " << lnum 
		<< std::setw(10) << "wire  = " << hitp->GetWire() << " "
		<< std::setw(20) << "WirePosition = "<<hitp->GetWirePosition() << " "
		<< std::setw(20) << "DriftLength = "<<hitp->GetDriftLength() << " "
		<< std::setw(20) << "hit position = "<<hitp->GetLocalHitPos()<< " " 
		<< std::endl;
#endif
    }
  }
  //std::cout<<"###here111"<<std::endl;
 
  std::size_t nn = z.size();
  double matrx[16], *mtp[4], fitp[4];
  mtp[0]=&matrx[0]; mtp[1]=&matrx[4]; mtp[2]=&matrx[8]; mtp[3]=&matrx[12];

  for( int i=0; i<4; ++i ){
    fitp[i]=0.0;
    for( int j=0; j<4; ++j ){
      mtp[i][j]=0.0;
    }
  }

  for( std::size_t i=0; i<nn; ++i ){
    double ww=w[i], zz=z[i], ss=s[i], ctt=ct[i], stt=st[i];
    mtp[0][0] += ww*ctt*ctt;
    mtp[0][1] += ww*zz*ctt*ctt;
    mtp[0][2] += ww*ctt*stt;
    mtp[0][3] += ww*zz*ctt*stt;
    mtp[1][1] += ww*zz*zz*ctt*ctt;
    mtp[1][2] += ww*zz*ctt*stt;
    mtp[1][3] += ww*zz*zz*ctt*stt;
    mtp[2][2] += ww*stt*stt;
    mtp[2][3] += ww*zz*stt*stt;
    mtp[3][3] += ww*zz*zz*stt*stt;

    fitp[0] += ww*ss*ctt;
    fitp[1] += ww*zz*ss*ctt;
    fitp[2] += ww*ss*stt;
    fitp[3] += ww*zz*ss*stt;
  }
  mtp[1][0]=mtp[0][1]; mtp[2][0]=mtp[0][2]; mtp[3][0]=mtp[0][3];
  mtp[2][1]=mtp[1][2]; mtp[3][1]=mtp[1][3]; mtp[3][2]=mtp[2][3];

  std::vector<int> indxc(nn), indxd(nn), ipiv(nn);

#if 0
  std::cout<<"             Vector:  Q_i               "<<std::endl;
  std::cout<<"   A1=  "<<fitp[0]<<"     A2=  "<<fitp[1]<<"     A3=  "<<fitp[2]
	   <<"     A4=  "<<fitp[3]<<std::endl;

  std::cout<<"             original matrix: M_ij        "<<std::endl;
  std::cout<<"    A11="<<mtp[0][0]<<"    A12="<<mtp[0][1]<<"    A13="<<mtp[0][2]
      	   <<"    A14="<<mtp[0][3]<<std::endl;
  std::cout<<"    A21="<<mtp[1][0]<<"    A22="<<mtp[1][1]<<"    A23="<<mtp[1][2]
	   <<"    A24="<<mtp[1][3]<<std::endl;
  std::cout<<"    A31="<<mtp[2][0]<<"    A32="<<mtp[2][1]<<"    A33="<<mtp[2][2]
	   <<"    A34="<<mtp[2][3]<<std::endl;
  std::cout<<"    A41="<<mtp[3][0]<<"    A42="<<mtp[3][1]<<"    A43="<<mtp[3][2]
	   <<"    A44="<<mtp[3][3]<<std::endl;
#endif

  double Org[4][4]={0},Red[4][4]={0},Final[4][4]={0};
  double Org_vec[4]={0}, Solution_vec[4]={0};
  for(int l=0; l<4;l++){
    for(int m=0; m<4; m++){
      Org[l][m]=mtp[l][m];
    }
  }

  for (int i=0; i<4; i++){
    Org_vec[i]=fitp[i];
  }

  if( MathTools::GaussJordan(mtp,4,fitp,&indxc[0],
			     &indxd[0],&ipiv[0])==false ){
    std::cerr << funcname << ": Fitting fails" << std::endl;
    return status_=false;
  }
  x0_=fitp[0]; y0_=fitp[2]; u0_=fitp[1]; v0_=fitp[3];

#if 0
  std::cout<<"             reduced matrix        "<<std::endl;
  std::cout<<"     A11="<<mtp[0][0]<<"    A12="<<mtp[0][1]<<"    A13="<<mtp[0][2]
	   <<"     A14="<<mtp[0][3]<<std::endl;
  std::cout<<"     A21="<<mtp[1][0]<<"    A22="<<mtp[1][1]<<"    A23="<<mtp[1][2]
	   <<"     A24="<<mtp[1][3]<<std::endl;
  std::cout<<"     A31="<<mtp[2][0]<<"    A32="<<mtp[2][1]<<"    A33="<<mtp[2][2]
	   <<"     A34="<<mtp[2][3]<<std::endl;
  std::cout<<"     A41="<<mtp[3][0]<<"    A42="<<mtp[3][1]<<"    A43="<<mtp[3][2]
	   <<"     A44="<<mtp[3][3]<<std::endl;
#endif

#if 0
  std::cout<<"             Solution from Gauss-Jordan             "<<std::endl;
  std::cout<<"     A1="<<fitp[0]<<"    A2="<<fitp[1]<<"    A3="<<fitp[2]
	   <<"     A4="<<fitp[3]<<std::endl;
#endif


  for(int l=0; l<4;l++){
    for(int m=0; m<4; m++){
      Red[l][m]=mtp[l][m];
    }
  }

  for(int i=0;i<4;i++){
    for(int j=0;j<4;j++){
      Solution_vec[i] +=Red[i][j]*Org_vec[j];
    }
  }

#if 0
  std::cout<<"             Solution from explicite calculation       "<<std::endl;
  std::cout<<"    A1="<<Solution_vec[0]<<"    A2="<<Solution_vec[1]<<"    A3="<<Solution_vec[2]
	   <<"    A4="<<Solution_vec[3]<<std::endl;
#endif

  for(int i=0; i<4; i++){
    for(int j=0; j<4; j++){
      for(int k=0; k<4; k++){
	Final[i][j] += Red[i][k]*Org[k][j];
      }
    }
  }
  for(int i=0; i<4; i++){
    for(int j=0; j<4; j++){
      if(Final[i][j]<1.0e-10) Final[i][j]=0.0;
    }
  }

#if 0
  std::cout<< ""<<std::endl;
  std::cout<<"            final matrix        "<<std::endl;
  std::cout<<"       A11="<<std::setw(10)<<Final[0][0]<<"       A12="<<std::setw(10)<<Final[0][1]
	   <<"       A13="<<std::setw(10)<<Final[0][2]<<"       A14="<<std::setw(10)<<Final[0][3]
	   <<std::endl;
  std::cout<<"       A21="<<std::setw(10)<<Final[1][0]<<"       A22="<<std::setw(10)<<Final[1][1]
	   <<"       A23="<<std::setw(10)<<Final[1][2]<<"       A24="<<std::setw(10)<<Final[1][3]
	   <<std::endl;
  std::cout<<"       A31="<<std::setw(10)<<Final[2][0]<<"       A32="<<std::setw(10)<<Final[2][1]
	   <<"       A33="<<std::setw(10)<<Final[2][2]<<"       A34="<<std::setw(10)<<Final[2][3]
	   <<std::endl;
  std::cout<<"       A41="<<std::setw(10)<<Final[3][0]<<"       A42="<<std::setw(10)<<Final[3][1]
	   <<"       A43="<<std::setw(10)<<Final[3][2]<<"       A44="<<std::setw(10)<<Final[3][3]
	   <<std::endl;

#endif

  double chisqr=0.0;
  for( std::size_t i=0; i<nn; ++i ){
    double ww=w[i], zz=z[i];
    double scal=GetX(zz)*ct[i]+GetY(zz)*st[i];
    chisqr += ww*(s[i]-scal)*(s[i]-scal);

#if 0
    if(1){
      //      std::cout<<" !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! "<<std::endl;
      //     std::cout<<"x coordinate = "<<(x0_+u0_*zz)<<std::endl;
      //     std::cout<<"y coordinate = "<<(y0_+v0_*zz)<<std::endl;
      std::cout<<std::setw(10)<<"layer = "<<i<<
	std::setw(10)<<"scal = "<<scal<<
	std::setw(10)<<"sdata = "<<s[i]<<std::endl;
      std::cout<<std::setw(10)<<"Res = "<<s[i]-scal<<std::endl;
      std::cout<<std::setw(10)<<"X = "<< GetX(zz)<<" Y = "<< GetY(zz)<<std::endl;
      std::cout<<std::setw(10)<<"chisqr = "<<chisqr<<std::endl;
    }
#endif
  }
  chisqr /= nn-4.;
  /*  
  if(chisqr<2){
    std::cout<<" !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! "<<std::endl;
    std::cout << "chisqr = " << chisqr << " nn-4 = " << nn-4 << std::endl;
    std::cout<<" !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! "<<std::endl;
  }
  */
  chisqr_=chisqr;
  for( std::size_t i=0; i<nn; ++i ){
    DCLTrackHit *hitp = hitArray[i];
    if( hitp ){
      int lnum = hitp->GetLayer();
      double zz = DCGeomMan::GetInstance().GetLocalZ( lnum );
      /*  
	  if(chisqr<2){
	  std::cout<<std::setw(10)<<"lnum = "<< lnum <<std::endl;
	  std::cout<<std::setw(10)<<"X = "<< GetX(zz)<<" Y = "<< GetY(zz)<<std::endl;
	  }
      */
      hitp->SetCalPosition( GetX(zz), GetY(zz) );
      
    }
  }

  // std::cout << "***********************************************************" << std::endl;

  return status_=true;
}

bool DCLocalTrack::DoFitBcSdc( void )
{
  const std::string funcname = "[DCLocalTrack::DoFitBcSdc()]";
  const DCGeomMan & geomMan=DCGeomMan::GetInstance();

  std::size_t n = hitArray.size();

  if( n < DCLocalMinNHits ) return status_ = false;

  std::vector <double> z, w, s, ct, st;
  z.reserve(n); w.reserve(n); s.reserve(n);
  ct.reserve(n); st.reserve(n);

  // BcOut local z -> SKS-frame z using the SURVEYED positions:
  // the target sits at GetLocalZ(130) in the BcOut frame and at
  // GetLocalZ(0) in the SKS frame (the previous code ignored the
  // latter, -145.3 mm, displacing the whole BcOut arm)
  double zK18tgt=geomMan.GetLocalZ( 130 )-geomMan.GetLocalZ( 0 );
  for( std::size_t i=0; i<n; ++i ){
    DCLTrackHit *hitp = hitArray[i];
    if( hitp ){
      int lnum = hitp->GetLayer();
      double ww = geomMan.GetResolution( lnum );
      double zz = geomMan.GetLocalZ( lnum );

      if (lnum>=113&&lnum<=124) //BcOut
        zz -= zK18tgt;

      double aa = hitp->GetTiltAngle()*Deg2Rad;

      z.push_back( zz ); w.push_back( 1./(ww*ww) ); 
      s.push_back( hitp->GetLocalHitPos() );
      ct.push_back( cos(aa) ); st.push_back( sin(aa) );

#if 0
      std::cout << std::setw(2) << i << " (" << lnum << ") " 
                << std::setw(3) << hitp->GetWire() << " "
                << hitp->GetLocalHitPos() << " "
                << zz << " " << aa 
                << std::endl;
#endif
    }
  }
  std::size_t nn = z.size();

  double matrx[16], *mtp[4], fitp[4];
  mtp[0]=&matrx[0]; mtp[1]=&matrx[4]; mtp[2]=&matrx[8]; mtp[3]=&matrx[12];

  for( int i=0; i<4; ++i ){
    fitp[i]=0.0;
    for( int j=0; j<4; ++j ){
      mtp[i][j]=0.0;
    }
  }

  for( std::size_t i=0; i<nn; ++i ){
    double ww=w[i], zz=z[i], ss=s[i], ctt=ct[i], stt=st[i];
    mtp[0][0] += ww*ctt*ctt;
    mtp[0][1] += ww*zz*ctt*ctt;
    mtp[0][2] += ww*ctt*stt;
    mtp[0][3] += ww*zz*ctt*stt;
    mtp[1][1] += ww*zz*zz*ctt*ctt;
    mtp[1][2] += ww*zz*ctt*stt;
    mtp[1][3] += ww*zz*zz*ctt*stt;
    mtp[2][2] += ww*stt*stt;
    mtp[2][3] += ww*zz*stt*stt;
    mtp[3][3] += ww*zz*zz*stt*stt;

    fitp[0] += ww*ss*ctt;
    fitp[1] += ww*zz*ss*ctt;
    fitp[2] += ww*ss*stt;
    fitp[3] += ww*zz*ss*stt;
  }
  mtp[1][0]=mtp[0][1]; mtp[2][0]=mtp[0][2]; mtp[3][0]=mtp[0][3];
  mtp[2][1]=mtp[1][2]; mtp[3][1]=mtp[1][3]; mtp[3][2]=mtp[2][3];

  std::vector<int> indxc(nn), indxd(nn), ipiv(nn);

  if( MathTools::GaussJordan(mtp,4,fitp,&indxc[0],
                             &indxd[0],&ipiv[0])==false ){
    std::cerr << funcname << ": Fitting fails" << std::endl;
    return status_=false;
  }
  x0_=fitp[0]; y0_=fitp[2]; u0_=fitp[1]; v0_=fitp[3];

  double chisqr=0.0;
  for( std::size_t i=0; i<nn; ++i ){
    double ww=w[i], zz=z[i];
    double scal=GetX(zz)*ct[i]+GetY(zz)*st[i];
    chisqr += ww*(s[i]-scal)*(s[i]-scal);
  }
  chisqr /= nn-4.;
  chisqr_=chisqr;
  for( std::size_t i=0; i<n; ++i ){
    DCLTrackHit *hitp = hitArray[i];
    if( hitp ){
      int lnum = hitp->GetLayer();
      double zz = DCGeomMan::GetInstance().GetLocalZ( lnum );
      if (lnum>=113&&lnum<=124) //BcOut
        zz -= zK18tgt;

      hitp->SetCalPosition( GetX(zz), GetY(zz) );
    }
  }

  return status_=true;
}

///////XUV Tracking
bool DCLocalTrack::DoFitVXU( void )
{
  const std::string funcname = "[DCLocalTrack::DoFitVXU()]";
  const DCGeomMan & geomMan=DCGeomMan::GetInstance();

  std::size_t n = hitArray.size();

  if(n < DCLocalMinNHitsVXU ) return status_ = false;

  double a=0,b=0;// <-Add !!
  double w[n+1],z[n+1],x[n+1];

  for( std::size_t i=0; i<n; ++i ){
    DCLTrackHit *hitp = hitArray[i];
    if( hitp ){
      int lnum = hitp->GetLayer();
      w[i] = geomMan.GetResolution( lnum );
      z[i] = geomMan.GetLocalZ( lnum );
      x[i] = hitp->GetLocalHitPos();

#if 0
      std::cout << "" << std::endl;
      std::cout << "**********" << std::endl;
      std::cout << std::setw(10) << "layer = " << lnum 
		<< std::setw(10) << "wire  = " << hitp->GetWire() << " "
		<< std::setw(20) << "WirePosition = "<<hitp->GetWirePosition() << " "
		<< std::setw(20) << "DriftLength = "<<hitp->GetDriftLength() << " "
		<< std::setw(20) << "hit position = "<<hitp->GetLocalHitPos()<< " " 
		<< std::endl;
      std::cout << "**********" << std::endl;
      std::cout << "" << std::endl;
#endif
    }
  }

  double A=0, B=0, C=0, D=0, E=0;// <-Add!! 
  for(int i=0; i<n; ++i){

    A += z[i]/(w[i]*w[i]);
    B += 1/(w[i]*w[i]);
    C += x[i]/(w[i]*w[i]);
    D += z[i]*z[i]/(w[i]*w[i]);
    E += x[i]*z[i]/(w[i]*w[i]);
  }

  a_ = (E*B-C*A)/(D*B-A*A);
  b_ = (D*C-E*A)/(D*B-A*A);

  double chisqr = 0.0;

  for(int i=0; i<n; ++i){
    chisqr += (x[i]-a_*z[i]-b_)*(x[i]-a_*z[i]-b_)/(w[i]*w[i]);
  }

  chisqr /= n-2.;
  chisqr_=chisqr;

  //std::cout << "chisqr = " << chisqr << std::endl;
  //std::cout << "end "  << std::endl;
  //std::cout << " "  << std::endl;
  
  for( std::size_t i=0; i<n; ++i ){
    DCLTrackHit *hitp = hitArray[i];
    if( hitp ){
      int lnum = hitp->GetLayer();
      double zz = DCGeomMan::GetInstance().GetLocalZ( lnum );
      /*  
	  if(chisqr<2){
	  std::cout<<std::setw(10)<<"lnum = "<< lnum <<std::endl;
	  std::cout<<std::setw(10)<<"X = "<< GetX(zz)<<" Y = "<< GetY(zz)<<std::endl;
	  }
      */
      hitp->SetLocalCalPosVXU( a_*zz+b_ );
      
    }
  }

  return status_=true;
}

bool DCLocalTrack::ReCalc( bool applyRecursively )
{
  static const std::string funcname = "[DCLocalTrack::ReCalc]";

  std::size_t n = hitArray.size();
  for( std::size_t i=0; i<n; ++i ){
    DCLTrackHit *hitp = hitArray[i];
    if( hitp ) hitp->ReCalc( applyRecursively );
  }
  
  bool ret=DoFit();
  if( !ret ){
    std::cerr << funcname << ": Recalculation fails" << std::endl;
  }
  return ret;
}
