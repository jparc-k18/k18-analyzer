/*
  SksTrack.cc

  2012/1/24
*/

#include "SksTrack.hh"

#include "DCLocalTrack.hh"
#include "TrackHit.hh"
#include "DCGeomMan.hh"
#include "MathTools.hh"
#include "DCAnalyzer.hh"

#include <iostream>
#include <iomanip>
#include <cstring>
#include <string>
#include <stdexcept>
#include <sstream>

#define WARNOUT 1

const double DefPini = 1.70;
const int MaxIteraction = 100;
const double InitialChiSqr = 1.E+10; 
const double MaxChiSqr = 10000.;
const double MinDeltaChiSqrR = 0.0002;

SksTrack::SksTrack( DCLocalTrack *In, DCLocalTrack *Out )
  : In_(In), Out_(Out), iniP_(DefPini), 
    status_(false), nIteration_(-1), chisqr_(InitialChiSqr),
    gfastatus_(true)
{
  fillHitArray();
}
 
SksTrack::~SksTrack()
{
  clearHitArray();
}

TrackHit *SksTrack::GetHit( std::size_t nth ) const
{
  if( nth<hitArray_.size() ){
    return hitArray_[nth];
  }
  else
    return 0;
}

TrackHit *SksTrack::GetHitOfLayerNumber( int lnum ) const
{
  for( std::size_t i=0; i<hitArray_.size(); ++i )
    if( hitArray_[i]->GetLayer()==lnum )
      return hitArray_[i];
  return 0;
}

void SksTrack::fillHitArray( void )
{
  std::size_t nIn = In_->GetNHit();
  std::size_t nOut = Out_->GetNHit();
  clearHitArray();
  hitArray_.reserve( nIn+nOut );

  for( std::size_t i=0; i<nIn; ++i ){
    DCLTrackHit *hit = In_->GetHit( i );
    TrackHit *thit = new TrackHit( hit );
    hitArray_.push_back( thit );
  }

  for( std::size_t i=0; i<nOut; ++i ){
    DCLTrackHit *hit = Out_->GetHit( i );
    TrackHit *thit = new TrackHit( hit );
    hitArray_.push_back( thit );
  }

}

void SksTrack::clearHitArray( void )
{
  int nh=hitArray_.size();
  for( int i=nh-1; i>=0; --i )
    delete hitArray_[i];

  hitArray_.clear();
}

bool SksTrack::ReCalc( bool applyRecursively )
{
  static const std::string funcname = "[SksTrack::ReCalc]";

  if( applyRecursively ){
    In_->ReCalc(applyRecursively);
    Out_->ReCalc(applyRecursively);

    //    int nh=hitArray_.size();
    //    for( int i=0; i<nh; ++i ){
    //      hitArray_[i]->ReCalc(applyRecursively);
    //    }

  }

  return doFitAgain(CPval_);
  //  return doFit();
}


bool SksTrack::doFit( void )
{
  static const std::string funcname = "[SksTrack::doFit]";

  //  clearHitArray();
  //  fillHitArray();

  int IdTof = DCGeomMan::GetInstance().GetTofId();
  double zTof = DCGeomMan::GetInstance().GetLocalZ( IdTof ); 
  double xTof = Out_->GetX( zTof );
  double yTof = Out_->GetY( zTof );
  double uTof = Out_->GetU0(); 
  double vTof = Out_->GetV0(); 

  ThreeVector gpTof = DCGeomMan::GetInstance().
    Local2GlobalPos( IdTof, ThreeVector( xTof, yTof, 0.0 ) );
  double pz = iniP_/sqrt(1.+uTof*uTof+vTof*vTof);
  ThreeVector gmTof = DCGeomMan::GetInstance().
    Local2GlobalDir( IdTof, ThreeVector( pz*uTof, pz*vTof, pz ) );

#if 0
  {
    std::ios::fmtflags oldFlags = std::cout.flags();
    std::size_t oldPrec = std::cout.precision();
    std::cout.setf( std::ios::fixed );
    std::cout.precision(2);
    std::cout << funcname << ": X (" 
	      << std::setw(9) << gpTof.x() << ","
	      << std::setw(9) << gpTof.y() << ","
	      << std::setw(9) << gpTof.z() << ")";
    std::cout.precision(5);
    std::cout << " P "
	      << std::setw(9) << gmTof.Mag() << " (" 
	      << std::setw(9) << gmTof.x() << ","
	      << std::setw(9) << gmTof.y() << ","
	      << std::setw(9) << gmTof.z() << ")"
	      << std::endl;

    std::cout.flags( oldFlags );
    std::cout.precision( oldPrec );
  }
#endif
  
  RKCordParameter iniCord( gpTof, gmTof );
  RKCordParameter prevCord;

  RKHitPointContainer preHPntCont;

  double chiSqr, prevChiSqr=InitialChiSqr;
  double lambdaCri=0.01, estDChisqr=InitialChiSqr;

  double dmp = 0.; 
  status_ = false;

  hPntCont_ = RKmakeHPContainer(); 
  RKHitPointContainer prevHPntCont;

  int iIter=0, iIterEf=1;

  while( ++iIter < MaxIteraction ){
    if( !RKtrace( iniCord, hPntCont_ ) ){
      // Error 
#ifdef WARNOUT 
      //std::cerr << funcname << ": Error in RKtrace. " << std::endl;
#endif
      break;
    }
    chiSqr = calcChiSqr( hPntCont_ );
    double dChiSqr = chiSqr-prevChiSqr;
    double dChiSqrR = dChiSqr/prevChiSqr;
    double Rchisqr = dChiSqr/estDChisqr;
#if 0
    { 
      std::ios::fmtflags oldFlags = std::cout.flags();
      std::size_t oldPrec = std::cout.precision();

      std::cout.setf( std::ios::scientific ); 
      std::cout.precision(3);
      std::cout << funcname << ": #"
		<< std::setw(3) << iIter << " ( " 
		<< std::setw(2) << iIterEf << " )"
		<< " chi=" << std::setw(10) << chiSqr;
      std::cout.precision(5);
      std::cout << " (" << std::fixed << std::setw(10) << dChiSqrR << " )"
		<< " [" << std::fixed << std::setw(10) << Rchisqr << " ]";
      std::cout.flags( oldFlags );
      std::cout.setf( std::ios::scientific ); 
      std::cout.precision(2); 
      std::cout << " df=" << std::setw(8) << dmp
		<< " (" << std::setw(8) << lambdaCri << ")" << std::endl;

      std::cout.flags( oldFlags );
      std::cout.precision( oldPrec );
    }
#endif

#if 0
    PrintCalcHits( hPntCont_, std::cout );
#endif
  
    if( fabs(dChiSqrR)<MinDeltaChiSqrR && 
	(chiSqr<MaxChiSqr || Rchisqr>1.) ){
      // Converged
      status_ = true; 
      if( dChiSqr>0. ){
	iniCord = prevCord;
	chiSqr = prevChiSqr;
	hPntCont_ = prevHPntCont;
      }
      break;
    }
  
    // Next Guess
    if( iIter==1 ){
      prevCord = iniCord;
      prevChiSqr = chiSqr;
      prevHPntCont = hPntCont_;
      ++iIterEf;
    }
    else if( dChiSqr <= 0.0 ){
      prevCord = iniCord;
      prevChiSqr = chiSqr;
      prevHPntCont = hPntCont_;
      ++iIterEf;
      if( Rchisqr>=0.75 ){
	dmp*=0.5;
	if( dmp < lambdaCri ) dmp=0.;
      }
      else if( Rchisqr>0.25 ){}
      else{
	if( dmp==0.0 ) dmp=lambdaCri; 
	else dmp*=2.;
      }
    }
    else {
      if( dmp==0.0 ) dmp=lambdaCri;
      else {
	double uf=2.;
	if( 2.-Rchisqr > 2. ) uf=2.-Rchisqr;
	dmp *= uf;
      }
      iniCord = prevCord;
      hPntCont_ = prevHPntCont;
    }
    if( !guessNextParameters( hPntCont_, iniCord, 
			      estDChisqr, lambdaCri, dmp ) )
      return status_=false;

#if 0
    { 
      int ii;
      std::cout << "####:" ;
      std::cin >> ii;
      if(ii<0) return status_;
    }
#endif

  }  /* End of Iteration */
  
  nIteration_ = iIter;
  chisqr_ = chiSqr;

  if( !RKtraceToLast( hPntCont_ ) )
    status_ = false;

  saveCalcPosition( hPntCont_ ); 

  if( !saveTrackParameters( iniCord ) )
    status_ = false;

#if 0
  { 
    std::ios::fmtflags oldFlags = std::cout.flags();
    std::size_t oldPrec = std::cout.precision();
    std::cout.setf( std::ios::fixed );
    std::cout.precision(5);

    std::cout << funcname << ": Status = " << status_
	      << " in " << std::setw(3) << iIter << " ( "
	      << std::setw(2) << iIterEf << " ) Iteractions "
	      << " with ChiSqr=" << std::setw(10) << chisqr_ << std::endl;
    std::cout << " ** TGTlocal ** X ( " << std::setprecision(2)
	      << std::setw(7) << priPos_.x() << ", "
	      << std::setw(7) << priPos_.y() << ", "
	      << std::setw(7) << priPos_.z() << " )"
	      << " P " << std::setprecision(5) 
	      << std::setw(7) << priMom_.Mag() << " ( "
	      << std::setw(7) << priMom_.x() << ", "
	      << std::setw(7) << priMom_.y() << ", "
	      << std::setw(7) << priMom_.z() << " )"
	      << " PL: " << std::setprecision(1) 
	      << std::setw(7) << pLenTtoT_ << " "
	      << std::setw(7) << pLenTot_ << std::endl;

    PrintCalcHits( hPntCont_, std::cout );

    std::cout.flags( oldFlags );
    std::cout.precision( oldPrec );
  }
#endif

    
  return status_;
}

bool SksTrack::doFitAgain( RKCordParameter iniCord )
{
  static const std::string funcname = "[SksTrack::doFitAgain]";

  //  clearHitArray();
  //  fillHitArray();

  RKCordParameter prevCord;

  RKHitPointContainer preHPntCont;

  double chiSqr, prevChiSqr=InitialChiSqr;
  double lambdaCri=0.01, estDChisqr=InitialChiSqr;

  double dmp = 0.; 
  status_ = false;

  hPntCont_ = RKmakeHPContainer(); 
  RKHitPointContainer prevHPntCont;

  int iIter=0, iIterEf=1;

  while( ++iIter < MaxIteraction ){
    if( !RKtrace( iniCord, hPntCont_ ) ){
      // Error 
#ifdef WARNOUT 
      //std::cerr << funcname << ": Error in RKtrace. " << std::endl;
#endif
      break;
    }
    chiSqr = calcChiSqr( hPntCont_ );
    double dChiSqr = chiSqr-prevChiSqr;
    double dChiSqrR = dChiSqr/prevChiSqr;
    double Rchisqr = dChiSqr/estDChisqr;
#if 0
    { 
      std::ios::fmtflags oldFlags = std::cout.flags();
      std::size_t oldPrec = std::cout.precision();

      std::cout.setf( std::ios::scientific ); 
      std::cout.precision(3);
      std::cout << funcname << ": #"
		<< std::setw(3) << iIter << " ( " 
		<< std::setw(2) << iIterEf << " )"
		<< " chi=" << std::setw(10) << chiSqr;
      std::cout.precision(5);
      std::cout << " (" << std::fixed << std::setw(10) << dChiSqrR << " )"
		<< " [" << std::fixed << std::setw(10) << Rchisqr << " ]";
      std::cout.flags( oldFlags );
      std::cout.setf( std::ios::scientific ); 
      std::cout.precision(2); 
      std::cout << " df=" << std::setw(8) << dmp
		<< " (" << std::setw(8) << lambdaCri << ")" << std::endl;

      std::cout.flags( oldFlags );
      std::cout.precision( oldPrec );
    }
#endif

#if 0
    PrintCalcHits( hPntCont_, std::cout );
#endif

    if( fabs(dChiSqrR)<MinDeltaChiSqrR && 
	(chiSqr<MaxChiSqr || Rchisqr>1.) ){
      // Converged
      status_ = true; 
      if( dChiSqr>0. ){
	iniCord = prevCord;
	chiSqr = prevChiSqr;
	hPntCont_ = prevHPntCont;
      }
      break;
    }

    // Next Guess
    if( iIter==1 ){
      prevCord = iniCord;
      prevChiSqr = chiSqr;
      prevHPntCont = hPntCont_;
      ++iIterEf;
    }
    else if( dChiSqr <= 0.0 ){
      prevCord = iniCord;
      prevChiSqr = chiSqr;
      prevHPntCont = hPntCont_;
      ++iIterEf;
      if( Rchisqr>=0.75){
	dmp*=0.5;
	if( dmp < lambdaCri ) dmp=0.;
      }
      else if( Rchisqr>0.25 ){}
      else{
	if( dmp==0.0 ) dmp=lambdaCri; 
	else dmp*=2.;
      }
    }
    else {
      if( dmp==0.0 ) dmp=lambdaCri;
      else {
	double uf=2.;
	if( 2.-Rchisqr > 2. ) uf=2.-Rchisqr;
	dmp *= uf;
      }
      iniCord = prevCord;
      hPntCont_ = prevHPntCont;
    }
    if( !guessNextParameters( hPntCont_, iniCord, 
			      estDChisqr, lambdaCri, dmp ) )
      return status_=false;

#if 0
    { 
      int ii;
      std::cout << "####:" ;
      std::cin >> ii;
      if(ii<0) return status_;
    }
#endif

  }  /* End of Iteration */
  
  nIteration_ = iIter;
  chisqr_ = chiSqr;

  if( !RKtraceToLast( hPntCont_ ) )
    status_ = false;

  saveCalcPosition( hPntCont_ ); 

  if( !saveTrackParameters( iniCord ) )
    status_ = false;

#if 0
  { 
    std::ios::fmtflags oldFlags = std::cout.flags();
    std::size_t oldPrec = std::cout.precision();
    std::cout.setf( std::ios::fixed );
    std::cout.precision(5);

    std::cout << funcname << ": Status = " << status_
	      << " in " << std::setw(3) << iIter << " ( "
	      << std::setw(2) << iIterEf << " ) Iteractions "
	      << " with ChiSqr=" << std::setw(10) << chisqr_ << std::endl;
    std::cout << " ** TGTlocal ** X ( " << std::setprecision(2)
	      << std::setw(7) << priPos_.x() << ", "
	      << std::setw(7) << priPos_.y() << ", "
	      << std::setw(7) << priPos_.z() << " )"
	      << " P " << std::setprecision(5) 
	      << std::setw(7) << priMom_.Mag() << " ( "
	      << std::setw(7) << priMom_.x() << ", "
	      << std::setw(7) << priMom_.y() << ", "
	      << std::setw(7) << priMom_.z() << " )"
	      << " PL: " << std::setprecision(1) 
	      << std::setw(7) << pLenTtoT_ << " "
	      << std::setw(7) << pLenTot_ << std::endl;

    PrintCalcHits( hPntCont_, std::cout );

    std::cout.flags( oldFlags );
    std::cout.precision( oldPrec );
  }
#endif

    
  return status_;
}


double SksTrack::calcChiSqr( const RKHitPointContainer &hpCont ) const
{
  std::size_t nh = hitArray_.size();

  double chisqr=0.0;
  int n=0;

  for( std::size_t i=0; i<nh; ++i ){
    TrackHit *thp = hitArray_[i];
    if(!thp) continue;

    int lnum = thp->GetLayer();
    const RKcalcHitPoint &calhp = hpCont.HitPointOfLayer( lnum );
    double hitpos = thp->GetLocalHitPos();
    double calpos = calhp.PositionInLocal();
    double resol = DCGeomMan::GetInstance().GetResolution( lnum ); 
    double weight = 1./(resol*resol);
    chisqr += weight*(hitpos-calpos)*(hitpos-calpos);
    ++n;
  }
  chisqr /= double(n-5);
  return chisqr;
}


bool SksTrack::
guessNextParameters( const RKHitPointContainer &hpCont,
		     RKCordParameter &Cord, double &estDeltaChisqr, 
		     double &lambdaCri,  double dmp ) const
{
  static const std::string funcname = "[SksTrack::geussNextParameters]";

  double *a2[10], a2c[10*5],  *v[5], vc[5*5];
  double *a3[5], a3c[5*5], *v3[5], v3c[5*5], w3[5]; 
  double dm[5];
  
  for( int i=0; i<10; ++i ){ a2[i]=&a2c[5*i]; }  
  for( int i=0; i<5; ++i ){  
    v[i]=&vc[5*i]; a3[i]=&a3c[5*i]; v3[i]=&v3c[5*i]; 
  }

  double cb2[10], wSvd[5], dcb[5];
  double wv[5];   // working space for SVD functions

  std::size_t nh = hitArray_.size();

  for( int i=0; i<10; ++i ){
    cb2[i]=0.0;
    for( int j=0; j<5; ++j ) a2[i][j]=0.0;
  }

  int nth=0;
  for( int i=0; i<nh; ++i ){
    TrackHit *thp = hitArray_[i];
    if(!thp) continue;

    int lnum = thp->GetLayer();
    const RKcalcHitPoint &calhp = hpCont.HitPointOfLayer( lnum );
    double cb = thp->GetLocalHitPos()-calhp.PositionInLocal();
    double wt = DCGeomMan::GetInstance().GetResolution( lnum );  
    wt = 1./(wt*wt);
    double cfy=calhp.coefY(), cfz=calhp.coefZ();
    double cfu=calhp.coefU(), cfv=calhp.coefV(), cfq=calhp.coefQ();
    ++nth;

    cb2[0] += 2.*cfy*wt*cb;  cb2[1] += 2.*cfz*wt*cb;  cb2[2] += 2.*cfu*wt*cb;
    cb2[3] += 2.*cfv*wt*cb;  cb2[4] += 2.*cfq*wt*cb;

    a2[0][0] += 2.*wt*(cfy*cfy - cb*calhp.coefYY());
    a2[0][1] += 2.*wt*(cfy*cfz - cb*calhp.coefYZ());
    a2[0][2] += 2.*wt*(cfy*cfu - cb*calhp.coefYU());
    a2[0][3] += 2.*wt*(cfy*cfv - cb*calhp.coefYV());
    a2[0][4] += 2.*wt*(cfy*cfq - cb*calhp.coefYQ());

    a2[1][0] += 2.*wt*(cfz*cfy - cb*calhp.coefZY());
    a2[1][1] += 2.*wt*(cfz*cfz - cb*calhp.coefZZ());
    a2[1][2] += 2.*wt*(cfz*cfu - cb*calhp.coefZU());
    a2[1][3] += 2.*wt*(cfz*cfv - cb*calhp.coefZV());
    a2[1][4] += 2.*wt*(cfz*cfq - cb*calhp.coefZQ());

    a2[2][0] += 2.*wt*(cfu*cfy - cb*calhp.coefUY());
    a2[2][1] += 2.*wt*(cfu*cfz - cb*calhp.coefUZ());
    a2[2][2] += 2.*wt*(cfu*cfu - cb*calhp.coefUU());
    a2[2][3] += 2.*wt*(cfu*cfv - cb*calhp.coefUV());
    a2[2][4] += 2.*wt*(cfu*cfq - cb*calhp.coefUQ());
    
    a2[3][0] += 2.*wt*(cfv*cfy - cb*calhp.coefVY());
    a2[3][1] += 2.*wt*(cfv*cfz - cb*calhp.coefVZ());
    a2[3][2] += 2.*wt*(cfv*cfu - cb*calhp.coefVU());
    a2[3][3] += 2.*wt*(cfv*cfv - cb*calhp.coefVV());
    a2[3][4] += 2.*wt*(cfv*cfq - cb*calhp.coefVQ());

    a2[4][0] += 2.*wt*(cfq*cfy - cb*calhp.coefQY());
    a2[4][1] += 2.*wt*(cfq*cfz - cb*calhp.coefQZ());
    a2[4][2] += 2.*wt*(cfq*cfu - cb*calhp.coefQU());
    a2[4][3] += 2.*wt*(cfq*cfv - cb*calhp.coefQV());
    a2[4][4] += 2.*wt*(cfq*cfq - cb*calhp.coefQQ());
  }

  for( int i=0; i<5; ++i )
    for( int j=0; j<5; ++j )
      a3[i][j]=a2[i][j];
      

  // Levenberg-Marqardt method
  double lambda=sqrt(dmp);
  //  a2[5][0]=a2[6][1]=a2[7][2]=a2[8][3]=a2[9][4]=lambda;

  for( int ii=0; ii<5; ++ii ){ 
    dm[ii]=a2[ii][ii]; a2[ii+5][ii]=lambda*sqrt(a2[ii][ii]);
  }
    

#if 0
  {
    std::ios::fmtflags oldFlags = std::cout.flags();
    std::size_t oldPrec = std::cout.precision();
    std::cout.setf( std::ios::scientific ); 
    std::cout.precision(3); 
    std::cout << funcname << ": A2 and CB2 before SVDcmp" 
	      <<  std::endl;
    for( int ii=0; ii<10; ++ii ) 
      std::cout << std::setw(12) << a2[ii][0] << ","
		<< std::setw(12) << a2[ii][1] << ","
		<< std::setw(12) << a2[ii][2] << ","
		<< std::setw(12) << a2[ii][3] << ","
		<< std::setw(12) << a2[ii][4] << "  "
		<< std::setw(12) << cb2[ii] << std::endl;

    std::cout.flags( oldFlags );
    std::cout.precision( oldPrec );
  }
#endif

  // Solve the Eq. with SVD (Singular Value Decomposition) Method
  if( !MathTools::SVDcmp( a2, 10, 5, wSvd, v, wv ) )
    return false;

#if 0
  {
    std::ios::fmtflags oldFlags = std::cout.flags();
    std::size_t oldPrec = std::cout.precision();
    std::cout.setf( std::ios::scientific ); 
    std::cout.precision(3); 
    std::cout << funcname << ": A2 after SVDcmp"
	      <<  std::endl;
    for( int ii=0; ii<10; ++ii ) 
      std::cout << std::setw(12) << a2[ii][0] << ","
		<< std::setw(12) << a2[ii][1] << ","
		<< std::setw(12) << a2[ii][2] << ","
		<< std::setw(12) << a2[ii][3] << ","
		<< std::setw(12) << a2[ii][4] << std::endl;
    std::cout.flags( oldFlags );
    std::cout.precision( oldPrec );
  }
#endif

#if 0
  // check orthogonality of decomposted matrics
  {
    std::ios::fmtflags oldFlags = std::cout.flags();
    std::size_t oldPrec = std::cout.precision();
    std::cout.setf( std::ios::fixed ); 
    std::cout.precision(5); 
    std::cout << funcname << ": Check V*~V" <<  std::endl;
    for( int i=0; i<5; ++i ){
      for( int j=0; j<5; ++j ){
	double f=0.0;
	for( int k=0; k<5; ++k )
	  f += v[i][k]*v[j][k];
	std::cout << std::setw(10) << f;
      }
      std::cout << std::endl;
    }
    std::cout << funcname << ": Check U*~U" <<  std::endl;
    for( int i=0; i<10; ++i ){
      for( int j=0; j<10; ++j ){
	double f=0.0;
	for( int k=0; k<5; ++k )
	  f += a2[i][k]*a2[j][k];
	std::cout << std::setw(10) << f;
      }
      std::cout << std::endl;
    }

    std::cout << funcname << ": Check ~U*U" <<  std::endl;
    for( int i=0; i<5; ++i ){
      for( int j=0; j<5; ++j ){
	double f=0.0;
	for( int k=0; k<10; ++k )
	  f += a2[k][i]*a2[k][j];
	std::cout << std::setw(10) << f;
      }
      std::cout << std::endl;
    }

    std::cout.flags( oldFlags );
    std::cout.precision( oldPrec );
  }
#endif

  double wmax=0.0;
  for( int i=0; i<5; ++i )
    if( wSvd[i]>wmax ) wmax=wSvd[i];

  double wmin=wmax*1.E-15;
  for( int i=0; i<5; ++i )
    if( wSvd[i]<wmin ) wSvd[i]=0.0;

#if 0
  {
    std::ios::fmtflags oldFlags = std::cout.flags();
    std::size_t oldPrec = std::cout.precision();
    std::cout.setf( std::ios::scientific ); 
    std::cout.precision(3); 
    std::cout << funcname << ": V and Wsvd after SVDcmp"
	      <<  std::endl;
    for( int ii=0; ii<5; ++ii ) 
      std::cout << std::setw(12) << v[ii][0] << ","
		<< std::setw(12) << v[ii][1] << ","
		<< std::setw(12) << v[ii][2] << ","
		<< std::setw(12) << v[ii][3] << ","
		<< std::setw(12) << v[ii][4] << "  "
		<< std::setw(12) << wSvd[ii] << std::endl;

    std::cout.flags( oldFlags );
    std::cout.precision( oldPrec );
  }
#endif

  MathTools::SVDksb( a2, wSvd, v, 10, 5, cb2, dcb, wv );

#if 0
  {
    std::ios::fmtflags oldFlags = std::cout.flags();
    std::size_t oldPrec = std::cout.precision();
    std::cout.setf( std::ios::fixed ); 
    std::cout.precision(5); 
    std::cout << funcname << ": " 
	      << std::setw(12) << Cord.X();
    std::cout.flags( oldFlags );
    std::cout.setf( std::ios::scientific );
    std::cout << "  Dumping Factor = " << std::setw(12) << dmp << std::endl;
    std::cout.flags( oldFlags );
    std::cout.setf( std::ios::fixed );
    std::cout << std::setw(12) << Cord.Y() << "  "
	      << std::setw(12) << dcb[0] << " ==>  "
	      << std::setw(12) << Cord.Y()+dcb[0] << std::endl;
    std::cout << std::setw(12) << Cord.Z() << "  "
	      << std::setw(12) << dcb[1] << " ==>  "
	      << std::setw(12) << Cord.Z()+dcb[1] << std::endl;
    std::cout << std::setw(12) << Cord.U() << "  "
	      << std::setw(12) << dcb[2] << " ==>  "
	      << std::setw(12) << Cord.U()+dcb[2] << std::endl;
    std::cout << std::setw(12) << Cord.V() << "  "
	      << std::setw(12) << dcb[3] << " ==>  "
	      << std::setw(12) << Cord.V()+dcb[3] << std::endl;
    std::cout << std::setw(12) << Cord.Q() << "  "
	      << std::setw(12) << dcb[4] << " ==>  "
	      << std::setw(12) << Cord.Q()+dcb[4] << std::endl;

    std::cout.flags( oldFlags );
    std::cout.precision( oldPrec );
  }
#endif

  Cord = RKCordParameter( Cord.X(),
			  Cord.Y()+dcb[0],
			  Cord.Z()+dcb[1],
			  Cord.U()+dcb[2],
			  Cord.V()+dcb[3],
			  Cord.Q()+dcb[4] );

  // calc. the critical dumping factor & est. delta-ChiSqr
  double s1=0., s2=0.;
  for(  int i=0; i<5; ++i ){
    s1 += dcb[i]*dcb[i]; s2 += dcb[i]*cb2[i];
  }
  estDeltaChisqr = (-s2-dmp*s1)/double(nth-5);

  if( !MathTools::SVDcmp( a3, 5, 5, w3, v3, wv ) )
    return false;

  double spur=0.;
  for( int i=0; i<5; ++i ){
    double s=0.;
    for( int j=0; j<5; ++j )
      s += v3[i][j]*a3[i][j];
    if( w3[i]!=0.0 )
      spur += s/w3[i]*dm[i];
  }

  lambdaCri = 1./spur;
 
  return true;
}

void SksTrack::saveCalcPosition( const RKHitPointContainer &hpCont )
{
  static const std::string funcname = "[SksTrack::saveCalcPosition]";

  std::size_t nh = hitArray_.size();

  for( int i=0; i<nh; ++i ){
    TrackHit *thp = hitArray_[i];
    if( !thp ) continue;
    int lnum = thp->GetLayer();
    const RKcalcHitPoint &calhp = hpCont.HitPointOfLayer( lnum );
    thp->SetCalGPos( calhp.PositionInGlobal() );
    thp->SetCalLPos( calhp.PositionInLocal() );
  }
}


void SksTrack::PrintCalcHits( const RKHitPointContainer &hpCont,
			      std::ostream &ost ) const
{
  std::size_t ndh = hitArray_.size();

  RKHitPointContainer::RKHpCIterator 
    itr=hpCont.begin(), end=hpCont.end();

  std::ios::fmtflags oldFlags = ost.flags();
  std::size_t oldPrec = ost.precision();

  ost.setf( std::ios::fixed );

  for( ; itr!=end; ++itr ){
    int lnum = itr->first;
    const RKcalcHitPoint &calhp = itr->second;
    ThreeVector pos = calhp.PositionInGlobal();
    ost.precision(2);
    ost << "PL#" << std::setw(2) << lnum << ":"
	<< " L: " << std::setw(9) << calhp.PathLength() 
	<< " X " << std::setw(7) << calhp.PositionInLocal() 
	<< " ( " << std::setw(8) << pos.x()
	<< ", " << std::setw(8) << pos.y()
	<< ", " << std::setw(8) << pos.z()
	<< " )" ;
    for( int i=0; i<ndh; ++i )
      if( hitArray_[i] && hitArray_[i]->GetLayer()==lnum ){
	TrackHit *thp=hitArray_[i];
	if( thp ){
	  ost << " " << std::setw(9) << thp->GetLocalHitPos()
	      << " -> " << std::setw(8) 
	      << thp->GetLocalHitPos()-calhp.PositionInLocal();
	}
      }
    ost << std::endl;
#if 0
    {
      ThreeVector mom = calhp.MomentumInGlobal();
      ost.precision(5);
      ost << "   P=" << std::setw(7) << mom.Mag()
	  << " (" << std::setw(9) << mom.x()
	  << ", " << std::setw(9) << mom.y()
	  << ", " << std::setw(9) << mom.z()
	  << " )" << std::endl;
    }
#endif
  }

  ost.flags( oldFlags );
  ost.precision( oldPrec );
}

bool SksTrack::saveTrackParameters( const RKCordParameter &cp )
{
  CPval_ = cp;
  const DCGeomMan &geomMan = DCGeomMan::GetInstance();
  int TOFid = geomMan.GetTofId();
  int TGTid = hPntCont_.begin()->first;

  const RKcalcHitPoint &hpTof  = hPntCont_.HitPointOfLayer( TOFid );
  const RKcalcHitPoint &hpTgt  = hPntCont_.begin()->second;
  const RKcalcHitPoint &hpLast = hPntCont_.rbegin()->second;

  const ThreeVector &pos = hpTgt.PositionInGlobal();
  const ThreeVector &mom = hpTgt.MomentumInGlobal();

  priPos_ = geomMan.Global2LocalPos( TGTid, pos );
  priMom_ = geomMan.Global2LocalDir( TGTid, mom );

  pLenTtoT_ = fabs( hpTgt.PathLength()-hpTof.PathLength() );
  pLenTot_  = fabs( hpTgt.PathLength()-hpLast.PathLength() );

  return true;
}

bool SksTrack::
GetTrajectoryLocalPosition( int layer, double & x, double & y ) const
{
  static const std::string funcname = 
    "[SksTrack:GetTrajectoryLocalPosition]";

  try {
    const RKcalcHitPoint &HP=hPntCont_.HitPointOfLayer( layer );
    const ThreeVector &gpos=HP.PositionInGlobal();
    ThreeVector lpos=DCGeomMan::GetInstance().Global2LocalPos(layer,gpos);
    x=lpos.x(); y=lpos.y();
    return true;
  }
  catch( std::out_of_range ) {
    return false;
  }
}
