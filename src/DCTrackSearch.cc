/*
  DCTrackSearch.cc

  2012/1/24
*/

#include "DCTrackSearch.hh"
#include "DCParameters.hh"
#include "DCLTrackHit.hh"
#include "DCPairHitCluster.hh"
#include "DCLocalTrack.hh"
#include "TemplateLib.hh"
#include "DetectorID.hh"
#include "DCGeomMan.hh"
#include "MWPCCluster.hh"
#include "TrackMaker.hh"

#include "Hodo2Hit.hh"
#include "Hodo1Hit.hh"

#include <string>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <cmath>

const double Deg2Rad = acos(-1.)/180.;
const double Rad2Deg = 180./acos(-1.);
const double MaxChisquare      = 30.;//Set to be More than 30
const double MaxChisquareSDC34 = 1000.;//Set to be More than 30

const double MaxNumberOfClusters = 20.;//Set to be Less than 30
const double MaxCombi = 1.0e6;//Set to be Less than 10^6

//SdcIn & BcOut for XUV Tracking routine
const double MaxChisquareVXU = 50.;//
const double ChisquareCutVXU = 50.;//

//SdcOut for XUV Tracking routine
const double MaxChisquareVXUSDC34 = 100.;
const double ChisquareCutVXUSDC34 = 100.;
int NhitGroup=0; 

/////////////////////////////////////
////SDC3&4 Tracking Rootine Selection

//ver.3 -> XUV Tracking
#define ver3  1 

//ver.2 -> Slow routine
#define ver2  0 

//ver.1 -> Old Routine
#define ver1  0 


///////////////////////////
////MWPC Tracking

//______________________________________________________________________________
// file local functions 
//______________________________________________________________________________
inline
void
clearCandidates(std::vector<std::vector<DCPairHitCluster*> >& candCont)
{
  for (std::vector<std::vector<DCPairHitCluster*> >::iterator 
	 i    = candCont.begin(), 
	 iEnd = candCont.end();
       i!=iEnd; ++i)
    {
      std::vector<DCPairHitCluster*>& c = *i;
      std::for_each(c.begin(), c.end(), DeleteObject());
    }
  return;
}

//______________________________________________________________________________
inline
void
clearFlags(std::vector<DCLocalTrack*>& trackCont)

{
  for(int i=0, n=trackCont.size(); i<n; ++i)
    {
      const DCLocalTrack* const tp=trackCont[i];
      if (!tp)
	continue;
      int nh=tp->GetNHit();
      for( int j=0; j<nh; ++j ) tp->GetHit(j)->clearFlags();
    }
  return;
}

//______________________________________________________________________________
inline
void
deleteDuplicatedTracks(std::vector<DCLocalTrack*>& trackCont)
		       
{
  // evaluate container size in every iteration
  for(std::size_t i=0; i<trackCont.size(); ++i){
    const DCLocalTrack* const tp=trackCont[i];
    if (!tp)
      continue;
    int nh=tp->GetNHit();
    for( int j=0; j<nh; ++j ) tp->GetHit(j)->setFlags();

    for(std::size_t i2=trackCont.size()-1; i2>i; --i2 ){
      const DCLocalTrack* tp2=trackCont[i2];
      int nh2=tp2->GetNHit(), flag=0;
      for( int j=0; j<nh2; ++j )
	if( tp2->GetHit(j)->showFlags() ) ++flag;

      if(flag>0){
	delete tp2;
	tp2 = 0;
	trackCont.erase(trackCont.begin()+i2);
      }
    }
  }

  return;
}

//______________________________________________________________________________
inline
void
calcIntersection(std::vector<DCLocalTrack*>& trackCont)
{
  for(int i=0, nn=trackCont.size(); i<nn; ++i ){
    const DCLocalTrack* const tp=trackCont[i];
    if (!tp)
      continue;
    int nh=tp->GetNHit();
    for( int j=0; j<nh; ++j ){
      int lnum = tp->GetHit(j)->GetLayer();
      double zz = DCGeomMan::GetInstance().GetLocalZ( lnum );
      tp->GetHit(j)->SetCalPosition(tp->GetX(zz), tp->GetY(zz));
    }
  }
  return;
}

//______________________________________________________________________________
void
debugPrint(const std::vector<int>& nCombi,
	   const std::string& funcname="",
	   const std::string& msg="")
{
  int n  =1;
  int nn =1;
  int sum=0;
  std::cout << "#D " << funcname << ":" ;
  for (std::vector<int>::const_iterator i=nCombi.begin(), iEnd= nCombi.end(); 
       i!=iEnd; ++i)
    {
      int val = *i;
      sum += val;
      nn *= (val+1);
      std::cout << " " << val;
      if (val!=0)
	{
	  n *= val;
	}
    }
  if (sum==0)
    n=0;
  std::cout << ": total = " << n << ", " << nn << ", " << std::endl;
  return;
}

//______________________________________________________________________________
void
debugPrint(const std::vector<DCLocalTrack*> trackCont,
	   const std::string& funcname="",
	   const std::string& msg="")
{
  int nn=trackCont.size();
  std::cout << funcname << msg
	    << nn << std::endl;
  
  for( int i=0; i<nn; ++i ){
    const DCLocalTrack * const track=trackCont[i];
    if (!track)
      continue;
    std::cout << std::setw(3) << i << " #Hits="
	      << std::setw(2) << track->GetNHit() 
	      << " ChiSqr=" << track->GetChiSquare()
	      << std::endl;
  }
  std::cout << std::endl;
  return;
}

//______________________________________________________________________________
void
debugPrint(const std::vector<int>& nCombi,
	   const std::vector<std::vector<DCPairHitCluster*> >& CandCont,
	   const std::string& funcname="",
	   const std::string& msg="")
{
  std::cout << funcname << ": #Hits of each group" << std::endl;
  int np = nCombi.size();
  for( int i=0; i<np; ++i ) std::cout << std::setw(4) << nCombi[i];
  std::cout << std::endl;
  for( int i=0; i<np; ++i ){
    int n=CandCont[i].size();
    std::cout << "[" << std::setw(3) << i << "]: "
	      << std::setw(3) << n << " ";
    for( int j=0; j<n; ++j ){
      std::cout << ((DCLTrackHit *)CandCont[i][j]->GetHit(0))->GetWire() << " ";
      //      std::cout << CandCont[i][j] << " ";
    }
    std::cout << std::endl;
  }
  return;
}

//______________________________________________________________________________
// file local functions (end)
//______________________________________________________________________________

//______________________________________________________________________________
int
findUnusedHits(int n,
	       const DCHitContainer* src,
	       std::vector<DCHitContainer>& dest)
{
  if (dest.size()<n)
    dest.resize(n);
  int ret=0;
  for (int iLayer=0; iLayer<n; ++iLayer)
    {
      const DCHitContainer& csrc = src[iLayer];
      DCHitContainer& cdest      = dest[iLayer];
      for (DCHitContainer::const_iterator 
	     it=csrc.begin(),
	     itEnd=csrc.end();
	   it!=itEnd; ++it)
	{
	  const DCHit* const h = *it;
	  if (h && !h->showFlags())
	    {
	      cdest.push_back(const_cast<DCHit*>(h));
	      ++ret;
	    }
	}
    }
  return ret;
}

///////////////////////////////////
//BC3&4, SDC1&2 Pair Plane tracking
///////////////////////////////////
//----------------------------------------------------------------------------
// Local Track Search with BH2Filter
//----------------------------------------------------------------------------
int LocalTrackSearch(const std::vector<std::vector<DCHitContainer> > &hcAssemble,
		     const DCPairPlaneInfo * PpInfo,
		     int npp, std::vector <DCLocalTrack *> &trackCont,
		     int MinNumOfHits ){
  
  static const std::string funcname = "[LocalTrackSearch (BH2Filter)]";
  for (std::vector<std::vector<DCHitContainer> >::const_iterator 
	 it    = hcAssemble.begin(),
	 itEnd = hcAssemble.end();
       it!=itEnd; ++it){
    const std::vector<DCHitContainer>& l = *it;
    std::vector<DCLocalTrack*> tc;
    LocalTrackSearch(&l[0], PpInfo, npp, tc, MinNumOfHits);
    trackCont.insert(trackCont.end(), tc.begin(), tc.end());
    //       std::cout << " tc " << tc.size()
    // 		<< " : " << trackCont.size() << std::endl;
  }
  clearFlags(trackCont);

  partial_sort( trackCont.begin(), trackCont.end(), 
		trackCont.end(), DCLTrackComp() );

  deleteDuplicatedTracks(trackCont);
//   std::cout << " trackCont.size() = " << trackCont.size() << std::endl;
  return trackCont.size();
}



//----------------------------------------------------------------------------
// Local Track Search without BH2Filter
//----------------------------------------------------------------------------
int LocalTrackSearch( const DCHitContainer * HC,
		      const DCPairPlaneInfo * PpInfo,
		      int npp, std::vector <DCLocalTrack *> &TrackCont,
		      int MinNumOfHits )
{
  static const std::string funcname = "[LocalTrackSearch]";

  std::vector < std::vector <DCPairHitCluster *> > CandCont;
  CandCont.resize(npp);

  for( int i=0; i<npp; ++i ){

    //std::cout << "npp = " << npp << ", i = " << i << std::endl;

    bool ppFlag=PpInfo[i].flag;
    int layer1=PpInfo[i].id1, layer2=PpInfo[i].id2;
    //std::cout << "layer1 = " << layer1 << ", layer2 = " << layer2 << std::endl;
    if(ppFlag) {
      //std::cout << "MakePairPlaneHitCluster" << std::endl;
      MakePairPlaneHitCluster( HC[layer1], HC[layer2], 
			       PpInfo[i].CellSize, CandCont[i] );
    }else if(layer1==layer2){
      MakeMWPCPairPlaneHitCluster( HC[layer1], CandCont[i]);
    }else{
      MakeUnPairPlaneHitCluster( HC[layer1], CandCont[i] );
    }
  }

  std::vector <int> nCombi(npp);
  for( int i=0; i<npp; ++i ){ 
    nCombi[i]=(CandCont[i]).size();
    if( nCombi[i]>MaxNumberOfClusters ) nCombi[i]=0;
  }
  //   std::vector <int> nCombi(npp);
  //   for( int i=0; i<npp; ++i ){ 
  //     nCombi[i]=(CandCont[i]).size();
  //     // If #Cluster>MaxNumberOfCluster,  error return
  //     if( nCombi[i]>MaxNumberOfClusters ){
  //       for( int i=0; i<npp; ++i )
  //   	for_each( CandCont[i].begin(), CandCont[i].end(), DeleteObject() );
  //       return 0;
  //     } 
  //   }

  //  if( nCombi[2]==0 ){
#if 0
  std::cout << funcname << ": #Hits of each group" << std::endl;
  for( int i=0; i<npp; ++i ) std::cout << std::setw(4) << nCombi[i];
  std::cout << std::endl;
  for( int i=0; i<npp; ++i ){
    int n=CandCont[i].size();
    std::cout << "[" << std::setw(3) << i << "]: "
	      << std::setw(3) << n << " ";
    for( int j=0; j<n; ++j ){
      std::cout << "layer: "<<((DCLTrackHit *)CandCont[i][j]->GetHit(0))->GetLayer() <<"  "<<"wire: " <<((DCLTrackHit *)CandCont[i][j]->GetHit(0))->GetWire()<<"     ";
      //      std::cout << CandCont[i][j] << " ";
    }
    std::cout << std::endl;
  }
#endif
  //       }

  std::vector < std::vector <int> > 
    CombiIndex = makeindex( npp, &nCombi[0] );
  int nnCombi=CombiIndex.size();

#if 0
  std::cout << " ===> " << nnCombi << " combinations will be checked.." 
	    << std::endl;
#endif
  if( nnCombi>MaxCombi ) return 0;


  for( int i=0; i<nnCombi; ++i ){
    DCLocalTrack *track = MakeTrack( CandCont, &((CombiIndex[i])[0]));
    if( !track ) continue;
    //std::cout<<"#########here-4-4-4-4-4-4-4"<<std::endl;
    if( 
       track->GetNHit()>=MinNumOfHits && 
       track->DoFit() && track->GetChiSquare()<MaxChisquare ){

#if 0      
      std::cout << "NHit = " << track->GetNHit() << " Chisqr = " << track->GetChiSquare() << std::endl;
      
      std::cout << "******** i = " << i << " *******************" << std::endl;
#endif

      TrackCont.push_back(track);
    }    
    else
      delete track;
  }

  // Clear Flags
  int nbefore=TrackCont.size();
  for( int i=0; i<nbefore; ++i ){
    DCLocalTrack *tp=TrackCont[i];
    int nh=tp->GetNHit();
    for( int j=0; j<nh; ++j ) tp->GetHit(j)->clearFlags();
  }

#if 0
  {
    int nn=TrackCont.size();
    std::cout << funcname << ": Before Sorting. #Tracks = " 
 	      << nn << std::endl;
    for( int i=0; i<nn; ++i ){
      DCLocalTrack *track=TrackCont[i];
      std::cout << std::setw(3) << i << " #Hits="
 		<< std::setw(2) << track->GetNHit() 
 		<< " ChiSqr=" << track->GetChiSquare()
 		<< std::endl;
    }
    std::cout << std::endl;
    
  }
#endif

  partial_sort( TrackCont.begin(), TrackCont.end(), 
		TrackCont.end(), DCLTrackComp() );
  
#if 0
  {
    int nn=TrackCont.size();
    std::cout << funcname << ": After Sorting. #Tracks = " 
 	      << nn << std::endl;
    
    for( int i=0; i<nn; ++i ){
      DCLocalTrack *track=TrackCont[i];
      std::cout << std::setw(3) << i << " #Hits="
 		<< std::setw(2) << track->GetNHit() 
 		<< " ChiSqr=" << track->GetChiSquare()
 		<< std::endl;
    }
    std::cout << std::endl;
    
  }
#endif
  

  // Delete Duplicated Tracks
  for( int i=0; i<int(TrackCont.size()); ++i ){
    DCLocalTrack *tp=TrackCont[i];
    int nh=tp->GetNHit();
    for( int j=0; j<nh; ++j ) tp->GetHit(j)->setFlags();
    
     for( int i2=TrackCont.size()-1; i2>i; --i2 ){
       DCLocalTrack *tp2=TrackCont[i2];
       int nh2=tp2->GetNHit(), flag=0;
       for( int j=0; j<nh2; ++j )
	 if( tp2->GetHit(j)->showFlags() ) ++flag;
       if(flag){
	 delete tp2;
	 TrackCont.erase(TrackCont.begin()+i2);
       }
     }      
   }

  {
    int nn=TrackCont.size();
    for(int i=0; i<nn; ++i ){
      DCLocalTrack *tp=TrackCont[i];
       int nh=tp->GetNHit();
       for( int j=0; j<nh; ++j ){
	 int lnum = tp->GetHit(j)->GetLayer();
	 double zz = DCGeomMan::GetInstance().GetLocalZ( lnum );
	 tp->GetHit(j)->SetCalPosition(tp->GetX(zz), tp->GetY(zz));
       }
    }
  }

#if 0
  //  if(TrackCont.size()>1){
    int nn=TrackCont.size();
     std::cout << funcname << ": After Deleting. #Tracks = " 
	       << nn << std::endl;
     
     for( int i=0; i<nn; ++i ){
       DCLocalTrack *track=TrackCont[i];
       std::cout << std::setw(3) << i << " #Hits="
		 << std::setw(2) << track->GetNHit() 
		 << " ChiSqr=" << track->GetChiSquare()
		 << std::endl;
       std::cout << std::endl;

       /*
       for( int j=0; j<(track->GetNHit()); ++j){
	 DCLTrackHit *hit = track->GetHit(j);
	 std::cout << "layer = " << hit->GetLayer()+1 << " Res = " << hit->GetResidual() << std::endl ;
	 std::cout << "hitp = " << hit->GetLocalHitPos() << " calp = " << hit->GetLocalCalPos() << std::endl ;
	 std::cout << "X = " << hit->GetXcal() << " Y = " << hit->GetYcal() << std::endl ;
	 std::cout << std::endl;
       }
       */
       std::cout << "*********************************************" << std::endl;
     }
     std::cout << std::endl;
     
     //  }
#endif
     
     for( int i=0; i<npp; ++i ){
       for_each( CandCont[i].begin(), CandCont[i].end(), DeleteObject() );
     }
     

  return TrackCont.size();
}


//----------------------------------------------------------------------------
// Local Track Search SdcIn w/Fiber   (added by yamamoto)
//----------------------------------------------------------------------------
int LocalTrackSearchSdcInFiber( const DCHitContainer * HC,
				const DCPairPlaneInfo * PpInfo,
				int npp, std::vector <DCLocalTrack *> &TrackCont,
				int MinNumOfHits )
{
  static const std::string funcname = "[LocalTrackSearch]";

  std::vector < std::vector <DCPairHitCluster *> > CandCont;
  CandCont.resize(npp);

  for( int i=0; i<npp; ++i ){

    bool ppFlag=PpInfo[i].flag;
    int layer1=PpInfo[i].id1, layer2=PpInfo[i].id2;
    //std::cout << "layer1 = " << layer1 << ", layer2 = " << layer2 << std::endl;
    if(ppFlag) {
      MakePairPlaneHitCluster( HC[layer1], HC[layer2], 
			       PpInfo[i].CellSize, CandCont[i] );
    }else if(layer1==layer2){
      MakeMWPCPairPlaneHitCluster( HC[layer1], CandCont[i]);
    }else{
      MakeUnPairPlaneHitCluster( HC[layer1], CandCont[i] );
    }
  }

  std::vector <int> nCombi(npp);
  for( int i=0; i<npp; ++i ){ 
    nCombi[i]=(CandCont[i]).size();
    if( nCombi[i]>MaxNumberOfClusters ) nCombi[i]=0;
  }

  std::vector < std::vector <int> > 
    CombiIndex = makeindex( npp, &nCombi[0] );
  int nnCombi=CombiIndex.size();

  if( nnCombi>MaxCombi ) return 0;


  for( int i=0; i<nnCombi; ++i ){
    DCLocalTrack *track = MakeTrack( CandCont, &((CombiIndex[i])[0]));
    if( !track ) continue;
    if( track->GetNHit()>=MinNumOfHits && track->DoFit() &&
	track->GetChiSquare()<MaxChisquare ){

      TrackCont.push_back(track);
    }    
    else
      delete track;
  }

  // Clear Flags
  int nbefore=TrackCont.size();
  for( int i=0; i<nbefore; ++i ){
    DCLocalTrack *tp=TrackCont[i];
    int nh=tp->GetNHit();
    for( int j=0; j<nh; ++j ) tp->GetHit(j)->clearFlags();
  }

  partial_sort( TrackCont.begin(), TrackCont.end(), 
		TrackCont.end(), DCLTrackCompSdcInFiber() );

  // Delete Duplicated Tracks
  for( int i=0; i<int(TrackCont.size()); ++i ){
    DCLocalTrack *tp=TrackCont[i];
    int nh=tp->GetNHit();
    for( int j=0; j<nh; ++j ) tp->GetHit(j)->setFlags();
    
     for( int i2=TrackCont.size()-1; i2>i; --i2 ){
       DCLocalTrack *tp2=TrackCont[i2];
       int nh2=tp2->GetNHit(), flag=0;
       for( int j=0; j<nh2; ++j )
	 if( tp2->GetHit(j)->showFlags() ) ++flag;
       if(flag){
	 delete tp2;
	 TrackCont.erase(TrackCont.begin()+i2);
       }
     }      
   }

  {
    int nn=TrackCont.size();
    for(int i=0; i<nn; ++i ){
      DCLocalTrack *tp=TrackCont[i];
       int nh=tp->GetNHit();
       for( int j=0; j<nh; ++j ){
	 int lnum = tp->GetHit(j)->GetLayer();
	 double zz = DCGeomMan::GetInstance().GetLocalZ( lnum );
	 tp->GetHit(j)->SetCalPosition(tp->GetX(zz), tp->GetY(zz));
       }
    }
  }
     
     for( int i=0; i<npp; ++i ){
       for_each( CandCont[i].begin(), CandCont[i].end(), DeleteObject() );
     }

  return TrackCont.size();
}



////////////////////////////
//BC3&4, SDC1&2 XUV Tracking
////////////////////////////
int LocalTrackSearchVUX( const DCHitContainer * HC,
			 const DCPairPlaneInfo * PpInfo,
			 int npp, std::vector <DCLocalTrack *> &TrackCont,
			 int MinNumOfHits )
{
  static const std::string funcname = "[LocalTrackSearchVUX]";
  
  std::vector <DCLocalTrack *>  TrackContV;
  std::vector <DCLocalTrack *>  TrackContX;
  std::vector <DCLocalTrack *>  TrackContU;

  std::vector < std::vector <DCPairHitCluster *> > CandCont;
  std::vector < std::vector <DCPairHitCluster *> > CandContV;
  std::vector < std::vector <DCPairHitCluster *> > CandContX;
  std::vector < std::vector <DCPairHitCluster *> > CandContU;

  int NumOfLayersDC_12 = 12 ; 

  CandCont.resize(npp);
  CandContV.resize(npp);
  CandContX.resize(npp);
  CandContU.resize(npp);

  int iV=0,iX=0,iU=0;
  int nV=0, nX=0, nU=0;

  for( int i=0; i<npp; ++i ){
    bool ppFlag=PpInfo[i].flag;
    int layer1=PpInfo[i].id1, layer2=PpInfo[i].id2;
    int nh1 = HC[layer1].size();
    int nh2 = HC[layer2].size();
    double TiltAngle ;
    if(nh1>0) TiltAngle = HC[layer1][0]->GetTiltAngle();
    if((ppFlag) && (nh1==0) && (nh2>0)) TiltAngle = HC[layer2][0]->GetTiltAngle();
    if((ppFlag) && (nh1>0) && (nh2>0)){
      if(TiltAngle<0){
	MakePairPlaneHitClusterVUX( HC[layer1], HC[layer2], 
				    PpInfo[i].CellSize, CandContV[iV] );
	++iV;
	nV = nV+2;
      }
      if(TiltAngle==0){
	MakePairPlaneHitClusterVUX( HC[layer1], HC[layer2], 
				    PpInfo[i].CellSize, CandContX[iX] );
	++iX;
	nX = nX+2;
      }
      if(TiltAngle>0){
	MakePairPlaneHitClusterVUX( HC[layer1], HC[layer2], 
				    PpInfo[i].CellSize, CandContU[iU] );
	++iU;
	nU = nU+2;
      }
    }
    if(!(ppFlag)){
      if(TiltAngle<0){
	MakeUnPairPlaneHitCluster( HC[layer1], CandContV[iV] );
	++nV;
	++iV;
      }
      if(TiltAngle==0){
	MakeUnPairPlaneHitCluster( HC[layer1], CandContX[iX] );
	++nX;
	++iX;
      }
      if(TiltAngle>0){
	MakeUnPairPlaneHitCluster( HC[layer1], CandContU[iU] );
	++nU;
	++iU;
      }
    }
    if((ppFlag) && (nh1==0) && (nh2>0)){

      if(TiltAngle<0){
	MakeUnPairPlaneHitCluster( HC[layer2], CandContV[iV] );
	++nV;
	++iV;
      }
      if(TiltAngle==0){
	MakeUnPairPlaneHitCluster( HC[layer2], CandContX[iX] );
	++nX;
	++iX;
      }
      if(TiltAngle>0){
	MakeUnPairPlaneHitCluster( HC[layer2], CandContU[iU] );
	++nU;
	++iU;
      }
    }
  }

  std::vector <int> nCombi(npp);
  std::vector <int> nCombiV(npp);
  std::vector <int> nCombiX(npp);
  std::vector <int> nCombiU(npp);

  for( int i=0; i<npp; ++i ){ 
    nCombiV[i]=(CandContV[i]).size();
    nCombiX[i]=(CandContX[i]).size();
    nCombiU[i]=(CandContU[i]).size();

    if( nCombiV[i]>MaxNumberOfClusters ) nCombiV[i]=0;
    if( nCombiX[i]>MaxNumberOfClusters ) nCombiX[i]=0;
    if( nCombiU[i]>MaxNumberOfClusters ) nCombiU[i]=0;
  }

#if 0
  std::cout << funcname << ": #Hits of each group V plane" << std::endl;
  for( int i=0; i<npp; ++i ) std::cout << std::setw(4) << nCombiV[i];
  std::cout << std::endl;
  for( int i=0; i<npp; ++i ){
    int n=CandContV[i].size();
    std::cout << "[" << std::setw(3) << i << "]: "
	      << std::setw(3) << n << " ";
    for( int j=0; j<n; ++j ){
      std::cout << ((DCLTrackHit *)CandContV[i][j]->GetHit(0))->GetWire() << " ";
    }
    std::cout << std::endl;
  }

  std::cout << funcname << ": #Hits of each group X plane" << std::endl;
  for( int i=0; i<npp; ++i ) std::cout << std::setw(4) << nCombiX[i];
  std::cout << std::endl;
  for( int i=0; i<npp; ++i ){
    int n=CandContX[i].size();
    std::cout << "[" << std::setw(3) << i << "]: "
	      << std::setw(3) << n << " ";
    for( int j=0; j<n; ++j ){
      std::cout << ((DCLTrackHit *)CandContX[i][j]->GetHit(0))->GetWire() << " ";
    }
    std::cout << std::endl;
  }

  std::cout << funcname << ": #Hits of each group U plane" << std::endl;
  for( int i=0; i<npp; ++i ) std::cout << std::setw(4) << nCombiU[i];
  std::cout << std::endl;
  for( int i=0; i<npp; ++i ){
    int n=CandContU[i].size();
    std::cout << "[" << std::setw(3) << i << "]: "
	      << std::setw(3) << n << " ";
    for( int j=0; j<n; ++j ){
      std::cout << ((DCLTrackHit *)CandContU[i][j]->GetHit(0))->GetWire() << " ";
    }
    std::cout << std::endl;
  }
#endif

  std::vector < std::vector <int> > 
    CombiIndexV = makeindex( npp, &nCombiV[0] );
  int nnCombiV=CombiIndexV.size();
  std::vector < std::vector <int> > 
    CombiIndexX = makeindex( npp, &nCombiX[0] );
  int nnCombiX=CombiIndexX.size();
  std::vector < std::vector <int> > 
    CombiIndexU = makeindex( npp, &nCombiU[0] );
  int nnCombiU=CombiIndexU.size();
  
#if 0
  std::cout << "V Plane  ===> " << nnCombiV << " combinations will be checked.." 
	    << std::endl;
  std::cout << "X Plane  ===> " << nnCombiX << " combinations will be checked.." 
	    << std::endl;
  std::cout << "U Plane  ===> " << nnCombiU << " combinations will be checked.." 
	    << std::endl;
#endif

  int NHitV=0,NHitX=0,NHitU=0;
  for( int i=0; i<nnCombiV; ++i ){
    DCLocalTrack *track = MakeTrack( CandContV, &((CombiIndexV[i])[0]) );
    if( !track ) continue;
    if( track->GetNHit()>=3 && track->DoFitVXU() &&
	track->GetChiSquare()<MaxChisquareVXU ){
      TrackContV.push_back(track);
      NHitV = track->GetNHit();
      double chisqr = track->GetChiSquare();
    }
    else{
      delete track;
    }
  }

  for( int i=0; i<nnCombiX; ++i ){
    DCLocalTrack *track = MakeTrack( CandContX, &((CombiIndexX[i])[0]) );
    if( !track ) continue;
    if( track->GetNHit()>=3 && track->DoFitVXU() &&
	track->GetChiSquare()<MaxChisquareVXU ){
      TrackContX.push_back(track);
      NHitX = track->GetNHit();
      double chisqr = track->GetChiSquare();
    }
    else{
      delete track;
    }
  }

  for( int i=0; i<nnCombiU; ++i ){
    DCLocalTrack *track = MakeTrack( CandContU, &((CombiIndexU[i])[0]) );
    if( !track ) continue;
    if( track->GetNHit()>=3 && track->DoFitVXU() &&
	track->GetChiSquare()<MaxChisquareVXU ){
      TrackContU.push_back(track);
      NHitU = track->GetNHit();
      double chisqr = track->GetChiSquare();
    }
    else{
      delete track;
    }
  }

  int nbeforeV=0, nbeforeX=0, nbeforeU=0;

    // Clear Flags
  if(nV>3){
    nbeforeV=TrackContV.size();
    for( int i=0; i<nbeforeV; ++i ){
      DCLocalTrack *tp=TrackContV[i];
      int nh=tp->GetNHit();
      for( int j=0; j<nh; ++j ) tp->GetHit(j)->clearFlags();
    }
  }

  if(nX>3){
    nbeforeX=TrackContX.size();
    for( int i=0; i<nbeforeX; ++i ){
      DCLocalTrack *tp=TrackContX[i];
      int nh=tp->GetNHit();
      for( int j=0; j<nh; ++j ) tp->GetHit(j)->clearFlags();
    }
  }

  if(nU>3){
    nbeforeU=TrackContU.size();
    for( int i=0; i<nbeforeU; ++i ){
      DCLocalTrack *tp=TrackContU[i];
      int nh=tp->GetNHit();
      for( int j=0; j<nh; ++j ) tp->GetHit(j)->clearFlags();
    }
  }

#if 0
  {
    /* V Plane  */
    {
      int nn=TrackContV.size();
      std::cout << funcname << ": [V Plane] Before Sorting. #Tracks = " 
		<< nn << " nV = " << nV << std::endl;
      
      for( int i=0; i<nn; ++i ){
	DCLocalTrack *track=TrackContV[i];
	std::cout << std::setw(3) << i << " #Hits="
		  << std::setw(2) << track->GetNHit() 
		  << " ChiSqr=" << track->GetChiSquare()
		  << std::endl;
      }
      std::cout << std::endl;
    }
    /* X Plane  */
    {
      int nn=TrackContX.size();
      std::cout << funcname << ": [X Plane] Before Sorting. #Tracks = " 
		<< nn << " nX = " << nX << std::endl;
      
      for( int i=0; i<nn; ++i ){
	DCLocalTrack *track=TrackContX[i];
	std::cout << std::setw(3) << i << " #Hits="
		  << std::setw(2) << track->GetNHit() 
		  << " ChiSqr=" << track->GetChiSquare()
		  << std::endl;
      }
      std::cout << std::endl;
    }
    /* U Plane  */
    {
      int nn=TrackContU.size();
      std::cout << funcname << ": [U Plane] Before Sorting. #Tracks = " 
		<< nn << " nU = " << nU << std::endl;
      
      for( int i=0; i<nn; ++i ){
	DCLocalTrack *track=TrackContU[i];
	std::cout << std::setw(3) << i << " #Hits="
		  << std::setw(2) << track->GetNHit() 
		  << " ChiSqr=" << track->GetChiSquare()
		  << std::endl;
      }
      std::cout << std::endl;
    }
  }
#endif

  //  if(nbeforeV>0)  
  partial_sort( TrackContV.begin(), TrackContV.end(), 
		TrackContV.end(), DCLTrackComp1() );
  //  if(nbeforeX>0)  
  partial_sort( TrackContX.begin(), TrackContX.end(), 
		TrackContX.end(), DCLTrackComp1() );
  //  if(nbeforeU>0)  
  partial_sort( TrackContU.begin(), TrackContU.end(), 
		TrackContU.end(), DCLTrackComp1() );

#if 0
  {
    /* V Plane  */
    {
      int nn=TrackContV.size();
      std::cout << funcname << ": [V Plane] After Sorting. #Tracks = " 
		<< nn << " nV = " << nV << std::endl;
      
      for( int i=0; i<nn; ++i ){
	DCLocalTrack *track=TrackContV[i];
	std::cout << std::setw(3) << i << " #Hits="
		  << std::setw(2) << track->GetNHit() 
		  << " ChiSqr=" << track->GetChiSquare()
		  << " dX/dZ(Av)=" << track->GetVXU_A()*cos(Deg2Rad*(-15.0))
		  << " dY/dZ(Av)=" << (-1)*track->GetVXU_A()*sin(Deg2Rad*(-15.0))
		  << std::endl;
      }
      std::cout << std::endl;
    }
    /* X Plane  */
    {
      int nn=TrackContX.size();
      std::cout << funcname << ": [X Plane] After Sorting. #Tracks = " 
		<< nn << " nX = " << nX << std::endl;
      
      for( int i=0; i<nn; ++i ){
	DCLocalTrack *track=TrackContX[i];
	std::cout << std::setw(3) << i << " #Hits="
		  << std::setw(2) << track->GetNHit() 
		  << " ChiSqr=" << track->GetChiSquare()
		  << " dX/dZ(Ax)=" << track->GetVXU_A()
		  << std::endl;
      }
      std::cout << std::endl;
    }
    /* U Plane  */
    {
      int nn=TrackContU.size();
      std::cout << funcname << ": [U Plane] After Sorting. #Tracks = " 
		<< nn << " nU = " << nU << std::endl;
      
      for( int i=0; i<nn; ++i ){
	DCLocalTrack *track=TrackContU[i];
	std::cout << std::setw(3) << i << " #Hits="
		  << std::setw(2) << track->GetNHit() 
		  << " ChiSqr=" << track->GetChiSquare()
		  << " dX/dZ(Au)=" << track->GetVXU_A()*cos(Deg2Rad*15.0)
		  << " dY/dZ(Au)=" << track->GetVXU_A()*sin(Deg2Rad*15.0)
		  << std::endl;
      }
      std::cout << std::endl;
    }
  }
#endif

  // Delete Duplicated Tracks (cut chisqr>100 & flag)
  double chiV=ChisquareCutVXU, chiX=ChisquareCutVXU, chiU=ChisquareCutVXU ;

#if 1
  {
    /* V Plane  */  
    for( int i=0; i<int(TrackContV.size()); ++i ){
      DCLocalTrack *tp=TrackContV[i];
      int nh=tp->GetNHit();
      if(chiV > tp->GetChiSquare()) chiV = tp->GetChiSquare();
      for( int j=0; j<nh; ++j ) tp->GetHit(j)->setFlags();
      
      for( int i2=TrackContV.size()-1; i2>i; --i2 ){
	DCLocalTrack *tp2=TrackContV[i2];
	int nh2=tp2->GetNHit(), flag=0;
	for( int j=0; j<nh2; ++j )
	  if( tp2->GetHit(j)->showFlags() ) ++flag;
	//if(tp2->GetChiSquare()>100){
	if((flag>2) && (tp2->GetChiSquare()>ChisquareCutVXU)){
	  //if((flag>2) || (tp2->GetChiSquare()>ChisquareCutVXU)){
	  delete tp2;
	  TrackContV.erase(TrackContV.begin()+i2);
	}
      }      
    }
    /* X Plane  */    
    for( int i=0; i<int(TrackContX.size()); ++i ){
      DCLocalTrack *tp=TrackContX[i];
      int nh=tp->GetNHit();
      if(chiX > tp->GetChiSquare()) chiX = tp->GetChiSquare();
      for( int j=0; j<nh; ++j ) tp->GetHit(j)->setFlags();
      
      for( int i2=TrackContX.size()-1; i2>i; --i2 ){
	DCLocalTrack *tp2=TrackContX[i2];
	int nh2=tp2->GetNHit(), flag=0;
	for( int j=0; j<nh2; ++j )
	  if( tp2->GetHit(j)->showFlags() ) ++flag;
	//if(tp2->GetChiSquare()>100){
	if((flag>2) && (tp2->GetChiSquare()>ChisquareCutVXU)){
	  //if((flag>2) || (tp2->GetChiSquare()>ChisquareCutVXU)){
	  delete tp2;
	  TrackContX.erase(TrackContX.begin()+i2);
	}
      }      
    }
    /* U Plane  */    
    for( int i=0; i<int(TrackContU.size()); ++i ){
      DCLocalTrack *tp=TrackContU[i];
      int nh=tp->GetNHit();
      if(chiU > tp->GetChiSquare()) chiU = tp->GetChiSquare();
      for( int j=0; j<nh; ++j ) tp->GetHit(j)->setFlags();
      
      for( int i2=TrackContU.size()-1; i2>i; --i2 ){
	DCLocalTrack *tp2=TrackContU[i2];
	int nh2=tp2->GetNHit(), flag=0;
	for( int j=0; j<nh2; ++j )
	  if( tp2->GetHit(j)->showFlags() ) ++flag;
	//if(tp2->GetChiSquare()>100){
	if((flag>2) && (tp2->GetChiSquare()>ChisquareCutVXU)){
	//if((flag>2) || (tp2->GetChiSquare()>ChisquareCutVXU)){
	  delete tp2;
	  TrackContU.erase(TrackContU.begin()+i2);
	}
      }      
    }
  }
#endif

  //  int cV=TrackContV.size(), cX=TrackContX.size(), cU=TrackContU.size() ;

  /*  
      int cV=0, cX=0, cU=0 ;
      if(chiV<ChisquareCutVXU) ++cV ;
      if(chiX<ChisquareCutVXU) ++cX ;
      if(chiU<ChisquareCutVXU) ++cU ;
  */

#if 0
  {
    /* V Plane  */
    {
      int nn=TrackContV.size();
      std::cout << funcname << ": [V Plane] After Delete. #Tracks = " 
		<< nn << " nV = " << nV << std::endl;
      
      for( int i=0; i<nn; ++i ){
	DCLocalTrack *track=TrackContV[i];
	std::cout << std::setw(3) << i << " #Hits="
		  << std::setw(2) << track->GetNHit() 
		  << " ChiSqr=" << track->GetChiSquare()
		  << std::endl;
      }
      std::cout << std::endl;
    }
    /* X Plane  */
    {
      int nn=TrackContX.size();
      std::cout << funcname << ": [X Plane] After Delete. #Tracks = " 
		<< nn << " nX = " << nX << std::endl;
      
      for( int i=0; i<nn; ++i ){
	DCLocalTrack *track=TrackContX[i];
	std::cout << std::setw(3) << i << " #Hits="
		  << std::setw(2) << track->GetNHit() 
		  << " ChiSqr=" << track->GetChiSquare()
		  << std::endl;
      }
      std::cout << std::endl;
    }
    /* U Plane  */
    {
      int nn=TrackContU.size();
      std::cout << funcname << ": [U Plane] After Delete. #Tracks = " 
		<< nn << " nU = " << nU << std::endl;
      
      for( int i=0; i<nn; ++i ){
	DCLocalTrack *track=TrackContU[i];
	std::cout << std::setw(3) << i << " #Hits="
		  << std::setw(2) << track->GetNHit() 
		  << " ChiSqr=" << track->GetChiSquare()
		  << std::endl;
      }
      std::cout << std::endl;
    }
  }
#endif

  int nnV=1,nnX=1,nnU=1;
  int nkV=0,nkX=0,nkU=0;
  int nnVT=1,nnXT=1,nnUT=1;
  int checkV=0,checkX=0,checkU=0;

  int cV=TrackContV.size();
  if(chiV>1.5 || cV<5 ) checkV++;
  int cX=TrackContX.size();
  if(chiX>1.5 || cX<5 ) checkX++;
  int cU=TrackContU.size();
  if(chiU>1.5 || cU<5 ) checkU++;

  std::vector < std::vector <int> > CombiIndexSV;
  std::vector < std::vector <int> > CombiIndexSX;
  std::vector < std::vector <int> > CombiIndexSU;

  {
    if((nV>=3) && (cV)){
      nnV = TrackContV.size();
      nnVT = TrackContV.size();
      ++nkV;
    }
    if((nV>0) && checkV){
      CombiIndexSV = makeindex_VXU( npp, 2, &(nCombiV[0]) );
      nnV = nnV +  CombiIndexSV.size();
    }
  }

  {
    if((nX>=3) && (cX)){
      nnX = TrackContX.size();
      nnXT = TrackContX.size();
      ++nkX;
    }
    if((nX>0) && checkX){
      CombiIndexSX = makeindex_VXU( npp, 2, &(nCombiX[0]) );
      nnX = nnX + CombiIndexSX.size();
    }
  }

  {
    if((nU>=3) && (cU)){
      nnU = TrackContU.size();
      nnUT = TrackContU.size();
      ++nkU;
    }
    if((nU>0) && checkU){
      CombiIndexSU = makeindex_VXU( npp, 2, &(nCombiU[0]) );
      nnU = nnU + CombiIndexSU.size();
    }
  }

  double DifVXU=0.0;
  double Av=0.0;
  double Ax=0.0;
  double Au=0.0;

  double chiv, chix, chiu;

#if 0
  for( int i=0; i<nnV; ++i){
    for( int j=0; j<nnX; ++j){
      for( int k=0; k<nnU; ++k){

	chiv=-1.0,chix=-1.0,chiu=-1.0;

	DCLocalTrack *track = new DCLocalTrack();
	
	int mV=0,mX=0,mU=0;
	int mmV=0,mmX=0,mmU=0;
	
	/* V Plane  */
	if(nkV){
	  DCLocalTrack *trackV=TrackContV[i];
	  //Av=trackV->GetVXU_A();
	  chiv=trackV->GetChiSquare();
	  for( int l=0; l<(trackV->GetNHit()); ++l){
	    DCLTrackHit *hitpV=trackV->GetHit(l);
	    if(hitpV){
	      track->AddHit( hitpV ) ;
	    }
	  }
	  //delete trackV;
	}
	if((!nkV) && (nV>0)){
	  DCLocalTrack *trackV = MakeTrack( CandContV, &((CombiIndexV[i])[0]) ); 
	  for( int l=0; l<(trackV->GetNHit()); ++l){
	    DCLTrackHit *hitpV=trackV->GetHit(l);
	    if(hitpV){
	      track->AddHit( hitpV ) ;
	    }
	  }
	  delete trackV;
	}
	
	/* X Plane  */
	if(nkX){
	  DCLocalTrack *trackX=TrackContX[j];
	  //Ax=trackX->GetVXU_A();
	  chix=trackX->GetChiSquare();
	  for( int l=0; l<(trackX->GetNHit()); ++l){
	    DCLTrackHit *hitpX=trackX->GetHit(l);
	    if(hitpX){
	      track->AddHit( hitpX ) ;
	    }
	  }
	  //delete trackX;
	}
	if((!nkX) && (nX>0)){
	  DCLocalTrack *trackX = MakeTrack( CandContX, &((CombiIndexX[j])[0]) ); 
	  for( int l=0; l<(trackX->GetNHit()); ++l){
	    DCLTrackHit *hitpX=trackX->GetHit(l);
	    if(hitpX){
	      track->AddHit( hitpX ) ;
	    }
	  }
	  delete trackX;
	}
	
	/* U Plane  */
	if(nkU){
	  DCLocalTrack *trackU=TrackContU[k];
	  //Au=trackU->GetVXU_A();
	  chiu=trackU->GetChiSquare();
	  for( int l=0; l<(trackU->GetNHit()); ++l){
	    DCLTrackHit *hitpU=trackU->GetHit(l);
	    if(hitpU){
	      track->AddHit( hitpU ) ;
	    }
	  }
	  //delete trackU;
	}
	if((!nkU) && (nU>0)){
	  DCLocalTrack *trackU = MakeTrack( CandContU, &((CombiIndexU[k])[0]) ); 
	  for( int l=0; l<(trackU->GetNHit()); ++l){
	    DCLTrackHit *hitpU=trackU->GetHit(l);
	    if(hitpU){
	      track->AddHit( hitpU ) ;
	    }
	  }
	  delete trackU;
	}
	//track->SetAv(Av);
	//track->SetAx(Ax);
	//track->SetAu(Au);
	//DifVXU = track->GetDifVXU();
	
	track->SetChiv(chiv);
	track->SetChix(chix);
	track->SetChiu(chiu);
	
	if(!track) continue;
	if( track->GetNHit()>=MinNumOfHits && track->DoFit() &&
	    track->GetChiSquare()<MaxChisquare ){
	  TrackCont.push_back(track);
	}    
	else{
	  delete track;
	}
      }
    }
  }
#endif

  for( int i=-1; i<nnV; ++i){
    for( int j=-1; j<nnX; ++j){
      for( int k=-1; k<nnU; ++k){
	if( ((i+j)==-2) || ((j+k)==-2) || ((k+i)==-2)) continue;

	chiv=-1.0,chix=-1.0,chiu=-1.0;

	DCLocalTrack *track = new DCLocalTrack();
	
	int mV=0,mX=0,mU=0;
	int mmV=0,mmX=0,mmU=0;
	
	/* V Plane  */
	if(i>-1){
	  if(nkV && i<nnVT){
	    DCLocalTrack *trackV=TrackContV[i];
	    Av=trackV->GetVXU_A();
	    chiv=trackV->GetChiSquare();
	    
	    for( int l=0; l<(trackV->GetNHit()); ++l){
	      DCLTrackHit *hitpV=trackV->GetHit(l);
	      if(hitpV){
		track->AddHit( hitpV ) ;
	      }
	    }
	  }
	  if((i>=nnVT) && (nV>0)){
	    DCLocalTrack *trackV = MakeTrack( CandContV, &((CombiIndexSV[i-nnVT])[0]) ); 
	    for( int l=0; l<(trackV->GetNHit()); ++l){
	      DCLTrackHit *hitpV=trackV->GetHit(l);
	      if(hitpV){
		track->AddHit( hitpV ) ;
	      }
	    }
	    delete trackV ;
	  }
	  int NHitV = track->GetNHit();
	}

	/* X Plane  */
	if(j>-1){
	  if(nkX && j<nnXT){
	    DCLocalTrack *trackX=TrackContX[j];
	    Ax=trackX->GetVXU_A();
	    chix=trackX->GetChiSquare();
	    for( int l=0; l<(trackX->GetNHit()); ++l){
	      DCLTrackHit *hitpX=trackX->GetHit(l);
	      if(hitpX){
		track->AddHit( hitpX ) ;
	      }
	    }
	  }
	  if((j>=nnXT) && (nX>0)){
	    DCLocalTrack *trackX = MakeTrack( CandContX, &((CombiIndexSX[j-nnXT])[0]) ); 
	    for( int l=0; l<(trackX->GetNHit()); ++l){
	      DCLTrackHit *hitpX=trackX->GetHit(l);
	      if(hitpX){
		track->AddHit( hitpX ) ;
	      }
	    }
	    delete trackX;
	  }


	  int NHitX = track->GetNHit();
	}
	/* U Plane  */
	if(k>-1){
	  if(nkU && k<nnUT){
	    DCLocalTrack *trackU=TrackContU[k];
	    Au=trackU->GetVXU_A();
	    chiu=trackU->GetChiSquare();
	    for( int l=0; l<(trackU->GetNHit()); ++l){
	      DCLTrackHit *hitpU=trackU->GetHit(l);
	      if(hitpU){
		track->AddHit( hitpU ) ;
	      }
	    }
	  }
	  if((k>=nnUT) && (nU>0)){
	    DCLocalTrack *trackU = MakeTrack( CandContU, &((CombiIndexSU[k-nnUT])[0]) ); 
	    for( int l=0; l<(trackU->GetNHit()); ++l){
	      DCLTrackHit *hitpU=trackU->GetHit(l);
	      if(hitpU){
		track->AddHit( hitpU ) ;
	      }
	    }
	    delete trackU;
	  }
	  int NHitU = track->GetNHit();
	}
	
	//track->SetAv(Av);
	//track->SetAx(Ax);
	//track->SetAu(Au);
	//DifVXU = track->GetDifVXU();
	
	track->SetChiv(chiv);
	track->SetChix(chix);
	track->SetChiu(chiu);
	
	if(!track) continue;
	if(track->GetNHit()>=MinNumOfHits && track->DoFit() && 
	   track->GetChiSquare()<MaxChisquare ){//MaXChisquare
	    TrackCont.push_back(track);
	}    
	else{
	  delete track;
	}
      }
    }
  }
  
  {

    for( int i=0; i<int(TrackContV.size()); ++i ){
      DCLocalTrack *tp=TrackContV[i];
      delete tp;
      TrackContV.erase(TrackContV.begin()+i);
    }
    for( int i=0; i<int(TrackContX.size()); ++i ){
      DCLocalTrack *tp=TrackContX[i];
      delete tp;
      TrackContX.erase(TrackContX.begin()+i);
    }
    for( int i=0; i<int(TrackContU.size()); ++i ){
      DCLocalTrack *tp=TrackContU[i];
      delete tp;
      TrackContU.erase(TrackContU.begin()+i);
    }

  }
  // Clear Flags
  {
    int nbefore=TrackCont.size();
    for( int i=0; i<nbefore; ++i ){
      DCLocalTrack *tp=TrackCont[i];
      int nh=tp->GetNHit();
      for( int j=0; j<nh; ++j ){
	tp->GetHit(j)->clearFlags();
      }
    }
  }
  
  partial_sort( TrackCont.begin(), TrackCont.end(), 
		TrackCont.end(), DCLTrackComp1() );    
  
#if 0
  {
    int nn=TrackCont.size();
    std::cout << funcname << ": Before Sorting. #Tracks = " 
	      << nn << std::endl;
    
    for( int i=0; i<nn; ++i ){
      //    for( int i=0; i<20; ++i ){
      DCLocalTrack *track=TrackCont[i];
      std::cout << std::setw(3) << i << " #Hits="
		<< std::setw(2) << track->GetNHit() 
		<< " ChiSqr=" << track->GetChiSquare()
	//		<< " DifVXU=" << track->GetDifVXU()
		<< std::endl;
    }
    std::cout << std::endl;
    
  }
#endif

  // Clear Flags
  {
    int nbefore=TrackCont.size();
    for( int i=0; i<nbefore; ++i ){
      DCLocalTrack *tp=TrackCont[i];
      int nh=tp->GetNHit();
      for( int j=0; j<nh; ++j ){
	tp->GetHit(j)->clearFlags();
      }
    }
  }
  
#if 1  
  // Delete Tracks about  (Nhit1 > Nhit2+1) (Nhit1 > Nhit2  && chi1 < chi2)
  for( int i=0; i<int(TrackCont.size()); ++i ){
    DCLocalTrack *tp=TrackCont[i];
    int nh=tp->GetNHit();
      double chi=tp->GetChiSquare();
    for( int j=0; j<nh; ++j ) tp->GetHit(j)->setFlags();
    
    for( int i2=TrackCont.size()-1; i2>i; --i2 ){
      DCLocalTrack *tp2=TrackCont[i2];
      int nh2=tp2->GetNHit(), flag=0;
      double chi2=tp2->GetChiSquare();
      for( int j=0; j<nh2; ++j )
	if( tp2->GetHit(j)->showFlags() ) ++flag;
      if((flag>=2) && ((nh==nh2) || ((nh>nh2) && (chi<chi2)))){
	//      if((flag) && ((nh>nh2+1) || ((nh==nh2) || (nh>nh2) && (chi<chi2)))){
	//if((nh>nh2) && (chi<chi2)){
        delete tp2;
	TrackCont.erase(TrackCont.begin()+i2);
      }
    }      
  }    
#endif  

#if 0
  {
    int nn=TrackCont.size();
    std::cout << funcname << ": Before Deleted Sorting. #Tracks = " 
	      << nn << std::endl;
    
    for( int i=0; i<nn; ++i ){
      DCLocalTrack *track=TrackCont[i];
      std::cout << std::setw(3) << i << " #Hits="
		<< std::setw(2) << track->GetNHit() 
		<< " ChiSqr=" << track->GetChiSquare()
		<< std::endl;
    }
    std::cout << std::endl;
    
  }
#endif
  
  partial_sort( TrackCont.begin(), TrackCont.end(), 
		TrackCont.end(), DCLTrackComp() );
  
#if 0
  {
    int nn=TrackCont.size();
    std::cout << funcname << ": After Sorting. #Tracks = " 
	      << nn << std::endl;
    
    for( int i=0; i<nn; ++i ){
      //    for( int i=0; i<20; ++i ){
      DCLocalTrack *track=TrackCont[i];
      std::cout << std::setw(3) << i << " #Hits="
		<< std::setw(2) << track->GetNHit() 
		<< " ChiSqr=" << track->GetChiSquare()
		<< std::endl;
    }
    std::cout << std::endl;
    
  }
#endif

  // Clear Flags
  {
    int nbefore=TrackCont.size();
    for( int i=0; i<nbefore; ++i ){
      DCLocalTrack *tp=TrackCont[i];
      int nh=tp->GetNHit();
      for( int j=0; j<nh; ++j ){
	tp->GetHit(j)->clearFlags();
      }
    }
  }

  // Delete Duplicated Tracks
  
  for( int i=0; i<int(TrackCont.size()); ++i ){
    DCLocalTrack *tp=TrackCont[i];
    int nh=tp->GetNHit();
    for( int j=0; j<nh; ++j ) tp->GetHit(j)->setFlags();
    
    for( int i2=TrackCont.size()-1; i2>i; --i2 ){
      DCLocalTrack *tp2=TrackCont[i2];
      int nh2=tp2->GetNHit(), flag=0;
      for( int j=0; j<nh2; ++j )
	if( tp2->GetHit(j)->showFlags() ) ++flag;
      if(flag){
	delete tp2;
	TrackCont.erase(TrackCont.begin()+i2);
      }
    }      
  }

  // Clear Flags
  {
    int nbefore=TrackCont.size();
    for( int i=0; i<nbefore; ++i ){
      DCLocalTrack *tp=TrackCont[i];
      int nh=tp->GetNHit();
      for( int j=0; j<nh; ++j ){
	tp->GetHit(j)->clearFlags();
      }
    }
  }  

  {
    int nn=TrackCont.size();
    for(int i=0; i<nn; ++i ){
      DCLocalTrack *tp=TrackCont[i];
      int nh=tp->GetNHit();
      for( int j=0; j<nh; ++j ){
	int lnum = tp->GetHit(j)->GetLayer();
	double zz = DCGeomMan::GetInstance().GetLocalZ( lnum );
	tp->GetHit(j)->SetCalPosition(tp->GetX(zz), tp->GetY(zz));
      }
    }
  }
  
#if 0
  if(TrackCont.size()>0){
    int nn=TrackCont.size();
     std::cout << funcname << ": After Deleting. #Tracks = " 
	       << nn << std::endl;
     
     for( int i=0; i<nn; ++i ){
       DCLocalTrack *track=TrackCont[i];
       std::cout << std::setw(3) << i << " #Hits="
		 << std::setw(2) << track->GetNHit() 
		 << " ChiSqr=" << track->GetChiSquare()
		 << std::endl;
       std::cout << std::endl;
       for( int j=0; j<(track->GetNHit()); ++j){
	 DCLTrackHit *hit = track->GetHit(j);
	 std::cout << "layer = " << hit->GetLayer()+1 << " Res = " << hit->GetResidual() << std::endl ;
	 std::cout << "hitp = " << hit->GetLocalHitPos() << " calp = " << hit->GetLocalCalPos() << std::endl ;
	 std::cout << "X = " << hit->GetXcal() << " Y = " << hit->GetYcal() << std::endl ;
	 std::cout << std::endl;
       }
       std::cout << "*********************************************" << std::endl;
     }
     std::cout << std::endl;
     
  }
#endif

  for( int i=0; i<npp; ++i ){
    for_each( CandContV[i].begin(), CandContV[i].end(), DeleteObject() );
    for_each( CandContX[i].begin(), CandContX[i].end(), DeleteObject() );
    for_each( CandContU[i].begin(), CandContU[i].end(), DeleteObject() );
  }
  for( int i=0; i<npp; ++i )
    for_each( CandCont[i].begin(), CandCont[i].end(), DeleteObject() );

  return TrackCont.size();    
}


///////////////////////////
//SDC3&4 XUV Tracking
///////////////////////////

#if ver3
int SdcOutLocalTrackSearch( const DCHitContainer * HC,
 			    std::vector <DCLocalTrack *> &TrackCont)
{
  static const std::string funcname = "[LocalTrackSearch]";
  
  std::vector <DCLocalTrack *>  TrackContV;
  std::vector <DCLocalTrack *>  TrackContX;
  std::vector <DCLocalTrack *>  TrackContU;

  std::vector < std::vector <DCPairHitCluster *> > CandCont;
  std::vector < std::vector <DCPairHitCluster *> > CandContV;
  std::vector < std::vector <DCPairHitCluster *> > CandContX;
  std::vector < std::vector <DCPairHitCluster *> > CandContU;

  CandCont.resize(NumOfLayersSdcOut);
  CandContV.resize(NumOfLayersSdcOut/3);
  CandContX.resize(NumOfLayersSdcOut/3);
  CandContU.resize(NumOfLayersSdcOut/3);

  for( int i=0; i<NumOfLayersSdcOut/3; ++i ){//NumOfLayersSdcOut=12
    MakeUnPairPlaneHitCluster( HC[3*i+1], CandContV[i] );
    MakeUnPairPlaneHitCluster( HC[3*i+2], CandContX[i] );
    MakeUnPairPlaneHitCluster( HC[3*i+3], CandContU[i] );
  }   

  std::vector <int> nCombi(NumOfLayersSdcOut);
  std::vector <int> nCombiV(NumOfLayersSdcOut/3);
  std::vector <int> nCombiX(NumOfLayersSdcOut/3);
  std::vector <int> nCombiU(NumOfLayersSdcOut/3);

  int nV=0, nX=0, nU=0;

  for( int i=0; i<NumOfLayersSdcOut/3; ++i ){ 
    nCombiV[i]=(CandContV[i]).size();
    nCombiX[i]=(CandContX[i]).size();
    nCombiU[i]=(CandContU[i]).size();
    if(nCombiV[i]>0) ++nV;
    if(nCombiX[i]>0) ++nX;
    if(nCombiU[i]>0) ++nU;

    if(nCombiV[i]>MaxNumberOfClusters || 
       nCombiX[i]>MaxNumberOfClusters || 
       nCombiU[i]>MaxNumberOfClusters){
      for( int ii=0; ii<NumOfLayersSdcOut; ++ii ){
	for_each( CandContV[ii].begin(), CandContV[ii].end(), DeleteObject() );
	for_each( CandContX[ii].begin(), CandContX[ii].end(), DeleteObject() );
	for_each( CandContU[ii].begin(), CandContU[ii].end(), DeleteObject() );
	return 0;
      }
    } 
  }

  //const int NhitGroup=0;
  for(int i=0; i<NumOfLayersSdcOut/3; ++i){
    if(nCombiV[i]>3 || nCombiX[i]>3 || nCombiU[i]>3 ){
      NhitGroup++;
    }
  }
#if 0
  //////////////////----------------  V Plane -------------------////////////////////////////
  std::cout << funcname << ": V plane #Hits of each group" << std::endl;
  for( int i=0; i<NumOfLayersSdcOut/3; ++i ) std::cout << std::setw(4) << nCombiV[i];
  std::cout << std::endl;
  for( int i=0; i<NumOfLayersSdcOut/3; ++i ){
    int n=CandContV[i].size();
    std::cout << "[" << std::setw(3) << i << "]: "
	      << std::setw(3) << n << " ";
    for( int j=0; j<n; ++j ){
      std::cout << ((DCLTrackHit *)CandContV[i][j]->GetHit(0))->GetWire() << " ";
    }
    std::cout << std::endl;
  }
  
  //////////////////----------------  X Plane -------------------////////////////////////////
  std::cout << funcname << ": X plane #Hits of each group" << std::endl;
  for( int i=0; i<NumOfLayersSdcOut/3; ++i ) std::cout << std::setw(4) << nCombiX[i];
  std::cout << std::endl;
  for( int i=0; i<NumOfLayersSdcOut/3; ++i ){
    int n=CandContX[i].size();
    std::cout << "[" << std::setw(3) << i << "]: "
	      << std::setw(3) << n << " ";
    for( int j=0; j<n; ++j ){
      std::cout << ((DCLTrackHit *)CandContX[i][j]->GetHit(0))->GetWire() << " ";
    }
    std::cout << std::endl;
  }
  
  //////////////////----------------  U Plane -------------------////////////////////////////
  std::cout << funcname << ": U plane #Hits of each group" << std::endl;
  for( int i=0; i<NumOfLayersSdcOut/3; ++i ) std::cout << std::setw(4) << nCombiU[i];
  std::cout << std::endl;
  for( int i=0; i<NumOfLayersSdcOut/3; ++i ){
    int n=CandContU[i].size();
    std::cout << "[" << std::setw(3) << i << "]: "
	      << std::setw(3) << n << " ";
    for( int j=0; j<n; ++j ){
      std::cout << ((DCLTrackHit *)CandContU[i][j]->GetHit(0))->GetWire() << " ";
    }
    std::cout << std::endl;
  }
#endif

  std::vector < std::vector <int> > CombiIndexV;
  if(nV>=3){
    CombiIndexV = makeindex_SdcOut( NumOfLayersSdcOut/3, 
				    3, 
				    NumOfLayersSdcOut/3, 
				    &(nCombiV[0]) );
  }
  std::vector < std::vector <int> > CombiIndexX;
  if(nX>=3){
    CombiIndexX = makeindex_SdcOut( NumOfLayersSdcOut/3, 
				    3, 
				    NumOfLayersSdcOut/3, 
				    &(nCombiX[0]) );
  }
  std::vector < std::vector <int> > CombiIndexU;
  if(nU>=3){
    CombiIndexU = makeindex_SdcOut( NumOfLayersSdcOut/3, 
				    3, 
				    NumOfLayersSdcOut/3, 
				    &(nCombiU[0]) );
  }

  int nnCombiV=0, nnCombiX=0, nnCombiU=0;
  if(nV>=3) nnCombiV=CombiIndexV.size();
  if(nX>=3) nnCombiX=CombiIndexX.size();
  if(nU>=3) nnCombiU=CombiIndexU.size();
  
#if 0
  std::cout << "V Plane  ===> " << nnCombiV << " combinations will be checked.." 
	    << std::endl;
  std::cout << "X Plane  ===> " << nnCombiX << " combinations will be checked.." 
	    << std::endl;
  std::cout << "U Plane  ===> " << nnCombiU << " combinations will be checked.." 
	    << std::endl;
#endif

  int NHitV=0,NHitX=0,NHitU=0,NResV=0,NResX=0,NResU=0;
  for( int i=0; i<nnCombiV; ++i ){
    DCLocalTrack *track = MakeTrack( CandContV, &((CombiIndexV[i])[0]) );
    
    if( !track ) continue;
    if( track->GetNHit()>=3 && track->DoFitVXU() &&
	track->GetChiSquare()<MaxChisquareVXUSDC34 ){
      TrackContV.push_back(track);
      NHitV = track->GetNHit();
      double chisqr = track->GetChiSquare();
    }
    else{
      delete track;
    }
  }

  for( int i=0; i<nnCombiX; ++i ){
    DCLocalTrack *track = MakeTrack( CandContX, &((CombiIndexX[i])[0]) );
    
    if( !track ){
      continue;
    }
    if( track->GetNHit()>=3 && track->DoFitVXU() &&
	track->GetChiSquare()<MaxChisquareVXUSDC34 ){
      TrackContX.push_back(track);
      NHitX = track->GetNHit();
      double chisqr = track->GetChiSquare();
    }
    else{
      delete track;
    }
  }

  for( int i=0; i<nnCombiU; ++i ){
    DCLocalTrack *track = MakeTrack( CandContU, &((CombiIndexU[i])[0]) );
    
    if( !track ) continue;
    if( track->GetNHit()>=3 && track->DoFitVXU() &&
	track->GetChiSquare()<MaxChisquareVXUSDC34 ){
      TrackContU.push_back(track);
      NHitU = track->GetNHit();
      double chisqr = track->GetChiSquare();
    }
    else{
      delete track;
    }
  }

  int nbeforeV=0, nbeforeX=0, nbeforeU=0;
  if(nV>=3){
    nbeforeV=TrackContV.size();
    for( int i=0; i<nbeforeV; ++i ){
      DCLocalTrack *tp=TrackContV[i];
      int nh=tp->GetNHit();
      for( int j=0; j<nh; ++j ) tp->GetHit(j)->clearFlags();
    }
  }

  if(nX>=3){
    nbeforeX=TrackContX.size();
    for( int i=0; i<nbeforeX; ++i ){
      DCLocalTrack *tp=TrackContX[i];
      int nh=tp->GetNHit();
      for( int j=0; j<nh; ++j ) tp->GetHit(j)->clearFlags();
    }
  }

  if(nU>=3){
    nbeforeU=TrackContU.size();
    for( int i=0; i<nbeforeU; ++i ){
      DCLocalTrack *tp=TrackContU[i];
      int nh=tp->GetNHit();
      for( int j=0; j<nh; ++j ) tp->GetHit(j)->clearFlags();
    }
  }

#if 0
  {
    /* V Plane  */
    {
      int nn=TrackContV.size();
      std::cout << funcname << ": [V Plane] Before Sorting. #Tracks = " 
		<< nn << " nV = " << nV << std::endl;
      
      for( int i=0; i<nn; ++i ){
	DCLocalTrack *track=TrackContV[i];
	std::cout << std::setw(3) << i << " #Hits="
		  << std::setw(2) << track->GetNHit() 
		  << " ChiSqr=" << track->GetChiSquare()
		  << std::endl;
      }
      std::cout << std::endl;
    }
    /* X Plane  */
    {
      int nn=TrackContX.size();
      std::cout << funcname << ": [X Plane] Before Sorting. #Tracks = " 
		<< nn << " nX = " << nX << std::endl;
      
      for( int i=0; i<nn; ++i ){
	DCLocalTrack *track=TrackContX[i];
	std::cout << std::setw(3) << i << " #Hits="
		  << std::setw(2) << track->GetNHit() 
		  << " ChiSqr=" << track->GetChiSquare()
		  << std::endl;
      }
      std::cout << std::endl;
    }
    /* U Plane  */
    {
      int nn=TrackContU.size();
      std::cout << funcname << ": [U Plane] Before Sorting. #Tracks = " 
		<< nn << " nU = " << nU << std::endl;
      
      for( int i=0; i<nn; ++i ){
	DCLocalTrack *track=TrackContU[i];
	std::cout << std::setw(3) << i << " #Hits="
		  << std::setw(2) << track->GetNHit() 
		  << " ChiSqr=" << track->GetChiSquare()
		  << std::endl;
      }
      std::cout << std::endl;
    }
  }
#endif

  //  if(nbeforeV>0)  
  partial_sort( TrackContV.begin(), TrackContV.end(), 
		TrackContV.end(), DCLTrackComp1() );
  //  if(nbeforeX>0)  
  partial_sort( TrackContX.begin(), TrackContX.end(), 
		TrackContX.end(), DCLTrackComp1() );
  //  if(nbeforeU>0)  
  partial_sort( TrackContU.begin(), TrackContU.end(), 
		TrackContU.end(), DCLTrackComp1() );

#if 0
  {
    /* V Plane  */
    {
      int nn=TrackContV.size();
      std::cout << funcname << ": [V Plane] After Sorting. #Tracks = " 
		<< nn << " nV = " << nV << std::endl;
      
      for( int i=0; i<nn; ++i ){
	DCLocalTrack *track=TrackContV[i];
	std::cout << std::setw(3) << i << " #Hits="
		  << std::setw(2) << track->GetNHit() 
		  << " ChiSqr=" << track->GetChiSquare()
		  << " dX/dZ(Av)=" << track->GetVXU_A()*cos(Deg2Rad*(-30.0))
		  << " dY/dZ(Av)=" << (-1)*track->GetVXU_A()*sin(Deg2Rad*(-30.0))
		  << std::endl;
      }
      std::cout << std::endl;
    }
    /* X Plane  */
    {
      int nn=TrackContX.size();
      std::cout << funcname << ": [X Plane] After Sorting. #Tracks = " 
		<< nn << " nX = " << nX << std::endl;
      
      for( int i=0; i<nn; ++i ){
	DCLocalTrack *track=TrackContX[i];
	std::cout << std::setw(3) << i << " #Hits="
		  << std::setw(2) << track->GetNHit() 
		  << " ChiSqr=" << track->GetChiSquare()
		  << " dX/dZ(Ax)=" << track->GetVXU_A()
		  << std::endl;
      }
      std::cout << std::endl;
    }
    /* U Plane  */
    {
      int nn=TrackContU.size();
      std::cout << funcname << ": [U Plane] After Sorting. #Tracks = " 
		<< nn << " nU = " << nU << std::endl;
      
      for( int i=0; i<nn; ++i ){
	DCLocalTrack *track=TrackContU[i];
	std::cout << std::setw(3) << i << " #Hits="
		  << std::setw(2) << track->GetNHit() 
		  << " ChiSqr=" << track->GetChiSquare()
		  << " dX/dZ(Au)=" << track->GetVXU_A()*cos(Deg2Rad*30.0)
		  << " dY/dZ(Au)=" << track->GetVXU_A()*sin(Deg2Rad*30.0)
		  << std::endl;
      }
      std::cout << std::endl;
    }
  }
#endif

  // Delete Duplicated Tracks (cut chisqr>100 & flag)
  double chiV=ChisquareCutVXUSDC34, chiX=ChisquareCutVXUSDC34, chiU=ChisquareCutVXUSDC34 ;

#if 1
  {
    /* V Plane  */  
    for( int i=0; i<int(TrackContV.size()); ++i ){
      DCLocalTrack *tp=TrackContV[i];
      int nh=tp->GetNHit();
      if(chiV > tp->GetChiSquare()) chiV = tp->GetChiSquare();
      //      if(nh==nV && chiV > tp->GetChiSquare()) chiV = tp->GetChiSquare();
      else{
	delete tp;
	TrackContV.erase(TrackContV.begin()+i);
	continue;
      }
      for( int j=0; j<nh; ++j ) tp->GetHit(j)->setFlags();
      
      for( int i2=TrackContV.size()-1; i2>i; --i2 ){
	DCLocalTrack *tp2=TrackContV[i2];
	int nh2=tp2->GetNHit(), flag=0;
	for( int j=0; j<nh2; ++j )
	  if( tp2->GetHit(j)->showFlags() ) ++flag;
	if((flag) && tp2->GetChiSquare()>ChisquareCutVXUSDC34){
	  delete tp2;
	  TrackContV.erase(TrackContV.begin()+i2);
	}
      }      
    }
    /* X Plane  */    
    for( int i=0; i<int(TrackContX.size()); ++i ){
      DCLocalTrack *tp=TrackContX[i];
      int nh=tp->GetNHit();
      if(chiX > tp->GetChiSquare()) chiX = tp->GetChiSquare();
      //      if(nh==nX && chiX > tp->GetChiSquare()) chiX = tp->GetChiSquare();
      else{
	delete tp;
	TrackContX.erase(TrackContX.begin()+i);
	continue;
      }
      for( int j=0; j<nh; ++j ) tp->GetHit(j)->setFlags();
      
      for( int i2=TrackContX.size()-1; i2>i; --i2 ){
	DCLocalTrack *tp2=TrackContX[i2];
	int nh2=tp2->GetNHit(), flag=0;
	for( int j=0; j<nh2; ++j )
	  if( tp2->GetHit(j)->showFlags() ) ++flag;
	if((flag) && tp2->GetChiSquare()>ChisquareCutVXUSDC34){
	  delete tp2;
	  TrackContX.erase(TrackContX.begin()+i2);
	}
      }      
    }
    /* U Plane  */    
    for( int i=0; i<int(TrackContU.size()); ++i ){
      DCLocalTrack *tp=TrackContU[i];
      int nh=tp->GetNHit();
      if(chiU > tp->GetChiSquare()) chiU = tp->GetChiSquare();
      //      if(nh==nU && chiU > tp->GetChiSquare()) chiU = tp->GetChiSquare();
      else{
	delete tp;
	TrackContU.erase(TrackContU.begin()+i);
	continue;
      }
      for( int j=0; j<nh; ++j ) tp->GetHit(j)->setFlags();
      
      for( int i2=TrackContU.size()-1; i2>i; --i2 ){
	DCLocalTrack *tp2=TrackContU[i2];
	int nh2=tp2->GetNHit(), flag=0;
	for( int j=0; j<nh2; ++j )
	  if( tp2->GetHit(j)->showFlags() ) ++flag;
	if((flag) && tp2->GetChiSquare()>ChisquareCutVXUSDC34){
	  delete tp2;
	  TrackContU.erase(TrackContU.begin()+i2);
	}
      }      
    }
  }
#endif


#if 0
  {
    /* V Plane  */
    {
      int nn=TrackContV.size();
      std::cout << funcname << ": [V Plane] After Delete. #Tracks = " 
		<< nn << " nV = " << nV << std::endl;
      
      for( int i=0; i<nn; ++i ){
	DCLocalTrack *track=TrackContV[i];
	std::cout << std::setw(3) << i << " #Hits="
		  << std::setw(2) << track->GetNHit() 
		  << " ChiSqr=" << track->GetChiSquare()
		  << std::endl;
      }
      std::cout << std::endl;
    }
    /* X Plane  */
    {
      int nn=TrackContX.size();
      std::cout << funcname << ": [X Plane] After Delete. #Tracks = " 
		<< nn << " nX = " << nX << std::endl;
      
      for( int i=0; i<nn; ++i ){
	DCLocalTrack *track=TrackContX[i];
	std::cout << std::setw(3) << i << " #Hits="
		  << std::setw(2) << track->GetNHit() 
		  << " ChiSqr=" << track->GetChiSquare()
		  << std::endl;

      }
      std::cout << std::endl;
    }
    /* U Plane  */
    {
      int nn=TrackContU.size();
      std::cout << funcname << ": [U Plane] After Delete. #Tracks = " 
		<< nn << " nU = " << nU << std::endl;
      
      for( int i=0; i<nn; ++i ){
	DCLocalTrack *track=TrackContU[i];
	std::cout << std::setw(3) << i << " #Hits="
		  << std::setw(2) << track->GetNHit() 
		  << " ChiSqr=" << track->GetChiSquare()
		  << std::endl;
      }
      std::cout << std::endl;
    }
  }
#endif

  int nnV=1,nnX=1,nnU=1;
  int nnVT=1,nnXT=1,nnUT=1;
  int checkV=0,checkX=0,checkU=0;

  int nkV=0,nkX=0,nkU=0;

  int cV=TrackContV.size();
  if(chiV>1.5 || cV<5 ) checkV++;
  int cX=TrackContX.size();
  if(chiX>1.5 || cX<5 ) checkX++;
  int cU=TrackContU.size();
  if(chiU>1.5 || cU<5 ) checkU++;

  std::vector < std::vector <int> > CombiIndexSV;
  std::vector < std::vector <int> > CombiIndexSX;
  std::vector < std::vector <int> > CombiIndexSU;

  {
    if((nV>=3) && (cV)){
      nnV = TrackContV.size();
      nnVT = TrackContV.size();
      ++nkV;
    }
    if((nV>0) && checkV){
      CombiIndexSV = makeindex_SdcOut_below( NumOfLayersSdcOut/3, 
					     2, 
					     NumOfLayersSdcOut/3, 
					     &(nCombiV[0]) );
      nnV = nnV +  CombiIndexSV.size();
    }
  }
  
  {
    if((nX>=3) && (cX)){
      nnX = TrackContX.size();
      nnXT = TrackContX.size();
      ++nkX;
    }
    if((nX>0) && checkX){
      CombiIndexSX = makeindex_SdcOut_below( NumOfLayersSdcOut/3, 
					     2, 
					     NumOfLayersSdcOut/3, 
					     &(nCombiX[0]) );
      nnX = nnX + CombiIndexSX.size();
    }
  }
  
  {
    if((nU>=3) && (cU)){
      nnU = TrackContU.size();
      nnUT = TrackContU.size();
      ++nkU;
    }
    if((nU>0) && checkU){
      CombiIndexSU = makeindex_SdcOut_below( NumOfLayersSdcOut/3, 
					     2, 
					     NumOfLayersSdcOut/3, 
					     &(nCombiU[0]) );
      nnU = nnU + CombiIndexSU.size();
    }
  }

  /*** Event Cut by using anular difference***/
  double DifVXU=0.0;
  double Av=0.0;
  double Ax=0.0;
  double Au=0.0;

  bool angularVcut[nnVT+1];
  bool angularXcut[nnXT+1];
  bool angularUcut[nnUT+1];

  for( int i=0; i<nnVT; ++i){
    for( int j=0; j<nnXT; ++j){
      for( int k=0; k<nnUT; ++k){

	angularVcut[i]=false;
	angularXcut[j]=false;
	angularUcut[k]=false;

	if(nkV){
	  DCLocalTrack *trackV=TrackContV[i];
	  Av=trackV->GetVXU_A();
	}
	if(nkX){
	  DCLocalTrack *trackX=TrackContX[j];
	  Ax=trackX->GetVXU_A();
	}
	if(nkU){
	  DCLocalTrack *trackU=TrackContU[k];
	  Au=trackU->GetVXU_A();
	}

	//Yoneeeeeeeeee!!!!
	DifVXU = (Av*cos(acos(-1.)/180.*(-30.0))-Ax)*(Av*cos(acos(-1.)/180.*(-30.0))-Ax)+(Ax-Au*cos(acos(-1.)/180.*(30.0)))*(Ax-Au*cos(acos(-1.)/180.*(30.0)))+(Au*cos(acos(-1.)/180.*(30.0))-Av*cos(acos(-1.)/180.*(-30.0)))*(Au*cos(acos(-1.)/180.*(30.0))-Av*cos(acos(-1.)/180.*(-30.0)))+ (Au*sin(acos(-1.)/180.*(30.0))-Av*sin(acos(-1.)/180.*(30.0)))*(Au*sin(acos(-1.)/180.*(30.0))-Av*sin(acos(-1.)/180.*(30.0)));

	if( (DifVXU>=0) ){
	  angularVcut[i]=true;
	  angularXcut[j]=true;
	  angularUcut[k]=true;
	}
      }
    }
  }

  DifVXU=0.0;
  Av=0.0;
  Ax=0.0;
  Au=0.0;

  /////////////////////////////
  //X,U,V Tracks -> Add Track// 
  /////////////////////////////

 double chiv, chix, chiu;
  for( int i=-1; i<nnV; ++i){
    for( int j=-1; j<nnX; ++j){
      for( int k=-1; k<nnU; ++k){
	
	chiv=-1.0;
	chix=-1.0;
	chiu=-1.0;

	DCLocalTrack *track = new DCLocalTrack();
	
	int mV=0,mX=0,mU=0;
	int mmV=0,mmX=0,mmU=0;
	
	/* V Plane  */
	if(i>-1){
	  if(nkV && i<nnVT){
	    DCLocalTrack *trackV=TrackContV[i];
	    Av=trackV->GetVXU_A();
	    chiv=trackV->GetChiSquare();
	    
	    for( int l=0; l<(trackV->GetNHit()); ++l){
	      DCLTrackHit *hitpV=trackV->GetHit(l);
	      if(hitpV){
		track->AddHit( hitpV ) ;
	      }
	    }
	  }
	  if((i>=nnVT) && (nV>0)){
	    DCLocalTrack *trackV = MakeTrack( CandContV, &((CombiIndexSV[i-nnVT])[0]) ); 
	    for( int l=0; l<(trackV->GetNHit()); ++l){
	      DCLTrackHit *hitpV=trackV->GetHit(l);
	      if(hitpV){
		track->AddHit( hitpV ) ;
	      }
	    }
	    delete trackV ;
	  }
	  int NHitV = track->GetNHit();
	  //std::cout << "NHit(V) = " << track->GetNHit() << std::endl ;
	}

	/* X Plane  */
	if(j>-1){
	  if(nkX && j<nnXT){
	    DCLocalTrack *trackX=TrackContX[j];
	    Ax=trackX->GetVXU_A();
	    chix=trackX->GetChiSquare();
	    for( int l=0; l<(trackX->GetNHit()); ++l){
	      DCLTrackHit *hitpX=trackX->GetHit(l);
	      if(hitpX){
		track->AddHit( hitpX ) ;
	      }
	    }
	  }
	  if((j>=nnXT) && (nX>0)){
	    DCLocalTrack *trackX = MakeTrack( CandContX, &((CombiIndexSX[j-nnXT])[0]) ); 
	    for( int l=0; l<(trackX->GetNHit()); ++l){
	      DCLTrackHit *hitpX=trackX->GetHit(l);
	      if(hitpX){
		track->AddHit( hitpX ) ;
	      }
	    }
	    delete trackX;
	  }


	  int NHitX = track->GetNHit();
	  //std::cout << "NHit(V+X) = " << track->GetNHit() << std::endl ;
	}
	/* U Plane  */
	if(k>-1){
	  if(nkU && k<nnUT){
	    DCLocalTrack *trackU=TrackContU[k];
	    Au=trackU->GetVXU_A();
	    chiu=trackU->GetChiSquare();
	    for( int l=0; l<(trackU->GetNHit()); ++l){
	      DCLTrackHit *hitpU=trackU->GetHit(l);
	      if(hitpU){
		track->AddHit( hitpU ) ;
	      }
	    }
	  }
	  if((k>=nnUT) && (nU>0)){
	    DCLocalTrack *trackU = MakeTrack( CandContU, &((CombiIndexSU[k-nnUT])[0]) ); 
	    
	    for( int l=0; l<(trackU->GetNHit()); ++l){
	      DCLTrackHit *hitpU=trackU->GetHit(l);
	      if(hitpU){
		track->AddHit( hitpU ) ;
	      }
	    }
	    delete trackU;
	  }
	  int NHitU = track->GetNHit();
	  //std::cout << "NHit(V+X+U) = " << track->GetNHit() << std::endl ;
	}
	
	track->SetAv(Av);
	track->SetAx(Ax);
	track->SetAu(Au);
	DifVXU = track->GetDifVXUSDC34();
	
	track->SetChiv(chiv);
	track->SetChix(chix);
	track->SetChiu(chiu);
	
	if(!track) continue;
	if(track->GetNHit()>=MinNumOfHitsSdcOut && track->DoFit() && 
	   track->GetChiSquare()<MaxChisquareSDC34 ){	  
	    TrackCont.push_back(track);
	}    
	else{
	  delete track;
	}
      }
    }
  }
  
  for( int i=0; i<NumOfLayersSdcOut/3; ++i ){
    for_each( CandContV[i].begin(), CandContV[i].end(), DeleteObject() );
    for_each( CandContX[i].begin(), CandContX[i].end(), DeleteObject() );
    for_each( CandContU[i].begin(), CandContU[i].end(), DeleteObject() );
  }

  {
    for_each( TrackContV.begin(), TrackContV.end(), DeleteObject() );
    for_each( TrackContX.begin(), TrackContX.end(), DeleteObject() );
    for_each( TrackContU.begin(), TrackContU.end(), DeleteObject() );
  }

  // Clear Flags
  {
    int nbefore=TrackCont.size();
    for( int i=0; i<nbefore; ++i ){
      DCLocalTrack *tp=TrackCont[i];
      int nh=tp->GetNHit();
      for( int j=0; j<nh; ++j ){
	tp->GetHit(j)->clearFlags();
      }
    }
  }
  
  partial_sort( TrackCont.begin(), TrackCont.end(), 
		TrackCont.end(), DCLTrackComp1() );    
  
#if 0
  {
    int nn=TrackCont.size();
    if(nn>0){
      std::cout << funcname << ": After Sorting. #Tracks = " 
		<< nn << std::endl;
      
      for( int i=0; i<nn; ++i ){
      //for( int i=0; i<100; ++i ){
	DCLocalTrack *track=TrackCont[i];
	std::cout << std::setw(3) << i << " #Hits="
		  << std::setw(2) << track->GetNHit() 
		  << " ChiSqr=" << track->GetChiSquare()
		  << " DifVXU=" << track->GetDifVXUSDC34()
		  << std::endl;
      }
      std::cout << std::endl;
    }
  }
#endif

#if 1  
  // Delete Tracks about (Nhit1 = Nhit2  && chi1 < chi2)
  for( int i=0; i<int(TrackCont.size()); ++i ){
    DCLocalTrack *tp=TrackCont[i];
    int nh=tp->GetNHit();
      double chi=tp->GetChiSquare();
    for( int j=0; j<nh; ++j ) tp->GetHit(j)->setFlags();
    
    for( int i2=TrackCont.size()-1; i2>i; --i2 ){
      DCLocalTrack *tp2=TrackCont[i2];
      int nh2=tp2->GetNHit(), flag=0;
      double chi2=tp2->GetChiSquare();
      for( int j=0; j<nh2; ++j )
	if( tp2->GetHit(j)->showFlags() ) ++flag;
      if((flag) && ((nh==nh2) && (chi<=chi2))){
        delete tp2;
	TrackCont.erase(TrackCont.begin()+i2);
      }
    }      
  }    
#endif  

#if 0
  {
    int nn=TrackCont.size();
    if(nn>0){
      std::cout << funcname << ": Before Deleted Sorting. #Tracks = " 
		<< nn << std::endl;
      
      //for( int i=0; i<100; ++i ){
      for( int i=0; i<nn; ++i ){
	DCLocalTrack *track=TrackCont[i];
	std::cout << std::setw(3) << i << " #Hits="
		  << std::setw(2) << track->GetNHit() 
		  << " ChiSqr=" << track->GetChiSquare()
		  << std::endl;
      }
      std::cout << std::endl;
    }
  }
#endif

  partial_sort( TrackCont.begin(), TrackCont.end(), 
  		TrackCont.end(), DCLTrackComp() );

  partial_sort( TrackCont.begin(), TrackCont.end(), 
		TrackCont.end(), DCLTrackComp() );

  // Clear Flags
  {
    int nbefore=TrackCont.size();
    for( int i=0; i<nbefore; ++i ){
      DCLocalTrack *tp=TrackCont[i];
      int nh=tp->GetNHit();
      for( int j=0; j<nh; ++j ){
	tp->GetHit(j)->clearFlags();
      }
    }
  }

#if 0
  {
    int nn=TrackCont.size();
    if(nn>0){
      std::cout << funcname << ": After Sorting. #Tracks = " 
      		<< nn << std::endl;
      
      for( int i=0; i<nn; ++i ){
      //for( int i=0; i<100; ++i ){
	DCLocalTrack *track=TrackCont[i];
	std::cout << std::setw(3) << i << " #Hits="
		  << std::setw(2) << track->GetNHit() 
		  << " ChiSqr=" << track->GetChiSquare()
		  << std::endl;
      }

      std::cout << std::endl;
    }
  }
#endif
  // Delete Duplicated Tracks
  for( int i=0; i<int(TrackCont.size()); ++i ){
    DCLocalTrack *tp=TrackCont[i];
    int nh=tp->GetNHit();
    for( int j=0; j<nh; ++j ) tp->GetHit(j)->setFlags();
    for( int i2=TrackCont.size()-1; i2>i; --i2 ){
      DCLocalTrack *tp2=TrackCont[i2];
      int nh2=tp2->GetNHit(), flag=0;  
      for( int j=0; j<nh2; ++j )
	if( tp2->GetHit(j)->showFlags() ) ++flag;
      if(flag){
	delete tp2;
	TrackCont.erase(TrackCont.begin()+i2);
      }
    }      
  }

  // Clear Flags
  {
    int nbefore=TrackCont.size();
    for( int i=0; i<nbefore; ++i ){
      DCLocalTrack *tp=TrackCont[i];
      int nh=tp->GetNHit();
      for( int j=0; j<nh; ++j ){
	tp->GetHit(j)->clearFlags();
      }
    }
  }  
  
  {
    int nn=TrackCont.size();
    for(int i=0; i<nn; ++i ){
      DCLocalTrack *tp=TrackCont[i];
      int nh=tp->GetNHit();
      for( int j=0; j<nh; ++j ){
	int lnum = tp->GetHit(j)->GetLayer();
	double zz = DCGeomMan::GetInstance().GetLocalZ( lnum );
	tp->GetHit(j)->SetCalPosition(tp->GetX(zz), tp->GetY(zz));
      }
    }
  }
  
#if 0
  {
    int nn=TrackCont.size();
    //    std::cout << funcname << ": After Deleting. #Tracks = " 
    //	      << nn << std::endl;
    
    for( int i=0; i<nn; ++i ){
      DCLocalTrack *track=TrackCont[i];
      if((track->GetChiSquare()>20)){
	std::cout << "" << std::endl;
	std::cout << std::setw(3) <<"ntrack = " << nn << " Track# =  " << i << " #Hits="
		  << std::setw(2) << track->GetNHit() 
		  << " ChiSqr=" << track->GetChiSquare()
		  << " DifVXU=" << track->GetDifVXU()
		  << std::endl;
	std::cout << "" << std::endl;

      }
    }
  }
#endif

  for( int i=0; i<NumOfLayersSdcOut; ++i )
    for_each( CandCont[i].begin(), CandCont[i].end(), DeleteObject() );


  return TrackCont.size();    
}
#endif

//For MWPC
//_____________________________________________________________________________
int MWPCLocalTrackSearch( const DCHitContainer * HC,
			  std::vector <DCLocalTrack *> &TrackCont)
		      
{
  static const std::string funcname = "[MWPCLocalTrackSearch]";

  std::vector < std::vector <DCPairHitCluster *> > CandCont;
  CandCont.resize(NumOfLayersBcIn);

  for( int i=0; i<NumOfLayersBcIn; ++i ){
    MakeMWPCPairPlaneHitCluster( HC[i], CandCont[i] );
  }

  
//   std::vector <int> nCombi(NumOfLayersBcIn);
//   for( int i=0; i<NumOfLayersBcIn; ++i ){ 
//     nCombi[i]=(CandCont[i]).size();

//     // If #Cluster>MaxNumerOfCluster,  error return

//     if(nCombi[i]>MaxNumberOfClusters){
//       std::cout << "#W  too many clusters " << funcname
// 		<< "  layer = " << i << " : " << nCombi[i] << std::endl;
//       clearCandidates(CandCont);
//       return 0;
//     } 
//   }

//   debugPrint(nCombi, funcname);

#if 0
  debugPrint(nCombi, CandCont, funcname);
#endif

  TrackMaker trackMaker(CandCont, MinNumOfHitsBcIn, MaxCombi, MaxChisquare);
  trackMaker.MakeTracks(TrackCont);

  if( TrackCont.empty() ) 
    {
      clearCandidates(CandCont);
      return 0;
    }
  
  // Clear Flags
  clearFlags(TrackCont);


#if 0
  debugPrint(TrackCont, funcname, ": Before Sorting. #Tracks = ");
#endif


  std::partial_sort( TrackCont.begin(), TrackCont.end(), 
		     TrackCont.end(), DCLTrackComp() );

#if 0
  debugPrint(TrackCont, funcname, ": After Sorting. #Tracks = ");
#endif

  // Delete Duplicated Tracks
  deleteDuplicatedTracks(TrackCont);

  calcIntersection(TrackCont);
#if 0
  debugPrint(TrackCont, funcname, ": After Deleting. #Tracks = ");
#endif

  clearCandidates(CandCont);

  return TrackCont.size();
}

//______________________________________________________________________________
int MWPCLocalTrackSearch( const std::vector<std::vector<DCHitContainer> >& hcList,
			  std::vector <DCLocalTrack *> &trackCont)
{
  static const std::string funcname = "[MWPCLocalTrackSearch]";
  for (std::vector<std::vector<DCHitContainer> >::const_iterator 
	 it    = hcList.begin(),
	 itEnd = hcList.end();
       it!=itEnd; ++it)
    {
      const std::vector<DCHitContainer>& l = *it;
      std::vector<DCLocalTrack*> tc;
      MWPCLocalTrackSearch(&l[0], tc);
      trackCont.insert(trackCont.end(), tc.begin(), tc.end());
//       std::cout << " tc " << tc.size()
// 		<< " : " << trackCont.size() << std::endl;
    }
  clearFlags(trackCont);
  deleteDuplicatedTracks(trackCont);
//   std::cout << " trackCont.size() = " << trackCont.size() << std::endl;
  return trackCont.size();
}

//______________________________________________________________________________
int LocalTrackSearchBcOutSdcIn( const DCHitContainer * BcHC,  
				const DCPairPlaneInfo * BcPpInfo,
				const DCHitContainer * SdcHC,  
				const DCPairPlaneInfo * SdcPpInfo,
				int BcNpp, int SdcNpp,
				std::vector <DCLocalTrack *> &TrackCont,
				int MinNumOfHits )
{
  static const std::string funcname = "[LocalTrackSearchBcSdc]";

  std::vector < std::vector <DCPairHitCluster *> > CandCont;
  CandCont.resize(BcNpp+SdcNpp);


  for( int i=0; i<BcNpp; ++i ){
    bool ppFlag=BcPpInfo[i].flag;
    int layer1=BcPpInfo[i].id1, layer2=BcPpInfo[i].id2;
    if(ppFlag) 
      MakePairPlaneHitCluster( BcHC[layer1], BcHC[layer2], 
                               BcPpInfo[i].CellSize, CandCont[i] );
    else
      MakeUnPairPlaneHitCluster( BcHC[layer1], CandCont[i] );
  }

  for( int i=0; i<SdcNpp; ++i ){
    bool ppFlag=SdcPpInfo[i].flag;
    int layer1=SdcPpInfo[i].id1, layer2=SdcPpInfo[i].id2;
    if(ppFlag) 
      MakePairPlaneHitCluster( SdcHC[layer1], SdcHC[layer2], 
                               SdcPpInfo[i].CellSize, CandCont[i+BcNpp] );
    else
      MakeUnPairPlaneHitCluster( SdcHC[layer1], CandCont[i+BcNpp] );
  }

  std::vector <int> nCombi(BcNpp+SdcNpp);

  // NOTE: the original code required exactly one cluster on every
  // pair plane, which discards every event when a plane is dead
  // (BC3 layers 1-2 in E05) or noisy. Enumerate the combinations
  // like the standard LocalTrackSearch instead, with a guard
  // against a combinatorial explosion.
  double totalCombi = 1.;
  for( int i=0; i<BcNpp+SdcNpp; ++i ){
    nCombi[i]=(CandCont[i]).size();
    if( nCombi[i]>MaxNumberOfClusters ) nCombi[i]=0;
    totalCombi *= double(nCombi[i]+1);
  }
  if( totalCombi > 1.e5 ){
    clearCandidates(CandCont);
    return 0;
  }

#if 0
  debugPrint(TrackCont, nCombi, CandCont, funcname);
#endif

  std::vector < std::vector <int> > CombiIndex
    = makeindex( BcNpp+SdcNpp, &nCombi[0] );
  int nnCombi=CombiIndex.size();

#if 0
  std::cout << " ===> " << nnCombi << " combinations will be checked.." 
            << std::endl;
  for( int i=0; i<nnCombi; ++i ){
    for (int j=0;j<BcNpp+SdcNpp;j++) {
      std::cout << CombiIndex[i][j] << " ";
    }
    std::cout << std::endl;
  }
#endif

  for( int i=0; i<nnCombi; ++i ){
    DCLocalTrack *track = MakeTrack( CandCont, &((CombiIndex[i])[0]) );
    if( !track ) continue;
    if( track->GetNHit()>=MinNumOfHits && track->DoFitBcSdc() &&
        track->GetChiSquare()<MaxChisquare ) 
      TrackCont.push_back(track);
    else
      delete track;
  }

  // Clear Flags
  clearFlags(TrackCont);


#if 0
  debugPrint(TrackCont, funcname, ": Before Sorting. #Tracks = ");
#endif

 std::partial_sort( TrackCont.begin(), TrackCont.end(), 
		    TrackCont.end(), DCLTrackComp() );

#if 0
 debugPrint(TrackCont, funcname, ": After Sorting. #Tracks = ");
#endif

 // Delete Duplicated Tracks
 deleteDuplicatedTracks(TrackCont);

 {
   // surveyed BcOut->SKS frame link (cf. DoFitBcSdc)
   double zK18tgt=DCGeomMan::GetInstance().GetLocalZ( 130 )
     -DCGeomMan::GetInstance().GetLocalZ( 0 );
   int nn=TrackCont.size();
   for(int i=0; i<nn; ++i ){
     DCLocalTrack *tp=TrackCont[i];
     int nh=tp->GetNHit();
     for( int j=0; j<nh; ++j ){
       int lnum = tp->GetHit(j)->GetLayer();
       double zz = DCGeomMan::GetInstance().GetLocalZ( lnum );
       if (lnum>=113&&lnum<=124)
	 zz -= zK18tgt;
       tp->GetHit(j)->SetCalPosition(tp->GetX(zz), tp->GetY(zz));
     }
   }
 }

#if 0
 debugPrint(TrackCont, funcname, ": After Deleting. #Tracks = ");
#endif

 clearCandidates(CandCont);
  
 return TrackCont.size();
}

////////////////////////////
bool MakePairPlaneHitCluster( const DCHitContainer & HC1,
			      const DCHitContainer & HC2,
			      double CellSize,
			      std::vector <DCPairHitCluster *> & Cont )
{
  int nh1=HC1.size(), nh2=HC2.size();
  std::vector <int> UsedFlag(nh2,0);   

  for( int i1=0; i1<nh1; ++i1 ){
    DCHit *hit1=HC1[i1];

    double wp1=hit1->GetWirePosition();
    bool flag=false;
    for( int i2=0; i2<nh2; ++i2 ){
      DCHit *hit2=HC2[i2];
      double wp2=hit2->GetWirePosition();
      if( fabs(wp1-wp2)<CellSize ){

	int multi1 = hit1->GetDriftLengthSize();
	int multi2 = hit2->GetDriftLengthSize();
	for (int m1=0; m1<multi1; m1++) {
	  if( !(hit1->rangecheck(m1)) ){
	    continue;
	  }
	  for (int m2=0; m2<multi2; m2++) {
	    if( !(hit2->rangecheck(m2)) ){
	      continue;
	    }
	    double x1,x2;
	    if( wp1<wp2 ){
	      x1=wp1+hit1->GetDriftLength(m1);
	      x2=wp2-hit2->GetDriftLength(m2);
	    }
	    else {
	      x1=wp1-hit1->GetDriftLength(m1);
	      x2=wp2+hit2->GetDriftLength(m2);
	    }
	    Cont.push_back( new DCPairHitCluster( new DCLTrackHit(hit1,x1,m1),
						  new DCLTrackHit(hit2,x2,m2) ) );
	    flag=true; ++UsedFlag[i2];
	  }
	}
      }
    }      
#if 1
    if(!flag){
      int multi1 = hit1->GetDriftLengthSize();
      for (int m1=0; m1<multi1; m1++) {
	if( !(hit1->rangecheck(m1)) ) continue;
	double dl=hit1->GetDriftLength(m1);
	Cont.push_back( new DCPairHitCluster( new DCLTrackHit(hit1,wp1+dl,m1) ) );
	Cont.push_back( new DCPairHitCluster( new DCLTrackHit(hit1,wp1-dl,m1) ) );
      }
    }
#endif
  }
#if 1
  for( int i2=0; i2<nh2; ++i2 ){
    if( UsedFlag[i2]==0 ) {
      DCHit *hit2=HC2[i2];
      int multi2 = hit2->GetDriftLengthSize();
      for (int m2=0; m2<multi2; m2++) {
	if( !(hit2->rangecheck(m2)) ) continue;
	
	double wp=hit2->GetWirePosition();
	double dl=hit2->GetDriftLength(m2);
	Cont.push_back( new DCPairHitCluster( new DCLTrackHit(hit2,wp+dl,m2) ) );
	Cont.push_back( new DCPairHitCluster( new DCLTrackHit(hit2,wp-dl,m2) ) );
      }
    }
  }
#endif
  return true;
}

bool MakeUnPairPlaneHitCluster( const DCHitContainer & HC,
				std::vector <DCPairHitCluster *> & Cont )
{
  int nh=HC.size();

  for( int i=0; i<nh; ++i ){
    DCHit *hit=HC[i];
    if( hit ){
      int multi = hit->GetDriftLengthSize();
      for (int m=0; m<multi; m++) {
	if( !(hit->rangecheck(m)) ) continue;

	double wp=hit->GetWirePosition();
	double dl=hit->GetDriftLength(m);

	Cont.push_back( new DCPairHitCluster( new DCLTrackHit(hit,wp+dl,m) ) );
	Cont.push_back( new DCPairHitCluster( new DCLTrackHit(hit,wp-dl,m) ) );
      }
    }
  }

  return true;
}

bool MakeMWPCPairPlaneHitCluster( const DCHitContainer & HC,
				  std::vector <DCPairHitCluster *> & Cont )
{
  int nh=HC.size();
  
  for( int i=0; i<nh; ++i ){
    DCHit *hit=HC[i];
    if( hit ){
      int multi = hit->GetDriftTimeSize();
      for (int m=0; m<multi; m++) {
	if( !(hit->rangecheck(m)) ) continue;
	double wp=hit->GetWirePosition();
	double dl=hit->GetDriftLength(m);
	Cont.push_back( new DCPairHitCluster( new DCLTrackHit(hit,wp,m) ) );
      }
    }
  }

  return true;
}

std::vector< std::vector<int> > makeindex( int ndim, const int *index1 )
{
  if(ndim==1){
    std::vector< std::vector<int> > index2;
    for( int i=-1; i<index1[0]; ++i ){
      std::vector <int> elem(1,i);
      index2.push_back(elem);
    }
    return index2;
  }

  std::vector< std::vector<int> > 
    index2=makeindex( ndim-1, index1+1 );

  std::vector< std::vector<int> > index;
  int n2=index2.size();
  for( int j=0; j<n2; ++j ){
    for( int i=-1; i<index1[0]; ++i ){
      std::vector <int> elem;
      int n3=index2[j].size();
      elem.reserve(n3+1);
      elem.push_back(i);
      for( int k=0; k<n3; ++k )
	elem.push_back(index2[j][k]);
      index.push_back(elem);
      int size1=index.size();
    }
  }

  return index;
}  

std::vector< std::vector<int> > 
makeindex_SdcOut( int ndim_org, int minimumHit, int ndim, const int *index1 )
{
  if(ndim==1){
    std::vector< std::vector<int> > index2;
    for( int i=-1; i<index1[0]; ++i ){
      std::vector <int> elem(1,i);
      index2.push_back(elem);
    }
    return index2;
  }

  std::vector< std::vector<int> > 
    index2=makeindex_SdcOut( ndim_org, minimumHit, ndim-1, index1+1 );
 
  std::vector< std::vector<int> > index;
  int n2=index2.size();
  for( int j=0; j<n2; ++j ){
    for( int i=-1; i<index1[0]; ++i ){
      std::vector <int> elem;
      int validHitNum=0;
      int n3=index2[j].size();
      elem.reserve(n3+1);
      elem.push_back(i);
      if (i != -1)
	validHitNum++;
      for( int k=0; k<n3; ++k ) {
        elem.push_back(index2[j][k]);
        if (index2[j][k] != -1)
          validHitNum++;
      }
      if (ndim==ndim_org) {
        if (validHitNum >= minimumHit)
          index.push_back(elem);
      } else {
        index.push_back(elem);
      }
      int size1=index.size();
    }
  }

  return index;
}

/////////////////////////
/////XUV Tracking routine
/////////////////////////
bool MakePairPlaneHitClusterVUX( const DCHitContainer & HC1,
				 const DCHitContainer & HC2,
				 double CellSize,
				 std::vector <DCPairHitCluster *> & Cont )
{
  int nh1=HC1.size(), nh2=HC2.size();
  std::vector <int> UsedFlag(nh2,0);   
  
  for( int i1=0; i1<nh1; ++i1 ){
    DCHit *hit1=HC1[i1];

    double wp1=hit1->GetWirePosition();
    bool flag=false;
    for( int i2=0; i2<nh2; ++i2 ){
      DCHit *hit2=HC2[i2];
      double wp2=hit2->GetWirePosition();
      if( fabs(wp1-wp2)<CellSize ){

	int multi1 = hit1->GetDriftLengthSize();
	int multi2 = hit2->GetDriftLengthSize();
	for (int m1=0; m1<multi1; m1++) {
	  if( !(hit1->rangecheck(m1)) ) continue;
	  for (int m2=0; m2<multi2; m2++) {
	    if( !(hit2->rangecheck(m2)) ) continue;
	    double dl1=hit1->GetDriftLength(m1);
	    double dl2=hit2->GetDriftLength(m2);

	    Cont.push_back( new DCPairHitCluster( new DCLTrackHit(hit1,wp1+dl1,m1),
						  new DCLTrackHit(hit2,wp2+dl2,m2) ) );
	    Cont.push_back( new DCPairHitCluster( new DCLTrackHit(hit1,wp1+dl1,m1),
						  new DCLTrackHit(hit2,wp2-dl2,m2) ) );
	    Cont.push_back( new DCPairHitCluster( new DCLTrackHit(hit1,wp1-dl1,m1),
						  new DCLTrackHit(hit2,wp2+dl2,m2) ) );
	    Cont.push_back( new DCPairHitCluster( new DCLTrackHit(hit1,wp1-dl1,m1),
						  new DCLTrackHit(hit2,wp2-dl2,m2) ) );

	    flag=true; ++UsedFlag[i2];
	  }
	}
      }
    }
    if(!flag){
      int multi1 = hit1->GetDriftLengthSize();
      for (int m1=0; m1<multi1; m1++) {
	if( !(hit1->rangecheck(m1)) ) continue;
	double dl=hit1->GetDriftLength(m1);
	Cont.push_back( new DCPairHitCluster( new DCLTrackHit(hit1,wp1+dl,m1) ) );
	Cont.push_back( new DCPairHitCluster( new DCLTrackHit(hit1,wp1-dl,m1) ) );
      }
    }      
  }
  for( int i2=0; i2<nh2; ++i2 ){
    if( UsedFlag[i2]==0 ) {
      DCHit *hit2=HC2[i2];
      int multi2 = hit2->GetDriftLengthSize();
      for (int m2=0; m2<multi2; m2++) {
	if( !(hit2->rangecheck(m2)) ) continue;
	
	double wp=hit2->GetWirePosition();
	double dl=hit2->GetDriftLength(m2);
	Cont.push_back( new DCPairHitCluster( new DCLTrackHit(hit2,wp+dl,m2) ) );
	Cont.push_back( new DCPairHitCluster( new DCLTrackHit(hit2,wp-dl,m2) ) );
      }
    }
  }

  return true;
}

std::vector< std::vector<int> > makeindex_VXU( int ndim,int maximumHit, const int *index1 )
{
  if(ndim==1){
    std::vector< std::vector<int> > index2;
    for( int i=-1; i<index1[0]; ++i ){
      std::vector <int> elem(1,i);
      index2.push_back(elem);
    }
    return index2;
  }

  std::vector< std::vector<int> > 
    index2=makeindex_VXU( ndim-1, maximumHit, index1+1 );

  std::vector< std::vector<int> > index;
  int n2=index2.size();
  for( int j=0; j<n2; ++j ){
    for( int i=-1; i<index1[0]; ++i ){
      std::vector <int> elem;
      int validHitNum=0;
      int n3=index2[j].size();
      elem.reserve(n3+1);
      elem.push_back(i);
      if (i != -1)
	validHitNum++;
      for( int k=0; k<n3; ++k ){
	elem.push_back(index2[j][k]);
        if (index2[j][k] != -1)
          validHitNum++;
      }
      if (validHitNum <= maximumHit)
	index.push_back(elem);
      int size1=index.size();
    }
  }

  return index;
}

std::vector< std::vector<int> > 
makeindex_SdcOut_below( int ndim_org, int maximumHit, int ndim, const int *index1 )
{
  if(ndim==1){
    std::vector< std::vector<int> > index2;
    for( int i=-1; i<index1[0]; ++i ){
      std::vector <int> elem(1,i);
      index2.push_back(elem);
    }
    return index2;
  }

  std::vector< std::vector<int> > 
    index2=makeindex_SdcOut_below( ndim_org, maximumHit, ndim-1, index1+1 );
 
  std::vector< std::vector<int> > index;
  int n2=index2.size();
  for( int j=0; j<n2; ++j ){
    for( int i=-1; i<index1[0]; ++i ){
      std::vector <int> elem;
      int validHitNum=0;
      int n3=index2[j].size();
      elem.reserve(n3+1);
      elem.push_back(i);
      if (i != -1)
	validHitNum++;
      for( int k=0; k<n3; ++k ) {
        elem.push_back(index2[j][k]);
        if (index2[j][k] != -1)
          validHitNum++;
      }
      if (ndim==ndim_org) {
        if ((validHitNum <= maximumHit) && (validHitNum>0))
          index.push_back(elem);
      } else {
        index.push_back(elem);
      }
      int size1=index.size();
    }
  }

  return index;
}

std::vector< std::vector<int> > 
makeindex_BcIn( int ndim_org, int minimumHit, int ndim, const int *index1 )
{
  if(ndim==1){
    std::vector< std::vector<int> > index2;
    for( int i=-1; i<index1[0]; ++i ){
      std::vector <int> elem(1,i);
      index2.push_back(elem);
    }
    return index2;
  }

  std::vector< std::vector<int> > 
    index2=makeindex_BcIn( ndim_org, minimumHit, ndim-1, index1+1 );
 
  std::vector< std::vector<int> > index;
  int n2=index2.size();
  for( int j=0; j<n2; ++j ){
    for( int i=-1; i<index1[0]; ++i ){
      std::vector <int> elem;
      int validHitNum=0;
      int n3=index2[j].size();
      elem.reserve(n3+1);
      elem.push_back(i);
      if (i != -1)
	validHitNum++;
      for( int k=0; k<n3; ++k ) {
        elem.push_back(index2[j][k]);
        if (index2[j][k] != -1)
          validHitNum++;
      }
      if (ndim==ndim_org) {
        if (validHitNum >= minimumHit)
          index.push_back(elem);
      } else {
        index.push_back(elem);
      }
      int size1=index.size();
    }
  }

  return index;
}  

DCLocalTrack *MakeTrack(  const std::vector < std::vector <DCPairHitCluster *> > &CandCont,
	      const int *combination )
{
  static const std::string funcname = "[MakeTrack]";

  DCLocalTrack *tp=new DCLocalTrack();
  //DCLocalTrack *tp ;

  if(!tp){
    std::cerr << funcname << ": new fail" << std::endl;
    return 0;
  }    

  int n=CandCont.size();

  for( int i=0; i<n; ++i ){
    int m=combination[i];
    DCPairHitCluster *cluster=0;
    if(m>=0) cluster=CandCont[i][m];
#if 0
    std::cout << funcname << ":" << std::setw(3)
	      << i << std::setw(3) << m  << " "
	      << CandCont[i][m] << std::endl; 
#endif


#if 0

    /*****************************************************************/
    // Only SdcIn & BcOut Local Tracking !!                          
    /*****************************************************************/
    int mm = 0 ;
    //std::cout << "-------------------------------------" << std::endl;
    if(cluster){
     mm = cluster->NumberOfHits();
     if(mm>2) std::cout << "mm = " << mm << std::endl;

     double DL[mm] ;
     double DT[mm] ;
     int Layer[mm] ;

      for(int j=0; j<mm; ++j ){
	DCLTrackHit *hitp=cluster->GetHit(j);
#if 0
	std::cout << "mm = " << mm 
	          << " Layer = " << hitp->GetLayer() 
		  << " Drift Length = " << hitp->GetDriftLength()
		  << std::endl;
#endif
	DL[j] = hitp->GetDriftLength();
	DT[j] = hitp->GetDriftTime();
	Layer[j] = hitp->GetLayer();

	if(mm<2 || !((Layer[j]>0 && Layer[j]<11) || (Layer[j]>112 && Layer[j]<125))){
	  if((Layer[j]>=1 && Layer[j]<=4) && (DT[j]>-10 && DT[j]<40))
	    if(hitp) tp->AddHit( hitp );
	  if((Layer[j]>=5 && Layer[j]<=10) && (DT[j]>-10 && DT[j]<60))
	    if(hitp) tp->AddHit( hitp );
	  if((Layer[j]>=113 && Layer[j]<=118) && (DT[j]>-10 && DT[j]<40))
	    if(hitp) tp->AddHit( hitp );
	  if((Layer[j]>=119 && Layer[j]<=124) && (DT[j]>-10 && DT[j]<60))
	    if(hitp) tp->AddHit( hitp );
	}
      }
      //std::cout << "***********************************" << std::endl;
      if(mm==2 && ((Layer[0]>0 && Layer[0]<11) || (Layer[0]>112 && Layer[0]<125))){
	int Ok = 0 ;
	/*	
	std::cout << "DL[ " << Layer[0] << "] = " << DL[0] 
		  << " DL[ " << Layer[1] << "] = " << DL[1] << std::endl; 
	*/
	//std::cout << "***********************************" << std::endl;
	if(Layer[0] ==1 || Layer[0]==3 || Layer[0]==113 || Layer[0]==115 || Layer[0]==117){
	  if(((DL[0]+DL[1])>0.7) && ((DL[0]+DL[1])<2.3)) Ok = 1 ;
	}  
	if(Layer[0] ==5 || Layer[0]==7 || Layer[0]==9 || Layer[0]==119 || Layer[0]==121 || Layer[0]==123){
	  if(((DL[0]+DL[1])>1.5) && ((DL[0]+DL[1])<3.5)) Ok = 1 ;
	}
	if(Ok){
	  for(int j=0; j<mm; ++j ){
	    DCLTrackHit *hitp=cluster->GetHit(j);
	    if(hitp) tp->AddHit( hitp );
	  }
	}
      }
    }
#endif


#if 1
    if(cluster){
     int mm=cluster->NumberOfHits();

      for(int j=0; j<mm; ++j ){
	DCLTrackHit *hitp=cluster->GetHit(j);

	  if(hitp) tp->AddHit( hitp );
      }
    }
#endif

  }
  return tp;
}


////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////
//SDC3&4 (**Calculate all hit combinations -> Too slow !!**)

#if ver2
int SdcOutLocalTrackSearch( const DCHitContainer * HC,
			    std::vector <DCLocalTrack *> &TrackCont)
		      
{
  static const std::string funcname = "[LocalTrackSearch]";
  
  std::vector < std::vector <DCPairHitCluster *> > CandCont;
  CandCont.resize(NumOfLayersSdcOut);



  //Reject X plane
  int Layer_A=-2 ,Layer_B=-2 , Layer_C=-2 ,Layer_D=-2 ,Layer_E=-2 , Layer_F=-2 ;

#if 1
  int W=0 ,NHIT[NumOfLayersSdcOut] ;

  for( int i=0; i<NumOfLayersSdcOut; ++i ){
    int nh ;
    nh = HC[i].size();
    if(nh==0 && W==5){
      Layer_F = i ;
      W++ ;
    }
    if(nh==0 && W==4){
      Layer_E = i ;
      W++ ;
    }

    if(nh==0 && W==3){
      Layer_D = i ;
      W++ ;
    }

    if(nh==0 && W==2){
      Layer_C = i ;
      W++ ;
    }

    if(nh==0 && W==1){
      Layer_B = i ;
      W++ ;
    }

    if(nh==0 && W==0){
      Layer_A = i ;
      W++ ;
    }

    NHIT[i] = nh ;

#if 0
    if(nh==0) std::cout << "#" << i+1 <<"  nh =  " << nh << " W = " << W << std::endl ;
#endif

  }

#if 0
  if(W>0 && W<4){
    for( int i=1;i<NumOfLayersSdcOut+1; ++i ){
      if((i!=Layer_A) && (i!=Layer_A +1) && (i!=Layer_A -1)){
	if((i!=Layer_B) && (i!=Layer_B +1) && (i!=Layer_B -1)){
	  if((i!=Layer_C) && (i!=Layer_C +1) && (i!=Layer_C -1)){
	    if((i!=Layer_D) && (i!=Layer_D +1) && (i!=Layer_D -1)){
	      if(Layer_A!=-2 && Layer_B!=-2 && Layer_C!=-2 && Layer_D!=-2)
		break ; 
	      if(Layer_A!=-2 && Layer_B!=-2 && Layer_C!=-2 && Layer_D==-2)
		Layer_D = i-1 ;
              if(Layer_A!=-2 && Layer_B!=-2 && Layer_C==-2)
		Layer_C = i-1 ;
              if(Layer_A!=-2 && Layer_B==-2)
		Layer_B = i-1 ;
	      if(Layer_A==-2)
		Layer_A = i-1 ;
	    }
	  }
	}
      }
    }
  }
#endif


  if(W==0){
    //Layer_F = 10 ;//SDC4 X2
    //Layer_E = 8  ;//SDC4 U2
    Layer_D = 9  ;//SDC4 V2
    Layer_C = 7  ;//SDC4 X1
    Layer_B = 4  ;//SDC3 X2
    Layer_A = 2  ;//SDC3 U1
  }

#if 1
  if(W>0 && W<4){

    /* Calcurate the average nh (not include nh=0 plane)   */

    int Sum_nh=0 , Num_nh=0 ;
    
    for(int i=0;i<NumOfLayersSdcOut; i++){
      
      int nh ;
      nh = HC[i].size();
      
      if(nh!=0){
	Sum_nh = nh ;
	Num_nh = Num_nh +1 ;
      }
    }

    double Ave_nh = Sum_nh/Num_nh ;



    for(int i=0;i<NumOfLayersSdcOut;i++){
      
      if((i!=Layer_A -1) && (i!=Layer_A) && (i!=Layer_A +1) && (i!=Layer_B -1) && (i!=Layer_B) && (i!=Layer_B +1) && (i!=Layer_C -1) && (i!=Layer_C) && (i!=Layer_C +1) && (i!=Layer_D -1) && (i!=Layer_D) && (i!=Layer_D +1) && (i!=Layer_E -1) && (i!=Layer_E) && (i!=Layer_E +1)){
      
	int nh ;
	nh = HC[i].size();

	/*
	if(nh>Ave_nh+2 && W==5){
	  Layer_F = i ;
	  W++ ;
	  //std::cout << "nh>Ave_nh+2 (fourth) = " << Layer_D << std::endl ;
	}
	if(nh>Ave_nh+2 && W==4){
	  Layer_E = i ;
	  W++ ;
	  //std::cout << "nh>Ave_nh+2 (fourth) = " << Layer_D << std::endl ;
	} 
	*/ 
	if(nh>Ave_nh+2 && W==3){
	  Layer_D = i ;
	  W++ ;
	  //std::cout << "nh>Ave_nh+2 (thrid) = " << Layer_C << std::endl ;
	}  
	if(nh>Ave_nh+2 && W==2){
	  Layer_C = i ;
	  W++ ;
	  //std::cout << "nh>Ave_nh+2 (second) = " << Layer_B << std::endl ;
	}  
	if(nh>Ave_nh+2 && W==1){
	  Layer_B = i ;
	  W++ ;
	  //std::cout << "nh>Ave_nh+2 (first) = " << Layer_A << std::endl ;
	}  
      }
    }



    for(int i=0;i<NumOfLayersSdcOut;i++){
      
      if((i!=Layer_A -1) && (i!=Layer_A) && (i!=Layer_A +1) && (i!=Layer_B -1) && (i!=Layer_B) && (i!=Layer_B +1) && (i!=Layer_C -1) && (i!=Layer_C) && (i!=Layer_C +1)  && (i!=Layer_D -1) && (i!=Layer_D) && (i!=Layer_D +1) && (i!=Layer_E -1) && (i!=Layer_E) && (i!=Layer_E +1)){

	/*
	if(W==5){
	  Layer_F = i ;
	  W++ ;
	}	
	if(W==4){
	  Layer_E = i ;
	  W++ ;
	} 
	*/ 
	if(W==3){
	  Layer_D = i ;
	  W++ ;
	} 
	if(W==2){
	  Layer_C = i ;
	  W++ ;
	}  
	if(W==1){
	  Layer_B = i ;
	  W++ ;
	}  
      }
    }
  }
#endif

  
#if 0
  if(W>=0){
  std::cout << "Layer_A = " << Layer_A+1 << " Layer_B = " << Layer_B+1 
	    << " Layer_C = " << Layer_C+1 << " Layer_D = " << Layer_D+1 
	    << " Layer_E = " << Layer_E+1 << " Layer_F = " << Layer_F+1 << std::endl;
  }
#endif 

#endif

 
  for( int i=0; i<NumOfLayersSdcOut; ++i ){//NumOfLayersSdcOut=12

    if((i!=Layer_A) && (i!=Layer_B) && (i!=Layer_C) && (i!=Layer_D) && (i!=Layer_E) && (i!=Layer_F))
      MakeUnPairPlaneHitCluster( HC[i], CandCont[i] );//Go to bottom
    
  }
  

  
  std::vector <int> nCombi(NumOfLayersSdcOut);
  for( int i=0; i<NumOfLayersSdcOut; ++i ){ 
    nCombi[i]=(CandCont[i]).size();//always nCombi[i]=2 (left & right by 1 hit)
    
#if 0
    std::cout << "i=  " << i << "   nCombi = " << nCombi[i] << std::endl;
#endif 
    
    if(nCombi[i]>MaxNumberOfClusters){//MaxNumberOfClusters=10
      for( int i=0; i<NumOfLayersSdcOut; ++i )
	for_each( CandCont[i].begin(), CandCont[i].end(), DeleteObject() );
      return 0;
    } 
  }
  
  
  /*--------HitWire & Hit Local position of each plane----*/
#if 0
  std::cout << funcname << ": #Hits of each group" << std::endl;
  for( int i=0; i<NumOfLayersSdcOut; ++i ) std::cout << std::setw(4) << nCombi[i];
  std::cout << std::endl;
  for( int i=0; i<NumOfLayersSdcOut; ++i ){
    int n=CandCont[i].size();
    std::cout << "[" << std::setw(3) << i << "]: "
	      << std::setw(3) << n << " ";
    for( int j=0; j<n; ++j ){
      std::cout << ((DCLTrackHit *)CandCont[i][j]->GetHit(0))->GetWire() << " "
                << ((DCLTrackHit *)CandCont[i][j]->GetHit(0))->GetLocalHitPos() << "  :  " ;
    }
    std::cout << std::endl;
  }
#endif
  
  std::vector < std::vector <int> > 
    CombiIndex = makeindex_SdcOut( NumOfLayersSdcOut, 
				   MinNumOfHitsSdcOut, 
				   NumOfLayersSdcOut, 
				   &(nCombi[0]) );//go to bottom
  int nnCombi=CombiIndex.size();
  
#if 0
  std::cout << " ===> " << nnCombi << " combinations will be checked.." 
	    << std::endl;
#endif
  
  for( int i=0; i<nnCombi; ++i ){
    DCLocalTrack *track = MakeTrack( CandCont, &((CombiIndex[i])[0]) );
    if( !track ) continue;
    
    if( track->GetNHit()>=MinNumOfHitsSdcOut && track->DoFit() &&
    	track->GetChiSquare()<MaxChisquare ){
      
      TrackCont.push_back(track);
      double chisqr = track->GetChiSquare();
      
    }
    else{
      delete track;
    }
  }
  
  // Clear Flags
  int nbefore=TrackCont.size();
  for( int i=0; i<nbefore; ++i ){
    DCLocalTrack *tp=TrackCont[i];
    int nh=tp->GetNHit();
    for( int j=0; j<nh; ++j ) tp->GetHit(j)->clearFlags();
  }
  
  /*-------------------------------------------------------*/


#if 1 //Add residual plane!!
  {
    int nn=TrackCont.size();    
    
    for( int ii=0; ii<nn; ++ii ){
      DCLocalTrack *track=TrackCont[ii];

#if 0
      std::cout << std::setw(3) << "nn = " << nn << " ii = " << ii << " #Hits="
		<< std::setw(2) << track->GetNHit() 
		<< " ChiSqr=" << track->GetChiSquare()
		<< std::endl;
#endif
      
      
      for( int i=0; i<NumOfLayersSdcOut; ++i ){
	
	if((i==Layer_A) || (i==Layer_B) || (i==Layer_C)  || (i==Layer_D) || (i==Layer_E) || (i==Layer_F)){
	  int nh ;
	  nh = HC[i].size();
	  double Del_A=0 ,Del_B=0 ,Del_1[nh],Del_2[nh] ,Wp_1[nh],Dl_1[nh],Wp_2[nh],Dl_2[nh] ;
	  
	  int I1[nh],I2[nh] ,IA=0 ,IB=0 ,II;

	  DCLTrackHit *hitA[nh] ;
	  DCLTrackHit *hitB[nh] ;

	  for(int j=0;j<nh;++j ){
 	    DCHit * hit = HC[i][j];

#if 0
	    std::cout << "nh = " << nh << " ------j(nh) = "<< j << " ------" << std::endl ;
#endif

 	    if((hit) && (nh!=0)){

	      int multi = hit->GetDriftLengthSize();
	      
	      double wp[nh][multi] , dl[nh][multi] ,theta[nh][multi],th[nh][multi] , th1[nh][multi] ;
	      double del_1[nh][multi] , del_2[nh][multi] ;
	      
	      for (int m=0; m<multi; m++) { 
		if( !(hit->rangecheck(m)) ) continue;
		
		wp[j][m]=hit->GetWirePosition();
		dl[j][m]=hit->GetDriftLength(m);
                theta[j][m]=hit->GetTiltAngle();
		
		double Z = DCGeomMan::GetInstance().GetLocalZ(31+i);
		
                th[j][m] = (theta[j][m]/180)*3.14159 ; 
		
#if 0 
		if(dl[j][m]<0.001){
		  std::cout << "i = " << i << " nh = " << nh << " j = " << j << " multi = " << multi 
			    << " m = " << m << std::endl;
		  std::cout << " wp = "  
			    << (wp[j][m]+dl[j][m])               
			    << " dl = " 
			    << dl[j][m]  << std::endl ;
		  std::cout << "(wp+dl) = "  
			    << (wp[j][m]+dl[j][m])               
			    << " (wp-dl) = " 
			    << (wp[j][m]-dl[j][m])  << std::endl ; 
		  std::cout << " |X - (wp+dl)| = "  
			    << fabs(track->GetX(Z)*cos(th[j][m])+track->GetY(Z)*sin(th[j][m]) - (wp[j][m]+dl[j][m]))               
			    << " |X - (wp-dl)| = " 
			    << fabs(track->GetX(Z)*cos(th[j][m])+track->GetY(Z)*sin(th[j][m]) - (wp[j][m]-dl[j][m]))
			    << std::endl ;
		}
#endif
		
		if((fabs(wp[j][m])>pow(10,-10))){

		del_1[j][m] = fabs(track->GetX(Z)*cos(th[j][m])+track->GetY(Z)*sin(th[j][m]) - (wp[j][m]+dl[j][m]));   
		del_2[j][m] = fabs(track->GetX(Z)*cos(th[j][m])+track->GetY(Z)*sin(th[j][m]) - (wp[j][m]-dl[j][m]));
	      }


	      else{

		del_1[j][m] = 0.0;   
		del_2[j][m] = 0.0;
	      }
	      } 
	      
	      
	      for(int m=0; m<multi; m++){
		if( (hit->rangecheck(m)) ) continue;

		wp[j][m] = 0.0 ; 
		dl[j][m] = 0.0 ;
		del_1[j][m] = 0.0 ;
		del_2[j][m] = 0.0 ;

#if 0
		std::cout << "hit->rangecheck(m)" << std::endl ;
#endif

	      }
	      

	      Del_1[j]=del_1[j][0],Del_2[j]=del_2[j][0] ;
	      
              if(multi==1){

		Dl_1[j] = dl[j][0]; 
		Wp_1[j] = wp[j][0]; 
		
		Dl_2[j] = dl[j][0]; 
		Wp_2[j] = wp[j][0]; 
		
		I1[j] = 0 ;
		I2[j] = 0 ;
		
#if 0
		std::cout << "dl[j][0] = " << dl[j][0] << std::endl ;
		std::cout << "wp[j][0] = " << wp[j][0] << std::endl ;
		std::cout << "del_1[j][0] = " << del_1[j][0] << std::endl ;
		std::cout << "del_2[j][0] = " << del_2[j][0] << std::endl ;
		std::cout << "Wp_1[j]+Dl_1[j] = " << Wp_1[j]+Dl_1[j] << std::endl ;
		std::cout << "Wp_2[j]-Dl_2[j] = " << Wp_2[j]-Dl_2[j] << std::endl ;
#endif
		
	      }
	      
	      if(multi>1){

		//Del_1[j]=del_1[j][0],Del_2[j]=del_2[j][0] ;

		for(int m=0; m<multi-1 ;m++){
		  for(int M=m+1 ; M<multi ;M++){
		    if((del_1[j][M]>pow(10,-10)) && (Del_1[j] > del_1[j][M])){
		      Del_1[j] = del_1[j][M] ;
		      I1[j] = M ; 
		    }
		  }
		}		
		Dl_1[j] = dl[j][I1[j]] ;
		Wp_1[j] = wp[j][I1[j]] ;
		
		for(int m=0; m<multi-1 ;m++){
		  for(int M=m+1 ; M<multi ;M++){
		    if((del_2[j][M]>pow(10,-10)) && (Del_2[j] > del_2[j][M])){
		      Del_2[j] = del_2[j][M] ;
		      I2[j] = M ;
		    }
		  }
		}

		Dl_2[j] = dl[j][I2[j]] ;
		Wp_2[j] = wp[j][I2[j]] ;

#if 0
		std::cout << " j = " << j << " Wp_1[j]+Dl_1[j] = " << Wp_1[j]+Dl_1[j] << std::endl ;	
		std::cout << " j = " << j << " Wp_2[j]-Dl_2[j] = " << Wp_2[j]+Dl_2[j] << std::endl ;
#endif

	      }
	    }

	    if(!((hit) && (nh!=0))){

	      Del_1[j] = 0.0 ;
	      Del_2[j] = 0.0 ; 
	      Wp_1[j] = 0.0 ;
	      Dl_1[j] = 0.0 ;
	      Wp_2[j] = 0.0 ;
	      Dl_2[j] = 0.0 ;

#if 0
	      std::cout << "!((hit) && (nh!=0))" << std::endl ;	      
#endif

	    }

	  }

	  if(nh>0){ 
	    for(int I=0; I<nh; I++){
	      if(Del_1[I]>pow(10,-10)){     
		Del_A = Del_1[I];
                IA = I ;
                break ;
	      }
	    }
	    for(int I=0; I<nh; I++){
	      if(Del_2[I]>pow(10,-10)){     
		Del_B = Del_2[I];
                IB = I ;
                break ;
	      }
	    }


#if 0
	    std::cout << "----------------------------" << std::endl ;

	    std::cout << "Del_A = " << Del_A << std::endl ;
	    std::cout << "Del_B = " << Del_B << std::endl ;
#endif


	    if(nh>1){

	      for(int m=0; m<nh-1 ;m++){
		for(int M=m+1 ; M<nh ;M++){
		  if((Del_1[M]>pow(10,-10)) && (Del_A > Del_1[M])){
		    Del_A = Del_1[M] ;
		    IA = M ; 
		  }
		}
	      }
	      
	      for(int m=0; m<nh-1 ;m++){
		for(int M=m+1 ; M<nh ;M++){
		  if((Del_2[M]>pow(10,-10)) && (Del_B > Del_2[M])){
		    Del_B = Del_2[M] ;
		    IB = M ;
		  }
		}
	      }
	    }

#if 0
	    std::cout << "********************************" << std::endl ;
	    std::cout << "IA = " << IA << " IB = " << IB << std::endl ;
#endif


            if((Del_A>pow(10,-10)) && (Del_B>pow(10,-10))){

	      if((Del_A<=Del_B) && (Del_A<2.0)){
		DCHit *Hit_A = HC[i][IA] ;

#if 0
		std::cout << "I1[IA] = " << I1[IA] << "Del_A = " << Del_A << std::endl ;
		std::cout << "Wp_1[IA]+Dl_1[IA] = " << Wp_1[IA]+Dl_1[IA] << std::endl;
		std::cout << "wp = " << Hit_A->GetWirePosition() << " dl = " << Hit_A->GetDriftLength(I1[IA]) << std::endl;
		std::cout << "wp+dl = " << Hit_A->GetWirePosition()+Hit_A->GetDriftLength(I1[IA]) << std::endl; 
#endif

		DCLTrackHit *hit_A = new DCLTrackHit(Hit_A,Wp_1[IA]+Dl_1[IA],I1[IA]);

		track->AddHit(hit_A);
	      }

	      if((Del_A>Del_B) && (Del_B<2.0)){ 
		DCHit *Hit_B = HC[i][IB] ; 

#if 0
		std::cout << "I2[IB] = " << I2[IB] << " Del_B = " << Del_B << std::endl ;
		std::cout << "Wp_2[IB]-Dl_2[IB] = " << Wp_2[IB]-Dl_2[IB] << std::endl ;
		std::cout << "wp = " << Hit_B->GetWirePosition() << " dl = " << Hit_B->GetDriftLength(I2[IB]) << std::endl;
		std::cout << "wp-dl = " << Hit_B->GetWirePosition()-Hit_B->GetDriftLength(I2[IB]) << std::endl;
#endif

		DCLTrackHit *hit_B = new DCLTrackHit(Hit_B,Wp_2[IB]-Dl_2[IB],I2[IB]);
		
		track->AddHit(hit_B);
	      }
	    
	    }
	  }
	}
      }
    }
  
    
    for(int i=0; i<nn; ++i ){
      DCLocalTrack *track=TrackCont[i];
   
      if( !track ) continue;
      track->DoFit();

#if 0
      if((track->GetChiSquare()>0)) std::cout << "i(nn) = " << i 
						<< " Chisquare = " << track->GetChiSquare() 
						<< " N Hit = " << track->GetNHit() << std::endl ;
#endif

    }
    
    {
      for(int i=0; i<nn; ++i ){
	DCLocalTrack *tp=TrackCont[i];
	int nh=tp->GetNHit();

#if 0
	if(nh<12) std::cout << "---------- N Hit = " << nh << " ----------"<< std::endl ;
#endif


	for( int j=0; j<nh; ++j ){
	  int lnum = tp->GetHit(j)->GetLayer();
	  double zz = DCGeomMan::GetInstance().GetLocalZ( lnum );
	  tp->GetHit(j)->SetCalPosition(tp->GetX(zz), tp->GetY(zz));
	}
      }
    }    
  }    
#endif


  /*-------------------------------------------------------*/



#if 0
  {
    int nn=TrackCont.size();
    std::cout << funcname << ": Before Sorting. #Tracks = " 
	      << nn << std::endl;
    if(nn>0){

      if(nn>20) nn = 20 ;
      for( int i=0; i<nn; ++i ){
	//for( int i=0; i<20; ++i ){
	DCLocalTrack *track=TrackCont[i];
	std::cout << std::setw(3) << i << " #Hits="
		  << std::setw(2) << track->GetNHit() 
		  << " ChiSqr=" << track->GetChiSquare()
		  << std::endl;
      }
      std::cout << std::endl;
    }
  }
#endif
  
  partial_sort( TrackCont.begin(), TrackCont.end(), 
		TrackCont.end(), DCLTrackComp() );
  
  
#if 0
  {
    int nn=TrackCont.size();
    std::cout << funcname << ": After Sorting. #Tracks = " 
	      << nn << std::endl;
    if(nn>0){
      if(nn>20) nn = 20 ;
      for( int i=0; i<nn; ++i ){
	DCLocalTrack *track=TrackCont[i];
	std::cout << std::setw(3) << i << " #Hits="
		  << std::setw(2) << track->GetNHit() 
		  << " ChiSqr=" << track->GetChiSquare()
		  << std::endl;
      }
      std::cout << std::endl;
    }
  }
#endif
  
  // Delete Duplicated Tracks
  
  int nTC = TrackCont.size(); 
  
#if 0
  std::cout << "nTC =  " << nTC << std::endl ;
#endif
  
  
  for( int i=0; i<int(TrackCont.size()); ++i ){
    DCLocalTrack *tp=TrackCont[i];
    int nh=tp->GetNHit();
    for( int j=0; j<nh; ++j ) tp->GetHit(j)->setFlags();
    
    for( int i2=TrackCont.size()-1; i2>i; --i2 ){
      DCLocalTrack *tp2=TrackCont[i2];
      int nh2=tp2->GetNHit(), flag=0;
      for( int j=0; j<nh2; ++j )
	if( tp2->GetHit(j)->showFlags()) ++flag;
      if(flag){
	delete tp2;
	TrackCont.erase(TrackCont.begin()+i2);
      }
    }      
  }
  
  {
    int nn=TrackCont.size();

#if 1
    if(nn>1) std::cout << "nn =  " << nn << std::endl ;
#endif

    for(int i=0; i<nn; ++i ){
      DCLocalTrack *tp=TrackCont[i];
      int nh=tp->GetNHit();

#if 1
      if((nn>1) && (tp->GetChiSquare()>0)) std::cout << "i(nn) = " << i 
					     << " Chisquare(after) = " << tp->GetChiSquare() 
					     << " N Hit = " << nh << std::endl ;
#endif

      for( int j=0; j<nh; ++j ){
	int lnum = tp->GetHit(j)->GetLayer();
	double zz = DCGeomMan::GetInstance().GetLocalZ( lnum );
	tp->GetHit(j)->SetCalPosition(tp->GetX(zz), tp->GetY(zz));

	DCLTrackHit *hit = tp->GetHit(j);
	//std::cout << "layer = " << hit->GetLayer()-30 << " Res = " << hit->GetResidual() << std::endl ;
	//std::cout << "hitp = " << hit->GetLocalHitPos() << " calp = " << hit->GetLocalCalPos() << std::endl ;
	//std::cout << "X = " << hit->GetXcal() << " Y = " << hit->GetYcal() << std::endl ;
      }
    }
  }
  
  
  for( int i=0; i<NumOfLayersSdcOut; ++i )
    for_each( CandCont[i].begin(), CandCont[i].end(), DeleteObject() );

  return TrackCont.size();
}
#endif

//////////////////////////////////////////////////////
//SDC3&4 Old routine (Only used for study)

#if ver1
int SdcOutLocalTrackSearch( const DCHitContainer * HC,
 			    std::vector <DCLocalTrack *> &TrackCont)
  
{
  static const std::string funcname = "[LocalTrackSearch]";
  
  std::vector < std::vector <DCPairHitCluster *> > CandCont;
  CandCont.resize(NumOfLayersSdcOut);
  
  
  for( int i=0; i<NumOfLayersSdcOut; ++i ){
    MakeUnPairPlaneHitCluster( HC[i+1], CandCont[i] );
  }
  
  std::vector <int> nCombi(NumOfLayersSdcOut);
  for( int i=0; i<NumOfLayersSdcOut; ++i ){ 
    nCombi[i]=(CandCont[i]).size();
    
    // If #Cluster>MaxNumerOfCluster,  error return
    
    if(nCombi[i]>MaxNumberOfClusters){
      for( int i=0; i<NumOfLayersSdcOut; ++i )
 	for_each( CandCont[i].begin(), CandCont[i].end(), DeleteObject() );
      return 0;
    } 
  }
  
#if 0
  std::cout << funcname << ": #Hits of each group" << std::endl;
  for( int i=0; i<NumOfLayersSdcOut; ++i ) std::cout << std::setw(4) << nCombi[i];
  std::cout << std::endl;
  for( int i=0; i<NumOfLayersSdcOut; ++i ){
    int n=CandCont[i].size();
    std::cout << "[" << std::setw(3) << i << "]: "
 	      << std::setw(3) << n << " ";
    for( int j=0; j<n; ++j ){
      std::cout << ((DCLTrackHit *)CandCont[i][j]->GetHit(0))->GetWire() << " ";
     }
    std::cout << std::endl;
  }
#endif
  
  std::vector < std::vector <int> > 
    CombiIndex = makeindex_SdcOut( NumOfLayersSdcOut, 
 				   MinNumOfHitsSdcOut, 
 				   NumOfLayersSdcOut, 
 				   &(nCombi[0]) );
  int nnCombi=CombiIndex.size();
  
#if 0
  std::cout << " ===> " << nnCombi << " combinations will be checked.." 
 	    << std::endl;
 #endif
  
  for( int i=0; i<nnCombi; ++i ){
    DCLocalTrack *track = MakeTrack( CandCont, &((CombiIndex[i])[0]) );
     if( !track ) continue;
     if( track->GetNHit()>=MinNumOfHitsSdcOut && track->DoFit() &&
     	 track->GetChiSquare()<MaxChisquare ){
     //     if( track->GetNHit()>=12 && track->DoFit() &&
     //	 track->GetChiSquare()<MaxChisquare ){
       TrackCont.push_back(track);
       double chisqr = track->GetChiSquare();
     }
     else{
       //      std::cout << "No tracks available" << std::endl;
       delete track;
     }
   }
  
  // Clear Flags
   int nbefore=TrackCont.size();
   for( int i=0; i<nbefore; ++i ){
     DCLocalTrack *tp=TrackCont[i];
     int nh=tp->GetNHit();
     for( int j=0; j<nh; ++j ) tp->GetHit(j)->clearFlags();
   }

   partial_sort( TrackCont.begin(), TrackCont.end(), 
		 TrackCont.end(), DCLTrackComp1() );


#if 0
   {
     int nn=TrackCont.size();
     if(nn>0){
       std::cout << funcname << ": Before Sorting. #Tracks = " 
		 << nn << std::endl;
       //for( int i=0; i<nn; ++i ){

       for( int i=0; i<50; ++i ){
	 DCLocalTrack *track=TrackCont[i];
	 std::cout << std::setw(3) << i << " #Hits="
		   << std::setw(2) << track->GetNHit() 
		   << " ChiSqr=" << track->GetChiSquare()
		   << std::endl;
       }
       std::cout << std::endl;
     }
   }
 #endif

   //   partial_sort( TrackCont.begin(), TrackCont.end(), 
   //		 TrackCont.end(), DCLTrackComp1() );

   partial_sort( TrackCont.begin(), TrackCont.end(), 
 		TrackCont.end(), DCLTrackComp() );
   
#if 1
   {
     int nn=TrackCont.size();
     if(nn>0){
       std::cout << funcname << ": After Sorting. #Tracks = " 
		 << nn << std::endl;
       int N9 =0;
       for( int i=0; i<nn; ++i ){
       //for( int i=0; i<10; ++i ){
	 DCLocalTrack *track=TrackCont[i];
	
	 std::cout << std::setw(3) << i << " #Hits="
		   << std::setw(2) << track->GetNHit() 
		   << " ChiSqr=" << track->GetChiSquare()
		   << std::endl;
	 N9++;
	 if(N9>100) break;
	 
	 /*
	   for( int j=0; j<(track->GetNHit()); ++j){
	     DCLTrackHit *hit = track->GetHit(j);
	     std::cout << "layer = " << hit->GetLayer()-30 << " Res = " << hit->GetResidual() << std::endl ;
	     std::cout << "hitp = " << hit->GetLocalHitPos() << " calp = " << hit->GetLocalCalPos() << std::endl ;
	     std::cout << "X = " << hit->GetXcal() << " Y = " << hit->GetYcal() << std::endl ;
	     
	   }
	   std::cout << "***************" << std::endl;
	   std::cout << std::endl;
	 */
       }
       std::cout << std::endl;
     }
   }
#endif
   
   // Delete Duplicated Tracks
   
   for( int i=0; i<int(TrackCont.size()); ++i ){
     DCLocalTrack *tp=TrackCont[i];
     int nh=tp->GetNHit();
     for( int j=0; j<nh; ++j ) tp->GetHit(j)->setFlags();
     for( int i2=TrackCont.size()-1; i2>i; --i2 ){
       DCLocalTrack *tp2=TrackCont[i2];
       int nh2=tp2->GetNHit(), flag=0;
       for( int j=0; j<nh2; ++j )
	 if( tp2->GetHit(j)->showFlags() ) ++flag;
       if(flag){
	 delete tp2;
	 TrackCont.erase(TrackCont.begin()+i2);
       }
     }      
   }

   {
     int nn=TrackCont.size();
     for(int i=0; i<nn; ++i ){
       DCLocalTrack *tp=TrackCont[i];
       int nh=tp->GetNHit();
       for( int j=0; j<nh; ++j ){
	 int lnum = tp->GetHit(j)->GetLayer();
	 double zz = DCGeomMan::GetInstance().GetLocalZ( lnum );
	 tp->GetHit(j)->SetCalPosition(tp->GetX(zz), tp->GetY(zz));
       }
     }
   }
   
#if 1
   {
     int nn=TrackCont.size();
     std::cout << funcname << ": After Deleting. #Tracks = " 
	       << nn << std::endl;
     
     for( int i=0; i<nn; ++i ){
       DCLocalTrack *track=TrackCont[i];
       std::cout << std::setw(3) << i << " #Hits="
		 << std::setw(2) << track->GetNHit() 
		 << " ChiSqr=" << track->GetChiSquare()
		 << std::endl;

       /*	 
       for( int j=0; j<(track->GetNHit()); ++j){
	 DCLTrackHit *hit = track->GetHit(j);
	 std::cout << "layer = " << hit->GetLayer()-30 << " Res = " << hit->GetResidual() << std::endl ;
	 std::cout << "hitp = " << hit->GetLocalHitPos() << " calp = " << hit->GetLocalCalPos() << std::endl ;
	 std::cout << "X = " << hit->GetXcal() << " Y = " << hit->GetYcal() << std::endl ;
	 
       }
       std::cout << "***************" << std::endl;
       std::cout << std::endl;
       */
     }
     std::cout << std::endl;
       
   }
#endif
   
   for( int i=0; i<NumOfLayersSdcOut; ++i )
     for_each( CandCont[i].begin(), CandCont[i].end(), DeleteObject() );
  
   return TrackCont.size();
}
#endif

// int LocalTrackSearchBcOutSdcIn( const DCHitContainer * BcHC,  
//                                  const DCPairPlaneInfo * BcPpInfo,
//                                  const DCHitContainer * SdcHC,  
//                                  const DCPairPlaneInfo * SdcPpInfo,
//                                  int BcNpp, int SdcNpp,
//                                  std::vector <DCLocalTrack *> &TrackCont,
//                                  int MinNumOfHits )
// {
//   static const std::string funcname = "[LocalTrackSearchBcSdc]";

//   std::vector < std::vector <DCPairHitCluster *> > CandCont;
//   CandCont.resize(BcNpp+SdcNpp);


//   for( int i=0; i<BcNpp; ++i ){
//     bool ppFlag=BcPpInfo[i].flag;
//     int layer1=BcPpInfo[i].id1, layer2=BcPpInfo[i].id2;
//     if(ppFlag) 
//       MakePairPlaneHitCluster( BcHC[layer1], BcHC[layer2], 
//                                BcPpInfo[i].CellSize, CandCont[i] );
//     else
//       MakeUnPairPlaneHitCluster( BcHC[layer1], CandCont[i] );
//   }

//   for( int i=0; i<SdcNpp; ++i ){
//     bool ppFlag=SdcPpInfo[i].flag;
//     int layer1=SdcPpInfo[i].id1, layer2=SdcPpInfo[i].id2;
//     if(ppFlag) 
//       MakePairPlaneHitCluster( SdcHC[layer1], SdcHC[layer2], 
//                                SdcPpInfo[i].CellSize, CandCont[i+BcNpp] );
//     else
//       MakeUnPairPlaneHitCluster( SdcHC[layer1], CandCont[i+BcNpp] );
//   }

//   std::vector <int> nCombi(BcNpp+SdcNpp);

//   for( int i=0; i<BcNpp+SdcNpp; ++i ){ 
//     nCombi[i]=(CandCont[i]).size();
//     int clNum;
//     clNum=1;
//     if(nCombi[i]!=clNum){
//       for( int i=0; i<BcNpp+SdcNpp; ++i )
//         for_each( CandCont[i].begin(), CandCont[i].end(), DeleteObject() );
//       return 0;
//     } 
//   }

// #if 0
//   std::cout << funcname << ": #Hits of each group" << std::endl;
//   for( int i=0; i<BcNpp+SdcNpp; ++i ) std::cout << std::setw(4) << nCombi[i];
//   std::cout << std::endl;
//   for( int i=0; i<BcNpp+SdcNpp; ++i ){
//     int n=CandCont[i].size();
//     std::cout << "[" << std::setw(3) << i << "]: "
//               << std::setw(3) << n << " ";
//     for( int j=0; j<n; ++j ){
//       std::cout << CandCont[i][j] << " ";
//     }
//     std::cout << std::endl;
//   }
// #endif

//   std::vector < std::vector <int> > CombiIndex;
//   for (int i=0;i<1;i++) {
//     std::vector <int> index;
//     index.reserve(BcNpp+SdcNpp);
//     for (int j=0;j<BcNpp+SdcNpp;j++) {
//       index.push_back(0);
//     }
//     CombiIndex.push_back(index);
//   }
//   int nnCombi=CombiIndex.size();

// #if 0
//   std::cout << " ===> " << nnCombi << " combinations will be checked.." 
//             << std::endl;
//   for( int i=0; i<nnCombi; ++i ){
//     for (int j=0;j<BcNpp+SdcNpp;j++) {
//       std::cout << CombiIndex[i][j] << " ";
//     }
//     std::cout << std::endl;
//   }
// #endif

//   for( int i=0; i<nnCombi; ++i ){
//     DCLocalTrack *track = MakeTrack( CandCont, &((CombiIndex[i])[0]) );
//     if( !track ) continue;
//     if( track->GetNHit()>=MinNumOfHits && track->DoFitBcSdc() &&
//         track->GetChiSquare()<MaxChisquare ) 
//       TrackCont.push_back(track);
//     else
//       delete track;
//   }

//   // Clear Flags
//   int nbefore=TrackCont.size();
//   for( int i=0; i<nbefore; ++i ){
//     DCLocalTrack *tp=TrackCont[i];
//     int nh=tp->GetNHit();
//     for( int j=0; j<nh; ++j ) tp->GetHit(j)->clearFlags();
//   }


// #if 0
//  {
//    int nn=TrackCont.size();
//    std::cout << funcname << ": Before Sorting. #Tracks = " 
// 	     << nn << std::endl;

//    for( int i=0; i<nn; ++i ){
//      DCLocalTrack *track=TrackCont[i];
//      std::cout << std::setw(3) << i << " #Hits="
// 	       << std::setw(2) << track->GetNHit() 
// 	       << " ChiSqr=" << track->GetChiSquare()
// 	       << std::endl;
//    }
//    std::cout << std::endl;

//  }
// #endif

//   partial_sort( TrackCont.begin(), TrackCont.end(), 
//                 TrackCont.end(), DCLTrackComp() );

// #if 0
//  {
//    int nn=TrackCont.size();
//    std::cout << funcname << ": After Sorting. #Tracks = " 
// 	     << nn << std::endl;
//    /*
//     for( int i=0; i<nn; ++i ){
//       DCLocalTrack *track=TrackCont[i];
//       std::cout << std::setw(3) << i << " #Hits="
//                 << std::setw(2) << track->GetNHit() 
//                 << " ChiSqr=" << track->GetChiSquare()
//                 << std::endl;
//     }
//     std::cout << std::endl;
//    */
//  }
// #endif

//  // Delete Duplicated Tracks

//  for( int i=0; i<int(TrackCont.size()); ++i ){
//    DCLocalTrack *tp=TrackCont[i];
//    int nh=tp->GetNHit();
//    for( int j=0; j<nh; ++j ) tp->GetHit(j)->setFlags();

//    for( int i2=TrackCont.size()-1; i2>i; --i2 ){
//      DCLocalTrack *tp2=TrackCont[i2];
//      int nh2=tp2->GetNHit(), flag=0;
//      for( int j=0; j<nh2; ++j )
//        if( tp2->GetHit(j)->showFlags() ) ++flag;
//      if(flag){
//        delete tp2;
//        TrackCont.erase(TrackCont.begin()+i2);
//      }
//    }      
//  }

//  {
//    double zK18tgt=DCGeomMan::GetInstance().GetLocalZ( 130 );
//    int nn=TrackCont.size();
//    for(int i=0; i<nn; ++i ){
//      DCLocalTrack *tp=TrackCont[i];
//      int nh=tp->GetNHit();
//      for( int j=0; j<nh; ++j ){
//        int lnum = tp->GetHit(j)->GetLayer();
//        double zz = DCGeomMan::GetInstance().GetLocalZ( lnum );
//        if (lnum>=113&&lnum<=124)
// 	 zz -= zK18tgt;
//        tp->GetHit(j)->SetCalPosition(tp->GetX(zz), tp->GetY(zz));
//      }
//    }
//  }

// #if 0
//  {
//    int nn=TrackCont.size();
//    std::cout << funcname << ": After Deleting. #Tracks = " 
// 	     << nn << std::endl;

//    for( int i=0; i<nn; ++i ){
//      DCLocalTrack *track=TrackCont[i];
//      std::cout << std::setw(3) << i << " #Hits="
// 	       << std::setw(2) << track->GetNHit() 
// 	       << " ChiSqr=" << track->GetChiSquare()
// 	       << std::endl;
//    }
//    std::cout << std::endl;

//  }
// #endif


//  for( int i=0; i<SdcNpp+BcNpp; ++i )
//    for_each( CandCont[i].begin(), CandCont[i].end(), DeleteObject() );
  
//  return TrackCont.size();
// }
