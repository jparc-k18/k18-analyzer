/*
  DCAnalyzer.cc

  2012/1/24
*/

#include <iostream>
#include <iomanip>
#include <string>
#include <algorithm>
#include <stdexcept>
#include <sstream>
#include <set>

#include "DCAnalyzer.hh"
#include "DCHit.hh"
#include "DCLocalTrack.hh"
#include "K18Track.hh"
#include "K18Track_BFT.hh"
#include "SksTrack.hh"
#include "RawData.hh"
#include "DCRawHit.hh"
#include "TemplateLib.hh"
#include "DCTrackSearch.hh"
#include "K18Parameters.hh"
#include "ConfMan.hh"
#include "MWPCCluster.hh"
#include "FiberCluster.hh"
#include "HodoAnalyzer.hh"
#include "DCDriftParamMan.hh"
#include "DCGeomMan.hh"
#include "Hodo1Hit.hh"

// #include "SimuData.hh"
// #include "SimuDetectorId.hh"

#define DefStatic
#include "DCParameters.hh"
#undef DefStatic

const double MaxChiSqrSksTrack = 10000.;
const double MaxTimeDifMWPC = 100.;

const double kMWPCClusteringWireExtension = 1.0;  // [mm]
const double kMWPCClusteringTimeExtension = 10.0; // [nsec]

const double Deg2Rad = acos(-1.)/180. ;

/////////////////////////Tracking routine selection

///////////////////////////////
////Loacal Tracking Rootine////
///////////////////////////////

/******** Bc Out Tracking  *****/
//(A) -> XUV Tracking (slow but accerate)
#define ver2_BcOut  0

//(B) -> Pair plane Tracking (fast but bad for large angle track)
#define ver1_BcOut  1 
/*******************************/


/******** Sdc In Tracking  *****/
//(A) -> XUV Tracking (slow but accerate)
#define ver2_SdcIn  0 

//(B) -> Pair plane Tracking (fast but bad for large angle track)
#define ver1_SdcIn  1 
/*******************************/


////////////////////////////
////SKS Tracking Rootine////
////////////////////////////
/*******************************/
//(A) -> New SksTracking Rootine (Roughly calculate initial momentum track by track)
#define ver2_Sks  0
//double MinP = 0.7 ;//400A Minimum Momentum
double MinP = 0.6 ;//300A Minimum Momentum

//Initial Momentum calculated from Kinematics
#define Pini_kine  0 

// Magnetic Field Strength of SKS
double B_Sks = 2.1 ;//(T)

//(B) -> Old SksTracking Rootine (Fixed initial momentum)
#define ver1_Sks  1 
double PIni = 1.20 ;//Initial Momentum (GeV/c)  
/*******************************/

#ifdef MemoryLeak
debug::Counter DCAnalyzer::sm_counter("DCAnalyzer");
#endif
  
//______________________________________________________________________________
// forward declaration
int
clusterizeMWPCHit(const DCHitContainer&,
 		  MWPCClusterContainer& clusters);

//______________________________________________________________________________

DCAnalyzer::DCAnalyzer()
{
  //std::cout << "aaa" << std::endl;
#ifdef MemoryLeak
  ++sm_counter;
#endif
}

DCAnalyzer::~DCAnalyzer()
{
  clearSksTracks();
  clearK18Tracks();
  clearK18Tracks_BFT();
  clearTracksSdcOut();
  clearTracksSdcIn();
  clearTracksBcOut();
  clearTracksBcIn();
  clearTracksBcOutSdcIn(); 
  clearDCHits();
  clearVtxHits();

#ifdef MemoryLeak
  --sm_counter;
#endif
}

/*
bool DCAnalyzer::DecodeFiberHits(FiberCluster * FiberCl, int layer )
{
  const std::string funcname = "[DCAnalyzer::DecodeFiberHits]";

  //Decode SFT for SdcIn
  double meanpos = FiberCl->MeanPosition();
  double meantime = FiberCl->CMeanTime();
    
  DCHit *dchit = new DCHit( layer+PlMinSdcIn-1 );
  dchit->SetTdcVal( ( int ) meantime );
  bool stReturn = true;
  if( dchit->CalcFiberObservables() ){
    dchit->SetWirePosition( meanpos );
    SdcInHC[layer].push_back(dchit);
  }
  else{
    delete dchit;
    stReturn = false;
  }

  return stReturn;
}
*/
bool DCAnalyzer::DecodeBcInHits( RawData *rawData )
{
  const std::string funcname = "[DCAnalyzer::DecodeBcInHits]";

  clearBcInHits();

  // BcIn
  for( int layer=1; layer<=NumOfLayersBcIn; ++layer ){
    const DCRHitContainer &cont=rawData->GetBcInRawHC(layer);
    int nh=cont.size();
    for( int i=0; i<nh; ++i ){
      DCRawHit *rhit=cont[i];

      DCHit *thit=new DCHit(rhit->PlaneId()+PlOffsBc, rhit->WireId());
      int nhtdc= rhit->GetTdcSize();
      for( int j=0; j<nhtdc; ++j ){
	thit->SetTdcVal( rhit->GetTdc(j) );
	thit->SetTdcTrailing( rhit->GetTrailing(j) );
      }
      if(!thit) continue;
    
      if(thit->CalcMWPCObservables())
	TempBcInHC[layer].push_back(thit);
      else
	delete thit;
    }

    //std::cout<<"*************************************"<<std::endl;
    int ncl = clusterizeMWPCHit( TempBcInHC[layer], MWPCClCont[layer]);
    //std::cout<<"numCl="<< ncl << std::endl;
     for( int i=0; i<ncl; ++i ){
      MWPCCluster *p=MWPCClCont[layer][i];
      if (!p) continue;

      const MWPCCluster::Statistics& mean  = p->GetMean();
      const MWPCCluster::Statistics& first = p->GetFirst();
      double mwire    = mean.m_wire;
      double mwirepos = mean.m_wpos;
      double mtime    = mean.m_leading;
      double mtrail   = mean.m_trailing;
      
      DCHit *hit=new DCHit( layer+PlOffsBc, mwire );
      hit->SetClusterSize( p->GetClusterSize() );
      hit->SetMWPCFlag( true );
      hit->SetWire( mwire );
      hit->SetMeanWire( mwire );
      hit->SetMeanWirePosition( mwirepos );
      hit->SetDriftTime( mtime );
      hit->SetTrailingTime( mtrail );
      hit->SetDriftLength( 0. );
      hit->SetTdcVal( 0 );
      hit->SetTdcTrailing( 0 );
      if(!hit) continue;

      if(hit->CalcMWPCObservables())
	BcInHC[layer].push_back(hit);
      else
	delete hit;
    }
    //std::cout << "nh="<< BcInHC[layer].size() <<std::endl;
  } 
  return true;
}


bool DCAnalyzer::DecodeBcOutHits( RawData *rawData )
{  
  const std::string funcname = "[DCAnalyzer::DecodeBcOutHits]";

  clearBcOutHits();

  // BcOut
  for( int layer=1; layer<=NumOfLayersBcOut; ++layer ){
    const DCRHitContainer &cont=rawData->GetBcOutRawHC(layer);
    int nh=cont.size();
    for( int i=0; i<nh; ++i ){
      DCRawHit *rhit=cont[i];
      
      DCHit *hit=new DCHit(rhit->PlaneId()+PlOffsBc, rhit->WireId());
       int nhtdc= rhit->GetTdcSize();
      for( int j=0; j<nhtdc; ++j ){
	hit->SetTdcVal( rhit->GetTdc(j) );
      }
      if(!hit) continue;
      
      if(hit->CalcDCObservables())
	BcOutHC[layer].push_back(hit);
      else
	delete hit;
    }
  }
  return true;
}
/*
bool DCAnalyzer::DecodeSsdHits( RawData *rawData )
{
  clearSsdHits();

  for( int layer=1; layer<NumOfLayersSsd+1;++layer)
    {
      const DCRHitContainer &cont=rawData->GetSsdRawHC(layer);
      int nh = cont.size();
      for(int i=0;i<nh;++i)
	{
	  DCRawHit *rhit = cont[i];
	  DCHit *hit = new DCHit(rhit->PlaneId()+PlOffsSsd,
				 rhit->WireId());
	  int nhadc = rhit->GetTdcSize();
	  for(int j=0;j<nhadc-1;++j)
	    hit->SetAdcVal( rhit->GetTdc(j) );
	  hit->SetTdcVal(rhit->GetTdc(nhadc-1));
	  if(!hit) continue;
	  //      SsdHC[layer].push_back(hit);
	  if(hit->CalcSsdObservables())
	    SsdHC[layer].push_back(hit);
	  else
	    delete hit;
	}
    }

  clusterizeSSDHit();

  return true;
}
*/

bool DCAnalyzer::DecodeSdcInHits( RawData *rawData )
{  
  const std::string funcname = "[DCAnalyzer::DecodeSdcInHits]";
  
  clearSdcInHits();
    for( int layer=1; layer<=NumOfLayersSdcIn; ++layer ){
      const DCRHitContainer &cont=rawData->GetSdcInRawHC(layer);
      int nh=cont.size();
      for( int i=0; i<nh; ++i ){
	DCRawHit *rhit=cont[i];
	
	DCHit *hit=new DCHit( rhit->PlaneId(), rhit->WireId() );
	int nhtdc= rhit->GetTdcSize();
	for( int j=0; j<nhtdc; ++j ){ hit->SetTdcVal( rhit->GetTdc(j) ); }
	if(!hit){ continue; }
	
	if(hit->CalcDCObservables()){ SdcInHC[layer].push_back(hit); }
	else{ delete hit; }
      }
    }
  
  return true;
}

bool DCAnalyzer::DecodeSdcOutHits( RawData *rawData )
{
  const std::string funcname = "[DCAnalyzer::DecodeSdcOutHits]";

  clearSdcOutHits();

  // SdcOut
  for( int layer=1; layer<=NumOfLayersSdcOut; ++layer ){
    const DCRHitContainer &cont=rawData->GetSdcOutRawHC(layer);
    int nh=cont.size();
    for( int i=0; i<nh; ++i ){
      DCRawHit *rhit=cont[i];

      DCHit *hit=new DCHit( rhit->PlaneId(), rhit->WireId() );
      int nhtdc= rhit->GetTdcSize();
      //std::cout<<"nhtdc=   "<<nhtdc<<std::endl;
      for( int j=0; j<nhtdc; ++j ){
	hit->SetTdcVal( rhit->GetTdc(j) );
	      }
      if(!hit) continue;

      if(hit->CalcDCObservables())
	SdcOutHC[layer].push_back(hit);
      else
	delete hit;
    }
  }
  return true;
}

bool DCAnalyzer::DecodeRawHits( RawData *rawData )
{
  const std::string funcname = "[DCAnalyzer::DecodeRawHits]";

  clearDCHits();

  DecodeBcInHits( rawData );
  DecodeBcOutHits( rawData );
  DecodeSdcInHits( rawData );
  DecodeSdcOutHits( rawData );
  //DecodeSsdHits( rawData );
  return true;
}

/*
bool DCAnalyzer::DecodeRawHits( RawData *rawData )
{
  const std::string funcname = "[DCAnalyzer::DecodeRawHits]";

  clearDCHits();

  // BcIn
  for( int layer=1; layer<=NumOfLayersBcIn; ++layer ){
    const DCRHitContainer &cont=rawData->GetBcInRawHC(layer);
    int nh=cont.size();
    for( int i=0; i<nh; ++i ){
      DCRawHit *rhit=cont[i];

      DCHit *thit=new DCHit(rhit->PlaneId()+PlOffsBc, rhit->WireId());
      int nhtdc= rhit->GetTdcSize();
      for( int j=0; j<nhtdc; ++j ){
	thit->SetTdcVal( rhit->GetTdc(j) );
	thit->SetTdcTrailing( rhit->GetTrailing(j) );
      }
      if(!thit) continue;
    
      if(thit->CalcMWPCObservables())
	TempBcInHC[layer].push_back(thit);
      else
	delete thit;
    }

    //std::cout<<"*************************************"<<std::endl;
    int ncl = clusterizeMWPCHit( TempBcInHC[layer], MWPCClCont[layer]);
    //std::cout<<"numCl="<< ncl << std::endl;
     for( int i=0; i<ncl; ++i ){
      MWPCCluster *p=MWPCClCont[layer][i];
      if (!p) continue;

      const MWPCCluster::Statistics& mean  = p->GetMean();
      const MWPCCluster::Statistics& first = p->GetFirst();
      double mwire    = mean.m_wire;
      double mwirepos = mean.m_wpos;
      double mtime    = mean.m_leading;
      double mtrail   = mean.m_trailing;
      
      DCHit *hit=new DCHit( layer+PlOffsBc, mwire );
      hit->SetClusterSize( p->GetClusterSize() );
      hit->SetMWPCFlag( true );
      hit->SetWire( mwire );
      hit->SetMeanWire( mwire );
      hit->SetMeanWirePosition( mwirepos );
      hit->SetDriftTime( mtime );
      hit->SetTrailingTime( mtrail );
      hit->SetDriftLength( 0. );
      hit->SetTdcVal( 0 );
      hit->SetTdcTrailing( 0 );
      if(!hit) continue;

      if(hit->CalcMWPCObservables())
	BcInHC[layer].push_back(hit);
      else
	delete hit;
    }
    //std::cout << "nh="<< BcInHC[layer].size() <<std::endl;
  } 
  
  // BcOut
  for( int layer=1; layer<=NumOfLayersBcOut; ++layer ){
    const DCRHitContainer &cont=rawData->GetBcOutRawHC(layer);
    int nh=cont.size();
    for( int i=0; i<nh; ++i ){
      DCRawHit *rhit=cont[i];
      
      DCHit *hit=new DCHit(rhit->PlaneId()+PlOffsBc, rhit->WireId());
       int nhtdc= rhit->GetTdcSize();
      for( int j=0; j<nhtdc; ++j ){
	hit->SetTdcVal( rhit->GetTdc(j) );
      }
      if(!hit) continue;
      
      if(hit->CalcDCObservables())
	BcOutHC[layer].push_back(hit);
      else
	delete hit;
    }
  }
      
  // SdcIn

  for( int layer=2; layer<=NumOfLayersSdcIn; ++layer ){

    if(layer>=2&&layer<2+NumOfLayersSFT){
      FiberCluster *Fcont = GetFiberCluster( layer );
      DecodeFiberHits( Fcont );
    }
    
    const DCRHitContainer &cont=rawData->GetSdcInRawHC(layer);
    int nh=cont.size();
    for( int i=0; i<nh; ++i ){
      DCRawHit *rhit=cont[i];

      DCHit *hit=new DCHit( rhit->PlaneId(), rhit->WireId() );
      int nhtdc= rhit->GetTdcSize();
      for( int j=0; j<nhtdc; ++j ){
	hit->SetTdcVal( rhit->GetTdc(j) );
      }
      if(!hit) continue; 

      if(hit->CalcDCObservables())
	SdcInHC[layer].push_back(hit);
      else
	delete hit;
    }
  }
      
  // SdcOut
  for( int layer=1; layer<=NumOfLayersSdcOut; ++layer ){
    const DCRHitContainer &cont=rawData->GetSdcOutRawHC(layer);
    int nh=cont.size();
    for( int i=0; i<nh; ++i ){
      DCRawHit *rhit=cont[i];

      DCHit *hit=new DCHit( rhit->PlaneId(), rhit->WireId() );
      int nhtdc= rhit->GetTdcSize();
      for( int j=0; j<nhtdc; ++j ){
	hit->SetTdcVal( rhit->GetTdc(j) );
      }
      if(!hit) continue;

      if(hit->CalcDCObservables())
	SdcOutHC[layer].push_back(hit);
      else
	delete hit;
    }
  }
      
  return true;
}

bool DCAnalyzer::TrackSearchBcIn( void )
{
  int ntrack =
    MWPCLocalTrackSearch(  &(BcInHC[1]), TrackBcInCol );

  return true;
}
*/

bool DCAnalyzer::TrackSearchBcIn( const std::vector<std::vector<DCHitContainer> >& hc)
{
  //   std::cout << __FILE__ << ":" << __func__ << std::endl;
  int ntracks = 0;
  MWPCLocalTrackSearch(hc, TrackBcInCol);

//   std::vector<DCHitContainer> unusedHits;
//   int toSearch=findUnusedHits(NumOfLayersBcIn, &(BcInHC[1]), unusedHits);
//   if (toSearch>0) {
// //     //       std::cout << "#D DCAnalyzer::" << __func__ 
// //     // 		<< " : " << ntracks << " tracks found. "
// //     // 		<< " trying to find tracks w/o filter" << std::endl;
//     std::vector<DCLocalTrack*> additionalTracks;
//     int toAdd = MWPCLocalTrackSearch(&(unusedHits[0]), additionalTracks);
//     if (toAdd>0){
// //     // 	std::cout << "#D " << __func__
// //     // 		  << " : " << toAdd << " tracks found" << std::endl;
//       TrackBcInCol.insert(TrackBcInCol.end(), 
// 			  additionalTracks.begin(), additionalTracks.end());
//       ntracks += toAdd;
//     }
//   }
// //   //   std::cout << "#D " << __FILE__ << ":" << __LINE__
// //   // 	    << " " << __func__ << " total tracks = "
// //   // 	    << TrackBcInCol.size() << std::endl;

  return true;
}


//Pair Plane Tracking Routine for BCOut
#if ver1_BcOut
bool DCAnalyzer::TrackSearchBcOut( void )
{
  int ntrack =
    LocalTrackSearch( BcOutHC, PPInfoBcOut, NPPInfoBcOut,
		      TrackBcOutCol, MinNumOfHitsBcOut );

  return true;
}

// Use with BH2Filter
bool DCAnalyzer::TrackSearchBcOut( const std::vector<std::vector<DCHitContainer> >& hc)
{
  int ntrack =
    LocalTrackSearch( hc, PPInfoBcOut, NPPInfoBcOut,
		      TrackBcOutCol, MinNumOfHitsBcOut );

//   std::vector<DCHitContainer> unusedHits;
//   int toSearch=findUnusedHits(NumOfLayersBcOut+1, &(BcOutHC[0]), unusedHits);
//   if (toSearch>0) {
//     //     std::cout << "#D DCAnalyzer::" << __func__ 
//     // 	      << " : " << ntrack << " tracks found. "
//     // 	      << " trying to find tracks w/o filter" << std::endl;
//     std::vector<DCLocalTrack*> additionalTracks;
//     int toAdd = LocalTrackSearch(&(unusedHits[0]), PPInfoBcOut, NPPInfoBcOut,
//   				 additionalTracks, MinNumOfHitsBcOut);
//     if (toAdd>0){
//     //       std::cout << "#D " << __func__
//     // 		<< " : " << toAdd << " tracks found" << std::endl;
//     //    }
//       TrackBcOutCol.insert(TrackBcOutCol.end(), 
//       			   additionalTracks.begin(), additionalTracks.end());
//       ntrack += toAdd;
//     }
//   }

  return true;
}
#endif

//XUV Tracking Routine for BCOut
#if ver2_BcOut         
bool DCAnalyzer::TrackSearchBcOut( void )
{
  int ntrack =
    LocalTrackSearchVUX( BcOutHC, PPInfoBcOut, NPPInfoBcOut,
			 TrackBcOutCol, MinNumOfHitsBcOut );
  
  return true;
}

// Use with BH2Filter
bool DCAnalyzer::TrackSearchBcOut( const std::vector<std::vector<DCHitContainer> > &hc)
{
  int ntrack =
    LocalTrackSearch( hc, PPInfoBcOut, NPPInfoBcOut,
		      TrackBcOutCol, MinNumOfHitsBcOut );

  return true;
}
#endif

//Pair Plane Tracking Routine for SDCIn
#if ver1_SdcIn
bool DCAnalyzer::TrackSearchSdcIn( void )
{
  int ntrack =
    LocalTrackSearch( SdcInHC, PPInfoSdcIn, NPPInfoSdcIn,
			   TrackSdcInCol, MinNumOfHitsSdcIn );

  return true;
}
/*
bool DCAnalyzer::TrackSearchSdcInFiber( void )
{
  int ntrack =
    LocalTrackSearchSdcInFiber( SdcInHC, PPInfoSdcIn, NPPInfoSdcIn,
				TrackSdcInCol, MinNumOfHitsSdcIn );

  return true;
}
*/
#endif

//XUV Tracking Routine for SDCIn
#if ver2_SdcIn         
bool DCAnalyzer::TrackSearchSdcIn( void )
{
  int ntrack =
    LocalTrackSearchVUX( SdcInHC, PPInfoSdcIn, NPPInfoSdcIn,
			 TrackSdcInCol, MinNumOfHitsSdcIn );
  
  //std::cout << "*************************" << std::endl;
  return true;
}
#endif

bool DCAnalyzer::TrackSearchSdcOut( void )
{

  //Tracking Routine for SDCOut (**Calculate all hit combinations -> Too slow !!**)
  //int ntrack =
  //SdcOutLocalTrackSearch( &(SdcOutHC[1]), TrackSdcOutCol, TrackSdcOutCol1, TrackSdcOutCol2);

  //XUV Tracking Routine for SDCOut  
  int ntrack =
    SdcOutLocalTrackSearch( SdcOutHC, TrackSdcOutCol );

  return true;
}

bool DCAnalyzer::TrackSearchBcOutSdcIn( void )
{

  // 16 of the 20 live layers (BC3 layers 1-2 are dead in E05)
  int ntrack =
    LocalTrackSearchBcOutSdcIn( BcOutHC, PPInfoBcOut,
				SdcInHC, PPInfoSdcIn,
				NPPInfoBcOut,NPPInfoSdcIn,
				TrackBcOutSdcInCol, 16 );

  return true;
}

bool DCAnalyzer::TrackSearchK18( void )
{
  static const std::string funcname = "[DCAnalyzer::TrackSearchK18]";

  clearK18Tracks();
  int nIn=TrackBcInCol.size(), nOut=TrackBcOutCol.size();
#if 0
  std::cout<<"**************************************"<<std::endl;
  std::cout << funcname << ": #TracksIn=" << std::setw(3) << nIn
	    << " #TracksOut=" << std::setw(3) << nOut << std::endl;
#endif

  if( nIn==0 || nOut==0 ) return true;

  double pK18=ConfMan::GetConfManager()->K18Momentum();

  for( int iIn=0; iIn<nIn; ++iIn ){
    DCLocalTrack *trIn=TrackBcInCol[iIn];
#if 0
    std::cout << "TrackIn  :" << std::setw(2) << iIn
	      << " X0=" << trIn->GetX0() << " Y0=" << trIn->GetY0()
	      << " U0=" << trIn->GetU0() << " V0=" << trIn->GetV0()
	      << std::endl;
#endif
    if( !trIn->GoodForTracking() ||
	trIn->GetX0()<MinK18InX || trIn->GetX0()>MaxK18InX ||
	trIn->GetY0()<MinK18InY || trIn->GetY0()>MaxK18InY ||
	trIn->GetU0()<MinK18InU || trIn->GetU0()>MaxK18InU ||
	trIn->GetV0()<MinK18InV || trIn->GetV0()>MaxK18InV ) continue;
    for( int iOut=0; iOut<nOut; ++iOut ){
      DCLocalTrack *trOut=TrackBcOutCol[iOut];
#if 0
      std::cout << "TrackOut :" << std::setw(2) << iOut
		<< " X0=" << trOut->GetX0() << " Y0=" << trOut->GetY0()
		<< " U0=" << trOut->GetU0() << " V0=" << trOut->GetV0()
		<< std::endl;
#endif
      if( !trOut->GoodForTracking() ||
	  trOut->GetX0()<MinK18OutX || trOut->GetX0()>MaxK18OutX ||
	  trOut->GetY0()<MinK18OutY || trOut->GetY0()>MaxK18OutY ||
	  trOut->GetU0()<MinK18OutU || trOut->GetU0()>MaxK18OutU ||
	  trOut->GetV0()<MinK18OutV || trOut->GetV0()>MaxK18OutV ) continue;

#if 0
      std::cout << funcname << ": In -> " << trIn->GetChiSquare() 
		<< " (" << std::setw(2) << trIn->GetNHit() << ") "
		<< "Out -> " << trOut->GetChiSquare()
		<< " (" << std::setw(2) << trOut->GetNHit() << ") "
		<< std::endl;
#endif

      K18Track *track=new K18Track( trIn, trOut, pK18 );
      if( track && track->doFit() )
	K18TrackCol.push_back(track);
      else
	delete track;
    }
  }

#if 0
  std::cout<<"********************"<<std::endl;
 {
   int nn=K18TrackCol.size();
   std::cout << funcname << ": Before sorting. #Track=" 
	     << nn << std::endl;
   for( int i=0; i<nn; ++i ){
     K18Track *tp=K18TrackCol[i];

     std::cout << std::setw(3) << i 
	       << " ChiSqr=" << tp->chisquare()
	       << " Delta=" << tp->Delta()
	       << " P=" << tp->P() << "\n";
//      std::cout<<"********************"<<std::endl;
//      std::cout << "In :"
// 	       << " X " << tp->Xin() << "(" << tp->TrackIn()->GetX0() << ")"
// 	       << " Y " << tp->Yin() << "(" << tp->TrackIn()->GetY0() << ")"
// 	       << " U " << tp->Uin() << "(" << tp->TrackIn()->GetU0() << ")"
// 	       << " V " << tp->Vin() << "(" << tp->TrackIn()->GetV0() << ")"
// 	       << "\n";
//      std::cout << "Out:"
// 	       << " X " << tp->Xout() << "(" << tp->TrackOut()->GetX0() << ")"
// 	       << " Y " << tp->Yout() << "(" << tp->TrackOut()->GetY0() << ")"
// 	       << " U " << tp->Uout() << "(" << tp->TrackOut()->GetU0() << ")"
// 	       << " V " << tp->Vout() << "(" << tp->TrackOut()->GetV0() << ")"
// 	       << std::endl;
   }
 }
#endif

  partial_sort( K18TrackCol.begin(), K18TrackCol.end(),
		K18TrackCol.end(), K18TrackComp() );

  return true;
}

bool DCAnalyzer::TrackSearchK18_BFT( const std::vector<double>& cont_x )
{
  static const std::string funcname = "[DCAnalyzer::TrackSearchK18_BFT]";

  clearK18Tracks_BFT();
  int nOut=TrackBcOutCol.size();

  if(0 == cont_x.size()){return false;}

  double pK18=ConfMan::GetConfManager()->K18Momentum();

  for( int iIn=0; iIn<cont_x.size(); ++iIn ){
    double LocalX = cont_x.at(iIn);
    
    for( int iOut=0; iOut<nOut; ++iOut ){
      DCLocalTrack *trOut=TrackBcOutCol[iOut];
#if 0
      std::cout << "TrackOut :" << std::setw(2) << iOut
		<< " X0=" << trOut->GetX0() << " Y0=" << trOut->GetY0()
		<< " U0=" << trOut->GetU0() << " V0=" << trOut->GetV0()
		<< std::endl;
#endif
      if( !trOut->GoodForTracking() ||
	  trOut->GetX0()<MinK18OutX || trOut->GetX0()>MaxK18OutX ||
	  trOut->GetY0()<MinK18OutY || trOut->GetY0()>MaxK18OutY ||
	  trOut->GetU0()<MinK18OutU || trOut->GetU0()>MaxK18OutU ||
	  trOut->GetV0()<MinK18OutV || trOut->GetV0()>MaxK18OutV ) continue;

#if 0
      std::cout << funcname 
		<< "Out -> " << trOut->GetChiSquare()
		<< " (" << std::setw(2) << trOut->GetNHit() << ") "
		<< std::endl;
#endif

      K18Track_BFT *track=new K18Track_BFT( LocalX, trOut, pK18 );
      if( track ) {
	track->CalcMomentumD2U();
	K18TrackBFTCol.push_back(track);
      } else
	delete track;
    }
  }

#if 0
  std::cout<<"********************"<<std::endl;
  {
   int nn=K18TrackBFTCol.size();
   std::cout << funcname << ": Before sorting. #Track=" 
	     << nn << std::endl;
   for( int i=0; i<nn; ++i ){
     K18Track *tp=K18TrackBFTCol[i];

     std::cout << std::setw(3) << i 
	       << " ChiSqr=" << tp->chisquare()
	       << " Delta=" << tp->Delta()
	       << " P=" << tp->P() << "\n";
//      std::cout<<"********************"<<std::endl;
//      std::cout << "In :"
// 	       << " X " << tp->Xin() << "(" << tp->TrackIn()->GetX0() << ")"
// 	       << " Y " << tp->Yin() << "(" << tp->TrackIn()->GetY0() << ")"
// 	       << " U " << tp->Uin() << "(" << tp->TrackIn()->GetU0() << ")"
// 	       << " V " << tp->Vin() << "(" << tp->TrackIn()->GetV0() << ")"
// 	       << "\n";
//      std::cout << "Out:"
// 	       << " X " << tp->Xout() << "(" << tp->TrackOut()->GetX0() << ")"
// 	       << " Y " << tp->Yout() << "(" << tp->TrackOut()->GetY0() << ")"
// 	       << " U " << tp->Uout() << "(" << tp->TrackOut()->GetU0() << ")"
// 	       << " V " << tp->Vout() << "(" << tp->TrackOut()->GetV0() << ")"
// 	       << std::endl;
   }
 }
#endif

  return true;
}

bool DCAnalyzer::TrackSearchSks( void )
{
  static const std::string funcname = "[DCAnalyzer::TrackSearchSks]";

  clearSksTracks();

  int nIn=TrackSdcInCol.size(), nOut=TrackSdcOutCol.size();
#if 0
  std::cout<<"*********************************************"<<std::endl;
  std::cout << funcname << ": #TracksIn=" << std::setw(3) << nIn
	    << " #TracksOut=" << std::setw(3) << nOut << std::endl;
#endif
  if( nIn==0 || nOut==0 ) return true;

  for( int iIn=0; iIn<nIn; ++iIn ){
    DCLocalTrack *trIn=TrackSdcInCol[iIn];
    if( !trIn->GoodForTracking() ) continue;
    for( int iOut=0; iOut<nOut; ++iOut ){
      DCLocalTrack *trOut=TrackSdcOutCol[iOut];
      if( !trOut->GoodForTracking() ) continue;


#if Pini_kine  
      /*************************************************************************/
      /***** Initial Momentum calculated from Kinematics *******/
      
      double X0In = trIn->GetX0(), Y0In = trIn->GetY0();
      double U0In = trIn->GetU0(), V0In = trIn->GetV0();
      double X0Out = trOut->GetX0(), Y0Out = trOut->GetY0();
      double U0Out = trOut->GetU0(), V0Out = trOut->GetV0();
      
      int SksVp1_Id = DCGeomMan::GetInstance().GetSksVp1Id();
      int SksVp5_Id = DCGeomMan::GetInstance().GetSksVp5Id();
      double zSksVp1 = DCGeomMan::GetInstance().GetLocalZ( SksVp1_Id );
      double zSksVp5 = DCGeomMan::GetInstance().GetLocalZ( SksVp5_Id );
  
      //Get Global Position at the SksVp1 & SksVp5
      ThreeVector gpVp1 = DCGeomMan::GetInstance().GetGlobalPosition(SksVp1_Id);
      ThreeVector gpVp5 = DCGeomMan::GetInstance().GetGlobalPosition(SksVp5_Id);
      
      double X0_Vp1 = gpVp1.x(); 
      double Y0_Vp1 = gpVp1.y(); 
      double X0_Vp5 = gpVp5.x(); 
      double Y0_Vp5 = gpVp5.y(); 
      
      //SdcIn Local Tracking -> X(Z = SksVp1) = X0 + U0*Z 
      double XIn_Vp1 = X0In + U0In*zSksVp1 ;
      double theta_In = atan(U0In);
      
      //triangle 
      double dX_Vp1 = fabs(XIn_Vp1);
      
      double L_Vp1 = dX_Vp1*sin((90-theta_In)*Deg2Rad)/sin((60+theta_In)*Deg2Rad);
      
      double X_Vp1,Y_Vp1;
      
      if(XIn_Vp1>=0){
	X_Vp1 = X0_Vp1 - L_Vp1*cos(70*Deg2Rad);
	Y_Vp1 = Y0_Vp1 - L_Vp1*sin(70*Deg2Rad);
      }
      if(XIn_Vp1<0){
	X_Vp1 = X0_Vp1 + L_Vp1*cos(70*Deg2Rad);
	Y_Vp1 = Y0_Vp1 + L_Vp1*sin(70*Deg2Rad);
      }
      
      //SdcOut Local Tracking -> X(Z = SksVp5) = X0 + U0*Z 
      double XOut_Vp5 = X0Out + U0Out*zSksVp5 ;
      double theta_Out = atan(U0Out);
      
      //triangle 
      double dX_Vp5 = fabs(XOut_Vp5);
      double L_Vp5 = dX_Vp5*sin((90+theta_Out)*Deg2Rad)/sin((60-theta_Out)*Deg2Rad);
      
      double X_Vp5,Y_Vp5;
      
      if(XOut_Vp5>=0){
	X_Vp5 = X0_Vp5 + L_Vp5*cos(70*Deg2Rad);
	Y_Vp5 = Y0_Vp5 - L_Vp5*sin(70*Deg2Rad);
      }
      if(XOut_Vp5<0){
	X_Vp5 = X0_Vp5 - L_Vp5*cos(70*Deg2Rad);
	Y_Vp5 = Y0_Vp5 + L_Vp5*sin(70*Deg2Rad);
      }

      double L_Sks = sqrt(pow( (X_Vp1-X_Vp5) ,2) 
			  + pow( (Y_Vp1-Y_Vp5),2))  ;
      double theta_Sks = 100 - theta_In + theta_Out ;
      double r_Sks = L_Sks/2/sin(theta_Sks*Deg2Rad/2);
      double P_xy = 0.3 * B_Sks * r_Sks *0.001;
      double P_Sks = P_xy*sqrt(1+V0In*V0In+U0In*U0In)/sqrt(1+U0In*U0In) ;

      //Momentum Correction (Momentum correction from the collelation between PSks and U0In) 
      double P_Sks_Cor = P_Sks+0.3204898*U0In + 0.12 ;
      
      //Initial momentum
      Pini = P_Sks_Cor ;
      /*************************************************************************/
#endif

//(A) -> New SksTracking Rootine 
//(Roughly calculate initial momentum track by track)
#if ver2_Sks
      
      SksTrack *pretrack = new SksTrack( trIn, trOut );      
      int    preK = 0, prekk = 2;
      double preiP, preiP_pre = PIni;
      int    preNiteratio n = 0;
      bool   prestatus;
      double chisqr2, prep_check, prechisqr2 = 10000., prechisqr3 = 1.0e+10, preP2 = PIni;
      double dchisqr2 = -1.;

      for(;;){
	if(preK==0) preiP = PIni;
	pretrack->SetInitialMomentum( preiP );
	pretrack->PreTracking( 1 );//2010-7-28
	pretrack->doFit();
	preNiteration++;
	//	preNiteration=pretrack->Niteration();
	
	prep_check=pretrack->PrimaryMomentum().Mag();
	prestatus=pretrack->Status();
	chisqr2=pretrack->chisqr();
	
	if((prechisqr2>chisqr2) && (prestatus)){
	  prechisqr2 = chisqr2 ;
	  preP2 = preiP ;
	}
	
#if 0
	  /*
	    std::cout << "Event # = " << Number << " : prestatus = " << prestatus << std::endl;
	    std::cout << "preK = " << preK << " : # of iteration = " << pretrack->Niteration() 
	    << " : preiP = " << preiP << " : P_g4 = " << gp << std::endl;
	    std::cout << "P = " << prep_check << " : chisqr2 = " << chisqr2  << std::endl;
	    std::cout << "preP2 = " << preP2 << " : prechisqr2 = " << prechisqr2  << std::endl;
	  */
	  // std::cout << pretrack->Niteration() << " " << chisqr2 << " " << preiP  << " " << PIni << std::endl;
#endif
	
	++preK;
	if((chisqr2>1.0e+05) || !(chisqr2>1.0e-100))
	  preiP=preiP-0.02;
	else if(chisqr2>1.0e+04)
	  preiP=preiP-0.01;
	else if(chisqr2>1.0e+03)
	  preiP=preiP-0.002;//0.001
	else if(chisqr2>1.0e+02)
	  preiP=preiP-0.001;//0.0005
	else 
	  preiP=preiP-0.0005;//0.0001

	prechisqr3 = chisqr2;

	//delete pretrack;
	if( ((log10(chisqr2)-log10(prechisqr2))>1) && (prechisqr2<1.0e+02) && (prechisqr2>0) ) break;
	if(preiP<MinP) break;
      }
      
      if((prechisqr2<1.0e+03) && (prechisqr2>0)){
	/***** Sks Traking Start!! ************/
	SksTrack *track = new SksTrack( trIn, trOut );
	int K=0,kk=2;
	double iP,iP_pre=preP2;
	int Niteration;
	bool status;
	double chisqr,p_check,prechisqr,preP;
	
	for(;;){
	  if(K==0) iP = iP_pre;
	  track->SetInitialMomentum( iP );
	  //track->SetInitialMomentum( PIni );
	  track->doFit();
	  Niteration=track->Niteration();
	  p_check=track->PrimaryMomentum().Mag();
	  status=track->Status();
	  chisqr=track->chisqr();
	  
	  if(K==0){
	    prechisqr = chisqr ;
	    preP = iP ;
	  }
	  if(K>0){
	    if(prechisqr>chisqr){
	      prechisqr = chisqr ;
	      preP = iP ;
	    }
	  }
	  
#if 0
	  std::cout << K << " " << track->Niteration() << " " << chisqr 
		    << " " << p_check << " " << preP2 << std::endl;
#endif
	  
	  if(status || (p_check>0 && p_check<=MinP) || (iP<=MinP)) break;
	  if((K==0) && !(chisqr<0) && (p_check>MinP)){
	    iP_pre=p_check-0.01;
	  }
	  if(K!=0 && !(chisqr<0) && (p_check>MinP)){
	    iP=p_check-0.01;
	    if((iP_pre-iP)>0.01){
	      iP_pre=p_check-0.01;
	      kk = 2;
	    }
	    if(((iP_pre-iP)<=0.01)){ 
	      if(iP>0) iP=p_check-0.01*kk ;
	      else iP=iP_pre-0.01*kk ;
	      ++kk;
	    }
	  }
	  if(chisqr<0 || !(p_check>0)){
	    iP_pre=iP-0.01;
	    iP=iP-0.01;
	  }
	  
	  ++K;
	  //  iP=iP-0.005;
	  if(iP<MinP) break;
	}
	
	//std::cout << "chisqr = " << chisqr << " : momentum = " << p_check << std::endl;
	
	//std::cout << "SKSTracking finish!! ( # iteration = " << K << " )" << std::endl;      
	
	if(!track) continue;
	if( track->doFit() && track->chisqr()<MaxChiSqrSksTrack ){
	  SksTrackCol.push_back(track);
	}
	else{
	  delete track;
	}
      }
      delete pretrack;

#endif
      
//(B) -> Old SksTracking Rootine
#if ver1_Sks
      SksTrack *tp = new SksTrack( trIn, trOut );

      if(!tp) continue;
      // initial momentum from the SdcIn/SdcOut slopes: the
      // empirical dispersion 1/p = c0 + c1 uIn + c2 uOut fitted
      // over the beam-through scan (1.1-1.8 GeV/c) and the CH2
      // production tracks reproduces p within +-3%, which is a
      // condition-independent seed (a fixed value is up to 50%
      // off across the settings)
      {
	static const double c0 =  0.70990;
	static const double c1 = -0.06302;
	static const double c2 =  0.86731;
	static const double MinPini = 0.60, MaxPini = 2.20;
	double invp = c0 + c1*(trIn->GetU0()) + c2*(trOut->GetU0());
	if( invp > 1./MaxPini && invp < 1./MinPini )
	  tp->SetInitialMomentum( 1./invp );
      }
      if( tp->doFit() && tp->chisqr()<MaxChiSqrSksTrack ){
	SksTrackCol.push_back(tp);
      }
      else{
	delete tp;
      }
#endif
    }
  }

  partial_sort( SksTrackCol.begin(), SksTrackCol.end(),
		SksTrackCol.end(), SksTrackComp() );

#if 0
  std::cout<<"********************"<<std::endl;
 {
   int nn=SksTrackCol.size();
   std::cout << funcname << ": Before Deleting. #Track="
	     << nn << std::endl;
   for( int i=0; i<nn; ++i ){
     SksTrack *tp=SksTrackCol[i];
     std::cout << std::setw(3) << i
	       << " Nitra=" << std::setw(3) << tp->Niteration()
	       << " ChiSqr=" << tp->chisqr()
	       << " P=" << tp->PrimaryMomentum().Mag()
	       << " PL(TOF)=" << tp->PathLengthToTOF()
	       << std::endl;
   }
 }
#endif

  return true;
}

// bool DCAnalyzer::DecodeSimuHits( SimuData *simuData )
// {
//   const std::string funcname = "[DCAnalyzer::DecodeSimuHits]";

//   clearDCHits();
// #if 0
//   // BcIn
//   for( int layer=1; layer<=NumOfLayersBcIn; ++layer ){
//     const DCRHitContainer &cont=rawData->GetBcInRawHC(layer);
//     int nh=cont.size();
//     for( int i=0; i<nh; ++i ){
//        DCRawHit *rhit=cont[i];
//        if(!rhit) continue;
//        DCHit *hit=new DCHit(rhit->PlaneId()+PlOffsBc,
// 			    rhit->WireId(),rhit->Tdc());
//        if(!hit) continue;

//        if(hit->CalcObservables())
// 	 BcInHC[layer].push_back(hit);
//        else
// 	 delete hit;
//     } /* for( int i... ) */
//   } /* for( int layer... ) */
//  #endif

// #if 0
//   // BcOut
//   for( int layer=SimuPlMinBdcOut+SimuPlOffsBdc; layer<=SimuPlMaxBdcOut+SimuPlOffsBdc; ++layer ){
//     int nh=simuData->GetDcMulti(layer);
//     for( int i=0; i<nh; ++i ){
//       //Add 2010/09/29
//       DCHit *hit=new DCHit(layer, simuData->GetDcHitWire(layer, i));
//       //      DCHit *hit=new DCHit(layer, simuData->GetDcHitWire(layer, i), -1);
//       if(!hit) continue;
      
//       if(hit->CalcObservablesSimulation(simuData->GetDcDrLength(layer, i)))
//  	BcOutHC[layer-(SimuPlMinBdcOut+SimuPlOffsBdc-1)].push_back(hit);
//        else
// 	 delete hit;
//     } /* for( int i... ) */
//   } /* for( int layer... ) */
// #endif  

//   // SdcIn
//   for( int layer=1; layer<=10; ++layer ){
//     int nh=simuData->GetDcMulti(layer);
//     for( int i=0; i<nh; ++i ){
//       //Add 2010/09/29
//       DCHit *hit=new DCHit(layer, simuData->GetDcHitWire(layer, i));
//       //      DCHit *hit=new DCHit(layer,simuData->GetDcHitWire(layer, i),-1);
//       if(!hit) continue;
//       hit->SetTdcVal(-1);

//       if(hit->CalcObservablesSimulation(simuData->GetDcDrLength(layer, i))){
//  	SdcInHC[layer].push_back(hit);
//       } 
//       else {
// 	delete hit;
//       }
//     } /* for( int i... ) */
//   } /* for( int layer... ) */
      

//   // SdcOut
//   for( int layer=31; layer<=42; ++layer ){
//      int nh=simuData->GetDcMulti(layer);
//      for( int i=0; i<nh; ++i ){
//       //Add 2010/09/29
//       DCHit *hit=new DCHit(layer, simuData->GetDcHitWire(layer, i));
//        //       DCHit *hit=new DCHit(layer,simuData->GetDcHitWire(layer, i),-1); 
//        if(!hit) continue;
//        hit->SetTdcVal(-1);
       
//        if(hit->CalcObservablesSimulation(simuData->GetDcDrLength(layer, i)))
// 	 SdcOutHC[layer-30].push_back(hit);
//        else
// 	 delete hit;
//      } /* for( int i... ) */
//   } /* for( int layer... ) */

//   return true;
// }


void DCAnalyzer::clearDCHits( void )
{
  clearBcInHits();
  clearBcOutHits();
  clearSdcInHits();
  clearSdcOutHits();
  clearSsdHits();
}

void DCAnalyzer::clearBcInHits(){
  for( int l=0; l<=NumOfLayersBcIn; ++l ){
    for_each( TempBcInHC[l].begin(),  TempBcInHC[l].end(),  DeleteObject() );
    TempBcInHC[l].clear();
  }

  for( int l=0; l<=NumOfLayersBcIn; ++l ){
    for_each( BcInHC[l].begin(),  BcInHC[l].end(),  DeleteObject() );
    BcInHC[l].clear();
  }

  for ( int l=0; l<=NumOfLayersBcIn; ++l ){
    std::for_each( MWPCClCont[l].begin(),  MWPCClCont[l].end(), DeleteObject() );
    MWPCClCont[l].clear();
  }
}

void DCAnalyzer::clearBcOutHits(){
  for( int l=0; l<=NumOfLayersBcOut; ++l ){
    for_each( BcOutHC[l].begin(), BcOutHC[l].end(), DeleteObject() );
    BcOutHC[l].clear();
  }  
}

void DCAnalyzer::clearSdcInHits(){
  for( int l=0; l<=NumOfLayersSdcIn; ++l ){
    for_each( SdcInHC[l].begin(),  SdcInHC[l].end(),  DeleteObject() );
    SdcInHC[l].clear();
  }
}

void DCAnalyzer::clearSdcOutHits(){
  for( int l=0; l<=NumOfLayersSdcOut; ++l ){
    for_each( SdcOutHC[l].begin(), SdcOutHC[l].end(), DeleteObject() );
    SdcOutHC[l].clear();
  }
}

void DCAnalyzer::clearSsdHits(){
  for( int l=0; l<=NumOfLayersSsd; ++l ){
    for_each( SsdHC[l].begin(),  SsdHC[l].end(),  DeleteObject() );
    SsdHC[l].clear();
  }

  for( int l=0; l<=NumOfLayersSsd; ++l ){
    for_each( TempSsdHC[l].begin(),  TempSsdHC[l].end(),  DeleteObject() );
    TempSsdHC[l].clear();
  }
}

void DCAnalyzer::clearVtxHits( void )
{
  for_each( VtxPoint.begin(),  VtxPoint.end(),  DeleteObject() );
  VtxPoint.clear();
}

void DCAnalyzer::clearTracksBcIn( void )
{
  for_each( TrackBcInCol.begin(), TrackBcInCol.end(), DeleteObject() );
  TrackBcInCol.clear();
}

void DCAnalyzer::clearTracksBcOut( void )
{
  for_each( TrackBcOutCol.begin(), TrackBcOutCol.end(), DeleteObject() );
  TrackBcOutCol.clear();
}

void DCAnalyzer::clearTracksSdcIn( void )
{
  for_each( TrackSdcInCol.begin(), TrackSdcInCol.end(), DeleteObject() );
  TrackSdcInCol.clear();
}

void DCAnalyzer::clearTracksSdcOut( void )
{
  for_each( TrackSdcOutCol.begin(), TrackSdcOutCol.end(), DeleteObject() );
  TrackSdcOutCol.clear();
}

void DCAnalyzer::clearK18Tracks( void )
{
  for_each( K18TrackCol.begin(), K18TrackCol.end(), DeleteObject() );
  K18TrackCol.clear();
}

void DCAnalyzer::clearK18Tracks_BFT( void )
{
  for_each( K18TrackBFTCol.begin(), K18TrackBFTCol.end(), DeleteObject() );
  K18TrackBFTCol.clear();
}

void DCAnalyzer::clearSksTracks( void )
{
  for_each( SksTrackCol.begin(), SksTrackCol.end(), DeleteObject() );
  SksTrackCol.clear();
}

void DCAnalyzer::clearTracksBcOutSdcIn( void )
{
  for_each( TrackBcOutSdcInCol.begin(), TrackBcOutSdcInCol.end(),
            DeleteObject() );
  TrackBcOutSdcInCol.clear();
}

bool DCAnalyzer::ReCalcDCHits( bool applyRecursively )
{
  for( int l=0; l<=NumOfLayersBcIn; ++l ){
    int n=TempBcInHC[l].size();
    for( int i=0; i<n; ++i ){
      DCHit *hit=(TempBcInHC[l])[i];
      if(hit) hit->ReCalcMWPC(applyRecursively);
    }
  }
  for( int l=0; l<=NumOfLayersBcIn; ++l ){
    int n=BcInHC[l].size();
    for( int i=0; i<n; ++i ){
      DCHit *hit=(BcInHC[l])[i];
      if(hit) hit->ReCalcMWPC(applyRecursively);
    }
  }
  for( int l=0; l<=NumOfLayersBcOut; ++l ){
    int n=BcOutHC[l].size();
    for( int i=0; i<n; ++i ){
      DCHit *hit=(BcOutHC[l])[i];
      if(hit) hit->ReCalcDC(applyRecursively);
    }
  }
  for( int l=0; l<=NumOfLayersSdcIn; ++l ){
    int n=SdcInHC[l].size();
    for( int i=0; i<n; ++i ){
      DCHit *hit=(SdcInHC[l])[i];
      if(hit) hit->ReCalcDC(applyRecursively);
    }
  }
  for( int l=0; l<=NumOfLayersSdcOut; ++l ){
    int n=SdcOutHC[l].size();
    for( int i=0; i<n; ++i ){
      DCHit *hit=(SdcOutHC[l])[i];
      if(hit) hit->ReCalcDC(applyRecursively);
    }
  }

  return true;
}

bool DCAnalyzer::ReCalcTrackBcIn( bool applyRecursively )
{
  int n=TrackBcInCol.size();
  for( int i=0; i<n; ++i ){
    DCLocalTrack *track=TrackBcInCol[i];
    if( track ) track->ReCalc( applyRecursively );
  }
  return true;
}

bool DCAnalyzer::ReCalcTrackBcOut( bool applyRecursively )
{
  int n=TrackBcOutCol.size();
  for( int i=0; i<n; ++i ){
    DCLocalTrack *track=TrackBcOutCol[i];
    if( track ) track->ReCalc( applyRecursively ); 
  }
  return true;
}

bool DCAnalyzer::ReCalcTrackSdcIn( bool applyRecursively )
{
  int n=TrackSdcInCol.size();
  for( int i=0; i<n; ++i ){
    DCLocalTrack *track=TrackSdcInCol[i];
    if( track ) track->ReCalc( applyRecursively );
  }
  return true;
}

bool DCAnalyzer::ReCalcTrackSdcOut( bool applyRecursively )
{
  int n=TrackSdcOutCol.size();
  for( int i=0; i<n; ++i ){
    DCLocalTrack *track=TrackSdcOutCol[i];
    if( track ) track->ReCalc( applyRecursively );
  }
  return true;
}

bool DCAnalyzer::ReCalcK18Track( bool applyRecursively )
{
  int n=K18TrackCol.size();
  for( int i=0; i<n; ++i ){
    K18Track *track=K18TrackCol[i];
    if( track ) track->ReCalc( applyRecursively );
  }
  return true;
}

bool DCAnalyzer::ReCalcSksTrack( bool applyRecursively )
{
  int n=SksTrackCol.size();
  for( int i=0; i<n; ++i ){
    SksTrack *track=SksTrackCol[i];
    if( track ) track->ReCalc( applyRecursively );
  }
  return true;
}

bool DCAnalyzer::ReCalcAll( void )
{
  ReCalcDCHits();

  ReCalcTrackBcIn();
  ReCalcTrackBcOut();
  ReCalcTrackSdcIn();
  ReCalcTrackSdcOut();

  ReCalcK18Track();
  ReCalcSksTrack();

  return true;
}

//______________________________________________________________________________
double
log2(double x)
{
  return std::log(x)/std::log(2.0);
}

//______________________________________________________________________________
bool
isConnectable(double wire1,
	      double leading1,
	      double trailing1,
	      double wire2,
	      double leading2,
	      double trailing2,
	      double wExt,
	      double tExt)
{
  double w1Min = wire1 - wExt;
  double w1Max = wire1 + wExt;
  double t1Min = leading1  - tExt;
  double t1Max = trailing1 + tExt;

  double w2Min = wire2 - wExt;
  double w2Max = wire2 + wExt;
  double t2Min = leading2  - tExt;
  double t2Max = trailing2 + tExt;

  bool isWireOk = !(w1Min>w2Max || w1Max<w2Min);
  bool isTimeOk = !(t1Min>t2Max || t1Max<t2Min);

  //   std::cout << " w1 = " << wire1
  // 	    << " le1 = " << leading1
  // 	    << " tr1 = " << trailing1 << "\n"
  // 	    << " w2 = " << wire2
  // 	    << " le2 = " << leading2
  // 	    << " tr2 = " << trailing2 << "\n"
  // 	    << " w1(" << w1Min << " -- " << w1Max << "), t1(" 
  // 	    << t1Min << " -- " << t1Max << ")\n" 
  // 	    << " w2(" << w2Min << " -- " << w2Max << "), t2(" 
  // 	    << t2Min << " -- " << t2Max << ")\n" 
  // 	    << " wire : " << isWireOk
  // 	    << ", time : " << isTimeOk
  // 	    << std::endl;
  
  return (isWireOk && isTimeOk);
}

//______________________________________________________________________________
void
printConnectionFlag(const std::vector<std::deque<bool> >& flag)
{
  for(int i=0, n=flag.size(); i<n; ++i){
    std::cout << "\n";
    for (int j=0, m=flag[i].size(); j<m; ++j){
      std::cout << " " << flag[i][j];
    }
  }
  std::cout << std::endl;

  return;
}

//______________________________________________________________________________
int
clusterizeMWPCHit(const DCHitContainer& hits,
		  MWPCClusterContainer& clusters)
{
  if (!clusters.empty()){
      std::for_each(clusters.begin(), clusters.end(), DeleteObject());
      clusters.clear();
  }
  
  const int nhits = hits.size();
  //   std::cout << "#D " << __func__ << " " << nhits << std::endl;
  if (nhits==0)
    return 0;

  int n = 0;
  for (int i=0; i<nhits; ++i){
    const DCHit* h = hits[i];
    if (!h)
      continue;
    n += h->GetTdcSize();
  }

  DCHitContainer singleHits;
  singleHits.reserve(n);
  for (int i=0; i<nhits; ++i){
    const DCHit* h = hits[i];
    if (!h)
      continue;
    int nn = h->GetTdcSize();
    for (int ii=0; ii<nn; ++ii){
      DCHit* htmp = new DCHit(h->GetLayer(), h->GetWire());
      htmp->SetTdcVal(h->GetTdcVal());
      htmp->SetTdcTrailing(h->GetTdcTrailing());
      htmp->SetTrailingTime(h->GetTrailingTime());
      htmp->SetDriftTime(h->GetDriftTime());
      htmp->SetDriftLength(h->GetDriftLength());
      htmp->SetTiltAngle(h->GetTiltAngle());
      htmp->SetWirePosition(h->GetWirePosition());
      htmp->setRangeCheckStatus(h->rangecheck(), 0);
      singleHits.push_back(htmp);
    }
  }
  
  std::vector<std::deque<bool> > flag(n, std::deque<bool>(n, false));
  n = singleHits.size();
  for (int i=0;  i<n; ++i){
    flag[i][i] = true;
    const DCHit* h1 = singleHits[i];
    //       h1->print("h1");
    for (int j=i+1; j<n; ++j){
      const DCHit* h2 = singleHits[j];
      // 	  h2->print("h2");
      // 	  std::cout << " (i,j) = (" << i << ", " << j << ")" << std::endl;
      bool val
	= isConnectable(h1->GetWirePosition(), 
			h1->GetDriftTime(),
			h1->GetTrailingTime(),
			h2->GetWirePosition(),
			h2->GetDriftTime(),
			h2->GetTrailingTime(),
			kMWPCClusteringWireExtension,
			kMWPCClusteringTimeExtension);
      // 	  std::cout << "#D val = " << val << std::endl;
      flag[i][j] = val;
      flag[j][i] = val;
    }
  }
  
  //   std::cout << "#D " << __func__ << "  before " << std::endl;
  //   printConnectionFlag(flag);
  
  const int maxLoop = static_cast<int>(log2(n))+1;
  for (int loop=0; loop<maxLoop; ++loop){
    std::vector<std::deque<bool> > tmp(n, std::deque<bool>(n, false));
    for (int i=0; i<n; ++i){
      for (int j=i; j<n; ++j){
	for (int k=0; k<n; ++k){
	  tmp[i][j] |= (flag[i][k] && flag[k][j]);
	  tmp[j][i] = tmp[i][j];
	}
      }
    }
    flag = tmp;
    //       std::cout << " n iteration = " << loop << std::endl;
    //       printConnectionFlag(flag);
  }
  
  //   std::cout << "#D " << __func__ << "  after " << std::endl;
  //   printConnectionFlag(flag);
  
  std::set<int> checked;
  for (int i=0; i<n; ++i){
    if (checked.find(i)!=checked.end())
      continue;
    MWPCCluster* c = 0;
    for (int j=i; j<n; ++j){
      if (flag[i][j]){
	checked.insert(j);
	if (!c) {
	  c = new MWPCCluster;
	  // 		  std::cout << " new cluster " << std::endl;
	}
	// 	      std::cout << " " << i << "---" << j << std::endl;
	c->Add(singleHits[j]);
      }
    }
    
    if (c){
      c->Calculate();
      clusters.push_back(c);
    }
  }
  
  //   std::cout << " end of " << __func__
  // 	    << " : n = " << n << ", " << checked.size()
  // 	    << std::endl;
  
  //   std::cout << __func__ << " n clusters = " << clusters.size() << std::endl;

  return clusters.size();
}

// 

// ChiSqrCutBcOut
void DCAnalyzer::ChiSqrCutBcOut(double chisqr){
  ChiSqrCut(TrackBcOutCol, chisqr);
}

// ChiSqrCut (implemention of ChiSqrCutXXX)
void DCAnalyzer::ChiSqrCut(std::vector<DCLocalTrack*>& cont, 
			   double chisqr)
{
  std::vector<DCLocalTrack*> DeleteCand;
  std::vector<DCLocalTrack*> ValidCand;
  int NofTrack = cont.size();
  for(int i = NofTrack-1; i>=0; --i){
    DCLocalTrack* tempTrack = cont.at(i);
    if(tempTrack->GetChiSquare() > chisqr){
      DeleteCand.push_back(tempTrack);
    }else{
      ValidCand.push_back(tempTrack);
    }
  }

  for_each(DeleteCand.begin(), DeleteCand.end(), DeleteObject());
  DeleteCand.clear();
  
  cont.clear();
  cont.resize(ValidCand.size());
  std::copy(ValidCand.begin(), ValidCand.end(), cont.begin());
  ValidCand.clear();
}
//________________________________________________________________________
/*
void DCAnalyzer::clusterizeSSDHit()
{
  
  for( int l=0; l<=NumOfLayersSsd; ++l )
    {
      std::for_each( TempSsdHC[l].begin(),  TempSsdHC[l].end(),  DeleteObject() );
      TempSsdHC[l].clear();
    }
  for(int layer=1;layer<NumOfLayersSsd+1;++layer)
    {
      int nhits = SsdHC[layer].size();
      //      std::cout<<"#D Before Clusterize SsdHC["<<layer<<"] size = "
      //	       <<nhits<<std::endl;
      for(int i=0;i<nhits;++i)
	{
	  DCHit *h1 = SsdHC[layer][i];
	  if(i+1<nhits)
	    {
	      DCHit *h2 = SsdHC[layer][i+1];
	      int wire1 = int(h1->GetWire());
	      int wire2 = int(h2->GetWire());
	      if(fabs(wire2-wire1)==1)
		{
		  DCHit *htmp = new DCHit(h1->GetLayer(),h1->GetWire());
		  for(int j=0;j<h1->GetAdcSize();++j)
		    htmp->SetAdcVal(h1->GetAdcVal(j));
		  for(int j=0;j<h1->GetTdcSize();++j)
		    htmp->SetTdcVal(h1->GetTdcVal(j));
		  TempSsdHC[layer].push_back(htmp);
		  ++i;
		}
	      else
		{
		  DCHit *htmp = new DCHit(h1->GetLayer(),h1->GetWire());
		  for(int j=0;j<h1->GetAdcSize();++j)
		    htmp->SetAdcVal(h1->GetAdcVal(j));
		  for(int j=0;j<h1->GetTdcSize();++j)
		    htmp->SetTdcVal(h1->GetTdcVal(j));
		  TempSsdHC[layer].push_back(htmp);
		}
	    }
	  else
	    {
	      DCHit *htmp = new DCHit(h1->GetLayer(),h1->GetWire());
	      for(int j=0;j<h1->GetAdcSize();++j)
		htmp->SetAdcVal(h1->GetAdcVal(j));
	      for(int j=0;j<h1->GetTdcSize();++j)
		htmp->SetTdcVal(h1->GetTdcVal(j));
	      TempSsdHC[layer].push_back(htmp);
	    }
	}
    }

  for( int l=0; l<=NumOfLayersSsd; ++l )
    {
      for_each( SsdHC[l].begin(),  SsdHC[l].end(),  DeleteObject() );
      SsdHC[l].clear();
    }
  // Re push_back SsdHC
  for(int layer=1;layer<NumOfLayersSsd+1;++layer)
    {
      for(unsigned int i=0;i<TempSsdHC[layer].size();++i)
	{
	  DCHit *h = TempSsdHC[layer][i];
	  DCHit *hnew = new DCHit(h->GetLayer(),h->GetWire());
	  for(int j=0;j<h->GetAdcSize();++j)
	    hnew->SetAdcVal(h->GetAdcVal(j));
	  for(int j=0;j<h->GetTdcSize();++j)
	    hnew->SetTdcVal(h->GetTdcVal(j));
	  if(hnew->CalcSsdObservables())
	    SsdHC[layer].push_back(hnew);
	  else
	    delete hnew;
	}
      //      std::cout<<"#D After Clusterize SsdHC["<<layer<<"] size = "
      //	       <<SsdHC[layer].size()<<std::endl;
    }
  return;
}
*/

bool DCAnalyzer::TDCrearrangeSDCout( DCAnalyzer *DCAna, HodoAnalyzer *HodoAna )
{
  DCAna->ClearDiffTrackSeg();
  std::cout<<"###########################TDC rearrange#########################"<<std::endl;
  //std::cout<<"###Size of check parameter before:   "<<DCAna->GetSizeDiffTrackSeg()<<std::endl;
  const std::string funcname = "[DCAnalyzer::TDCrearrangeSDCout]";
  int nt = DCAna -> GetNtracksSdcOut();
  int nttof = HodoAna -> GetNHitsTOFMT();
  std::cout<<"###Size of NTsdcout:   "<<nt<<std::endl;
  std::cout<<"###Size of NHtof:   "<<nttof<<std::endl;
  int rearrange_flag = 0;
  //DCAna->ResizeDiffTrackSeg(nt);
  if (nt==0) 
    {
      std::cout<<"#there is no sdc out tracking"<<std::endl;
      //DCAna->SetDiffTrackSeg(-10000.);
      return false;
    }
 
  else if (nt>0 && nttof==0)
    {
      std::cout<<"#sdc out tracking but no tofmt"<<std::endl;
      for( int l=0 ; l<nt; l++)
	{
	  DCAna->SetDiffTrackSeg(-2000.);
	}
     }
  else if (nt>0 && nttof>0)
    {
      if (nttof==1) std::cout<<"#only one tof hit"<<std::endl;
      //std::cout<<"#####check that TDC stop come from real segment"<<std::endl;
      ConfMan *confMan =ConfMan::GetConfManager();
      DCDriftParamMan *driftMan=confMan->GetDCDriftParamManager();
      double tofz = DCGeomMan::GetInstance().GetLocalZ(51);
      for ( int i = 0; i<nt ; i++)
	{
  	  DCLocalTrack *cont1 = DCAna->GetTrackSdcOut(i);
  	  //TOF segment and tracking X match
	  double x = cont1->GetX(tofz) ;
	  double getsegtof = DCAna->GetSegTof(x);
	  int  good_seg = -10;
	  double good_diff = -500.;
	  int good_hitnum=0;
	  for ( int k = 0; k<nttof; k++)
	    {
	      Hodo1Hit *cont_tof1 = HodoAna->GetHitTOFMT(k);
	      int SegTof = cont_tof1->SegmentId();
	      double diff = getsegtof-SegTof;
	      if(fabs(good_diff)>diff) 
		{
		  good_diff=diff;
		  good_seg=SegTof;
		  good_hitnum = k;
		}
	      //std::cout<<"##TofSeg from Tracking:   "<<getsegtof<<std::endl;
	      //std::cout<<"##TofSeg hit:   "<<SegTof<<std::endl;
	      //std::cout<<"##Good TofSeg hit:   "<<good_seg<<std::endl;
	    }
	  //std::cout<<"##Good Diff:   "<<good_diff<<std::endl;
 	  DCAna->SetDiffTrackSeg(good_diff);
	  Hodo1Hit *cont_tof2 = HodoAna->GetHitTOFMT(good_hitnum);
	  double diff_dt = cont_tof2 -> GetT();
	  //std::cout<<"###dt difference:   "<<diff_dt<<std::endl;
	  //std::cout<<"###check find seg:   "<<cont_tof2->SegmentId()<<std::endl;

 	  if (fabs(good_diff)>2.0 && diff_dt>2. &&diff_dt<100. && fabs(diff_dt) < .5 )
		{
		  std::cout<<"#####There is no good TOF seg"<<std::endl;
		  continue;
		  
		}
 	  int nh = cont1->GetNHit();
	  rearrange_flag++;
	  //std::cout<<"###rearrange_flag:     "<<rearrange_flag<<std::endl;
	  //std::cout << "##nh of track: "<< nh <<std::endl;
	  for( int j = 0; j <nh ; j++)
	    {
	      DCLTrackHit *cont2 = cont1->GetHit(j);
	      int layer = cont2->GetLayer();
	      double wire = cont2->GetWire();
	      double dt1 = cont2->GetDriftTime();
	      //std::cout << "TDC value before reset:           "<<dt1<<std::endl;
	      cont2->ClearDriftTime();
	      cont2->SetDriftTime(dt1+diff_dt);
	      double dt = cont2->GetDriftTime();
	      double dl;
	      bool status=driftMan->calcDrift_sdc34( layer, wire, dt, dl);
	      if(!status)
		{
		  //std::cout <<"####SDC34 TDC rearrange driftMan error"<<std::endl;
		  continue;
		}
	      else{
		cont2->ClearDriftLength();
		cont2->SetDriftLength(dl);
		//TDC = cont2->GetDriftTime();
		//std::cout << "####TDC value after reset:           "<<dt<<std::endl;
	      }
	    }
	}
      //std::cout<<"###Siz of check parameter after   "<<DCAna->GetSizeDiffTrackSeg()<<std::endl;
      if(rearrange_flag>0) return true;
      else return false;
    }
   else {
     //std::cout<<"###Siz of check parameter after   "<<DCAna->GetSizeDiffTrackSeg()<<std::endl;
    return false;}
  return false;
}

double DCAnalyzer::GetSegTof( double x )
{
  double p0=16.83;
  double p1=0.0137949;
  return p0+x*p1;
}
