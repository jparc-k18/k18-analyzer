/*
  HodoAnalyzer.cc

  2012/1/24
*/

#include "HodoAnalyzer.hh"
#include "RawData.hh"
#include "Hodo2Hit.hh"
#include "Hodo1Hit.hh"
#include "BH2Hit.hh"
#include "FiberHit.hh"
#include "HodoCluster.hh"
#include "BH2Cluster.hh"
#include "FLHit.hh"
#include "FiberCluster.hh"

#include "TemplateLib.hh"
#include <iostream>
#include <iomanip>
#include <string>
#include <cmath>
#include <algorithm>
#include <cstdlib>

static const std::string MyName = "HodoAnalyzer";
HodoAnalyzer *HodoAnalyzer::MySelf_ = NULL;

const double MaxTimeDifBh1 = 2.0;
const double MaxTimeDifBh2 = 2.0;
const double MaxTimeDifTof = 3.5;
const double MaxTimeDifLc  = 5.5;
const double MaxTimeDifBFT = 8.0;
const double MaxTimeDifSFT = 8.0;

#define Cluster 1
#define preSelection_BFT 0
#define MultiHit_BH      0

HodoAnalyzer::HodoAnalyzer()
//  :fl_decoded_SFT(false)
{
  if(MySelf_){
    std::cerr << "#E : [HodoAnalyzer::Constructor] instance is already created" << std::endl;
    exit(-1);
  }else{
    MySelf_ = this;
  }
}

HodoAnalyzer::~HodoAnalyzer()
{
  clearBFTHits();
  clearLCHits();
  clearTOFHits();
  clearBH2Hits();
  clearBH1Hits();

  clearBMWHits();
  clearBACHits();
  clearSFVHits();
  clearTOFMTHits();

  MySelf_ = NULL;
}


void HodoAnalyzer::clearBH1Hits()
{
  for_each(BH1Cont.begin(),BH1Cont.end(),DeleteObject());
  for_each(BH1ClCont.begin(),BH1ClCont.end(),DeleteObject());
}

void HodoAnalyzer::clearBH2Hits()
{
  for_each(BH2Cont.begin(),BH2Cont.end(),DeleteObject());
  for_each(BH2ClCont.begin(),BH2ClCont.end(),DeleteObject());
}

void HodoAnalyzer::clearBFTHits()
{
  for_each(BFTCont.begin(), BFTCont.end(), DeleteObject());
  for_each(BFTClCont.begin(), BFTClCont.end(), DeleteObject());
}


void HodoAnalyzer::clearTOFHits( void )
{
  for_each( TOFCont.begin(), TOFCont.end(), DeleteObject() );
  for_each( TOFClCont.begin(), TOFClCont.end(), DeleteObject() );
}

void HodoAnalyzer::clearLCHits( void )
{
  for_each( LCCont.begin(), LCCont.end(), DeleteObject() );
  for_each( LCClCont.begin(), LCClCont.end(), DeleteObject() );
}

void HodoAnalyzer::clearBMWHits()
{
  for_each(BMWCont.begin(),BMWCont.end(),DeleteObject());
}
void HodoAnalyzer::clearBACHits()
{
  for_each(BACCont.begin(),BACCont.end(),DeleteObject());
}
void HodoAnalyzer::clearSFVHits()
{
  for_each(SFVCont.begin(),SFVCont.end(),DeleteObject());
}
void HodoAnalyzer::clearTOFMTHits()
{
  for_each(TOFMTCont.begin(),TOFMTCont.end(),DeleteObject());
}


bool HodoAnalyzer::DecodeRawHits( RawData *rawData )
{
  //DecodeGCHits( rawData );
  DecodeBH1Hits( rawData );
  DecodeBH2Hits( rawData );
  DecodeBFTHits( rawData );
  //DecodeSFTHits( rawData );
  DecodeTOFHits( rawData );
  //DecodeACHits( rawData );
  DecodeLCHits( rawData );

  DecodeBMWHits( rawData );
  DecodeBACHits( rawData );
  DecodeSFVHits( rawData );
  DecodeTOFMTHits( rawData );

  return true;
}

bool HodoAnalyzer::DecodeBH1Hits( RawData *rawData )
{
  clearBH1Hits();
  
  const HodoRHitContainer &cont=rawData->GetBH1RawHC();
  int nh=cont.size();
  for( int i=0; i<nh; ++i ){
    HodoRawHit *hit=cont[i];

    int NofMHit = 1;
#if multihit_BH
    int NofU = hit->GetSizeTdcUp(); int NofD = hit->GetSizeTdcDown();
    if(NofU == NofD){ NofMHit = NofU; }
    else{NofMHit = (NofU < NofD) ? NofU : NofD;}
#endif

    if( !hit ) continue;
    for(int m = 0; m<NofMHit; ++m){
      int Tu=hit->GetTdcUp(m), Td=hit->GetTdcDown(m);
      if( Tu>0 && Td>0 ){
	Hodo2Hit *hp=new Hodo2Hit( hit, m );
	if( !hp ) continue;
	if( hp->calculate() )
	  BH1Cont.push_back(hp);
	else
	  delete hp;
      }    
    }
  }
  
#if Cluster
  MakeUpClusters( BH1Cont, BH1ClCont, MaxTimeDifBh1 );
#endif
  
  return true;
}

bool HodoAnalyzer::DecodeBFTHits(RawData* rawData){

  clearBFTHits();

#if preSelection_BFT
  static const double min = 725;
  static const double max = 750;
  
  static const double a = -5.641e-3;
  static const double b = 7.396e-1;
#endif

  for(int p = 0; p<NumOfPlaneBFT; ++p){
    const HodoRHitContainer &cont = rawData->GetBFTRawHC(p);
    int nh=cont.size();
    for( int i=0; i<nh; ++i ){
      HodoRawHit *hit=cont[i];
      if( !hit ) continue;

#if preSelection_BFT
      // Find valid signal from multi-hit
      int mhit = hit->GetSizeTdcUp();
      int bestTdcPair[2] = {-1, -1};
      for(int m = 0; m<mhit; ++m){
	int leading  = hit->GetTdcUp(m); // Up means leading edge in this case.
	int trailing = hit->GetTdcDown(m); // Down means trailing.
	double width = leading - trailing;
	if(true 
	   && (double)leading < (a*width*width+b*width+max)
	   && (double)leading > (a*width*width+b*width+min)
	   && leading > bestTdcPair[0]
	   ){
	  bestTdcPair[0] = leading;
	  bestTdcPair[1] = trailing;
	}
      }

      // There is no valid signal
      if(-1 == bestTdcPair[0]){continue;}

      // Reset with valid signal (other hits are removed)
      hit->clear();
      hit->SetTdcUp(bestTdcPair[0]);
      hit->SetTdcDown(bestTdcPair[1]);
#endif
      
      FiberHit *hp = new FiberHit(hit, "BFT");
      if(!hp) continue;
      if(hp->calculate()){
	BFTCont.push_back(hp);
      }
      else{
	delete hp;
	hp = NULL;
      }
    }
  }
  std::sort(BFTCont.begin(), BFTCont.end(), FiberHit::CompFiberHit);

#if Cluster
  MakeUpClusters( BFTCont, BFTClCont, MaxTimeDifBFT );
#endif

  return true;
}

bool HodoAnalyzer::DecodeBH2Hits( RawData *rawData )
{
  clearBH2Hits();
  
  const HodoRHitContainer &cont=rawData->GetBH2RawHC();
  int nh=cont.size();
  for( int i=0; i<nh; ++i ){
    HodoRawHit *hit=cont[i];

    if( !hit ) continue;
    int T=hit->GetTdcUp(0);
    if( T>0 ){
      BH2Hit *hp=new BH2Hit( hit );
      if( !hp ) continue;
      if( hp->calculate() )
	BH2Cont.push_back(hp);
      else
	delete hp;
    }
  }
  
#if Cluster
  MakeUpClusters( BH2Cont, BH2ClCont, MaxTimeDifBh2 );
#endif
    
  return true;
}

bool HodoAnalyzer::DecodeTOFHits( RawData *rawData )
{
  clearTOFHits();

  const HodoRHitContainer &cont=rawData->GetTOFRawHC();
  int nh=cont.size();
  for( int i=0; i<nh; ++i ){
    HodoRawHit *hit=cont[i];
    if( !hit ) continue;
    int Tu=hit->GetTdcUp(), Td=hit->GetTdcDown();
    if( Tu>0 && Td>0 ){
      Hodo2Hit *hp=new Hodo2Hit( hit );
      if( !hp ) continue;
      if( hp->calculate() )
	TOFCont.push_back(hp);
      else
	delete hp;
    }
  }

#if Cluster
  MakeUpClusters( TOFCont, TOFClCont, MaxTimeDifTof );
#endif

  return true;
}

bool HodoAnalyzer::DecodeLCHits( RawData *rawData )
{
  clearLCHits();

  const HodoRHitContainer &cont=rawData->GetLCRawHC();
  int nh=cont.size();
  for( int i=0; i<nh; ++i ){
    HodoRawHit *hit=cont[i];
    if( !hit ) continue;
    int Tu=hit->GetTdcUp(), Td=hit->GetTdcDown();
    if( Tu>0 && Td>0 ){
      Hodo2Hit *hp=new Hodo2Hit( hit );
      if( !hp ) continue;
      if( hp->calculate() )
	LCCont.push_back(hp);
      else
	delete hp;
    }
  }

#if Cluster
  MakeUpClusters( LCCont, LCClCont, MaxTimeDifLc );
#endif

  return true;
}


bool HodoAnalyzer::DecodeBMWHits( RawData *rawData )
{
  clearBMWHits();
  const HodoRHitContainer &cont=rawData->GetBMWRawHC();
  int nh=cont.size();
  for( int i=0; i<nh; ++i ){
    HodoRawHit *hit=cont[i];
    if( !hit ) continue;
    Hodo1Hit *hp=new Hodo1Hit( hit );
    if( !hp ) continue;
    if( hp->calculate() )
      BMWCont.push_back(hp);
    else
      delete hp;
  }
  return true;
}
bool HodoAnalyzer::DecodeBACHits( RawData *rawData )
{
  clearBACHits();
  const HodoRHitContainer &cont=rawData->GetBACRawHC();
  int nh=cont.size();
  for( int i=0; i<nh; ++i ){
    HodoRawHit *hit=cont[i];
    if( !hit ) continue;
    Hodo1Hit *hp=new Hodo1Hit( hit );
    if( !hp ) continue;
    if( hp->calculate() )
      BACCont.push_back(hp);
    else
      delete hp;
  }
  return true;
}
bool HodoAnalyzer::DecodeSFVHits( RawData *rawData )
{
  clearSFVHits();
  const HodoRHitContainer &cont=rawData->GetSFVRawHC();
  int nh=cont.size();
  for( int i=0; i<nh; ++i ){
    HodoRawHit *hit=cont[i];
    if( !hit ) continue;
    Hodo1Hit *hp=new Hodo1Hit( hit );
    if( !hp ) continue;
    if( hp->calculate() )
      SFVCont.push_back(hp);
    else
      delete hp;
  }
  return true;
}
bool HodoAnalyzer::DecodeTOFMTHits( RawData *rawData )
{
  clearTOFMTHits();
  const HodoRHitContainer &cont=rawData->GetTOFMTRawHC();
  int nh=cont.size();
  for( int i=0; i<nh; ++i ){
    HodoRawHit *hit=cont[i];
    if( !hit ) continue;
    Hodo1Hit *hp=new Hodo1Hit( hit );
    if( !hp ) continue;
    if( hp->calculate() )
      TOFMTCont.push_back(hp);
    else
      delete hp;
  }
  return true;
}


int HodoAnalyzer::
MakeUpClusters( const Hodo2HitContainer & HitCont,
		HodoClusterContainer & ClusterCont, double maxTimeDif )
{
  static const std::string funcname = "HodoAnalyzer::MakeUpClusters";
  
  if( !ClusterCont.empty() )
    for_each( ClusterCont.begin(), ClusterCont.end(), DeleteObject() );
  
  int nh=HitCont.size();
  
  std::vector <int> flag(nh,0);
  
  for(int i=0; i<nh; ++i ){
    if( flag[i] ) continue;
    Hodo2Hit *hitA=HitCont[i];
    int segA=hitA->SegmentId();
    double cmtA=hitA->CMeanTime();
    Hodo2Hit *hitB=0;
    int iB=-1;
    double cmtB;
    int segB;
    for( int j=i+1; j<nh; ++j ){
      Hodo2Hit *hit=HitCont[j];
      int seg=hit->SegmentId();
      double cmt=hit->CMeanTime();
      if( abs(seg-segA)==1 && fabs(cmt-cmtA)<maxTimeDif ){
	hitB=hit; ++flag[j]; iB=j; segB=seg; cmtB=cmt; break;
      }
    }
    if(hitB){
      Hodo2Hit *hitC=0;
      for( int j=i+1; j<nh; ++j ){
	if( j==iB ) continue;
	Hodo2Hit *hit=HitCont[j];
	int seg=hit->SegmentId();
	double cmt=hit->CMeanTime();
	if( (abs(seg-segA)==1 && fabs(cmt-cmtA)<maxTimeDif) ||
	    (abs(seg-segB)==1 && fabs(cmt-cmtB)<maxTimeDif) ){
	  hitC=hit; ++flag[j]; break;
	}
      }
      if(hitC){
	HodoCluster *cluster=new HodoCluster(hitA,hitB,hitC);
	if( cluster ) ClusterCont.push_back(cluster);
      }
      else{
	HodoCluster *cluster=new HodoCluster(hitA,hitB);
	if( cluster ) ClusterCont.push_back(cluster);
      }
    }
    else{
      HodoCluster *cluster=new HodoCluster(hitA);
      if( cluster ) ClusterCont.push_back(cluster);
    }
  }

  return ClusterCont.size(); 
}				  

int HodoAnalyzer::
MakeUpClusters( const BH2HitContainer & HitCont,
		BH2ClusterContainer & ClusterCont, double maxTimeDif )
{
  static const std::string funcname = "HodoAnalyzer::MakeUpClusters";

  if( !ClusterCont.empty() )
    for_each( ClusterCont.begin(), ClusterCont.end(), DeleteObject() );

  int nh=HitCont.size();

  std::vector <int> flag(nh,0);

  for(int i=0; i<nh; ++i ){
    if( flag[i] ) continue;
    BH2Hit *hitA=HitCont[i];
    int segA=hitA->SegmentId();
    double cmtA=hitA->CTime();
    BH2Hit *hitB=0;
    int iB=-1;
    double cmtB;
    int segB;
    for( int j=i+1; j<nh; ++j ){
      BH2Hit *hit=HitCont[j];
      int seg=hit->SegmentId();
      double cmt=hit->CTime();
      if( abs(seg-segA)==1 && fabs(cmt-cmtA)<maxTimeDif ){
	hitB=hit; ++flag[j]; iB=j; segB=seg; cmtB=cmt; break;
      }
    }
    if(hitB){
      BH2Hit *hitC=0;
      for( int j=i+1; j<nh; ++j ){
        if( j==iB ) continue;
        BH2Hit *hit=HitCont[j];
        int seg=hit->SegmentId();
        double cmt=hit->CTime();
        if( (abs(seg-segA)==1 && fabs(cmt-cmtA)<maxTimeDif) ||
            (abs(seg-segB)==1 && fabs(cmt-cmtB)<maxTimeDif) ){
          hitC=hit; ++flag[j]; break;
        }
      }
      if(hitC){
        BH2Cluster *cluster=new BH2Cluster(hitA,hitB,hitC);
        if( cluster ) ClusterCont.push_back(cluster);
      }
      else{
	BH2Cluster *cluster=new BH2Cluster(hitA,hitB);
	if( cluster ) ClusterCont.push_back(cluster);
      }
    }
    else{
      BH2Cluster *cluster=new BH2Cluster(hitA);
      if( cluster ) ClusterCont.push_back(cluster);
    }
  }

  return ClusterCont.size(); 
}				  

int HodoAnalyzer::MakeUpClusters(const FiberHitContainer& cont,
				 FiberClusterContainer& ClusterCont,
				 double maxTimeDif)
{
  if( !ClusterCont.empty() ){
    for_each( ClusterCont.begin(), ClusterCont.end(), DeleteObject() );
  }

  int NofSeg = cont.size();
  for(int seg = 0; seg<NofSeg; ++seg){
    FiberHit* HitA = cont.at(seg);

    bool fl_ClCandA  = false;
    if(seg != (NofSeg -1)){
      if(3 > (cont.at(seg+1)->PairId() - HitA->PairId())){fl_ClCandA = true;}
    }

    int NofHitA = HitA->GetNumOfHit();
    for(int mhitA = 0; mhitA<NofHitA; ++mhitA){
      if(HitA->Joined(mhitA)){continue;}

      FiberCluster *cluster = new FiberCluster();
      cluster->push_back(new FLHit(HitA, mhitA));
      
      if(!fl_ClCandA){
	// there is no more candidates
	if(cluster->calculate()){ ClusterCont.push_back(cluster); }
	continue;
      }
      
      // Start Search HitB
      double cmtA    = (double)HitA->GetCTime(mhitA);
      int    NofHitB = cont.at(seg+1)->GetNumOfHit();
      bool   fl_HitB = false;
      double cmtB    = -1;
      int    CurrentPair = HitA->PairId();
      for(int mhitB = 0; mhitB<NofHitB; ++mhitB){
	if(cont.at(seg+1)->Joined(mhitB)){continue;}

	FiberHit* HitB = cont.at(seg+1);
	
	cmtB = (double)HitB->GetCTime(mhitB);
	if(fabs(cmtB-cmtA)<maxTimeDif){
	  cluster->push_back(new FLHit(HitB, mhitB));
	  CurrentPair = HitB->PairId();
	  fl_HitB = true;
	  break;
	}
      }
      
      bool fl_ClCandB  = false;
      if((seg+1) != (NofSeg -1)){
	if(3 > (cont.at(seg+2)->PairId() - CurrentPair)){fl_ClCandB = true;}
      }

      if(!fl_ClCandB){
	// there is no more candidates
	if(cluster->calculate()){ ClusterCont.push_back(cluster); }
	continue;
      }

      // Start Search HitC
      int    NofHitC = cont.at(seg+2)->GetNumOfHit();
      bool   fl_HitC = false;
      double cmtC    = -1;
      for(int mhitC = 0; mhitC<NofHitC; ++mhitC){
	if(cont.at(seg+2)->Joined(mhitC)){continue;}

	FiberHit* HitC = cont.at(seg+2);
	
	cmtC = (double)HitC->GetCTime(mhitC);
	if(true
	   && fabs(cmtC-cmtA)<maxTimeDif
	   && !(fl_HitB && (fabs(cmtC-cmtB)>maxTimeDif))
	   ){
	  cluster->push_back(new FLHit(HitC, mhitC));
	  CurrentPair = HitC->PairId();
	  fl_HitC = true;
	  break;
	}
      }

      bool fl_ClCandC  = false;
      if((seg+2) != (NofSeg -1)){
	if(3 > (cont.at(seg+3)->PairId() - CurrentPair)){fl_ClCandC = true;}
      }

      if(!fl_ClCandC){
	// there is no more candidates
	if(cluster->calculate()){ ClusterCont.push_back(cluster); }
	continue;
      }

      // Start Search HitD
      int    NofHitD = cont.at(seg+3)->GetNumOfHit();
      double cmtD    = -1;
      for(int mhitD = 0; mhitD<NofHitD; ++mhitD){
	if(cont.at(seg+3)->Joined(mhitD)){continue;}

	FiberHit* HitD = cont.at(seg+3);
	
	cmtD = (double)HitD->GetCTime(mhitD);
	if(true
	   && fabs(cmtD-cmtA)<maxTimeDif
	   && !(fl_HitB && (fabs(cmtD-cmtB)>maxTimeDif))
	   && !(fl_HitC && (fabs(cmtD-cmtC)>maxTimeDif))
	   ){
	  cluster->push_back(new FLHit(HitD, mhitD));
	  break;
	}
      }
      
      // Finish
      if(cluster->calculate()){ ClusterCont.push_back(cluster); }
    }
  }

  return ClusterCont.size();
}

// bool HodoAnalyzer::ReCalcGCHits( bool applyRecursively )
// {
//   int n=GCCont.size();
//   for( int i=0; i<n; ++i ){
//     Hodo1Hit *hit=GCCont[i];
//     if(hit) hit->ReCalc(applyRecursively);
//   }
//   return true;
// }

bool HodoAnalyzer::ReCalcBH1Hits( bool applyRecursively )
{
  int n=BH1Cont.size();
  for( int i=0; i<n; ++i ){
    Hodo2Hit *hit=BH1Cont[i];
    if(hit) hit->ReCalc(applyRecursively);
  }
  return true;
}

bool HodoAnalyzer::ReCalcBH2Hits( bool applyRecursively )
{
  int n=BH2Cont.size();
  for( int i=0; i<n; ++i ){
    BH2Hit *hit=BH2Cont[i];
    if(hit) hit->ReCalc(applyRecursively);
  }
  return true;
}

bool HodoAnalyzer::ReCalcTOFHits( bool applyRecursively )
{
  int n=TOFCont.size();
  for( int i=0; i<n; ++i ){
    Hodo2Hit *hit=TOFCont[i];
    if(hit) hit->ReCalc(applyRecursively);
  }
  return true;

}

// bool HodoAnalyzer::ReCalcACHits( bool applyRecursively )
// {
//   int n=ACCont.size();
//   for( int i=0; i<n; ++i ){
//     Hodo1Hit *hit=ACCont[i];
//     if(hit) hit->ReCalc(applyRecursively);
//   }
//   return true;
// }

bool HodoAnalyzer::ReCalcLCHits( bool applyRecursively )
{
  int n=LCCont.size();
  for( int i=0; i<n; ++i ){
    Hodo2Hit *hit=LCCont[i];
    if(hit) hit->ReCalc(applyRecursively);
  }
  return true;
}

bool HodoAnalyzer::ReCalcBMWHits( bool applyRecursively )
{
  int n=BMWCont.size();
  for( int i=0; i<n; ++i ){
    Hodo1Hit *hit=BMWCont[i];
    if(hit) hit->ReCalc(applyRecursively);
  }
  return true;
}
bool HodoAnalyzer::ReCalcBACHits( bool applyRecursively )
{
  int n=BACCont.size();
  for( int i=0; i<n; ++i ){
    Hodo1Hit *hit=BACCont[i];
    if(hit) hit->ReCalc(applyRecursively);
  }
  return true;
}
bool HodoAnalyzer::ReCalcSFVHits( bool applyRecursively )
{
  int n=SFVCont.size();
  for( int i=0; i<n; ++i ){
    Hodo1Hit *hit=SFVCont[i];
    if(hit) hit->ReCalc(applyRecursively);
  }
  return true;
}
bool HodoAnalyzer::ReCalcTOFMTHits( bool applyRecursively )
{
  int n=TOFMTCont.size();
  for( int i=0; i<n; ++i ){
    Hodo1Hit *hit=TOFMTCont[i];
    if(hit) hit->ReCalc(applyRecursively);
  }
  return true;
}


bool HodoAnalyzer::ReCalcBH1Clusters( bool applyRecursively )
{
  int n=BH1ClCont.size();
  for( int i=0; i<n; ++i ){
    HodoCluster *cl=BH1ClCont[i];
    if(cl) cl->ReCalc(applyRecursively);
  }
  return true;
}

bool HodoAnalyzer::ReCalcBH2Clusters( bool applyRecursively )
{
  int n=BH2ClCont.size();
  for( int i=0; i<n; ++i ){
    BH2Cluster *cl=BH2ClCont[i];
    if(cl) cl->ReCalc(applyRecursively);
  }
  return true;
}

bool HodoAnalyzer::ReCalcTOFClusters( bool applyRecursively )
{
  int n=TOFClCont.size();
  for( int i=0; i<n; ++i ){
    HodoCluster *cl=TOFClCont[i];
    if(cl) cl->ReCalc(applyRecursively);
  }
  return true;
}

bool HodoAnalyzer::ReCalcLCClusters( bool applyRecursively )
{
  int n=LCClCont.size();
  for( int i=0; i<n; ++i ){
    HodoCluster *cl=LCClCont[i];
    if(cl) cl->ReCalc(applyRecursively);
  }
  return true;
}

bool HodoAnalyzer::ReCalcAll( void )
{
  //ReCalcGCHits();
  ReCalcBH1Hits();
  ReCalcBH2Hits();
  ReCalcTOFHits();
  //ReCalcACHits();
  ReCalcLCHits();

  ReCalcBMWHits();
  ReCalcBACHits();
  ReCalcSFVHits();
  ReCalcTOFMTHits();

  ReCalcBH1Clusters();
  ReCalcBH2Clusters();
  ReCalcTOFClusters();
  ReCalcLCClusters();

  return true;
}

// Time cut for the BH1 ClusterContainer
void HodoAnalyzer::TimeCutBH1(double tmin, double tmax){
  TimeCut(BH1ClCont, tmin, tmax);
}

// Time cut for the BH2 ClusterContainer
void HodoAnalyzer::TimeCutBH2(double tmin, double tmax){
  TimeCut(BH2ClCont, tmin, tmax);
}

// Time cut for the BFTClusterContainer
void HodoAnalyzer::TimeCutBFT(double tmin, double tmax){
  TimeCut(BFTClCont, tmin, tmax);
}

// Time cut for the SFTClusterContainer
// void HodoAnalyzer::TimeCutSFT(int layer, double tmin, double tmax){
//   TimeCut(SFTClCont[layer], tmin, tmax);
// }

//Implementation of Time cut for the cluster container
template <typename TypeCluster>
void HodoAnalyzer::TimeCut(std::vector<TypeCluster>& cont,
			   double tmin, double tmax)
{
  std::vector<TypeCluster> DeleteCand;
  std::vector<TypeCluster> ValidCand;
  int NofCl = cont.size();
  for(int i = 0; i<NofCl; ++i){
    double ctime = cont.at(i)->CMeanTime();
    if(tmin < ctime && ctime < tmax){
      ValidCand.push_back(cont.at(i));
    }else{
      DeleteCand.push_back(cont.at(i));
    }
  }

  for_each(DeleteCand.begin(), DeleteCand.end(), DeleteObject());
  DeleteCand.clear();
  
  cont.clear();
  cont.resize(ValidCand.size());
  std::copy(ValidCand.begin(), ValidCand.end(), cont.begin());
  ValidCand.clear();
}
