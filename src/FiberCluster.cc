#include"FiberHit.hh"
#include"FLHit.hh"
#include"FiberCluster.hh"

#include<string>
#include<math.h>

#define DEBUG 0

static const std::string MyName = "FiberCluster";

// Constructor
FiberCluster::FiberCluster()
  : csize_(0)
{
  Initializer();
}

// Destructor
FiberCluster::~FiberCluster(){
  CleanUp();
}

// Initializer
void FiberCluster::Initializer(){
  for(int i = 0; i<sizeFlagsFiber; ++i){
    flag_[i] = false;
  }
}

// CleanUp
void FiberCluster::CleanUp(){
  HitAssemblage_.clear();
}

// GetHit
FLHit *FiberCluster::GetHit( int i ) const
{
  return HitAssemblage_.at(i);
}

// Calculate
bool FiberCluster::calculate(){
  static const std::string MyFunc = "::calculate ";
  if(flag_[Initialized]){
    std::cerr << "#E " << MyName << MyFunc
	      << "Already initialied" << std::endl;
    return false;
  }

  if(0 == (csize_ = HitAssemblage_.size())){
    std::cerr << "#E " << MyName << MyFunc
	      << "No FiberHit in local container" << std::endl;
    return false;
  }

  int width = (int)HitAssemblage_.at(0)->GetWidth();
  if(-1 == width){
    calc_wo_width();
  }else{
    calc_w_width();
  }

  flag_[Initialized] = true;

#if DEBUG
  Debug();
#endif

  return true;
}

// calc_wo_width
// Time closest to 0 is selected as mean time
void FiberCluster::calc_wo_width(){
  double mx = 0., bt = 0.;
  
  mx = HitAssemblage_.at(0)->GetPosition();
  bt = HitAssemblage_.at(0)->GetCTime();

  for(int i = 1; i<csize_; ++i){
    mx += HitAssemblage_.at(i)->GetPosition();
    double t = HitAssemblage_.at(i)->GetCTime();
    if(fabs(t) < fabs(bt)){
      bt = t;
    }    
  }

  mx /= double(csize_);
  
  MeanTime_=bt; bestWidth_=-1; MeanPos_=mx;
}

// calc_w_width
// Time closest to 0 is selected as mean time
void FiberCluster::calc_w_width(){
  double mx = 0., mt = 0., bw = 0.;

  mx = HitAssemblage_.at(0)->GetPosition();
  mt = HitAssemblage_.at(0)->GetCTime();
  bw = HitAssemblage_.at(0)->GetWidth();
  
  for(int i = 1; i<csize_; ++i){
    mx += HitAssemblage_.at(i)->GetPosition();
    double t = HitAssemblage_.at(i)->GetCTime();
    if(fabs(t) < fabs(mt)){
      mt = t;
    }    

    if(bw < HitAssemblage_.at(i)->GetWidth()){
      bw = HitAssemblage_.at(i)->GetWidth();
    }
  }

  mx /= double(csize_);
  
  MeanTime_=mt; bestWidth_=bw; MeanPos_=mx;
}

// Debug
void FiberCluster::Debug(){
  std::cout << "Used hit\n";
  for(int i = 0; i<csize_; ++i){
    HitAssemblage_.at(i)->Dump();
  }
}
