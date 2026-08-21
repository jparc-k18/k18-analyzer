#include"FiberHit.hh"
#include"FLHit.hh"

#include"ConfMan.hh"
#include"HodoParamMan.hh"
#include"HodoPHCMan.hh"
#include"DCGeomMan.hh"
#include"RawData.hh"

#include<iostream>
#include<iomanip>
#include<string>
#include<cstring>
#include<sstream>

static const std::string MyName = "FiberHit";

// Constructor
FiberHit::FiberHit(HodoRawHit *object, const char* name):
  Hodo1Hit(object),
  DetectorName_(name),
  position_(-999.), offset_(0),
  pair_id_(0),
  Status_(false)
{

}

// Destructor
FiberHit::~FiberHit(){
  CleanUp();
}

// CleanUp
void FiberHit::CleanUp(){
  int NofPtr = PtrCont_.size();
  for(int i = 0; i<NofPtr; ++i){
    delete PtrCont_.at(i);
    PtrCont_.at(i) = NULL;
  }
  PtrCont_.clear();
}

// calcurate
bool FiberHit::calculate(){
  static const std::string MyFunc = "::calculate ";

  if(Status_){
    std::cerr << "#E " << MyName << MyFunc
	      << "Already calculated" << std::endl;
    return false;
  }

  if("" == DetectorName_){
    std::cerr << "#E " << MyName << MyFunc
	      << "DetectorName is NULL" << std::endl;
    return Status_;
  }

  // Detector information
  int cid  = raw_->DetectorId();
  int plid = raw_->PlaneId();
  int seg  = raw_->SegmentId();

  // Geometry calibration
  if("BFT" == DetectorName_ || "SFT_X" == DetectorName_){
    // case of BFT and SFT X layers
    // They have up and down planes in 1 layer.
    // We treat these 2 planes as 1 dimentional plane.
    if(1 == raw_->PlaneId()){
      // case of down plane
      offset_ = 0.5;
      pair_id_  = 1;
    }
    pair_id_ += 2*raw_->SegmentId();
  }else{
    // case of SFT UV layers
    // They have only 1 plane in 1 layer.
    pair_id_ = raw_->SegmentId();
  }
  
  DCGeomMan& geom = DCGeomMan::GetInstance();
  int DetectorId = geom.GetDetectorId(DetectorName_);
  position_      = geom.calcWirePosition(DetectorId, seg);

  ConfMan *confMan=ConfMan::GetConfManager();
  if(!confMan){
    std::cerr << "#E " << MyName << MyFunc
	      << "Cannot get confManager"
	      << std::endl; 
    return Status_;
  }
  HodoParamMan *hodoMan = confMan->GetHodoParamManager();
  HodoPHCMan   *phcMan  = confMan->GetHodoPHCManager();
  if(!hodoMan){
    std::cerr << "#E " << MyName << MyFunc
	      << "Cannot get HodoParamManager"
	      << std::endl; 
    return Status_;
  }

  // hit information
  multi_hit = raw_->SizeTdc1();
  int UorD=0; // BFT is allways U
  int fl_EASIROC = (-1 != raw_->GetTdc2());

  for(int m = 0; m<multi_hit; ++m){
    int leading  = raw_->GetTdc1(m); // leading
    int width    = -1;
    if(fl_EASIROC){
      int trailing = raw_->GetTdc2(m); // trailing
      width = leading - trailing; // width
    }
    a_.push_back(width);

    // Time calibration
    if( leading<0 ){
      std::cerr << "#E " << MyName << MyFunc
		<< "No valid TDC\n" 
		<< "TDC = " << leading
		<< std::endl;
      return Status_;
    }
  
    // time is overwriten
    double time = 0;
    if( !hodoMan->GetTime(cid, plid, seg, UorD, leading, time) ){
      std::cerr << "#E " << MyName << MyFunc
		<< "Something wrong at GetTime"
		<< std::endl; 
      std::cerr << cid << " "
		<< plid << " "
		<< seg << " "
		<< UorD << " "
		<< leading << " "
		<< time << std::endl;
      
      return Status_;
    }
    t_.push_back(time);

    // ctime is overwriten
    double ctime = time;
    if( phcMan ){
      phcMan->doCorrection(cid, plid, seg, UorD, time, width, ctime );
    }
    ct_.push_back(ctime);

    flJoin_.push_back(false);
  }

  return Status_ = true;
}
