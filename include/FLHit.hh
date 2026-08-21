//----------------------------------------------------------------------------
// FLHit
//  This class corresponds to the DCLTrackHit
//----------------------------------------------------------------------------

#ifndef FLHit_H
#define FLHit_H

#include"FiberHit.hh"

class FLHit{
private:
  FiberHit* Hit_;
  int       nth_Hit_;

public:
  FLHit(FiberHit* ptr, int index):
    Hit_(ptr), nth_Hit_(index)
  {
    Hit_->RegisterHits(this);
    Hit_->SetJoined(index);
  }

  // leading and trailing are raw value (not 0)
  double GetLeading() const {return Hit_->GetLeading(nth_Hit_);}
  double GetTrailing() const {return Hit_->GetTrailing(nth_Hit_);}

  // t_ is time for PID (should be 0 for pion)
  double GetTime() const {return Hit_->GetTime(nth_Hit_);}

  // ct_ is corrected time after PHC
  double GetCTime() const {return Hit_->GetCTime(nth_Hit_);}

  // In this class, a_ is width = leading - trailing
  // If there is no trailing data, -1 is inserted into width
  double GetWidth() const {return Hit_->GetWidth(nth_Hit_);}

  // Return local position.
  // Note : BFT is a position detector, there is no z information.
  //        Z position should be written in the geometry file.
  double GetPosition() const {return Hit_->GetPosition();}

  // Treat BFT as 1 dimentional detector
  int    PairId() const {return Hit_->PairId();}

  // Join in Cluster
  bool   Joined() const {return Hit_->Joined(nth_Hit_);}

  // Dump
  void   Dump(){Debug();}

  friend class FiberHit;

private:
  ~FLHit(){;};

  FLHit();
  FLHit(const FLHit& object);
  FLHit& operator =(const FLHit& object);

  void Debug(){
    std::cout << "plid " << Hit_->PairId() << " ";
    std::cout << "Pos " << GetPosition() << " ";
    std::cout << "Time " << GetCTime() << std::endl;
  }
};

#endif
