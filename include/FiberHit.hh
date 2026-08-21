//-----------------------------------------------------------------------------
// Class for BFT & SFT
//-----------------------------------------------------------------------------

#ifndef FIBERHIT
#define FIBERHIT

#include"Hodo1Hit.hh"
#include<string>

class FLHit;

class FiberHit : public Hodo1Hit{
protected:
  std::string DetectorName_;
  double      position_;
  double      offset_;
  int         pair_id_;
  bool        Status_;
  std::vector<bool>   flJoin_;
  std::vector<FLHit*> PtrCont_;

public:
  void   SetDetectorName(const char* name){DetectorName_ = name;};
  bool   calculate();

  // leading and trailing are raw value (not 0)
  double GetLeading() const {return raw_->GetTdc1(0);}
  double GetLeading(int n) const {return raw_->GetTdc1(n);}
  double GetTrailing() const {return raw_->GetTdc2(0);}
  double GetTrailing(int n) const {return raw_->GetTdc2(n);}

  // t_ is time for PID (should be 0 for pion)
  double GetTime() const {return Hodo1Hit::GetT();}
  double GetTime(int n) const {return Hodo1Hit::GetT(n);}

  // ct_ is corrected time after PHC
  double GetCTime() const {return Hodo1Hit::GetCT();}
  double GetCTime(int n) const {return Hodo1Hit::GetCT(n);}

  // In this class, a_ is width = leading - trailing
  // If there is no trailing data, -1 is inserted into width
  double GetWidth() const {return double(Hodo1Hit::GetA());}
  double GetWidth(int n) const {return double(Hodo1Hit::GetA(n));}

  // Return local position
  // Note : BFT is a position detector, there is no z information.
  //        Z position should be written in the geometry file.
  double GetPosition() const {return position_ + offset_;}

  // Treat BFT & SFT as 1 dimentional detector
  int    PairId() const {return pair_id_;}

  // make flag_[joined]
  void SetJoined(int m){flJoin_.at(m) = true;}
  bool Joined(int m) const {return flJoin_.at(m);}

  // Memorize the pointer of FLHit
  void RegisterHits(FLHit* ptr){PtrCont_.push_back(ptr);}

  virtual bool ReCalc(bool allpyRecursively=false)
  {return FiberHit::calculate();}

  explicit FiberHit(HodoRawHit *object, const char* name);
  virtual ~FiberHit();

  static bool CompFiberHit(const FiberHit* rLeft,
			   const FiberHit* rRight);

private:
  FiberHit();
  FiberHit(const FiberHit& object);
  FiberHit& operator =(const FiberHit& object);

  void CleanUp();
};

inline bool FiberHit::CompFiberHit(const FiberHit* rLeft,
				   const FiberHit* rRight)
{
  return rLeft->PairId() < rRight->PairId();
}

#endif
