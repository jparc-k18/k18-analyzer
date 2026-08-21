//-----------------------------------------------------------------------------
// Class for BFT & SFT cluster
//-----------------------------------------------------------------------------

#ifndef FIBERCLUSTER
#define FIBERCLUSTER

#include<vector>

class FLHit;

class FiberCluster{
protected:
  typedef std::vector<FLHit*> HitContainer;
  
  HitContainer HitAssemblage_;
  int          csize_;
  double       MeanTime_;
  double       bestWidth_;
  double       MeanPos_;
  enum         FlagsFiber{Initialized, gfastatus,
			  sizeFlagsFiber};
  bool         flag_[sizeFlagsFiber];

public:
  bool      calculate();
  void      push_back(FLHit* ptr){ HitAssemblage_.push_back(ptr); };

  // Return # of used hit to create cluster
  int       ClusterSize() const {return csize_;};
  
  // Return most closest value to 0
  double    CMeanTime() const {return MeanTime_;};

  // Return best width
  double    Width() const {return bestWidth_;};

  // Return mean local position
  double    MeanPosition() const {return MeanPos_;};
  
  FLHit*    GetHit(int i) const;
  bool      GoodForAnalysis() const {return flag_[gfastatus];};
  bool      GoodForAnalysis(bool status)
  {bool ret = flag_[gfastatus]; flag_[gfastatus]=status; return ret;};
  
  bool      ReCalc( bool applyRecusively=false );

  FiberCluster();
  virtual ~FiberCluster();

private:
  void calc_wo_width();
  void calc_w_width();

  void Initializer();
  void CleanUp();
  void Debug();

  FiberCluster(const FiberCluster& object);
  FiberCluster& operator =(const FiberCluster& object);

};

#endif
