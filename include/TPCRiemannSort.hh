#ifndef TPC_RIEMANN_SORT_HH
#define TPC_RIEMANN_SORT_HH

class TPCRiemannSort {
public:
  enum Criteria {
    kSortX = 0,
    kSortY = 1,
    kSortReverseY = -1,
    kSortZ = 2,
    kSortReverseZ = -2,
    kSortR = 3,
    kSortReverseR = -3,
    kSortDistance = 4,
    kSortReverseDistance = -4,
    kSortAlpha = 5,
    kSortReverseAlpha = -5,
    kSortLengthInc = 6,
    kSortLengthDec = -6,
    kSortCharge = 7,
    kSortReverseCharge = -7,
    kSortPhi = 8,
    kSortReversePhi = -8,
  };
};

#endif
