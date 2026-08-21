/*
  HodoRawHit.hh

  2012/6/10
*/

#ifndef HodoRawHit_h 
#define HodoRawHit_h

#include <cstddef>
#include <vector>
#include <iostream>

typedef std::vector<int> dType;

class HodoRawHit
{

private:
  int DetId_, PlId_, SegId_;
  dType Adc1_, Adc2_;
  dType Tdc1_, Tdc2_;
  int NhitsTdc_;

public:
  HodoRawHit( int detid, int plid, int segid )
    : DetId_(detid), PlId_(plid), SegId_(segid),
      Adc1_(1, -1), Adc2_(1, -1), Tdc1_(1, -1), Tdc2_(1, -1), NhitsTdc_(0)
  {

  };

  ~HodoRawHit() {};

public:
  void SetAdc1(int adc) { 
    if(-1 == Adc1_.at(0)){
      Adc1_.at(0) = adc;
    }else{
      Adc1_.push_back(adc); 
    };
  }

  void SetAdc2(int adc) {
    if(-1 == Adc2_.at(0)){
      Adc2_.at(0) = adc;
    }else{
      Adc2_.push_back(adc);
    };
  }

  void SetTdc1(int tdc) {
    if(-1 == Tdc1_.at(0)){
      Tdc1_.at(0) = tdc;
      ++NhitsTdc_;
    }else{
      Tdc1_.push_back(tdc);
    }
  };

  void SetTdc2(int tdc) {
    if(-1 == Tdc2_.at(0)){
      Tdc2_.at(0) = tdc;
      ++NhitsTdc_;
    }else{
      Tdc2_.push_back(tdc);
    }
  };

  void SetAdcUp( int adc) { SetAdc1(adc); };
  void SetAdcLeft( int adc) { SetAdc1(adc); };
  void SetAdcDown( int adc) { SetAdc2(adc); };
  void SetAdcRight( int adc) { SetAdc2(adc); };
  void SetTdcUp( int tdc) { SetTdc1(tdc); };
  void SetTdcLeft( int tdc) { SetTdc1(tdc); };
  void SetTdcDown( int tdc) { SetTdc2(tdc); };
  void SetTdcRight( int tdc) { SetTdc2(tdc); };

  int DetectorId( void ) const { return DetId_; };
  int PlaneId( void ) const { return PlId_; };
  int SegmentId( void ) const { return SegId_; };

  int GetAdc1( void ) const { return Adc1_.at(0); };
  int GetAdc2( void ) const { return Adc2_.at(0); };
  int GetTdc1( void ) const { return Tdc1_.at(0); };
  int GetTdc2( void ) const { return Tdc2_.at(0); };

  int GetNumOfTdcHits( void ) const { return NhitsTdc_; };

  int GetAdcUp( void ) const { return GetAdc1(); };
  int GetAdcLeft( void ) const { return GetAdc1(); };
  int GetAdcDown( void ) const { return GetAdc2(); };
  int GetAdcRight( void ) const { return GetAdc2(); };
  int GetTdcUp( void ) const { return GetTdc1(); };
  int GetTdcLeft( void ) const { return GetTdc1(); };
  int GetTdcDown( void ) const { return GetTdc2(); };
  int GetTdcRight( void ) const { return GetTdc2(); };

  // for Multi-hit method
  int GetAdc1( int i ) const { return Adc1_.at(i); };
  int GetAdc2( int i ) const { return Adc2_.at(i); };
  int GetTdc1( int i ) const { return Tdc1_.at(i); };
  int GetTdc2( int i ) const { return Tdc2_.at(i); };

  int GetAdcUp( int i ) const { return GetAdc1(i); };
  int GetAdcLeft( int i ) const { return GetAdc1(i); };
  int GetAdcDown( int i ) const { return GetAdc2(i); };
  int GetAdcRight( int i ) const { return GetAdc2(i); };
  int GetTdcUp( int i ) const { return GetTdc1(i); };
  int GetTdcLeft( int i ) const { return GetTdc1(i); };
  int GetTdcDown( int i ) const { return GetTdc2(i); };
  int GetTdcRight( int i ) const { return GetTdc2(i); };

  int SizeAdc1(void) const {
    if(-1 == Adc1_.at(0)){
      return 0;
    }else{
      return Adc1_.size();
    }
  }

  int SizeAdc2(void) const {
    if(-1 == Adc2_.at(0)){
      return 0;
    }else{
      return Adc2_.size();
    }
  }

  int SizeTdc1(void) const {
    if(-1 == Tdc1_.at(0)){
      return 0;
    }else{
      return Tdc1_.size();
    }
  }

  int SizeTdc2(void) const {
    if(-1 == Tdc2_.at(0)){
      return 0;
    }else{
      return Tdc2_.size();
    }
  }

  int GetSizeAdcUp( void ) const { return SizeAdc1(); };
  int GetSizeAdcLeft( void ) const { return SizeAdc1(); };
  int GetSizeAdcDown( void ) const { return SizeAdc2(); };
  int GetSizeAdcRight( void ) const { return SizeAdc2(); };
  int GetSizeTdcUp( void ) const { return SizeTdc1(); };
  int GetSizeTdcLeft( void ) const { return SizeTdc1(); };
  int GetSizeTdcDown( void ) const { return SizeTdc2(); };
  int GetSizeTdcRight( void ) const { return SizeTdc2(); };

  void clear(){
    NhitsTdc_ = 0;
    Adc1_.clear(); Adc2_.clear(); Tdc1_.clear(); Tdc2_.clear();
    Adc1_.push_back(-1); Adc2_.push_back(-1);
    Tdc1_.push_back(-1); Tdc2_.push_back(-1);
  }
};
#endif
