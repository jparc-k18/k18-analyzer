/*
  BH1Filter.hh

  2012/1/24
*/

// -*- C++ -*-

#ifndef BH1Filter_h
#define BH1Filter_h

#include <string>
#include <vector>
#include <set>

class DCAnalyzer;
class HodoAnalyzer;

class BH1Filter
{

public:
  struct Param
  {
    // [plane]
    std::vector<double> m_xmin;
    std::vector<double> m_xmax;

    Param();
    ~Param();
    void Print(const std::string& arg="") const;
  };

  enum EParam
    {
      kBH1Segment,
      kLayerID,
      kXMin,
      kXMax,
      kNParam
    };


private:
  // [segment]
  std::vector<Param>  m_param;
  const DCAnalyzer*   m_dc;
  const HodoAnalyzer* m_hodo;

  void   BuildCandidates(std::set<int> &seg,
			 std::vector<std::vector<std::vector<DCHit*> > >& candidates);
  
public:
  static BH1Filter& GetInstance();
  ~BH1Filter();

  void   Apply(const HodoAnalyzer& hodo,
	       const DCAnalyzer& dc,
	       std::vector<std::vector<std::vector<DCHit*> > >& candidates);
  void   Apply_withT0Seg(int T0Seg,
			    const DCAnalyzer& dc,
			    std::vector<std::vector<std::vector<DCHit*> > >& candidates);
  const std::vector<double>& GetXmax(int seg) const;
  const std::vector<double>& GetXmin(int seg) const;
  void   Initialize(const std::string& filename);
  void   Print(const std::string& arg="") const;

private:
  BH1Filter();
  BH1Filter(const BH1Filter&);
  BH1Filter& operator=(const BH1Filter&);

};

#endif
