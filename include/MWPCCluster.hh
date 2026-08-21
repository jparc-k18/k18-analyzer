/*
  MWPCCluster.hh

  2012/1/24
*/

#ifndef MWPCCluster_h
#define MWPCCluster_h

#include <string>
#include <vector>

class DCHit;

#ifdef MemoryLeak
#include "DebugCounter.hh"
#endif

class MWPCCluster
{

public:
  struct Statistics
  {
    double m_wire;
    double m_wpos;
    double m_leading;
    double m_trailing;
    double m_length;
    double m_totalLength;
    int    m_clusterSize;

    Statistics();
    ~Statistics();

    void Print(const std::string& arg="") const;
  };

private:
  std::vector<DCHit*> m_hits;
  Statistics m_mean;
  Statistics m_first;
  bool       m_status;

#ifdef MemoryLeak
  static debug::Counter sm_counter;
#endif

public:
  MWPCCluster();
  ~MWPCCluster();
  
  void   Add(DCHit* h);
  void   Calculate();
  int    GetClusterSize() const;
  const Statistics& GetFirst() const;
  const std::vector<DCHit*>& GetHits() const;
  const Statistics& GetMean() const;
  double GetMeanTime() const;
  double GetMeanWire() const;
  double GetMeanWirePos() const;
  int    GetNumOfHits() const;
  bool   IsGoodForAnalysis() const;
  void   Print(const std::string& arg="") const;
  void   SetStatus(bool status);

private:
  MWPCCluster(const MWPCCluster&);
  MWPCCluster& operator=(const MWPCCluster&);

};

#endif
