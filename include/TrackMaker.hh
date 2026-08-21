/*
  TrackMaker.hh

  2012/1/24
*/

// -*- C++ -*-

#ifndef TrackMaker_h
#define TrackMaker_h

#include <string>
#include <vector>
#include <deque>


class DCPairHitCluster;
class DCLocalTrack;

#define KeepOrderAsMakeIndex

class TrackMaker
{

public:
  typedef std::vector<DCPairHitCluster*> ClusterList;
#ifdef KeepOrderAsMakeIndex
  typedef std::vector<ClusterList>::const_reverse_iterator  PlaneIterator;
#else
  typedef std::vector<ClusterList>::const_iterator          PlaneIterator;
#endif  

protected:
  const std::vector<ClusterList>& m_candidates;
  PlaneIterator                 m_itr;
  std::deque<DCPairHitCluster*> m_trackCandidate;
  std::vector<DCLocalTrack*>*   m_trackList;
  int    m_nCombi;
  int    m_nValidCombi;
  double m_maxChiSquare;
  int    m_maxCombi;
  int    m_minNumOfHits;

public:
  TrackMaker(const std::vector<ClusterList>& candCont,
	     int minNumOfHits,
	     int maxCombi,
	     double maxChiSquare);
  virtual ~TrackMaker();

  inline double GetMaxChiSquare() const;
  inline int    GetMaxCombi() const;
  inline int    GetMinNumOfHits() const;
  inline int    GetNumOfCombi() const;
  inline int    GetNumOfValidCombi() const;
  virtual void  MakeTracks(std::vector<DCLocalTrack*>& trackList);
  inline void   SetMaxChiSquare(double maxChiSquare);
  inline void   SetMaxCombi(int maxCombi);
  inline void   SetMinNumOfHits(int minNumOfHits);

protected:
  virtual bool IsGood(const DCLocalTrack* track) const;

private:
  TrackMaker(const TrackMaker&);
  TrackMaker& operator=(const TrackMaker&);

  void MakeTrack();
  void Next();

};

//______________________________________________________________________________
// inline accessors of TrackMaker
//______________________________________________________________________________
inline
double
TrackMaker::GetMaxChiSquare() const
{
  return m_maxChiSquare;
}

//______________________________________________________________________________
inline
int
TrackMaker::GetMaxCombi() const
{
  return m_maxCombi;
}

//______________________________________________________________________________
inline
int
TrackMaker::GetMinNumOfHits() const
{
  return m_minNumOfHits;
}

//______________________________________________________________________________
inline
int
TrackMaker::GetNumOfCombi() const
{
  return m_nCombi;
}

//______________________________________________________________________________
inline
int
TrackMaker::GetNumOfValidCombi() const
{
  return m_nValidCombi;
}

//______________________________________________________________________________
inline
void 
TrackMaker::SetMaxChiSquare(double maxChiSquare)
{
  m_maxChiSquare = maxChiSquare;
  return;
}

//______________________________________________________________________________
inline
void 
TrackMaker::SetMaxCombi(int maxCombi)
{
  m_maxCombi = maxCombi;
  return;
}

//______________________________________________________________________________
inline
void
TrackMaker::SetMinNumOfHits(int minNumOfHits)
{
  m_minNumOfHits = minNumOfHits;
  return;
}

#endif
