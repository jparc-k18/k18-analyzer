/*
  TrackMaker.cc

  2012/1/24
*/

// -*- C++ -*-

#include <cstdlib>
#include <iostream>

#include "DCPairHitCluster.hh"
#include "DCLocalTrack.hh"
#include "DCLTrackHit.hh"
#include "TrackMaker.hh"

const std::string ks_className("TrackMaker");

//______________________________________________________________________________
// template functions to resolve iterator-type automatically at compilation
//______________________________________________________________________________
template <typename IteratorType>
inline
void
initIterator(const std::vector<TrackMaker::ClusterList>& container,
	     IteratorType& itr)
{
  itr = container.begin();
  return;
}

//______________________________________________________________________________
template <>
inline
void
initIterator(const std::vector<TrackMaker::ClusterList>& container,
	     std::vector<TrackMaker::ClusterList>::const_reverse_iterator& itr)
{
  itr = container.rbegin();
  return;
}

//______________________________________________________________________________
template <typename IteratorType>
inline
bool
isEnd(const std::vector<TrackMaker::ClusterList>& container,
      IteratorType& itr)
{
  return itr==container.end();
}

//______________________________________________________________________________
template <>
inline
bool
isEnd(const std::vector<TrackMaker::ClusterList>& container,
      std::vector<TrackMaker::ClusterList>::const_reverse_iterator& itr)
{
  return itr==container.rend();
}

//______________________________________________________________________________
template <typename IteratorType>
inline
void
addCluster(const IteratorType& itr,
	   const DCPairHitCluster* cluster,
	   std::deque<DCPairHitCluster*>& container)
{
  container.push_back(const_cast<DCPairHitCluster*>(cluster));
  return;
}

//______________________________________________________________________________
template <>
inline
void
addCluster(const std::vector<TrackMaker::ClusterList>::const_reverse_iterator& itr,
	   const DCPairHitCluster* cluster,
	   std::deque<DCPairHitCluster*>& container)
{
  container.push_front(const_cast<DCPairHitCluster*>(cluster));
  return;
}

//______________________________________________________________________________
template <typename IteratorType>
inline
void
removeCluster(const IteratorType& itr,
	      std::deque<DCPairHitCluster*>& container)
{
  container.pop_back();
  return;
}

//______________________________________________________________________________
template <>
inline
void
removeCluster(const std::vector<TrackMaker::ClusterList>::const_reverse_iterator& itr,
	      std::deque<DCPairHitCluster*>& container)
{
  container.pop_front();
  return;
}

//______________________________________________________________________________
// class TrackMaker
//______________________________________________________________________________
TrackMaker::TrackMaker(const std::vector<ClusterList>& candCont,
		       int minNumOfHits,
		       int maxCombi,
		       double maxChiSquare)
  : m_candidates(candCont),
    m_itr(),
    m_trackCandidate(0),
    m_trackList(0),
    m_nCombi(0),
    m_nValidCombi(0),
    m_maxChiSquare(maxChiSquare),
    m_maxCombi(maxCombi),
    m_minNumOfHits(minNumOfHits)
{
}

//______________________________________________________________________________
TrackMaker::~TrackMaker()
{
}

//______________________________________________________________________________
bool
TrackMaker::IsGood(const DCLocalTrack* track) const
{
  return 
    (static_cast<int>(track->GetNHit())>=m_minNumOfHits)
    && 
    (const_cast<DCLocalTrack*>(track)->DoFit())
    &&
    (track->GetChiSquare()<m_maxChiSquare);
}

//______________________________________________________________________________
void
TrackMaker::MakeTrack()
{
  static const std::string funcname("["+ks_className+"::"+__func__+"]");
  DCLocalTrack* track = 0;
  int nhit = 0;
  
  if (!m_trackList)
    ++m_nCombi;
  else
    track = new DCLocalTrack;

  for (std::deque<DCPairHitCluster*>::iterator 
	 itr    = m_trackCandidate.begin(),
	 itrEnd = m_trackCandidate.end();
       itr!=itrEnd; ++itr)
    {
      const DCPairHitCluster* cluster = *itr;
      if (!cluster)
	continue;
      for (int i=0, n=cluster->NumberOfHits(); i<n; ++i)
	{
	  DCLTrackHit* hit = cluster->GetHit(i);
	  if (!hit)
	    continue;

	  if (!m_trackList)
	    ++nhit;
	  else if (track)
	    track->AddHit(hit);
	}
    }

  if (!m_trackList)
    {
      if (nhit>=m_minNumOfHits)
	++m_nValidCombi;
      return;
    }

  if (!IsGood(track))
    {
      delete track;
      track = 0;
      return;
    }

  if (track)
    m_trackList->push_back(track);
  return;
}

//______________________________________________________________________________
void
TrackMaker::MakeTracks(std::vector<DCLocalTrack*>& trackList)
{
  static const std::string funcname("["+ks_className+"::"+__func__+"]");
  
  // estimate the number of tracks, not creating tracks
  m_trackList = 0;
  initIterator(m_candidates, m_itr);
  Next();

//   std::cout << "#D " << funcname <<  " n combi = " << m_nValidCombi 
// 	    << std::endl;

  if (m_nValidCombi>m_maxCombi)
    {
      std::cout << "#W " << funcname << " : too many combinations."
		<< " n valid combi = " << m_nValidCombi << std::endl;
      return;
    }

  // create new tracks
  m_trackList = &(trackList);
  initIterator(m_candidates, m_itr);
  Next();
  return;
}

//______________________________________________________________________________
void
TrackMaker::Next()
{
  // decrement iterator on every return from this function

  if (isEnd(m_candidates, m_itr))
    {
      MakeTrack();

      --m_itr;
      return;
    }

  addCluster(m_itr, 0, m_trackCandidate);
  ++m_itr;
  Next();
  removeCluster(m_itr, m_trackCandidate);

  const ClusterList& clusters = *m_itr;
  for (ClusterList::const_iterator i=clusters.begin(), iEnd=clusters.end(); 
       i!=iEnd; ++i)
    {
      addCluster(m_itr, *i, m_trackCandidate);
      ++m_itr;
      Next();
      removeCluster(m_itr, m_trackCandidate);
    }

  --m_itr;
  return;
}
