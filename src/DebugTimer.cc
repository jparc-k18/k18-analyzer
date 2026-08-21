/*
  DebugTimer.hh

  2012/1/24
*/

// -*- C++ -*-

#include <ctime>
#include <iostream>

#include "DebugTimer.hh"

namespace debug
{
//______________________________________________________________________________
Timer::Timer(const std::string& msg)
  : m_start(new ::timespec),
    m_stop(0),
    m_msg(msg)
{
  ::clock_gettime(CLOCK_REALTIME, m_start);
}

//______________________________________________________________________________
Timer::~Timer()
{
  if (!m_stop)
    {
      stop();
      print();
    }
  
  if (m_stop)
    {
      delete m_stop;
      m_stop = 0;
    }
  if (m_start)
    {
      delete m_start;
      m_start = 0;
    }
}

//______________________________________________________________________________
void
Timer::print() const
{
  const double sec  = m_stop->tv_sec  - m_start->tv_sec;
  const double nsec = m_stop->tv_nsec - m_start->tv_nsec;
  std::cout << "#DTimer " << m_msg
	    << " : " << (sec*.1e9 + nsec) << " nsec ("
	    << (sec + nsec*1.e-9) << " sec)" << std::endl;
  return;
}

//______________________________________________________________________________
void
Timer::stop()
{
  m_stop = new ::timespec;
  ::clock_gettime(CLOCK_REALTIME, m_stop);
  return;
}

}
