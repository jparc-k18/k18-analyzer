/*
  DebugCounter.cc

  2012/1/24
*/

// -*- C++ -*-

#include <iostream>

#include "DebugCounter.hh"

namespace debug
{

//______________________________________________________________________________
Counter::Counter(const std::string& name)
  : m_count(0),
    m_name(name),
    m_verbose(false)
{
#ifdef MemoryLeakVerbose
  m_verbose = true;
#endif
}

//_____________________________________________________________________________
Counter::~Counter()
{
}

//_____________________________________________________________________________
void
Counter::operator++()
{
  if (m_verbose || m_count==0)
    std::cout << "#D(+) " << m_name << " : " << m_count << std::endl;
  ++m_count;
  return;
}

//_____________________________________________________________________________
void
Counter::operator--()
{
  --m_count;
  if (m_verbose || m_count==0)
    std::cout << "#D(-) " << m_name << " : " << m_count << std::endl;
  return;
}

//_____________________________________________________________________________
void
Counter::print(const std::string& arg) const
{
  std::cout << "#D (debug) "
	    << m_name << " : " << m_count << std::endl;
  return;
}

}
