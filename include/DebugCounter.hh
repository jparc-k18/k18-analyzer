/*
  DebugCounter.hh

  2012/1/24
*/

// -*- C++ -*-
#ifndef DebugCounter_h
#define DebugCounter_h

#include <string>

namespace debug
{
  class Counter
  {

  private:
    int m_count;
    std::string m_name;
    bool m_verbose;

  public:
    explicit Counter(const std::string& name);
    ~Counter();

    void operator++();
    void operator--();
    void print(const std::string& arg="") const;

  private:
    Counter(const Counter&);
    Counter& operator=(const Counter&);

  };
}

#endif
