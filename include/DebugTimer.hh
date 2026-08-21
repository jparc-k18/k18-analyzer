/*
  DebugTimer.hh

  2012/1/24
*/

// -*- C++ -*-

#ifndef DebugTimer_h
#define DebugTimer_h

#include <string>

namespace debug
{
  class Timer
  {

  public:
    struct ::timespec;

  private:
    ::timespec* m_start;
    ::timespec* m_stop;
    std::string m_msg;

  public:
    explicit Timer(const std::string& msg="");
    ~Timer();

    void stop();
    void print() const;

  private:
    Timer(const Timer&);
    Timer& operator=(const Timer&);

  };

}
#endif
