/*
  VEvent.hh

  2012/1/24
*/

#ifndef VEVENT_H
#define VEVENT_H

class VEvent
{

public:
  VEvent();
  virtual ~VEvent() = 0;

  virtual bool ProcessingBegin() = 0;
  virtual bool ProcessingEnd() = 0;
  virtual bool ProcessingNormal() = 0;

};

#endif
