#include "VEvent.hh"

#include <iostream>
#include <sstream>
#include <cmath>
#include <string>

#include "UnpackerManager.hh"
#include "ConfMan.hh"
#include "HistHelper.hh"
#include "DetectorID.hh"

using namespace hddaq::unpacker;
using namespace hddaq;

 VEvent::VEvent()
 {
 }

 VEvent::~VEvent()
 {
 }

class EventScaler
  : public VEvent
{

private:

public:
  EventScaler();
  ~EventScaler();

  bool ProcessingBegin();
  bool ProcessingEnd();
  bool ProcessingNormal();
  bool InitializeHistograms();
  void InitializeEvent();
};

struct Event{
  enum name1{BH1OR, BH2OR, sizename1};
  enum name2{BH2Seg1, BH2Seg2, BH2Seg3, BH2Seg4,
	     BH2Seg5, BH2Seg6, BH2Seg7,
	     BH1Seg1, BH1Seg2, BH1Seg3, BH1Seg4, BH1Seg5, BH1Seg6,
	     BH1Seg7, BH1Seg8, BH1Seg9, BH1Seg10, BH1Seg11,
	     sizename2};

  int scaler1[sizename1];
  int scaler2[sizename2];
};
static Event event;

VEvent *ConfMan::EventAllocator()
{
  return new EventScaler;
}

EventScaler::EventScaler()
  : VEvent()
{
}

EventScaler::~EventScaler()
{
}

bool EventScaler::ProcessingBegin()
{
 return true;
}

bool EventScaler::ProcessingNormal()
{
  static UnpackerManager& g_unpacker = GUnpacker::get_instance();
  static const int index1[] = {0, 1};
  static const int index2[] = { 2, 3, 4, 5, 6, 7, 8,
			        17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27
  };

  std::string ch_name_v820[2];
  ch_name_v820[0]="BH1OR   :";//1
  ch_name_v820[1]="BH2OR   :";//2

  std::string ch_name_v830[18];
  ch_name_v830[0]="BH2Seg1   :";//1
  ch_name_v830[1]="BH2Seg2   :";//2
  ch_name_v830[2]="BH2Seg3   :";//3
  ch_name_v830[3]="BH2Seg4   :";//4
  ch_name_v830[4]="BH2Seg5   :";//5
  ch_name_v830[5]="BH2Seg6   :";//6
  ch_name_v830[6]="BH2Seg7   :";//7
  ch_name_v830[7]="BH1Seg1   :";//16
  ch_name_v830[8]="BH1Seg2   :";//17
  ch_name_v830[9]="BH1Seg3   :";//18
  ch_name_v830[10]="BH1Seg4   :";//19
  ch_name_v830[11]="BH1Seg5   :";//20
  ch_name_v830[12]="BH1Seg6   :";//21
  ch_name_v830[13]="BH1Seg7   :";//22
  ch_name_v830[14]="BH1Seg8   :";//23
  ch_name_v830[15]="BH1Seg9   :";//24
  ch_name_v830[16]="BH1Seg10  :";//25
  ch_name_v830[17]="BH1Seg11  :";//26
  
  //check
  static int prev = 0;
  int current     = 0;
  if(0 != g_unpacker.get_entries(210, 0, 0, 0, 0)){
    current = g_unpacker.get(210, 0, 0, 0, 0);
  }
  
  if(prev > current){
    for(int i = 0; i<Event::sizename1; ++i){
      std::cout << ch_name_v820[i] << " "
		<< event.scaler1[i] << std::endl;
    }
    // sumup
    int validBH2OR = 0;
    int validBH1OR = 0;
    for(int i = (int)Event::BH2Seg1; i<(int)Event::BH2Seg7+1; ++i){
      validBH2OR += event.scaler2[i];
    }

    for(int i = (int)Event::BH1Seg1; i<(int)Event::BH1Seg11+1; ++i){
      validBH1OR += event.scaler2[i];
    }
    
    std::cout << "ValidBH1OR " << validBH1OR << std::endl;
    std::cout << "ValidBH2OR " << validBH2OR << std::endl;
    std::cout << std::endl;

    for(int i = 0; i<Event::sizename2; ++i){
      std::cout << ch_name_v830[i] << " "
		<< event.scaler2[i] << std::endl;
    } 

    std::cout << "#D End of spill" << std::endl;
    std::cout << std::endl;
  }

  for(int i = 0; i<Event::sizename1; ++i){
    int nh = g_unpacker.get_entries(210, 0, 0, index1[i], 0);
    if(nh == 0){continue;};
    event.scaler1[i] = g_unpacker.get(210, 0, 0, index1[i], 0);
  }

  for(int i = 0; i<Event::sizename2; ++i){
    int nh = g_unpacker.get_entries(210, 1, 0, index2[i], 0);
    if(nh == 0){continue;};
    event.scaler2[i] = g_unpacker.get(210, 1, 0, index2[i], 0);
  }

  prev = current;

  return true;
}

void EventScaler::InitializeEvent( void )
{

}

bool EventScaler::ProcessingEnd()
{
  return true;
}

bool ConfMan:: InitializeHistograms()
{
  return true;
}

bool ConfMan::InitializeParameterFiles( void )
{

  return true;
}
