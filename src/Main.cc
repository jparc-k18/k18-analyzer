/*
  Main.cc

  2012/1/24
*/

#include <iostream>

#include <string>
#include <vector>
#include <signal.h>

#include "ConfMan.hh"
#include "UnpackerManager.hh"
#include "VEvent.hh"

#include <TROOT.h>
#include <TFile.h>

enum eArg
  {
    kArgProcessName,
    kArgConfFile,
    kArgInFile,
    kArgOutRootFile,
    kArgc
  };


void closeFile(int sig)
{
  if(gFile)
    {
      gFile->Write();
      gFile->Close();
    }
}

bool UserStop = false;
void CatchSignal(int sig){
  UserStop = true;
  std::cout << "#D Catch Signal : Goto Exit process" << std::endl;
}


TROOT theROOT("sample-ana", "test");

int main(int argc, char* argv[])
{
  if (argc!=kArgc)
    {
      std::cout <<" #D Usage: "
		<<" [analyzer config file]"
		<<" [data input stream] "
		<<" [output root file]"
		<< std::endl;
      return 0;
    }

  std::vector<std::string> arg(argv, argv + argc);

  const std::string& confFile = arg[kArgConfFile];
  const std::string& inFile   = arg[kArgInFile];
  const std::string& rootFile = arg[kArgOutRootFile];

  //  ConfMan& gConfManager = ConfMan::GetInstance();
  ConfMan* gconfManager =new ConfMan(confFile);

  signal(SIGINT,CatchSignal);

  TFile f(rootFile.c_str(), "recreate");
  std::cout << "#D recreate root file : "
	    << rootFile << std::endl;

  if (!gconfManager->Initialize()){
    std::cerr << "#E gconfManager->Initialize failed ! in Main.cc" << std::endl;
    return 0;
  }
  gconfManager->WriteParameters();

  hddaq::unpacker::UnpackerManager& gUnpacker 
    = hddaq::unpacker::GUnpacker::get_instance();

  if (gconfManager->GetEvDispFlag()) {
    std::cout << "Create Event Display" << std::endl;
    gconfManager->InitializeEvDisp();
  }

  gUnpacker.set_istream(inFile);
  gUnpacker.initialize();
  //gUnpacker.set_skip_event(12324);
  //gUnpacker.set_skip_event(12325);
  //gUnpacker.set_skip_event(12326);
  //gUnpacker.set_skip_event(12327);
  //gUnpacker.set_skip_event(12323);

  //gUnpacker.set_skip_event(12331);
  //gUnpacker.set_skip_event(21182);
  //gUnpacker.set_skip_event(21183);
  //gUnpacker.set_skip_event(21184);
  //gUnpacker.set_skip_event(21185);
  //gUnpacker.set_skip_event(21186);
  //gUnpacker.set_skip_event(21187);

  int count = 0;
  for(;!gUnpacker.eof() && !UserStop;++gUnpacker){
    //  for(;!gUnpacker.eof();++gUnpacker){
    count++;
    //std::cout<<"Counter Number:    "<<count<<std::endl;
    VEvent* event = gconfManager->EventAllocator();
    event->ProcessingBegin();
    event->ProcessingNormal();
    event->ProcessingEnd();
    delete event;
  }

  gFile->Write();
  gFile->Close();
  
  std::cout << "#D : End of Analyzer" << std::endl;
  
  return 0;
}

