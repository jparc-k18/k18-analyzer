#include"FileManager.hh"
#include"ParameterManager.hh"
#include"DeleteObject.hh"
#include<iostream>
#include<sstream>
#include<algorithm>

using namespace h_Utility;

static const std::string MyName  = "FileManager";

//Functions -----------------------------------------------------------------

// GetParamMan (char*)
ParamMan& FileManager::GetParamMan(const char* ruName){
  const std::string ruTempName = ruName;
  return GetParamMan(ruTempName);
}

// GetParamMan (std::string)
ParamMan& FileManager::GetParamMan(const std::string& ruName){
  static const std::string MyFunc = "::GetParamMan ";
  ParamMan* ptrPMReturn = NULL;
  
  ParamContainor::iterator itrPMA = PM_Assemblage.find(ruName);
  if(itrPMA == PM_Assemblage.end()){
    std::cerr << "#E " << MyName << MyFunc
	      << ruName << " is Not Found" << std::endl;
    static ParamMan dummy; // Dummy Instance
    ptrPMReturn = &dummy;
  }else{
    ptrPMReturn = itrPMA->second;
  }

  return *ptrPMReturn;
}

// GetStrMan (char*)
StrMan& FileManager::GetStrMan(const char* ruName){
  const std::string ruTempName = ruName;
  return GetStrMan(ruTempName);
}

// GetStrMan (std::string)
StrMan& FileManager::GetStrMan(const std::string& ruName){
  static const std::string MyFunc = "::GetStrMan ";
  StrMan* ptrSMReturn = NULL;
  
  StringContainor::iterator itrSMA = SM_Assemblage.find(ruName);
  if(itrSMA == SM_Assemblage.end()){
    std::cerr << "#E " << MyName << MyFunc
	      << ruName << " in Not Found" << std::endl;
    static StrMan dummy; // Dummy Instance
    ptrSMReturn = &dummy;
  }else{
    ptrSMReturn = itrSMA->second;
  }

  return *ptrSMReturn;
}

//Controller ----------------------------------------------------------------
// Constructor
FileManager::FileManager(){
  
}

// Destructor
FileManager::~FileManager(){
  CleanUp();
}

// stInitialize (char*)
int FileManager::stInitialize(const char* ruName){
  std::string ruTempName = ruName;
  return stInitialize(ruTempName);
}

// stInitialize (std::string)
// Load file and make container of ParameterManager
int FileManager::stInitialize(const std::string& ruName){
  static const std::string MyFunc = "::stInitialize ";
  int status = 0;
  
  Initializer();

  std::ifstream rFile(ruName.c_str()); // system from user
  if(rFile.fail()){
    std::cerr << "#E : " << MyName << MyFunc
	      << "Auntocofig file not exists" << std::endl;
    status = -1;
  }

  std::string uLineBuffer;
  while(getline(rFile, uLineBuffer)){
    std::istringstream uLine_ToU(uLineBuffer);

    std::string uWordBuffer;
    while(uLine_ToU >> uWordBuffer){
      if(0 == uWordBuffer.size()){
	// empty
	break;
      }

      // What is this line ? Comment, Param or String
      if('#' == uWordBuffer.at(0)){
	// comment
	break;
      }

      if(uWordBuffer == "NameParam"){
	// ParameterManager
	uLine_ToU >> uWordBuffer;
	std::string rParamManName = uWordBuffer; // system from user

	std::cout << MyName << " Loading " << rParamManName << std::endl;
	getline(rFile, uLineBuffer);
	std::ifstream rParamFile(uLineBuffer.c_str());
	if(!rParamFile.is_open()){
	  std::cerr << "#E " << MyName << MyFunc
		    << uLineBuffer << " is Not Found"
		    << std::endl;
	  flag_[ParamError] = true;
	}else{
	  ParamMan *OneParamMan = new ParamMan(rParamFile);
	  PM_Assemblage.insert(ParamPair(rParamManName, OneParamMan));
	}
	
	break;
      }

      if(uWordBuffer == "NameString"){
	// StringManager
	uLine_ToU >> uWordBuffer;
	std::string rStringManName = uWordBuffer; // system from user
	
	std::cout << MyName << " Loading " << rStringManName << std::endl;
	getline(rFile, uLineBuffer);
	std::ifstream rStringFile(uLineBuffer.c_str());
	if(!rStringFile.is_open()){
	  std::cerr << "#E "<<MyName << MyFunc
		    << uLineBuffer << " is Not Found"
		    << std::endl;
	  flag_[StrError] = true;
	}else{
	  StrMan *OneStringMan  = new StrMan(rStringFile);
	  SM_Assemblage.insert(StringPair(rStringManName, OneStringMan));
	}

	break;
      }
    }
  }

  if(!flag_[ParamError] && !flag_[StrError]){
    sizePM_Assemblage  = PM_Assemblage.size();
    sizeSM_Assemblage = SM_Assemblage.size();
    flag_[Me] = true;
  }else{
    CleanUp();
  }

  return status;
}

// Private -------------------------------------------------------------------
// _Initializer
void FileManager::Initializer(){
  for(int i = 0; i < (int)size_Flags; ++i){
    flag_[i] = false;
  }
  return;
}

// _CleanUp
void FileManager::CleanUp(){
  ParamContainor::iterator itrPM     = PM_Assemblage.begin();
  ParamContainor::iterator itrPM_end = PM_Assemblage.end();
  for_each(itrPM, itrPM_end, h_Utility::DeleteObject());
  PM_Assemblage.clear();

  StringContainor::iterator itrSM     = SM_Assemblage.begin();
  StringContainor::iterator itrSM_end = SM_Assemblage.end();
  for_each(itrSM, itrSM_end, h_Utility::DeleteObject());
  SM_Assemblage.clear();

  return;
}

// _Debug
void FileManager::Debug(){
  for(int i = 0; i < (int)size_Flags; ++i){
    if(flag_[i]){
      std::cout << "!";
    }else{
      std::cout << ".";
    }
  }
  std::cout << std::endl;
  return;
}
