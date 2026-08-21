#ifndef PARAMMAN_H
#define PARAMMAN_H

#include<iostream>
#include<fstream>
#include<vector>
#include<string>
#include<map>

//---------------------------------------------------------------------------
// class definition of ParameterManager
//---------------------------------------------------------------------------

/*                         Parameter file structure                        */
/*

  # Paremeter file hoge hoge        <= Line starting with "#" is ignored
  StartDefinition SetName           <= Write "StartDefinition" then "name"
   NameA 1 2 3 4 ....               <= Write "ParameName" then values
   NameB 100 200 ...                
  EndDefinition                     <= End of one parameter set
  
  StartDefinition ....              <= Start of next parameter set
  ...

*/
/*                                                                         */

namespace h_Utility{

  template<typename pType>
  class ParameterManager{
  public:
    //Functions ---------------------------------------------------------------
    pType GetParameter(int ruIndex, const char* ruPramName);
    pType GetParameter(int ruIndex, const std::string& ruPramName);
    pType GetParameter(int ruIndex,
			const char* ruPramName,
			const char* ruSetName);
    pType GetParameter(int ruIndex,
			const std::string& ruPramName,
			const std::string& ruSetName);
    
    int    stChangeParamSet(const char* ruName);
    int    stChangeParamSet(const std::string& ruName);

    int    sizeData_ROOT();
    int    sizeData(const char* ruParamName);          // looking CurrentData
    int    sizeData(const std::string& ruParamName);   // looking CurrentData
    int    sizeData(const char* ruParamName,
		    const char* ruSetName);
    int    sizeData(const std::string& ruParamName,
		    const std::string& ruSetName);

    int    MyStatus();

    //Controller --------------------------------------------------------------
     ParameterManager(); // Dummy Instance (Do nothing)
     ParameterManager(std::ifstream& rFile);
    ~ParameterManager();

  private:
    typedef std::vector<pType>                   ValContainor;
    typedef std::map<std::string, ValContainor>   DataContainor;
    typedef std::map<std::string, DataContainor>  ROOTContainor;

    typedef std::pair<std::string, ValContainor>  DataPair;
    typedef std::pair<std::string, DataContainor> ROOTPair;

    static  std::string MyName;

    enum    Flags{Me, StartDef,        // Success reports
		  UndefData, NoEndDef, // Error reports
		  size_Flags};

    bool          flag_[size_Flags];
    ROOTContainor ROOTData;
    DataContainor CurrentData;
    int           sizeROOTData;

    std::string   NameData(std::ifstream& rFile);
    int           stSetOneData(std::ifstream& rFile, const std::string& rName);
    void          Initializer();
    void          CleanUp();
    void          Debug();

    //forbidden ---------------------------------------------------------------

    ParameterManager(const ParameterManager& object);
    ParameterManager& operator = (const ParameterManager& object);
  };
}

using namespace h_Utility;

template<typename pType>
std::string ParameterManager<pType>::MyName = "ParameterManager";

// Implementation of methods -------------------------------------------------
// sizeData_ROOT
template<typename pType>
inline
int ParameterManager<pType>::sizeData_ROOT(){
  return ROOTData.size();
}

// sizeData (char*)
template<typename pType>
inline
int ParameterManager<pType>::sizeData(const char* ruName){
  const std::string ruTempName = ruName;
  return sizeData(ruTempName);
}

// sizeData (std::string)
template<typename pType>
inline
int ParameterManager<pType>::sizeData(const std::string& ruName){
  static const std::string MyFunc = "::sizeData ";

  int sizeReturn = -1;
  typename DataContainor::iterator itrData = CurrentData.find(ruName);
  if(itrData == CurrentData.end()){
    std::cerr << "#E "<< MyName << MyFunc
	      << ruName << " is Not Found" << std::endl;
  }else{
    // Process normal
    sizeReturn = itrData->second.size();
  }

  return sizeReturn;
}

// sizeData 2 arg (char*) 
template<typename pType>
inline
int ParameterManager<pType>::sizeData(const char* ruParamName,
				      const char* ruSetName)
{
  const std::string ruTempParamName = ruParamName;
  const std::string ruTempSetName   = ruSetName;
  return sizeData(ruTempParamName, ruTempSetName);
}

// sizeData 2 arg (std::string) 
template<typename pType>
inline
int ParameterManager<pType>::sizeData(const std::string& ruParamName,
				      const std::string& ruSetName)
{
  static const std::string MyFunc     = "::sizeData ";
               int         sizeReturn = -1;

  typename ROOTContainor::iterator itrROOT = ROOTData.find(ruSetName);
  typename DataContainor::iterator itrData = itrROOT->second.find(ruParamName);
  if(false 
     || itrROOT == ROOTData.end()
     || itrData == CurrentData.end()
     ){
    std::cerr << "#E "<< MyName << MyFunc
	      << "Not Found" << std::endl;
  }else{
    // Process normal
    sizeReturn = itrData->second.size();
  }

  return sizeReturn;
}

template<typename pType>
inline
int ParameterManager<pType>::MyStatus(){
  int stReturn = -1;
  if(flag_[Me]){stReturn = 0;};
  return stReturn;
}

#include<sstream>

// Functions -----------------------------------------------------------------
// GetParameter (Same parameter set) (char*)
template<typename pType>
pType ParameterManager<pType>::GetParameter(int ruIndex,
					    const char* ruName)
{
  const std::string ruTempName = ruName;
  return GetParameter(ruIndex, ruTempName);
}

// GetParameter (Same parameter set) (std::string)
template<typename pType>
pType ParameterManager<pType>::GetParameter(int ruIndex,
					    const std::string& ruName)
{
  static const std::string MyFunc     = "::GetParameter ";
  static pType pInitialize;
  pType upValReturn = pInitialize;

  if(0 == CurrentData.size()){
    std::cerr << "#E " << MyName << MyFunc
	      << "Data Set " << ruName << " is empty" << std::endl;
  }

  typename DataContainor::iterator itrData = CurrentData.find(ruName);
  if(itrData == CurrentData.end()){
    std::cerr << "#E " << MyName << MyFunc 
	      << ruName << " is Not Found" << std::endl;
  }else{
    if(ruIndex+1 > (int)itrData->second.size()){
      std::cerr << "#E " << MyName << MyFunc
		<< "Index:" << ruIndex << " is too big for " 
		<< ruName << std::endl;
    }else{
      // Process normal
      int Index = ruIndex; // System From User
      upValReturn = itrData->second.at(Index);
    }
  }

  return upValReturn;;
}

// GetParameter (Chose parameter set) (char*)
template<typename pType>
pType ParameterManager<pType>::GetParameter(int ruIndex,
					    const char* ruParamName,
					    const char* ruSetName)
{
  const std::string ruTempParamName = ruParamName;
  const std::string ruTempSetName   = ruSetName;
  return GetParameter(ruIndex, ruTempParamName, ruTempSetName);
}

// GetParameter (Chose parameter set) (std::string)
template<typename pType>
pType ParameterManager<pType>::GetParameter(int ruIndex,
					    const std::string& ruParamName,
					    const std::string& ruSetName)
{
  static const std::string MyFunc     = "::GetParameter ";
  static pType pInitialize;
  pType upValReturn = pInitialize;

  typename ROOTContainor::iterator itrROOT = ROOTData.find(ruSetName);
  if(itrROOT == ROOTData.end()){
    std::cerr << "#E " << MyName << MyFunc
	      << ruSetName << " is Not Found in ROOT" << std::endl;
  }else{
    typename DataContainor::iterator itrData = itrROOT->second.find(ruParamName);
    if(itrData == itrROOT->second.end()){
      std::cerr << "#E " << MyName << MyFunc
		<< ruParamName << " is Not Found" << std::endl;
    }else{
      if(ruIndex+1 > (int)itrData->second.size()){
	std::cerr << "#E " << MyName << MyFunc
		  << "Index:" << ruIndex << " is too big for " 
		  << ruParamName << std::endl;
      }else{
	// Process normal
	int Index = ruIndex; // System From User
	upValReturn = itrData->second.at(Index);
      }
    }
  }
  
  return upValReturn;;
}

// SetParamSet (char*)
template<typename pType>
int ParameterManager<pType>::stChangeParamSet(const char* ruSetName)
{
  const std::string ruTempSetName = ruSetName;
  return stChangeParamSet(ruTempSetName);
}

// SetParamSet (std::string)
template<typename pType>
int ParameterManager<pType>::stChangeParamSet(const std::string& ruSetName)
{
  static const std::string MyFunc   = "::stChangeParamSet ";
               int         stReturn = -1;
  
  typename ROOTContainor::iterator itrROOT = ROOTData.find(ruSetName);
  if(itrROOT == ROOTData.end()){
    std::cerr << "#E " << MyName << MyFunc
	      << ruSetName << " is Not Found in ROOT" << std::endl;
  }else{
    // Process normal
    CurrentData = itrROOT->second;
    stReturn = 0;
  }

  return stReturn;
}

// Controller ----------------------------------------------------------------
// Constructor
template<typename pType>
ParameterManager<pType>::ParameterManager(std::ifstream& rFile){
  static const std::string MyFunc = "::Constructor ";
  
  Initializer();
  
  while(!rFile.eof()){
    std::string rName = NameData(rFile);
    // No start definition
    if("E R R O R" == rName){
      // EOF or UndefData (look flag_[UndefData])
      continue;
    }else{
      // Process normal
      if(-1 == stSetOneData(rFile, rName)){
	// No end definition
	flag_[NoEndDef] = true;
	std::cerr << "#E " << MyName << MyFunc
		  << "No EndDefinition in " << rName
		  << std::endl;
      }
    }
  }

  if(!flag_[UndefData] && !flag_[NoEndDef]){
    flag_[Me] = true;
  }else{
    CleanUp();
  }

  return;
}

// Dummy constructor
template<typename pType>
ParameterManager<pType>::ParameterManager(){
  static const std::string MyFunc = "::Constructor ";
  Initializer();
  CleanUp();
  std::cerr << "#E " << MyName << MyFunc
	    << "Dummy Instance is created" << std::endl;
}

// Destructor
template<typename pType>
ParameterManager<pType>::~ParameterManager(){
  CleanUp();
}

// Private Functions ---------------------------------------------------------
// NameData
template<typename pType>
std::string ParameterManager<pType>::NameData(std::ifstream& rFile){
  static const std::string MyFunc     = "::Constructor ";
               std::string NameReturn = "E R R O R";

  std::string uLineBuffer;
  while(getline(rFile, uLineBuffer)){
    // No chars
    if(0 == uLineBuffer.size()){
      continue;
    }
    
    // Comment
    if('#' == uLineBuffer.at(0)){
      continue;
    }
    
    std::istringstream uLine_ToU(uLineBuffer);
    std::string uWordBuffer;
    uLine_ToU >> uWordBuffer;
    // Normal process
    if(uWordBuffer == "StartDefinition"){
      flag_[StartDef] = true;
      if(uLine_ToU >> uWordBuffer){
	NameReturn = uWordBuffer; // System From User
	break;
      }else{
	flag_[UndefData] = true;
	std::cerr << "#E " << MyName << MyFunc
		  << "Undefined data" << std::endl;
      }
    }
  }
  
  return NameReturn;
}

// stSetOneData
template<typename pType>
int ParameterManager<pType>::stSetOneData(std::ifstream& rFile,
					  const std::string& rName)
{
  std::string   uLineBuffer;
  std::string   WordBuffer;
  DataContainor OneData;
  int           stReturn = -1;

  while(getline(rFile, uLineBuffer)){
    // No chars
    if(0 == uLineBuffer.size()){
      continue;
    }

    std::istringstream uLine_To(uLineBuffer);
    uLine_To >> WordBuffer;

    // Comment
    if(WordBuffer.at(0) == '#'){
      continue;
    }else if(WordBuffer == "EndDefinition"){
      // End of one data set
      stReturn = 0;
      break;
    }
    
    // Now WordBuffer is "Name"
    std::string& rParamName = WordBuffer;
    ValContainor OneValContainor;
    pType       OneParam;
    while(uLine_To >> OneParam){
      OneValContainor.push_back(OneParam);
    }
    
    // One data set was made successfully
    OneData.insert(DataPair(rParamName, OneValContainor));
  }

  if(0 == stReturn){
    // Process normal
    ROOTData.insert(ROOTPair(rName, OneData));
    CurrentData = OneData;
  }

  return stReturn;
}

// Initializer
template<typename pType>
void ParameterManager<pType>::Initializer(){
  for(int i = 0; i < (int)size_Flags; ++i){
    flag_[i] = false;
  }
  return;
}

// CleanUp
template<typename pType>
void ParameterManager<pType>::CleanUp(){
  ROOTData.clear();
  return;
}

// Debug
template<typename pType>
void ParameterManager<pType>::Debug(){
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
#endif

