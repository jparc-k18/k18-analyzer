#ifndef FILEMAN_H
#define FILEMAN_H

#include<fstream>
#include<string>
#include<map>

//---------------------------------------------------------------------------
// class definition of FileManager
//---------------------------------------------------------------------------

/*                        Autoconfig file structure                        */
/*

  # Autoconfig file hoge hoge       <= Line starting with "#" is ignored
  NameParam ManagerName             <= Write "ParamName" then "name"
  ../../param/hoge.param            <= Write Path to param file

  NameString ManagerName            <= Write "StringName" then "name"
  ../../param/hoge.string           <= Write Path to string file
  
  ...

  ParamName  is Declaration of ParameterManager.
  StringName is Declaration of StringManager.

  ManagerName is used to call GetParamMan or GetStrMan in your program.
*/
/*                                                                         */

namespace h_Utility{

  template<typename pType>
  class ParameterManager;

  typedef ParameterManager<double>      ParamMan;
  typedef ParameterManager<std::string> StrMan;

  class FileManager{
  public:
    //Functions ----------------------------------------------------------------
    ParamMan& GetParamMan(const char* ruName);
    ParamMan& GetParamMan(const std::string& ruName);
    StrMan&   GetStrMan(const char* Name);
    StrMan&   GetStrMan(const std::string& Name);
    int       MyStatus();
    int       stInitialize(const char* ruName);
    int       stInitialize(const std::string& ruName);

    static FileManager& GetInstance();

    //Controller -------------------------------------------------------------
    FileManager();
      ~FileManager();

  private:
    typedef std::map<std::string, ParamMan*>  ParamContainor;
    typedef std::pair<std::string, ParamMan*> ParamPair;

    typedef std::map<std::string, StrMan*>    StringContainor;
    typedef std::pair<std::string, StrMan*>   StringPair;

    enum            Flags{Me,                     // success reports
			  ParamError, StrError,   // error reports
			  size_Flags};
    bool            flag_[size_Flags];

    ParamContainor  PM_Assemblage;
    unsigned int    sizePM_Assemblage;

    StringContainor SM_Assemblage;
    unsigned int    sizeSM_Assemblage;

    void            Initializer();
    void            CleanUp();
    void            Debug();

    //forbidden --------------------------------------------------------------
    FileManager(const FileManager& object);
    FileManager& operator =(const FileManager& object);
  };


  // Get Instance
  inline FileManager& FileManager::GetInstance(){
    static FileManager object;
    return object;
  }

  inline int FileManager::MyStatus(){
    int stReturn = -1;
    if(flag_[Me]){stReturn = 0;}
    return stReturn;
  }
}
#endif
