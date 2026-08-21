#ifndef DECODE_EASIROC
#define DECODE_EASIROC


#include<iostream>
#include<fstream>
#include<cstring>

#include"gzfilter"
#include"DetectorID.hh"

namespace h_Utility{
  
  class decodeEASIROC{
  public:
    int  stDecodeOneEvent();

    int  GetMaxSize();
    int  GetSize(int plane, int seg);
    int  GetTDC(int plane, int seg, int mhit, int data_type); // 0 leading, 1 traling;
    void WriteContainerSize(int plane, int* size);
    void WriteContainer(int plane, int seg,
			int* dContainer_L, int* dContainer_T);
    
    static decodeEASIROC& GetInstance();
    int stLoadDataFile(std::string& runNo);
    int MyStatus();
    ~decodeEASIROC();

  private:
    enum Flag{Me, size_Flag};
    bool flag_[size_Flag];

    static const int NofFile = 10;
    static const int NofMHit = 8;

    int MaxSize;
    int size_L[NumOfPlaneBFT][NumOfSegBFT];
    int size_T[NumOfPlaneBFT][NumOfSegBFT];
    int dataContainer_L[NumOfPlaneBFT][NumOfSegBFT][NofMHit];
    int dataContainer_T[NumOfPlaneBFT][NumOfSegBFT][NofMHit];

    std::ifstream        *file_[NofFile];
    h_Utility::igzfilter *file[NofFile];

    void fill_Upstream_0246(int FileNum, int NofWord);
    void fill_Downstream_1357(int FileNum, int NofWord);
    void fill_shita_8(int FileNum, int NofWord);
    void fill_ue_9(int FileNum, int NofWord);

    void Initializer();
    void CleanUp();

    decodeEASIROC();
    decodeEASIROC(const decodeEASIROC& object);
    decodeEASIROC& operator =(const decodeEASIROC& object);
  };

  inline decodeEASIROC& decodeEASIROC::GetInstance(){
    static decodeEASIROC object;
    return object;
  }

  inline int decodeEASIROC::MyStatus(){
    int stReturn = -1;
    if(flag_[Me]){stReturn = 0;};
    return stReturn;
  }

  inline int decodeEASIROC::GetMaxSize(){
    return MaxSize;
  }

  inline int decodeEASIROC::GetSize(int plane, int seg){
    int dReturn = 0;
    if(plane < 0 || plane > (NumOfPlaneBFT-1) ){
      std::cerr << "#E : GetSize invalid plane" 
		<< plane << std::endl;
    }else{
      if(seg < 0 || seg > (NumOfSegBFT-1) ){
	std::cerr << "#E : GetSize invalid seg" << std::endl;
      }else{
	dReturn = size_L[plane][seg];
      }
    }
    return dReturn;
  }

  inline int decodeEASIROC::GetTDC(int plane, int seg,
				   int mhit, int data_type)
  {
    int dReturn = 0;
    if(plane < 0 || plane > (NumOfPlaneBFT-1) ){
      std::cerr << "#E : GetTDC invalid plane" 
		<< plane << std::endl;
    }else{
      if(seg < 0 || seg > (NumOfSegBFT-1) ){
	std::cerr << "#E : GetTDC invalid seg" << std::endl;
      }else{
	if(0 == data_type){
	  dReturn = dataContainer_L[plane][seg][mhit];
	}else if(1 == data_type){
	  dReturn = dataContainer_T[plane][seg][mhit];
	}
      }
    }
    return dReturn;
  }
    
  inline void decodeEASIROC::WriteContainerSize(int plane,int* size)
  {
    if(plane < 0 || plane > (NumOfPlaneBFT-1) ){
      std::cerr << "#E : WriteContainerSize invalid plane" 
		<< plane << std::endl;
    }else{
      memcpy(size, size_L[plane], NumOfSegBFT*sizeof(int));
    }

    return;
  }

  inline void decodeEASIROC::WriteContainer(int plane,
					    int seg,
					    int* dContainer_L,
					    int* dContainer_T)
  {
    if(plane < 0 || plane > (NumOfPlaneBFT-1) ){
      std::cerr << "#E : WriteContainerSize invalid plane" 
		<< plane << std::endl;
    }else{
      if(seg < 0 || seg > (NumOfSegBFT-1) ){
	std::cerr << "#E : WriteContainer invalid seg" << std::endl;
      }else{
	memcpy(dContainer_L, dataContainer_L[plane][seg], 8*sizeof(int));
	memcpy(dContainer_T, dataContainer_T[plane][seg], 8*sizeof(int));
      }
  }
    return;
  }
}
#endif
