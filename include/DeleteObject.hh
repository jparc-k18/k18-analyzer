#ifndef OBJDELETE
#define OBJDELETE

#include<utility>
// This is tomonori san no pakuri

namespace h_Utility{

  struct DeleteObject{
    
    template<typename Type>
    void operator ()(Type*& object) const {
      if(object){
	delete object;
	object = NULL;
      }
      return;
    }
    
    template<typename Type1, typename Type2>
    void operator ()(std::pair<Type1, Type2*>& object) const {
      operator ()(object.second);
      return;
    }
  };
}

#endif
