#include <iomanip>
#include <iostream>
#include "ParameterList.h"

#define COLUMN_WIDTH 10

namespace VMC {
namespace Parameters {
  
void ParameterList::push(ParamUniquePtr param_ptr) {
  ParamSharedPtr shared_param_ptr=std::move(param_ptr); 
  switch((param_ptr)->get_type()) {
    case JASTROW:
      //_jas_params.push_back(shared_param_ptr);
      _njas+=1;
      break;
    case AUXILIARY:
      //_aux_params.push_back(shared_param_ptr);
      _naux+=1;
      break;
  }
  _params.push_back(shared_param_ptr);
}

void ParameterList::report_jas_params() {
  for(auto it=_params.begin(); it!=_params.end(); it++) {
    switch((*it)->get_type()) {
      case JASTROW:
        std::cout << "JASTROW TYPE" << std::endl;
        break;
      case AUXILIARY:
        std::cout << "AUXILIARY TYPE" << std::endl;
        break;
    } 
    std::cout << '\t' << std::setw(COLUMN_WIDTH) 
              << std::left << "Name:" 
              << std::setw(COLUMN_WIDTH) 
              << (*it)->name << std::endl; 
    std::cout << '\t' << std::setw(COLUMN_WIDTH) 
              << std::left << "Value:" 
              << std::setw(COLUMN_WIDTH) 
              << (*it)->val << std::endl; 
    std::cout << '\t' << std::setw(COLUMN_WIDTH) 
              << std::left << "site_i:" 
              << std::setw(COLUMN_WIDTH) 
              << (*it)->site_i << std::endl; 
    std::cout << '\t' << std::setw(COLUMN_WIDTH) 
              << std::left << "site_j:" 
              << std::setw(COLUMN_WIDTH) 
              << (*it)->site_j << std::endl; 
  }
}

void ParameterList::report_aux_params() {
  for(auto it=_params.begin(); it!=_params.end(); it++) {
    switch((*it)->get_type()) {
      case JASTROW:
        std::cout << "JASTROW TYPE" << std::endl;
        break;
      case AUXILIARY:
        std::cout << "AUXILIARY TYPE" << std::endl;
        break;
    } 
    std::cout << '\t' << std::setw(COLUMN_WIDTH) 
              << std::left << "Name:" 
              << std::setw(COLUMN_WIDTH) 
              << (*it)->name << std::endl; 
    std::cout << '\t' << std::setw(COLUMN_WIDTH) 
              << std::left << "Value:" 
              << std::setw(COLUMN_WIDTH) 
              << (*it)->val << std::endl; 
    std::cout << '\t' << std::setw(COLUMN_WIDTH) 
              << std::left << "site_i:" 
              << std::setw(COLUMN_WIDTH) 
              << (*it)->site_i << std::endl; 
    std::cout << '\t' << std::setw(COLUMN_WIDTH) 
              << std::left << "site_j:" 
              << std::setw(COLUMN_WIDTH) 
              << (*it)->site_j << std::endl; 
  }
}

} // Parameters namespace
} // VMC namespace
