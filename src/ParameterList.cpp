#include <iomanip>
#include <iostream>
#include <memory>
#include "ParameterList.h"

#define COLUMN_WIDTH 10

namespace VMC {
namespace Parameters {
  
void ParameterList::push_aux_param(AuxParamUPtr auxptr) {
  //_aux_params.push_back(std::move(auxptr));
  _naux+=1;
}

void ParameterList::push_jas_param(JasParamUPtr jasptr) {
  //_jas_params.push_back(std::move(jasptr));
  _njas+=1;
}

void ParameterList::report_jas_params() {
  for(auto it=_jas_params.begin(); it!=_jas_params.end(); it++) {
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
  for(auto it=_aux_params.begin(); it!=_aux_params.end(); it++) {
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
