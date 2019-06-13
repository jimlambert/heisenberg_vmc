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
  _aux_params.push_back(std::move(auxptr));
}

void ParameterList::push_jas_param(JasParamUPtr jasptr) {
  //_jas_params.push_back(std::move(jasptr));
  _njas+=1;
  _jas_params.push_back(std::move(jasptr));
}

void ParameterList::build_aux_param(
  ParameterSubtype   apt, // auxiliary parameter type
  const std::string& n,   // parameter name 
  const double&      v,   // value of parameter
  const size_t&      s,   // parameter site 
  const size_t&      d,   // change in position
  const bool&        ti,  // translation invariant (true)
  const size_t&      bs   // binsize for local_meas 
) {
  AuxParamUPtr auxptr=std::make_unique<AuxiliaryParameter>
                      (apt, n, v, s, d, ti, bs);
  _aux_params.push_back(std::move(auxptr)); 
}

//void ParameterList::build_jas_param(
//  ParameterSubtype   jpt, // subtype for Jastrow parameter
//  const std::string& n,   // parameter name 
//  const double&      v,   // value of parameter
//  const size_t&      s,   // parameter site i
//  const size_t&      d,   // parameter site j
//  const bool&        ti,  // translation invariant (true)
//  const size_t&      bs   // binsize for local_meas   
//) {
//  JasParamUPtr jasptr=std::make_unique<JastrowParameter>
//                      (jpt, n, v, s, d, ti, bs);
//  _jas_params.push_back(std::move(jasptr));
//}

void ParameterList::report_aux_params() {
  std::cout << "====================" << std::endl;
  std::cout << "Auxiliary Parameters" << std::endl;
  std::cout << "====================" << std::endl << std::endl;
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
              << std::left << "site:" 
              << std::setw(COLUMN_WIDTH) 
              << (*it)->site << std::endl; 
    std::cout << '\t' << std::setw(COLUMN_WIDTH) 
              << std::left << "dr:" 
              << std::setw(COLUMN_WIDTH) 
              << (*it)->dr << std::endl; 
  }
}

void ParameterList::report_jas_params() {
  std::cout << "==================" << std::endl;
  std::cout << "Jastrow Parameters" << std::endl;
  std::cout << "==================" << std::endl << std::endl;
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
              << std::left << "site:" 
              << std::setw(COLUMN_WIDTH) 
              << (*it)->site << std::endl; 
    std::cout << '\t' << std::setw(COLUMN_WIDTH) 
              << std::left << "dr:" 
              << std::setw(COLUMN_WIDTH) 
              << (*it)->dr << std::endl; 
  }
}

} // Parameters namespace
} // VMC namespace
