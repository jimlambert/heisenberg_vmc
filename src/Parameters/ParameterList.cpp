#include <iomanip>
#include <iostream>
#include <memory>
#include "AuxiliaryParameter.h"
#include "SpinSpin.h"
#include "ParameterList.h"
#include "iomacros.h"

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
  _naux+=1;
}

void ParameterList::build_jas_param(
  ParameterSubtype   jpt, // subtype for Jastrow parameter
  const std::string& n,   // parameter name 
  const double&      v,   // value of parameter
  const size_t&      s,   // parameter site i
  const size_t&      d,   // parameter site j
  const bool&        ti,  // translation invariant (true)
  const size_t&      bs   // binsize for local_meas   
) {
  switch(jpt) {
    case SPIN:
    {
      JasParamUPtr jasptr=std::make_unique<SpinSpin>
                          (n, v, s, d, ti, bs);
      _jas_params.push_back(std::move(jasptr));
      _njas+=1;
    }
    break;
    default:
      std::cout << "Unknown parameters in build_jas_param" << std::endl;
      exit(1);
    break;  
  }
}

Parameter& ParameterList::operator[] (const size_t& i) {
  if(i<_naux) return *_aux_params[i];
  else if(i<(_njas+_naux)) return *_jas_params[i-_naux];
  else {
    std::cout << "Parameter list access is out of bounds at index: " 
              << i 
              << std::endl;
    exit(1);
  }
}

void ParameterList::report_aux_params() {
  std::cout << "Auxiliary Parameters" << std::endl;
  std::cout << std::setw(INDENT_WIDTH+2*TERM_COLUMN_WIDTH)
            << std::setfill('=') << "=" 
            << std::endl;
  for(auto it=_aux_params.begin(); it!=_aux_params.end(); it++) {
    std::cout << std::setfill(' ')
              << std::setw(INDENT_WIDTH) << " "
              << std::setw(TERM_COLUMN_WIDTH) 
              << std::left << "Type:" 
              << (*it)->subtype << std::endl; 
    std::cout << std::setw(INDENT_WIDTH) << " "
              << std::setw(TERM_COLUMN_WIDTH) 
              << std::left << "Name:" 
              << (*it)->name << std::endl; 
    std::cout << std::setw(INDENT_WIDTH) << " "
              << std::setw(TERM_COLUMN_WIDTH) 
              << std::left << "Value:" 
              << (*it)->val << std::endl; 
    std::cout << std::setw(INDENT_WIDTH) << " " 
              << std::setw(TERM_COLUMN_WIDTH) 
              << std::left << "Site:" 
              << (*it)->site << std::endl; 
    std::cout << std::setw(INDENT_WIDTH) << " " 
              << std::setw(TERM_COLUMN_WIDTH) 
              << std::left << "dr:" 
              << (*it)->dr << std::endl; 
    std::cout << std::setw(INDENT_WIDTH) << " " 
              << std::setw(TERM_COLUMN_WIDTH) 
              << std::left << "trans. inv.:" 
              << (*it)->trans_inv << std::endl; 
    std::cout << std::setw(INDENT_WIDTH+2*TERM_COLUMN_WIDTH)
              << std::setfill('-') << "-" 
              << std::endl;
  }
  std::cout << std::setfill(' ');
}

void ParameterList::report_jas_params() {
  std::cout << "Jastrow Parameters" << std::endl;
  std::cout << std::setw(INDENT_WIDTH+2*TERM_COLUMN_WIDTH)
            << std::setfill('=') << "=" 
            << std::endl;
  for(auto it=_jas_params.begin(); it!=_jas_params.end(); it++) {
    std::cout << std::setfill(' ')
              << std::setw(INDENT_WIDTH) << " "
              << std::setw(TERM_COLUMN_WIDTH) 
              << std::left << "Type:" 
              << (*it)->subtype << std::endl; 
    std::cout << std::setw(INDENT_WIDTH) << " "
              << std::setw(TERM_COLUMN_WIDTH) 
              << std::left << "Name:" 
              << (*it)->name << std::endl; 
    std::cout << std::setw(INDENT_WIDTH) << " "
              << std::setw(TERM_COLUMN_WIDTH) 
              << std::left << "Value:" 
              << (*it)->val << std::endl; 
    std::cout << std::setw(INDENT_WIDTH) << " " 
              << std::setw(TERM_COLUMN_WIDTH) 
              << std::left << "Site:" 
              << (*it)->site << std::endl; 
    std::cout << std::setw(INDENT_WIDTH) << " " 
              << std::setw(TERM_COLUMN_WIDTH) 
              << std::left << "dr:" 
              << (*it)->dr << std::endl; 
    std::cout << std::setw(INDENT_WIDTH) << " " 
              << std::setw(TERM_COLUMN_WIDTH) 
              << std::left << "trans. inv.:" 
              << (*it)->trans_inv << std::endl; 
    std::cout << std::setw(INDENT_WIDTH+2*TERM_COLUMN_WIDTH)
              << std::setfill('-') << "-" 
              << std::endl;
  }
  std::cout << std::setfill(' ');
}

} // Parameters namespace
} // VMC namespace
