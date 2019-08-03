#include <cmath>
#include <iostream>
#include <iomanip>
#include <vector>
#include <algorithm>
#include "iomacros.h"
#include "SpinState.h"

namespace VMC {
namespace Utils {

void SpinState::random_state() {
  
  _spins_list.resize(_phys_size, 0);
  _opers_list.resize(_extd_size, 0);

  for(u_int phys_pos=0; phys_pos<_phys_size; phys_pos++) {
    double spin=_hilbert_space_ptr->random_spin();
    u_int extd_pos=get_extd_position(phys_pos, spin);
    _spins_list[phys_pos]=spin;

    // operators are ordered by their physical position, i.e. the first spin,
    // wherever it is in extended space, is the first operator.
    _opers_list[extd_pos]=phys_pos+1;    
  } 
}


u_int SpinState::find(const u_int& lindex) const {
  auto start=_opers_list.begin();
  auto end=_opers_list.end();
  return std::distance(start, std::find(start, end, lindex));
}


u_int SpinState::get_extd_position(const u_int& phys_pos) 
{ return get_extd_position(phys_pos, _spins_list[phys_pos]); }


u_int SpinState::get_extd_position(const u_int& phys_pos, const double& spin) {
  u_int index=spin+_hilbert_space_ptr->max_spin();  
  return phys_pos + (index*_phys_size);
}


void SpinState::report_spin_state() const {
  std::cout << "Spin state:" << std::endl;
  std::cout << "-----------" << std::endl;
  for(u_int i=0; i<_phys_size; i++) 
    std::cout << std::setw(8) << std::left << std::setfill(' ')
              << i+1;
  std::cout << std::endl;
  for(u_int i=0; i<_phys_size; i++) 
    std::cout << std::setw(8) << std::left << std::setfill(' ')
              << _spins_list[i];
  std::cout << std::endl;
}


void SpinState::report_oper_state() const {
  std::cout << "Operator list:" << std::endl;
  std::cout << "--------------" << std::endl;
  for(u_int i=0; i<_extd_size; i++) 
    std::cout << std::setw(4) << std::left << std::setfill(' ')
              << i+1; 
  std::cout << std::endl;
  for(u_int i=0; i<_extd_size; i++) 
    std::cout << std::setw(4) << std::left << std::setfill(' ')
              << _opers_list[i];
  std::cout << std::endl;
}


} // namespace Utils
} // namespace VMC
