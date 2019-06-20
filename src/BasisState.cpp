#include "BasisState.h"

namespace VMC {


BasisState::BasisState(const size_t& size) : _size(size) {
  _spinstate.resize(_size,0); 
  _operslist.resize(2*_size,0); 
}

size_t BasisState::find(const size_t& lindex) const {
  auto start=_operslist.begin();
  auto end=_operslist.end();
  return std::distance(start, std::find(start, end, lindex)); 
}  

void BasisState::refresh() {
  std::fill(_spinstate.begin(), _spinstate.end(), 0);
  std::fill(_operslist.begin(), _operslist.end(), 0);
}

}
