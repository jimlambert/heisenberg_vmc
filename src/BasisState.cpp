#include "BasisState.h"

namespace VMC {

size_t BasisState::find(const size_t& lindex) {
  auto start=_operslist.begin();
  auto end=_operslist.end();
  return std::distance(start, std::find(start, end, lindex+1)); 
}  

}
