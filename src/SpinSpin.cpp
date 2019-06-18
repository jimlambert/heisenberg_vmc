#include "SpinSpin.h"

namespace VMC {
namespace Parameters {

void SpinSpin::operator() (BasisState& state) {
  size_t size=state.size();  
  std::complex<double> total=0.0;
  if(trans_inv) {
   for(size_t r=0; r<size; r++) total+=0.25*state[r]*state[(r+dr)%size];
   total=total/(double)size;
  }
  else total+=state[site]*state[(site+dr)%size];
  local_meas.push(total);
}

} // namespace Parameters 
} // namespace VMC 
