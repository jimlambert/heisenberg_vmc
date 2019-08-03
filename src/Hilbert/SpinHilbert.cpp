#include <iostream>
#include "SpinHilbert.h"

namespace VMC {
namespace Hilbert {


SpinHilbert::SpinHilbert(const double& max_spin) : _max_spin(max_spin) {
  _nstates=(2*_max_spin) + 1;
  for(size_t i=0; i<_nstates; i++) _states.push_back(-max_spin+i);
}

double SpinHilbert::random_spin() {
  // generate random number and partition unity
  double rnum=_random_double(_mteng);
  double part=1.0/(double)_nstates;
  double spin;
  for(size_t i=1; i<_nstates+1; i++) {
    if(rnum<(double)(i*part)) {
      spin=_states[i-1];
      break;
    }
    else continue;
  }
  return spin;
}

} // namespace Hilbert
} // namespace VMC
