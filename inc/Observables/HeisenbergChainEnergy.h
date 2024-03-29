#ifndef HEISENBERG_ENERGY_H
#define HEISENBERG_ENERGY_H

#include "AbstractObservable.h"

namespace VMC {
namespace Observables {

// Heisenberg chain energy functor
// =============================================================================

struct HeisenbergChainEnergy : Observable {
  HeisenbergChainEnergy(
    const std::string &n, // name of observable
    const size_t& bs,     // binsize for LocalObservable
    const double& jz,     // j(i)=spin exchange along ith axis
    const double& jx,     // .
    const double& jy      // .
  ) : Observable(n, bs) {
    coeffs.push_back(jz);
    coeffs.push_back(jx);
    coeffs.push_back(jy);
  }
  void operator()(
    const BasisState&, 
    const Eigen::MatrixXcd&,
    const JasParamUVec&
  );
};

// =============================================================================

} // namespace Observables
} // namespace VMC

#endif // HEISENBERG_ENERGY_H
