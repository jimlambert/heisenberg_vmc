#ifndef ISING_CHAIN_ENERGY_H
#define ISING_CHAIN_ENERGY_H

#include "AbstractObservable.h"

namespace VMC {
namespace Observables {

// Ising chain energy functor
// =============================================================================

struct IsingChainEnergy : Observable {
  IsingChainEnergy(
    const std::string& n, // observable name 
    const size_t& bs,     // binsize for LocalObservable
    const double& j,      // zz exchange
    const double& h       // transverse magnetic field
  ) : Observable(n, bs) {
    coeffs.push_back(j);
    coeffs.push_back(h);
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

#endif // ISING_CHAIN_ENERGY_H
