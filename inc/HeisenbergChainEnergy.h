#ifndef HEISENBERG_ENERGY_H
#define HEISENBERG_ENERGY_H

#include "AbstractObservable.h"

namespace VMC {
namespace Observables {

struct HeisenbergChainEnergy : Observable {
  HeisenbergChainEnergy(
    const std::string &n, 
    const size_t& bs,
    const double& jz,
    const double& jx,
    const double& jy
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

} // namespace Observables
} // namespace VMC

#endif // HEISENBERG_ENERGY_H
