#ifndef KJG_LADDER_ENERGY_H
#define KJG_LADDER_ENERGY_H

#include "AbstractObservable.h"

namespace VMC {
namespace Observables {

// KJG ladder energy functor
// =============================================================================

struct KjgLadderEnergy : Observable {
  KjgLadderEnergy(
    const std::string &n, // name of observable
    const size_t& bs,     // binsize 
    const double& j,      // heisenberg exchange
    const double& k,      // kitaev exchange
    const double& g       // symmetric off-diagonal exchange
  ) : Observable(n, bs) {
    coeffs.push_back(j);
    coeffs.push_back(k);
    coeffs.push_back(g);
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

#endif // KJG_LADDER_ENERGY_H
