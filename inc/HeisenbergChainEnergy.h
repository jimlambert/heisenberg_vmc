#ifndef HEISENBERG_ENERGY_H
#define HEISENBERG_ENERGY_H

#include "AbstractObservable.h"

namespace VMC {
namespace Observables {

struct HeisenbergChainEnergy : Oberservable {
  void operator(const BasisState&, const Eigen::MatrixXcd&);
};

} // namespace Observables
} // namespace VMC

#endif // HEISENBERG_ENERGY_H
