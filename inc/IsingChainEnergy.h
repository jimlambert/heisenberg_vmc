#ifndef ISING_CHAIN_ENERGY_H
#define ISING_CHAIN_ENERGY_H

#include "AbstractObservable.h"

namespace VMC {
namespace Observables {

struct IsingChainEnergy : Observable {
  IsingChainEnergy(
    const std::string& n, 
    const size_t& bs, 
    const double& j
  ) : Observable(n, bs) {coeffs.push_back(j);}
  void operator()(
    const BasisState&, 
    const Eigen::MatrixXcd&,
    const JasParamUVec&
  );
};

}
}

#endif // ISING_CHAIN_ENERGY_H
