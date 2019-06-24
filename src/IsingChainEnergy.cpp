#include "IsingChainEnergy.h"

namespace VMC {
namespace Observables {

void IsingChainEnergy::operator()(
  const BasisState&       state, 
  const Eigen::MatrixXcd& gmat,
  const JasParamUVec&     jas_vec
) {
  std::complex<double> total=0.0;
  size_t l=state.size();
  for(size_t i=0; i<l; i++) {
    total+=coeffs[0]*0.25*state[i]*state[(i+1)%l]; 
  }
  local_meas.push(total);
}

}
}
