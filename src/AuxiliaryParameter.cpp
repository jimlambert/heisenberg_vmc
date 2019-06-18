#include "AuxiliaryParameter.h"

namespace VMC {
namespace Parameters {

void AuxiliaryParameter::operator() (
  BasisState& state, 
  Eigen::MatrixXcd& gmat
) {
  size_t size=state.size();
  Eigen::MatrixXcd redmat(size, 2*size); // N_e X 2L matrix
  for(size_t l=0; l<size; l++) {
    size_t lpos=state.find(l+1);
    for(size_t r=0; r<2*size; r++) redmat(l,r)=mmat(lpos,r);
  }
  local_meas.push((redmat*gmat).trace()); 
}

} // namespace Parameters
} // namespace VMC
