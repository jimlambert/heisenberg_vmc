#include "utls.h"
#include "HeisenbergChainEnergy.h"

namespace VMC {
namespace Observables {

void HeisenbergChainEnergy::operator()(
  const BasisState&       state,
  const Eigen::MatrixXcd& gmat,
  const JasParamUVec&     jas_vec
) {
  std::complex<double> total=0.0;
  size_t size=state.size();
  for(size_t i=0; i<size; i++) {
    size_t j=(i+1)%size;
    total+=0.25*coeffs[0]*state[i]*state[j];
    if(state[i]==state[j]) continue;
    else if (state[i]==1){
      size_t iexpos1=i;
      size_t nexpos1=i+size;
      size_t iexpos2=j+size;
      size_t nexpos2=j;
      size_t lindex1=state(iexpos1)-1;
      size_t lindex2=state(iexpos2)-1;
      std::complex<double> det=gmat(nexpos1,lindex1)*gmat(nexpos2,lindex2)
        - gmat(nexpos2,lindex1)*(gmat(nexpos1, lindex2));
      double dj_sum=Utls::compute_dj_exchange(state, jas_vec, i, j, -1, 1); 
      total-=0.5*det*std::exp(dj_sum);
    }
    else {
      size_t iexpos1=i+size;
      size_t nexpos1=i;
      size_t iexpos2=j;
      size_t nexpos2=j+size;
      size_t lindex1=state(iexpos1)-1;
      size_t lindex2=state(iexpos2)-1;
      std::complex<double> det=gmat(nexpos1,lindex1)*gmat(nexpos2,lindex2)
        - gmat(nexpos2,lindex1)*gmat(nexpos1, lindex2);
      double dj_sum=Utls::compute_dj_exchange(state, jas_vec, i, j, 1, -1); 
      total-=0.5*det*std::exp(dj_sum);
    }
  }
  //std::cout << total << std::endl;
  local_meas.push(total);
}


} // namespace Observables
} // namespace VMC

