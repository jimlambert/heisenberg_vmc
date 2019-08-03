#include "utils.h"
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
    if(state[i]==1) {
      size_t iexpos=i;
      size_t nexpos=i+l;
      size_t lindex=state(iexpos)-1;
      std::complex<double> amp=gmat(nexpos, lindex);
      double dj=Utils::compute_dj_flip(state, jas_vec, i, -2);
      total+=0.5*coeffs[1]*exp(dj)*amp;
      //std::cout << "I:" << '\t' << iexpos << std::endl;
      //std::cout << "K:" << '\t' << nexpos << std::endl;
      //std::cout << "l:" << '\t' << lindex+1 << std::endl;
      //std::cout << "amp:" << '\t' << amp << std::endl;
    }
    else {
      size_t iexpos=i+l;
      size_t nexpos=i;
      size_t lindex=state(iexpos)-1;
      std::complex<double> amp=gmat(nexpos, lindex);
      double dj=Utils::compute_dj_flip(state, jas_vec, i, 2);
      total+=0.5*coeffs[1]*exp(dj)*amp;
      //std::cout << "I:" << '\t' << iexpos << std::endl;
      //std::cout << "K:" << '\t' << nexpos << std::endl;
      //std::cout << "l:" << '\t' << lindex+1 << std::endl;
      //std::cout << "amp:" << '\t' << amp << std::endl;
    }
  }
  //std::cout << "Energy" << total << std::endl;
  local_meas.push(total);
}

}
}
