#include "utils.h"
#include "KjgLadderEnergy.h"

namespace VMC {
namespace Observables {

void KjgLadderEnergy::operator()(
  const BasisState&       state, 
  const Eigen::MatrixXcd& gmat,
  const JasParamUVec&     jas_vec
) {
  double Jl=coeffs[0];
  double Jr=coeffs[1];
  double K=coeffs[2];
  double g=coeffs[3];
  std::complex<double> total=0.0;
  size_t size=state.size();
  size_t nrungs=size/2;
  for(size_t s=0; s<nrungs; s++) {
    for(size_t c=0; c<2; c++) { // intrachain coupling
      size_t i=Utils::ladder_site(nrungs, c, s); 
      size_t j=Utils::ladder_site(nrungs, c, (s+1)%nrungs);
      int    x=std::pow(-1, (1+c+s));
      total+=0.25*Jl*state[i]*state[j];
      if(state[i]==state[j]) total+=0.25*x*K;
      else if (state[i]==1){
        size_t iexpos1=i;
        size_t nexpos1=i+size;
        size_t iexpos2=j+size;
        size_t nexpos2=j;
        size_t lindex1=state(iexpos1)-1;
        size_t lindex2=state(iexpos2)-1;
        std::complex<double> det=gmat(nexpos1,lindex1)*gmat(nexpos2,lindex2)
          - gmat(nexpos2,lindex1)*(gmat(nexpos1, lindex2));
        double dj_sum=Utils::compute_dj_exchange(state, jas_vec, i, j, -1, 1); 
        total-=(Jl*0.5 + K*0.25)*det*std::exp(dj_sum);
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
        double dj_sum=Utils::compute_dj_exchange(state, jas_vec, i, j, 1, -1); 
        total-=(Jl*0.5 + K*0.25)*det*std::exp(dj_sum);
      }
    }
    size_t i=Utils::ladder_site(nrungs, 0, s);
    size_t j=Utils::ladder_site(nrungs, 1, s);
    total+=0.25*(Jr+K)*state[i]*state[j];
  }
  //std::cout << total << std::endl;
  local_meas.push(total);
}

} // namespace Observables
} // namespace VMC
