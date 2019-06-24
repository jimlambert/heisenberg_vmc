#include <cmath>
#include "utls.h"

namespace VMC {
namespace Utls {

double compute_dj_exchange(
  const BasisState&   state, 
  const JasParamUVec& jas_vec, 
  const size_t&       k,
  const size_t&       l,
  const int&          dk,
  const int&          dl
) {
  double sum=0.0;
  size_t size=state.size();
  for(size_t i=0; i<jas_vec.size(); i++) {
    std::complex<double> v=jas_vec[i]->val;
    size_t dr=jas_vec[i]->dr;
    size_t kl=(k-dr+size)%size;
    size_t kr=(k+dr)%size;
    size_t ll=(l-dr+size)%size;
    size_t lr=(l+dr)%size;
    if(kl!=l) sum+=0.5*(double)dk*v.real()*(double)state[kl];
    if(kr!=l) sum+=0.5*(double)dk*v.real()*(double)state[kr];
    if(ll!=k) sum+=0.5*(double)dl*v.real()*(double)state[ll];
    if(lr!=k) sum+=0.5*(double)dl*v.real()*(double)state[lr];
  }
  return std::exp(sum);
}

} // namespace VMC
} // namespace Utls
