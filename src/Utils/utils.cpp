#include <cmath>
#include "utils.h"

namespace VMC {
namespace Utils {

double compute_dj_exchange(
  const BasisState&   state, 
  const JasParamUVec& jas_vec, 
  const size_t&       k,
  const size_t&       l,
  const int&          dk,
  const int&          dl
) {
  double total=0.0;
  size_t size=state.size();
  for(size_t i=0; i<jas_vec.size(); i++) {
    std::complex<double> v=jas_vec[i]->val;
    size_t dr=jas_vec[i]->dr;
    size_t kl=(k-dr+size)%size;
    size_t kr=(k+dr)%size;
    size_t ll=(l-dr+size)%size;
    size_t lr=(l+dr)%size;
    if(kl!=l) total+=0.5*(double)dk*v.real()*(double)state[kl];
    if(kr!=l) total+=0.5*(double)dk*v.real()*(double)state[kr];
    if(ll!=k) total+=0.5*(double)dl*v.real()*(double)state[ll];
    if(lr!=k) total+=0.5*(double)dl*v.real()*(double)state[lr];
  }
  return total;
}

double compute_dj_flip(
  const BasisState&   state, 
  const JasParamUVec& jas_vec, 
  const size_t&       k,
  const int&          dk 
) {
  double total=0.0;
  size_t size=state.size();
  for(size_t i=0; i<jas_vec.size(); i++) {
    size_t dr=jas_vec[i]->dr;
    std::complex<double> v=jas_vec[i]->val;
    size_t rl=(k-dr+size)%size;
    size_t rr=(k+dr)%size;
    total+=0.50*v.real()*(double)(state[rl]+state[rr])*dk;
  }
  return total;
}


double ladder_site(const size_t& nr, const size_t& ci, const size_t& si) {
  return (ci*nr) + si;
}

} // namespace VMC
} // namespace Utls
