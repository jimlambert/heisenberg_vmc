#ifndef ZZ_CORRELATOR_H
#define ZZ_CORRELATOR_H

namespace VMC {
namespace Observables {

struct ZzCorrelator : Observable {
  bool trans_inv;  
  ZzCorrelator(
    const std::string& n, 
    const size_t& bs,
    const bool& ti
  ) : Observable(n, bs), trans_inv(ti) {}
  void operator()(
    const BasisState&,
    const Eigen::MatrixXcd&,
    const JasParamUVec&
  )
}

}
}

#endif // ZZ_CORRELATOR_H
