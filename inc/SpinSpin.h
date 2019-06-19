#ifndef SPIN_SPIN_H
#define SPIN_SPIN_H

#include "JastrowParameter.h"

namespace VMC {
namespace Parameters {

// struct SpinSpin
// =============================================================================

struct SpinSpin : JastrowParameter {
  SpinSpin(
    const std::string& n,  // parameter name 
    const double&      v,  // value of parameter
    const size_t&      s,  // parameter site
    const size_t&      d,  // distance to next parameter
    const bool&        ti, // translation invariant (true)
    const size_t&      bs  // binsize for local_meas 
  ) : JastrowParameter(SPIN, n, v, s, d, ti, bs) {}  
  void operator() (BasisState&, Eigen::MatrixXcd&);
};

// =============================================================================

} // namespace Parameters
} // namespace VMC

#endif // SPIN_SPIN_H
