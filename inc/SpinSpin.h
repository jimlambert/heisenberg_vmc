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
    const size_t&      si, // parameter site 1
    const size_t&      sj, // parameter site 2
    const bool&        ti, // translation invariant (true)
    const size_t&      bs  // binsize for local_meas 
  ) : JastrowParameter(SPIN, n, v, si, sj, ti, bs) {}  
  void operator() (BasisState&);
};

// =============================================================================

} // Parameters namespace

} // VMC namespace

#endif // SPIN_SPIN_H
