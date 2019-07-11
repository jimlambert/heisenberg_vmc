#ifndef UTILS_H
#define UTILS_H

#include "BasisState.h"
#include "ParameterList.h"

namespace VMC {
namespace Utils {

// Utility functions
// =============================================================================
// Assortment of functions that are useful in the course of the VMC algorithm
// but which don't fit conveniently into a particular class
// =============================================================================

double compute_dj_exchange(
  const BasisState&, 
  const JasParamUVec&, 
  const size_t&,
  const size_t&,
  const int&,
  const int&
);

double compute_dj_flip(
  const BasisState&, 
  const JasParamUVec&, 
  const size_t&,
  const int&
);

// =============================================================================

} // namespace Utils
} // namespace VMC

#endif
