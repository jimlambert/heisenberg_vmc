#ifndef ABSTRACT_OBSERVABLE_H
#define ABSTRACT_OBSERVABLE_H

#include <memory>
#include <Eigen/Dense>
#include <vector>
#include <string>
#include "LocalMeasurement.h"
#include "BasisState.h"

namespace VMC {
namespace Observables {

// Observable parent struct
// =============================================================================
// All observables including the energy inherit from this class. For the energy 
// in particular the coeffs vector is used to store value of J used in the
// calculation
// =============================================================================

struct Observable {
  std::string name;
  std::vector<double> coeffs; // stores coefficient used in computing energies
  MeasCd local_meas;
  Observable(const std::string& n, const size_t& bs) 
    : name(n), local_meas(bs) {} 
  virtual void operator()(const BasisState&, const Eigen::MatrixXcd&)=0;
};

// =============================================================================

} // namespace Observables

typedef std::unique_ptr<Observables::Observable> ObsUPtr;
typedef std::shared_ptr<Observables::Observable> ObsSPtr;
typedef std::vector<ObsUPtr>                     ObsUVec;
typedef std::vector<ObsSPtr>                     ObsSVec;

} // namespace VMC

#endif // ABSTRACT_OBSERVABLE_H
