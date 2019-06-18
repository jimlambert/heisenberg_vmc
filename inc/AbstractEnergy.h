#ifndef ABSTRACT_ENERGY_H
#define ABSTRACT_ENERGY_H

#include <memory>
#include <complex>
#include "BasisState.h"

namespace VMC {
namespace Energies {

struct Energy {
  virtual std::complex<double> operator(BasisState&)=0;
};

} // namespace Energies

typedef std::unique_ptr<Energy> EnergyUPtr;

} // namespace VMC

#endif // ABSTRACT_ENERGY_H
