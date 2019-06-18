#ifndef ISING_ENERGY_H
#define ISING_ENERGY_H

#include "AbstractEnergy.h"

namespace VMC {
namespace Energies {

struct IsingEnergy : Energy {
  std::complex<double> operator(BasisState&);
};

}
}

#endif // ISING_ENERGY_H
