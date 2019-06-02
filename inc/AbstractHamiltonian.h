#ifndef ABSTRACT_HAMILTONIAN_H
#define ABSTRACT_HAMILTONIAN_H

#include "AuxiliaryParameter.h"

namespace VMC {
namespace AuxiliaryHamiltonians {

// Parent class for auxiliary Hamiltonians
// ============================================================================= 

class AuxiliaryHamiltonian {
  public:
    virtual void solve() = 0;
    virtual void init(const AuxParamsSVec&) = 0;
    virtual Eigen::MatrixXcd get_reduced_matrix(const size_t&)=0;
};

// ============================================================================= 

} // namespace AuxiliaryHamiltonian
} // namespace VMC

#endif // ABSTRACT_HAMILTONIAN_H
