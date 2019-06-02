#ifndef HOPPING_HAMILTONIAN_H
#define HOPPING_HAMILTONIAN_H

#include "AbstractHamiltonian.h"

namespace VMC{
namespace AuxiliaryHamiltonians {

// Hopping Hamiltonian class 
// ============================================================================= 

class HoppingHamiltonian : AuxiliaryHamiltonian {
  private:
    bool                                            _trans_inv;
    size_t                                          _size;
    Eigen::MatrixXcd                                _hopping_matrix;
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> _solver;
  public:
    HoppingHamiltonian(
        const bool&,         // translation invariant
        const size_t&,       // size of expanded basis state 
        const AuxParamsSVec& // variational parameters 
    );
    void solve() {_solver.compute(_hopping_matrix);}
    void init(const AuxParamsSVec&);
    void set_vmats(const AuxParamsSVec&);
    void set_mmats(const AuxParamsSVec&);
    Eigen::MatrixXcd get_reduced_matrix(const size_t&);
};

// ============================================================================= 

} // namespace AuxiliaryHamiltonian
} // namespace VMC


#endif // HOPPING_HAMILTONIAN_H
