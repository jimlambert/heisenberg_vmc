#ifndef HOPPING_CHAIN_HAMILTONIAN_H
#define HOPPING_CHAIN_HAMILTONIAN_H

#include "AbstractHamiltonian.h"

namespace VMC{
namespace AuxiliaryHamiltonians {

// Hopping Chain Hamiltonian class 
// ============================================================================= 

class HoppingChainHamiltonian : AuxiliaryHamiltonian {
  private:
    bool                                            _bc; // true=PBC, false=APBC
    size_t                                          _size;
    Eigen::MatrixXcd                                _hopping_matrix;
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> _solver;

    void _check_init(const AuxParamUPtr&);
  public:
    HoppingChainHamiltonian(
        const bool&,        // boundary condition, PBC or APBC
        const size_t&,      // size of expanded basis state 
        const AuxParamUVec& // variational parameters 
    );
    void solve() {_solver.compute(_hopping_matrix);}
    void set_vmats(AuxParamUVec&);
    void set_mmats(AuxParamUVec&);
    void init(const AuxParamUVec&);
    Eigen::MatrixXcd get_reduced_matrix(const size_t&);
};

// ============================================================================= 

} // namespace AuxiliaryHamiltonian
} // namespace VMC


#endif // HOPPING_CHAIN_HAMILTONIAN_H
