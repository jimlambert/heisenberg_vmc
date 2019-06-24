#ifndef HOPPING_CHAIN_HAMILTONIAN_H
#define HOPPING_CHAIN_HAMILTONIAN_H

#include "AbstractHamiltonian.h"

namespace VMC{
namespace AuxiliaryHamiltonians {

// Hopping Chain Hamiltonian class 
// ============================================================================= 

class HoppingChainHamiltonian : public AuxiliaryHamiltonian {
  private:
    bool                                            _bc; // true=PBC, false=APBC
    size_t                                          _size;
    Eigen::MatrixXcd                                _hopping_matrix;
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> _solver;

    void _check_vmat_init(const AuxParamUPtr&);
  public:
    HoppingChainHamiltonian(
        const bool&,   // boundary condition, PBC or APBC
        const size_t&, // size of expanded basis state 
        AuxParamUVec&  // variational parameters 
    );
    void solve() {_solver.compute(_hopping_matrix);}
    void set_vmats(AuxParamUVec&);
    void set_mmats(AuxParamUVec&);
    void init(AuxParamUVec&);
    Eigen::MatrixXcd  get_reduced_matrix(const size_t&);
    Eigen::MatrixXcd  get_matrix(){return _hopping_matrix;}
    Eigen::MatrixXcd  get_eigenvectors(){return _solver.eigenvectors();}
    Eigen::VectorXd   get_eigenvalues(){return _solver.eigenvalues();}
};

// ============================================================================= 

} // namespace AuxiliaryHamiltonians

typedef AuxiliaryHamiltonians::HoppingChainHamiltonian HopChainHam;

} // namespace VMC


#endif // HOPPING_CHAIN_HAMILTONIAN_H
