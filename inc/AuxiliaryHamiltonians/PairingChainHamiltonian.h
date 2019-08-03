#ifndef PAIRING_CHAIN_HAMILTONIAN_H
#define PAIRING_CHAIN_HAMILTONIAN_H

#include "AbstractHamiltonian.h"

namespace VMC {
namespace AuxiliaryHamiltonians {

// Pairing Chain Hamiltonian class 
// ============================================================================= 

class PairingChainHamiltonian : public AbstractHamiltonian {

  private:
    bool                                            _bc; // true=PBC, false=APBC
    size_t                                          _size;
    double                                          _noise_coeff;
    Eigen::MatrixXcd                                _hopping_matrix;
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> _solver;

    AuxParamUVec _noise_vec;
    void _check_vmat_init(const AuxParamUPtr&);
  public:
    PairingChainHamiltonian(
        const bool&,   // boundary condition, PBC or APBC
        const size_t&, // size of expanded basis state 
        AuxParamUVec&, // variational parameters 
        const double&,  // noise parameter
        AuxParamUVec&  // noise vector 
    );
    void solve() {_solver.compute(_hopping_matrix);}
    void set_vmats(AuxParamUVec&);
    void set_mmats(AuxParamUVec&);
    void init(AuxParamUVec&);
    double& noise() { return _noise_coeff; }
    const double& noise() const { return _noise_coeff; }
    Eigen::MatrixXcd  get_reduced_matrix(const size_t&);
    Eigen::MatrixXcd  get_matrix(){return _hopping_matrix;}
    Eigen::MatrixXcd  get_eigenvectors(){return _solver.eigenvectors();}
    Eigen::VectorXd   get_eigenvalues(){return _solver.eigenvalues();}
};

// ============================================================================= 

} // namespace AuxiliaryHamiltonians

typedef AuxiliaryHamiltonians::PairingChainHamiltonian BCSChainHam;

} // namespace VMC


#endif // PAIRING_CHAIN_HAMILTONIAN_H
