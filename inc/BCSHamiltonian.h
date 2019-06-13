// -----------------------------------------------------------------------------
// CLASS: BCSChainHamiltonian
// -----------------------------------------------------------------------------
// Class for a BCSChainHamiltonian constructed in the particle-hole basis. The
// lattice of the original Hamiltonian would be of length L, with the particle-
// hole transformed Hamiltonian having length 2L. 
//
// --------------------------
//  PUBLIC MEMBER FUNCTIONS:
// --------------------------
// 
//  --> void init(ParamList_t&)
//  ----------------------------
//  accepts a parameter list (defined in VarParam.h) and constructs a BCS 
//  Hamiltonian with corresponding operators. Includes onsite operators, 
//  hopping operators, and pairing operators. 
// -----------------------------------------------------------------------------

#ifndef BCS_HAMILTONIAN
#define BCS_HAMILTONIAN

#include <vector>
#include <Eigen/Dense>
#include <iostream>
#include "VarParam.h"

namespace VMC {

class BCSChainHamiltonian {
  private:
    size_t _L;  
    // not actually Hamiltonian, but tunneling matrix T.
    Eigen::MatrixXcd _bcsmatrix;
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> _solver;   
  public:
    // accepts size of auxiliary Hamiltonian and parameter list.
    BCSChainHamiltonian(const size_t&, ParamList_t&);
    void solve(){_solver.compute(_bcsmatrix);} 
    void init(ParamList_t&); 
    void setopers(ParamList_t&); 
    void print_matrix(){std::cout << _bcsmatrix << std::endl;} 
    Eigen::MatrixXcd get_hamiltonian(){return _bcsmatrix;}
    Eigen::MatrixXcd get_reduced_matrix(const size_t&); 
    Eigen::MatrixXcd get_eigenvecs(){return _solver.eigenvectors();}
    Eigen::MatrixXd get_eigenvals(){return _solver.eigenvalues();}
};

} // namespace VMC

#endif // BCS_HAMILTONIAN
