// -----------------------------------------------------------------------------
// CLASS: BCSChainHamiltonian
// -----------------------------------------------------------------------------
// Class for a BCSChainHamiltonian constructed in the particle-hole basis. The
// lattice of the original Hamiltonian would be of length L, with the
// transformed Hamiltonian having length 2L. 
//
//  METHODS:
//
//  init  -   accepts a list of variatonal parameters, constructs BCS
//            hamiltonian, calls eigen solver, and creates correct operator
//            representations for each variational parameter.
// -----------------------------------------------------------------------------

#ifndef BCS_HAMILTONIAN
#define BCS_HAMILTONIAN

#include <vector>
#include <Eigen/Dense>
#include <iostream>
#include "VarParam.h"

namespace VMC{

class BCSChainHamiltonian{
  private:
    size_t _size;
    Eigen::MatrixXd _bcsmatrix;
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> _solver;  
  public:
    BCSChainHamiltonian(const size_t&, ParamList_t&);
    void solve(){_solver.compute(_bcsmatrix);}
    void init(ParamList_t&);
    void setopers(ParamList_t&);
    void print_matrix(){std::cout << _bcsmatrix << std::endl;}
    Eigen::MatrixXd get_reduced_matrix(const size_t&);
    Eigen::MatrixXd get_eigenvecs(){return _solver.eigenvectors();}
    Eigen::MatrixXd get_eigenvals(){return _solver.eigenvalues();}
};
}

#endif // BCS_HAMILTONIAN
