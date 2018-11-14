#ifndef BCS_HAMILTONIAN
#define BCS_HAMILTONIAN

#include <vector>
#include <Eigen/Dense>
#include "VarParam.h"

namespace VMC{

class BCSChainHamiltonian{
  private:
    size_t _size;
    size_t _nele;
    Eigen::MatrixXd _bcsmatrix;
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> _solver;  
  public:
    BCSChainHamiltonian(const size_t&,    // _size 
                        const size_t&,    // _nele
                        const ParamList_t&  // list of parameters as vector
                       );
    void solve(){_solver.compute(_bcsmatrix);}
    void reinit(const ParamList_t&);
    Eigen::MatrixXd get_reduced_matrix();
    Eigen::MatrixXd get_eigenvecs(){return _solver.eigenvectors();}
    Eigen::MatrixXd get_eigenvals(){return _solver.eigenvalues();}
};
}

#endif // BCS_HAMILTONIAN
