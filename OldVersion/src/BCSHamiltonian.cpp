#include <Eigen/Dense>
#include <vector>
#include <iostream>
#include "VarParam.h"
#include "BCSHamiltonian.h"

namespace VMC{

BCSChainHamiltonian::BCSChainHamiltonian
(const size_t& L, ParamList_t& params) : _L(L) {init(params);}

void BCSChainHamiltonian::init(ParamList_t& params) {
  _bcsmatrix = Eigen::MatrixXcd::Zero(_L, _L);
  for(auto it=params.begin(); it!=params.end(); it++) {
    _bcsmatrix+=(it->val)*(it->vmat); 
  }
  _solver.compute(_bcsmatrix);
  setopers(params);
  std::cout << "Operators initialized" << std::endl;
}

void BCSChainHamiltonian::setopers(ParamList_t& params) {
  Eigen::MatrixXcd U; // unitary transformation
  Eigen::VectorXd e; // single particle energies
  U=_solver.eigenvectors();
  e=_solver.eigenvalues();
  for(auto it=params.begin(); it!=params.end(); it++) {
    Eigen::MatrixXcd UVU(_L, _L);
    Eigen::MatrixXcd Q(_L, _L);
    Q = Eigen::MatrixXcd::Zero(_L, _L);
    UVU = U.adjoint()*(it->vmat)*U;
    for(size_t i=0; i<_L; i++)
    for(size_t j=0; j<_L; j++)
      if((i>((_L/2)-1))&&(j<=((_L/2)-1))) Q(i,j)=UVU(i,j)/(e(j)-e(i));
    it->mmat=U*Q*(U.adjoint());
    //std::cout << it->name << std::endl;
    //std::cout << "Q" << std::endl;
    //std::cout << Q << std::endl;
    //std::cout << "UVU" << std::endl;
    //std::cout << UVU << std::endl;
  }
  //std::cout << "====================" << std::endl; 
}

Eigen::MatrixXcd BCSChainHamiltonian::get_reduced_matrix
(const size_t& nele) {
  Eigen::MatrixXcd _evecs = get_eigenvecs();
  Eigen::MatrixXcd redmat(_L, nele);
  for(size_t i=0; i<_L; i++)
  for(size_t j=0; j<nele; j++)
    redmat(i, j)=_evecs(i, j);
  return redmat;
}
}

