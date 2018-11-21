#include <Eigen/Dense>
#include <vector>
#include <iostream>
#include "VarParam.h"
#include "BCSHamiltonian.h"

namespace VMC{

BCSChainHamiltonian::BCSChainHamiltonian
(const size_t& L, ParamList_t& params) : _size(L) {init(params);}

void BCSChainHamiltonian::init(ParamList_t& params) {
  _bcsmatrix = Eigen::MatrixXd::Zero(_size, _size);
  for(auto it=params.begin(); it!=params.end(); it++)
  switch(it->type) { 
    case Onsite:
      for(size_t i=0; i<_size; i++) {
        if(i<(_size/2)) _bcsmatrix(i, i)=-1.0*it->val;
        else _bcsmatrix(i, i)=it->val;
      }
      break;
    case Hopping:
      for(size_t i=0; i<_size; i++){
        if(i<(_size/2)) {
          _bcsmatrix(i, (i+it->space)%(_size/2))+=it->val;
          _bcsmatrix((i+it->space)%(_size/2), i)+=it->val;
        } 
        else {
          _bcsmatrix(i, (i+it->space)%(_size/2)+_size/2)+=-1.0*it->val;
          _bcsmatrix((i+it->space)%(_size/2)+_size/2, i)+=-1.0*it->val;
        } 
      }
      break;
    case Pairing:
      for(size_t i=0; i<_size/2; i++) {
        _bcsmatrix(i, (i+it->space)%(_size/2)+_size/2)+=it->val;
        _bcsmatrix((i+it->space)%(_size/2)+_size/2, i)+=it->val;
      }
      break;
  }
  _solver.compute(_bcsmatrix);
  setopers(params);
}

void BCSChainHamiltonian::setopers(ParamList_t& params) {
  Eigen::MatrixXd U; // unitary transformation
  Eigen::VectorXd e; // single particle energies
  U=_solver.eigenvectors();
  e=_solver.eigenvalues();
  for(auto it=params.begin(); it!=params.end(); it++) {
    Eigen::MatrixXd UVU(_size, _size);
    Eigen::MatrixXd Q(_size, _size);
    Q = Eigen::MatrixXd::Zero(_size, _size);
    std::cout << "Structure matrix: " << std::endl;
    std::cout << it->vmat << std::endl;
    std::cout << "-----" << std::endl;
    UVU = U.adjoint()*(it->vmat)*U;
    for(size_t i=0; i<_size; i++)
    for(size_t j=0; j<_size; j++)
    if((i>(_size/2))&&(j<=(_size/2))) Q(i,j)=UVU(i,j)/(e(i)-e(j));
    it->mmat=U*Q*U.adjoint();
  } 
}

Eigen::MatrixXd BCSChainHamiltonian::get_reduced_matrix
(const size_t& nele) {
  Eigen::MatrixXd _evecs = get_eigenvecs();
  Eigen::MatrixXd redmat(_size, nele);
  for(size_t i=0; i<_size; i++)
  for(size_t j=0; j<nele; j++)
    redmat(i, j)=_evecs(i, j);
  return redmat;
}
}

