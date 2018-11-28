#include <Eigen/Dense>
#include <vector>
#include <iostream>
#include "VarParam.h"
#include "BCSHamiltonian.h"

namespace VMC{

BCSChainHamiltonian::BCSChainHamiltonian
(const size_t& L, ParamList_t& params) : _L(L) {init(params);}

void BCSChainHamiltonian::init(ParamList_t& params) {
  _bcsmatrix = Eigen::MatrixXd::Zero(_L, _L);
  for(auto it=params.begin(); it!=params.end(); it++)
  switch(it->type) { 
    case Onsite:
      for(size_t i=0; i<_L; i++) {
        if(i<(_L/2)) _bcsmatrix(i, i)=-1.0*it->val;
        else _bcsmatrix(i, i)=it->val;
      }
      break;
    case Hopping:
      for(size_t i=0; i<_L; i++){
        if(i<(_L/2)) {
          _bcsmatrix(i, (i+it->space)%(_L/2))+=it->val;
          _bcsmatrix((i+it->space)%(_L/2), i)+=it->val;
        } 
        else {
          _bcsmatrix(i, (i+it->space)%(_L/2)+_L/2)+=-1.0*it->val;
          _bcsmatrix((i+it->space)%(_L/2)+_L/2, i)+=-1.0*it->val;
        } 
      }
      break;
    case Pairing:
      for(size_t i=0; i<_L/2; i++) {
        _bcsmatrix(i, (i+it->space)%(_L/2)+_L/2)+=it->val;
        _bcsmatrix((i+it->space)%(_L/2)+_L/2, i)+=it->val;
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
  std::cout << U << std::endl;
  for(auto it=params.begin(); it!=params.end(); it++) {
    Eigen::MatrixXd UVU(_L, _L);
    Eigen::MatrixXd Q(_L, _L);
    Q = Eigen::MatrixXd::Zero(_L, _L);
    UVU = U.adjoint()*(it->vmat)*U;
    for(size_t i=0; i<_L; i++)
    for(size_t j=0; j<_L; j++)
    if((i>(_L/2))&&(j<=(_L/2))) Q(i,j)=UVU(i,j)/(e(i)-e(j));
    it->mmat=U*Q*(U.adjoint());
  } 
}

Eigen::MatrixXd BCSChainHamiltonian::get_reduced_matrix
(const size_t& nele) {
  Eigen::MatrixXd _evecs = get_eigenvecs();
  Eigen::MatrixXd redmat(_L, nele);
  for(size_t i=0; i<_L; i++)
  for(size_t j=0; j<nele; j++)
    redmat(i, j)=_evecs(i, j);
  return redmat;
}
}

