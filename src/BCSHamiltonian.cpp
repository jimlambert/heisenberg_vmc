#include <Eigen/Dense>
#include <vector>
#include <iostream>
#include "VarParam.h"
#include "BCSHamiltonian.h"

namespace VMC{

BCSChainHamiltonian::BCSChainHamiltonian
(const size_t& L, const ParamList_t& params) : _size(L) {
  reinit(params);
  _solver.compute(_bcsmatrix);
}

void BCSChainHamiltonian::reinit(const ParamList_t& params) {
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
}

void BCSChainHamiltonian::setopers(const ParamList_t&) {
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

