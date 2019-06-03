#include <iostream>
#include <Eigen/Dense>
#include "HoppingHamiltonian.h"

namespace VMC {
namespace AuxiliaryHamiltonians {

HoppingHamiltonian::HoppingHamiltonian(
  const bool&          ti,
  const size_t&        n,
  const AuxParamSVec&  params_vec
) : _trans_inv(ti), _size(n) {
  init(params_vec);
  solve();
}

void HoppingHamiltonian::init(const AuxParamSVec& params_vec) {
  _hopping_matrix = Eigen::MatrixXcd::Zero(_size,_size);
}

void HoppingHamiltonian::set_vmats(const AuxParamSVec& params_vec) {
}

void HoppingHamiltonian::set_mmats(const AuxParamSVec& params_vec) {
  Eigen::MatrixXcd U=_solver.eigenvectors();
  Eigen::VectorXd  e=_solver.eigenvalues();
  for(auto it=params_vec.begin(); it!=params_vec.end(); it++) {
    Eigen::MatrixXcd UVU(_size,_size);
    Eigen::MatrixXcd Q=Eigen::MatrixXcd::Zero(_size,_size);
    UVU=U.adjoint()*((*it)->vmat)*U;
    for(size_t i=0; i<_size; i++)
    for(size_t j=0; j<_size; j++) 
      if((i>((_size/2)-1))&&(j<=(_size/2)-1)) Q(i,j)=UVU(i,j)/(e(j)-e(i));
    (*it)->mmat=U*Q*(U.adjoint());
  }
}

Eigen::MatrixXcd HoppingHamiltonian::get_reduced_matrix(const size_t& nele) {
  Eigen::MatrixXcd _evecs = _solver.eigenvectors();
  Eigen::MatrixXcd redmat(_size, nele);
  for(size_t i=0; i<_size; i++)
  for(size_t j=0; j<nele; j++)
    redmat(i, j)=_evecs(i, j);
  return redmat;
}

} // namespace AuxiliaryHamiltonians
} // namespace VMC
