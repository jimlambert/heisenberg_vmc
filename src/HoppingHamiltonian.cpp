#include <iostream>
#include <Eigen/Dense>
#include "HoppingHamiltonian.h"

namespace VMC {
namespace AuxiliaryHamiltonians {

HoppingHamiltonian::HoppingHamiltonian(
  const bool&          ti,
  const size_t&        n,
  const AuxParamsSVec& params_vec
) : _trans_inv(ti), _size(n) {
  init(params_vec);
  solve();
}

void HoppingHamiltonian::init(const AuxParamsSVec& params_vec) {
  _hopping_matrix = Eigen::MatrixXcd::Zero(_size,_size);
  for(auto it=params_vec.begin(); it!=params_vec.end(); it++) 
  switch((*it)->get_type()) {
    case Parameters::JASTROW: // ignore jastrow parameters
      continue; 
      break;
    case Parameters::AUXILIARY: 
      switch((*it)->get_subtype()) {
        case Parameters::ONSITE:
          if(_trans_inv) 
            for(size_t r=0; r<_size; r++) _hopping_matrix(r,r)+=(*it)->val;
          else {
            size_t r=(*it)->site_i;
            _hopping_matrix(r,r)+=(*it)->val;
          }
          break; // case Parameters::ONSITE
        case Parameters::HOPPING:  
          if(_trans_inv) {
            size_t dr=((*it)->site_i)-((*it)->site_j);
            for(size_t r=0; r<_size; r++) {
              _hopping_matrix(r,(r+dr)%_size)+=(*it)->val;
              _hopping_matrix((r+dr)%_size, r)=_hopping_matrix(r,(r+dr)%_size);
            }
          } 
          else {
            size_t ri=(*it)->site_i;
            size_t rj=(*it)->site_j;
            _hopping_matrix(ri,rj)=(*it)->val;
          }
          break; // case Parameters::HOPPING
        default:
          std::cout << "Invalid variational parameter." << std::endl;
          std::cout << (*it)->name 
                    << " is incompatible with hopping Hamiltonian"
                    << std::endl;
          exit(1);
          break; // default
      } // inner switch
      break; // case Parameters::AUXILIARY
  } // outer switch
}

void HoppingHamiltonian::set_vmats(const AuxParamsSVec& params_vec) {
  for(auto it=params_vec.begin(); it!=params_vec.end(); it++) {
    (*it)->vmat=Eigen::MatrixXd::Zero(_size,_size);
    switch((*it)->get_type()) {
      case Parameters::JASTROW: // ignore jastrow parameters
        continue; 
        break;
      case Parameters::AUXILIARY: 
        switch((*it)->get_subtype()) {
          case Parameters::ONSITE:
            if(_trans_inv) 
            for(size_t r=0; r<_size; r++) {
              (*it)->vmat(r,r)=1;
            }
            else {
              size_t r=(*it)->site_i;
              (*it)->vmat(r,r)=1;
            }
            break; // case Parameters::ONSITE
          case Parameters::HOPPING:  
            if(_trans_inv) {
              size_t dr=((*it)->site_i)-((*it)->site_j);
              for(size_t r=0; r<_size; r++) {
                (*it)->vmat(r,(r+dr)%_size)=1;
                (*it)->vmat((r+dr)%_size, r)=(*it)->vmat(r,(r+dr)%_size);
              }
            } 
            else {
              size_t ri=(*it)->site_i;
              size_t rj=(*it)->site_j;
              (*it)->vmat(ri,rj)=1;
              (*it)->vmat(rj,ri)=(*it)->vmat(ri,rj);
            }
            break; // case Parameters::HOPPING
          default:
            std::cout << "Invalid variational parameter." << std::endl;
            std::cout << (*it)->name 
                      << " is incompatible with hopping Hamiltonian"
                      << std::endl;
            exit(1);
            break; // default
        } // inner switch
        break; // case Parameters::AUXILIARY
    } // outer switch
  } // parameter loop
}

void HoppingHamiltonian::set_mmats(const AuxParamsSVec& params_vec) {
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
