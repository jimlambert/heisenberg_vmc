#include <iostream>
#include <Eigen/Dense>
#include "AuxiliaryParameter.h"
#include "PairingLadderHamiltonian.h"

namespace VMC {
namespace AuxiliaryHamiltonians {

// Private Members
// ============================================================================= 
 


void PairingLadderHamiltonian::_check_vmat_init(const AuxParamUPtr& aux_ptr) {
  if(!(aux_ptr->vinit)) {
    std::cout << "Parameter " << aux_ptr->name << " has unitialized vmatrix" 
              << std::endl;
    exit(1);
  }
}


void PairingLadderHamiltonian::_set_onsite(const AuxParamUPtr&) {


}


void PairingLadderHamiltonian::_set_hopping(const AuxParamUPtr&) {


}


void PairingLadderHamiltonian::_set_pairing(const AuxParamUPtr&) {
  

}

// ============================================================================= 
// Public Members
// ============================================================================= 

PairingLadderHamiltonian::PairingLadderHamiltonian(
  const bool&          b,
  const size_t&        nr,
  AuxParamUVec&  params_vec
) : _bc(b), _nrungs(nr), _size(2*nr) {
  set_vmats(params_vec);
  init(params_vec);
}

void PairingLadderHamiltonian::init(AuxParamUVec& params_vec) {
  _hopping_matrix = Eigen::MatrixXcd::Zero(_size,_size);
  for(auto it=params_vec.begin(); it!=params_vec.end(); it++) {
    _check_vmat_init((*it));  
    _hopping_matrix+=((*it)->val)*((*it)->vmat);
  }
  solve();
  set_mmats(params_vec);
  Eigen::VectorXd e_vals=_solver.eigenvalues();
  if(fabs(e_vals(_size/2-1)-e_vals(_size/2))<1e-3) {
    std::cout << "Error: open shell configuration, choose different parameters"
              << std::endl;
    //std::cout << e_vals << std::endl;
    //exit(1);
  }
}


void PairingLadderHamiltonian::set_vmats(AuxParamUVec& params_vec) {
  for(auto it=params_vec.begin(); it!=params_vec.end(); it++) {
    bool ti=(*it)->trans_inv;
    (*it)->vmat=Eigen::MatrixXd::Zero(_size,_size);
    switch((*it)->subtype) {
      case Parameters::ONSITE:
      {
        _set_onsite(*it);
      }
      break; // case ONSITE
      case Parameters::HOPPING:
      {
        _set_hopping(*it);
      }
      break; // case Hopping
      case Parameters::PAIRING:
      {
        _set_pairing(*it);
      }
      break; // case PAIRING
      default:
      {
        std::cout << "Forbidden parameter type for " << (*it)->name 
                  << std::endl;
        exit(1);
      }
      break; // case default
    } // subtype switch
    (*it)->vinit=true;
  }
}

void PairingLadderHamiltonian::set_mmats(AuxParamUVec& params_vec) {
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
    (*it)->minit=true;
  }
}

Eigen::MatrixXcd PairingLadderHamiltonian::get_reduced_matrix(const size_t& nele) {
  Eigen::MatrixXcd _evecs = _solver.eigenvectors();
  Eigen::MatrixXcd redmat(_size, nele);
  for(size_t i=0; i<_size; i++)
  for(size_t j=0; j<nele; j++)
    redmat(i,j)=_evecs(i, j);
  return redmat;
}

// ============================================================================= 

} // namespace AuxiliaryHamiltonians
} // namespace VMC
