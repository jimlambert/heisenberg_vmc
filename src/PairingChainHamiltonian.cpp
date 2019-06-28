
#include <iostream>
#include <Eigen/Dense>
#include "AuxiliaryParameter.h"
#include "PairingChainHamiltonian.h"

namespace VMC {
namespace AuxiliaryHamiltonians {

// Private Members
// ============================================================================= 
  
void PairingChainHamiltonian::_check_vmat_init(const AuxParamUPtr& aux_ptr) {
  if(!(aux_ptr->vinit)) {
    std::cout << "Parameter " << aux_ptr->name << " has unitialized vmatrix" 
              << std::endl;
    exit(1);
  }
}

// ============================================================================= 
// Public Members
// ============================================================================= 

PairingChainHamiltonian::PairingChainHamiltonian(
  const bool&          b,
  const size_t&        n,
  AuxParamUVec&  params_vec
) : _bc(b), _size(n) {
  set_vmats(params_vec);
  init(params_vec);
}

void PairingChainHamiltonian::init(AuxParamUVec& params_vec) {
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


void PairingChainHamiltonian::set_vmats(AuxParamUVec& params_vec) {
  for(auto it=params_vec.begin(); it!=params_vec.end(); it++) {
    bool ti=(*it)->trans_inv;
    (*it)->vmat=Eigen::MatrixXd::Zero(_size,_size);
    switch((*it)->subtype) {
      case Parameters::ONSITE: 
      {
        if(ti) for(size_t i=0; i<_size; i++) {
          if(i<(_size/2)) (*it)->vmat(i,i)+=-1;
          else            (*it)->vmat(i,i)+=1;
        }
        else {
          size_t r=(*it)->site;
          if(r<(_size/2)) (*it)->vmat(r,r)+=-1;
          else        (*it)->vmat(r,r)+=1;
        }
      }
      break; // case ONSITE
      case Parameters::HOPPING:
      {
        size_t r=(*it)->site;
        size_t dr=(*it)->dr;
        if(ti) {
          for(size_t i=0; i<_size/2; i++) {
            double bfactor=1.0;
            size_t j=(i+dr)%(_size/2);
            if(((i+dr)>=(_size/2)) && !_bc) bfactor=-1.0; 
            (*it)->vmat(i,j)+=bfactor;
            (*it)->vmat(j,i)+=bfactor;
          }
          for(size_t i=(_size/2); i<_size; i++) {
            double bfactor=1.0;
            size_t j=((i+dr)%_size)+(_size/2)*(size_t)((i+dr)/(_size));
            if(((i+dr)>=_size) && !_bc) bfactor=-1.0; 
            (*it)->vmat(i,j)+=-1*bfactor;
            (*it)->vmat(j,i)+=-1*bfactor;
          }
        }
        else {
          double bfactor=1.0;
          if(r<_size/2) {
            size_t l=(r+dr)%_size/2;
            if(((r+dr)>=_size/2) && !_bc) bfactor=-1.0;
            (*it)->vmat(r,l)+=bfactor;
            (*it)->vmat(l,r)+=bfactor;
          }
          else {
            size_t l=(r+dr)%_size+(_size/2)*(size_t)((r+dr)/(_size));
            if(((r+dr)>=_size) && !_bc) bfactor=-1.0;
            (*it)->vmat(r,l)+=-1*bfactor;
            (*it)->vmat(l,r)+=-1*bfactor;
          } 
        }
      }
      break; // case Hopping
      case Parameters::PAIRING:
      {
        size_t r=(*it)->site;
        size_t dr=(*it)->dr;
        if(ti) {
          for(size_t i=0; i<_size/2; i++) {
            double bfactor=1.0;
            size_t j=((i+dr)%(_size/2))+_size/2;
            if(((i+dr)>=(_size/2)) && !_bc) bfactor=-1;
            (*it)->vmat(i,j)+=1*bfactor;
            (*it)->vmat(j,i)+=1*bfactor;
          }
        }
        else {
          double bfactor=1.0;
          size_t l=((r+dr)%(_size/2))+_size/2;
          if(((r+dr)>=(_size/2)) && !_bc) bfactor=-1;
          (*it)->vmat(r,l)+=1*bfactor;
          (*it)->vmat(l,r)+=1*bfactor;
        } 
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

void PairingChainHamiltonian::set_mmats(AuxParamUVec& params_vec) {
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

Eigen::MatrixXcd PairingChainHamiltonian::get_reduced_matrix(const size_t& nele) {
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
