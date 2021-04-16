#include <iostream>
#include <Eigen/Dense>
#include <cmath>
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


void PairingLadderHamiltonian::_set_onsite_vmat(const AuxParamUPtr& it) {
  bool   ti=(it)->trans_inv;
  if(ti) for(size_t i=0; i<_size; i++) {
    if(i<(_size/2)) (it)->vmat(i,i)+=-1;
    else            (it)->vmat(i,i)+=1;
  }
  else {
    size_t r=_pair2index((it)->site1[0], (it)->site1[1]);
    if(r<(_size/2)) (it)->vmat(r,r)+=-1;
    else            (it)->vmat(r,r)+=1;
  }
}


void PairingLadderHamiltonian::_set_hopping_vmat(const AuxParamUPtr& it) {
  size_t c1=(it)->site1[0];
  size_t s1=(it)->site1[1];
  size_t c2=(it)->site2[0];
  size_t s2=(it)->site2[1];
  //int ds=abs(s1-s2);
  int ds;
  if(s1>s2) ds=s1-s2;
  else ds=s2-s1;
  
  bool   ti=(it)->trans_inv;
  if(ti) {  // translation invariant case
    // Up spin loop
   for(size_t i=0; i<_nrungs/2; i++) {
      double bfactor=1.0;
      size_t j=(i+ds)%(_nrungs/2);
      if(((i+ds)>=(_nrungs/2)) && !_bc) bfactor=-1.0;
      if(c1!=c2) bfactor=1.0;
      size_t n=_pair2index(c1, i);
      size_t m=_pair2index(c2, j);
      (it)->vmat(n,m)+=bfactor;
      (it)->vmat(m,n)+=bfactor;
    } 
    // Down spin loop
    for(size_t i=(_nrungs/2); i<_nrungs; i++) {
      double bfactor=1.0;
      size_t j=((i+ds)%_nrungs)+(_nrungs/2)*(size_t)((i+ds)/_nrungs);
      if(((i+ds)>=_nrungs) && !_bc) bfactor=-1.0;
      if(c1!=c2) bfactor=1.0;
      size_t n=_pair2index(c1, i);
      size_t m=_pair2index(c2, j);
      (it)->vmat(n,m)+=-1*bfactor;
      (it)->vmat(m,n)+=-1*bfactor; 
    }
  }  
  else { // non-translation invariant case
    double bfactor=1.0;
    if(s1<_nrungs/2) { // up spin
      if((((int)s2-(int)s1)<0) && !_bc) bfactor=-1.0;
      if(c1!=c2) bfactor=1.0;
      size_t n=_pair2index(c1, s1);
      size_t m=_pair2index(c2, s2);
      (it)->vmat(n,m)+=bfactor;
      (it)->vmat(m,n)+=bfactor;  
    }
    else { // down spin
      if((((int)s2-(int)s1)<0) && !_bc) bfactor=-1.0;
      if(c1!=c2) bfactor=1.0;
      size_t n=_pair2index(c1, s1);
      size_t m=_pair2index(c2, s2);
      (it)->vmat(n,m)+=-1*bfactor;
      (it)->vmat(m,n)+=-1*bfactor;  
    } 
  }
}


void PairingLadderHamiltonian::_set_pairing_vmat(const AuxParamUPtr& it) {
  size_t c1=(it)->site1[0];
  size_t s1=(it)->site1[1];
  size_t c2=(it)->site2[0];
  size_t s2=(it)->site2[1];
  //int ds=abs(s1-s2);
  int ds;
  if(s1>s2) ds=s1-s2;
  else ds=s2-s1;
  bool   ti=(it)->trans_inv;
  if(ti) { // translation invariant
    for(size_t i=0; i<(_nrungs/2); i++) {
      double bfactor=1.0;
      size_t j=((i+ds)%(_nrungs/2))+_nrungs/2;
      if(((i+ds)>=(_nrungs/2)) && !_bc) bfactor=-1;
      if(c1==c2) bfactor=1;
      size_t n=_pair2index(c1, i);
      size_t m=_pair2index(c2, j);
      (it)->vmat(n,m)+=bfactor;
      (it)->vmat(m,n)+=bfactor;
    }
  }  
  else { // non-translation invariant
    double bfactor=1.0;
    if((((int)s2-(int)s1)<0) && !_bc) bfactor=-1.0;
    if(c1==c2) bfactor=1;
    size_t n=_pair2index(c1, s1);
    size_t m=_pair2index(c2, s2);
    (it)->vmat(n,m)+=bfactor;
    (it)->vmat(m,n)+=bfactor;
  }
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
    (*it)->vmat=Eigen::MatrixXd::Zero(_size,_size);
    switch((*it)->subtype) {
      case Parameters::ONSITE:
      {_set_onsite_vmat(*it);}
      break; // case ONSITE
      case Parameters::HOPPING:
      {_set_hopping_vmat(*it);}
      break; // case Hopping
      case Parameters::PAIRING:
      {_set_pairing_vmat(*it);}
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
