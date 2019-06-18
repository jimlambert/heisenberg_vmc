#include "Wavefunction.h"

namespace VMC {
namespace Wavefunctions {

SpinWavefunction::SpinWavefunction (
  const size_t&     l,   // system size
  const double&     df,  // delta in the variational update
  ParamListUPtr&    plp, // unique pointer to parameter list
  AuxHamUPtr&       ahp, // unique pointer to auxiliary Hamiltonian
  ObsUPtr&          erg, // unique pointer to energy functor
  const std::string& vn, // name of file for variational parameters
  const std::string& on  // name of file for observables
) : 
  _size(l), 
  _delta(df),
  _jas_sum(0.0),
  _par_lst_ptr(std::move(plp)),
  _aux_ham_ptr(std::move(ahp)),
  _state(l),
  _var_param_file_name(vn),
  _obs_file_name(on) {
  
  // add energy to observable vector
  _obsvec.push_back(std::move(erg));
  
  // allocate reduced matrix and create initial state
  Eigen::MatrixXcd redmat(2*_size,_size);
  redmat=_aux_ham_ptr->get_reduced_matrix(_size);
  Eigen::MatrixXcd projmat(_size,_size);
  do {
   _init_state();
   for(size_t i=0; i<_size; i++) {
     int pos;
     if(_state[i]==1) pos=i;
     else pos=i+_size;
     for(size_t j=0; j<_size; j++) {
       projmat(i, j) = redmat(pos, j);
     }
   }
  } while(fabs(projmat.determinant())<1e-12);

  // build Green's function matrix
  _gmat=redmat*(projmat.inverse());
  _compute_jas_sum();
}

// setup functions -------------------------------------------------------------

void SpinWavefunction::_init_state() {
  _state.refresh();  
  for(size_t i=0; i<_size; i++) {
    double r=_rnum(_mteng);
    if(r<0.5) {
      _state[i]=1;
      _state(i)=i+1;
    }
    else {
      _state[i]=-1;
      _state(i+_size)=i+1;
    }
  } 
}

void SpinWavefunction::_compute_jas_sum() {
  // compute initial value of Jastrow sum
  for(size_t i=0; i<_par_lst_ptr->njas(); i++) {
    size_t dr=_par_lst_ptr->jas(i).dr;
    std::complex<double> v=_par_lst_ptr->jas(i).val;
    for(size_t r=0; r<_size; r++) {
      _jas_sum+=v.real()*_state[r]*_state[(r+dr)%_size]; 
    }
  }
  std::cout << _jas_sum << std::endl;
  print_state(); 
}

// update functions ------------------------------------------------------------

void SpinWavefunction::_flipspin(const size_t& r) {
  size_t iexpos; // initial position in extended basis
  size_t nexpos; // new position in extended basis
  size_t lindex; // position of fermion operator
  int    ds;     // change in spin
  if(_state[r]==1) {
    iexpos=r;
    nexpos=iexpos+_size;
    ds=-2;
  }
  else {
    iexpos=r+_size;
    nexpos=iexpos-_size;
    ds=2;
  }
}

// output functions ------------------------------------------------------------

void SpinWavefunction::print_state() {
  std::cout << "Spin state:" << std::endl;
  for(size_t i=0; i<_size; i++) 
   std::cout << std::setw(STATE_WIDTH) << std::left << std::setfill(' ')
             << _state[i];
  std::cout << std::endl; 
  std::cout << "Operator list:" << std::endl;
  for(size_t i=0; i<2*_size; i++) 
    std::cout << std::setw(STATE_WIDTH) << std::left << std::setfill(' ')
              << _state(i);
  std::cout << std::endl;

}

} // namespace VMC
} // namespace Wavefunction
