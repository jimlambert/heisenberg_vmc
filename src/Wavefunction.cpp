#include "Wavefunction.h"

namespace VMC {
namespace Wavefunctions {

// Constructor
// =============================================================================
SpinWavefunction::SpinWavefunction (
  const size_t&     l,   // system size
  ParamListUPtr&    plp, // unique pointer to parameter list
  AuxHamUPtr&       ahp, // unique pointer to auxiliary Hamiltonian
  ObsUPtr&          erg, // unique pointer to energy functor
  const std::string& vn, // name of file for variational parameters
  const std::string& on  // name of file for observables
) : 
  _size(l), 
  _jas_sum(0.0),
  _par_lst_ptr(std::move(plp)),
  _aux_ham_ptr(std::move(ahp)),
  _state(l),
  _varfile_name(vn),
  _obsfile_name(on) {
 
  // initialize random site generator
  _rand_pos=new std::uniform_int_distribution<>(0, _size-1);

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
  
  // setup output file for variational parameters
  _setup_varfile();

  print_state();
  _flipspin(0);
  print_state();
}

// =============================================================================
// Private member functions
// =============================================================================

// output streams --------------------------------------------------------------
void SpinWavefunction::_setup_varfile() {
  _varfile_output.open(_varfile_name, std::ios::trunc);  
  for(size_t i=0; i<_par_lst_ptr->npar(); i++) {
    _varfile_output << std::setw(FILE_COLUMN_WIDTH) 
                    << std::left
                    << std::setfill(' ')
                    << (*_par_lst_ptr)[i].name;
  } 
  _varfile_output << '\n';
  for(size_t i=0; i<_par_lst_ptr->npar(); i++) {
    _varfile_output << std::setw(FILE_COLUMN_WIDTH) 
                    << std::left
                    << std::setfill(' ')
                    << std::setprecision(FILE_PRECISION)
                    << (*_par_lst_ptr)[i].val;
  }
  _varfile_output.close();
}

// setup functions -------------------------------------------------------------
void SpinWavefunction::_init_state() {
  _state.refresh();  
  for(size_t i=0; i<_size; i++) {
    double rand_num=_rand_num(_mteng);
    if(rand_num<0.5) {
      _state[i]=1;
      _state(i)=i+1;
    }
    else {
      _state[i]=-1;
      _state(i+_size)=i+1;
    }
  } 
}


void SpinWavefunction::_reinit_gmat() {
  Eigen::MatrixXcd projmat(_size, _size);
  Eigen::MatrixXcd redmat(2*_size, _size);
  redmat=_aux_ham_ptr->get_reduced_matrix(_size);
  for(size_t ri=0; ri<_size; ri++) {
    int expos;
    if(_state[ri]==1) expos=ri;
    else expos=ri+_size;
    for(size_t rj=0; rj<_size; rj++) {
      projmat(ri, rj) = redmat(expos, rj);
    }
  }
  _gmat = redmat*projmat.inverse();
}


void SpinWavefunction::_clear_obs() {
  for(auto it=_obsvec.begin(); it!=_obsvec.end(); it++) {
    (*it)->local_meas.clear();
  }
}


void SpinWavefunction::_compute_jas_sum() {
  // compute initial value of Jastrow sum
  for(size_t i=0; i<_par_lst_ptr->njas(); i++) {
    size_t dr=_par_lst_ptr->jas(i).dr;
    std::complex<double> v=_par_lst_ptr->jas(i).val;
    for(size_t r=0; r<_size; r++) {
      _jas_sum+=0.25*v.real()*_state[r]*_state[(r+dr)%_size]; 
    }
  }
}


double SpinWavefunction::_dj_sum(const size_t& r, const int& ds) {
  double total=0.0;
  for(size_t i=0; i<_par_lst_ptr->njas(); i++) {
    size_t dr=_par_lst_ptr->jas(i).dr;
    std::complex<double> v=_par_lst_ptr->jas(i).val;
    size_t rl=(r-dr+_size)%_size;
    size_t rr=(r+dr)%_size;
    total+=0.25*v.real()*(double)(_state[rl]+_state[rr])*ds;
  }
  return total;
}

// update functions ------------------------------------------------------------
void SpinWavefunction::_flipspin(const int& r) {
  size_t iexpos, nexpos, lindex;
  int ds;
  if(_state[r]==1) {
    iexpos=r;
    nexpos=iexpos+_size;
    ds=-2;
  }
  else {
    iexpos=r+_size;
    nexpos=r;
    ds=2;
  }
  lindex=_state(iexpos)-1;
  double dj=_dj_sum(r, ds);
  double amp=fabs(std::exp(dj)*_gmat(nexpos, lindex));
  double rand_num=_rand_num(_mteng);
  if(rand_num<amp*amp) {
    _state(nexpos)=_state(iexpos);
    _state(iexpos)=0;
    _state[r]+=ds;
    Eigen::MatrixXcd newgmat(2*_size, _size);
    for(size_t i=0; i<2*_size; i++)
    for(size_t j=0; j<_size; j++) {
      double d=0;
      if(j==lindex) d=0;
      std::complex<double> den=_gmat(nexpos, lindex);
      std::complex<double> num=_gmat(i, lindex);
      newgmat(i,j)=_gmat(i,j)-(num/den)*(_gmat(nexpos,j)-d);
    }
    _gmat=newgmat;
    _jas_sum=_jas_sum+dj;
  }
}


void SpinWavefunction::_sweep() {
  for(size_t i=0; i<2*_size; i++) {
    int rand_pos=(*_rand_pos)(_mteng);
    _flipspin(rand_pos);   
  }
  // we recompute these values every sweep just to avoid roundoff errors. This
  // may be excessive and it might instead be best to compute these fewer times.
  _reinit_gmat();
  _compute_jas_sum(); 
}


void SpinWavefunction::_update_params(const double&) {

}

// =============================================================================
// Public members functions
// =============================================================================

void SpinWavefunction::optimize(
  const size_t& equil_steps, 
  const size_t& simul_steps, 
  const size_t& opt_steps,
  const double& delta
) {
  for(size_t opt=0; opt<opt_steps; opt++) {
    _clear_obs();
    for(size_t equil=0; equil<equil_steps; equil++) _sweep();  
    for(size_t simul=0; simul<simul_steps; simul++) {
      _sweep(); 
      // measure expectation values for variational parameters
      for(size_t i=0; i<_par_lst_ptr->npar(); i++) 
        (*_par_lst_ptr)[i](_state, _gmat);
      // measure local  expectation value of the energy
      (*_obsvec[0])(_state, _gmat);   
    } // simulation loop
  } // optimization loop
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

} // namespace Wavefunction
} // namespace VMC
