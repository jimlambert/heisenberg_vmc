#include <cmath>
#include "Wavefunction.h"

namespace VMC {
namespace Wavefunctions {

// Constructor
// =============================================================================
SpinWavefunction::SpinWavefunction (
  const size_t&     l,   // system size
  ParamListSPtr&    plp, // unique pointer to parameter list
  AuxHamUPtr&       ahp, // unique pointer to auxiliary Hamiltonian
  ObsUPtr&          erg, // unique pointer to energy functor
  const std::string& vn, // name of file for variational parameters
  const std::string& on  // name of file for observables
) : 
  _sfac(1),
  _size(l), 
  _jas_sum(0.0),
  _par_lst_ptr(plp),
  _aux_ham_ptr(std::move(ahp)),
  _state(l),
  _varfile_name(vn),
  _obsfile_name(on) {
 
  // initialize random site generator
  _rand_pos=new std::uniform_int_distribution<>(0, _size-1);

  // add energy to observable vector
  _obsvec.push_back(std::move(erg));
  
  
  // setup output file for variational parameters
  _setup_varfile();
  
  //Eigen::MatrixXcd redmat(2*_size,_size);
  //redmat=_aux_ham_ptr->get_reduced_matrix(_size);
  //Eigen::MatrixXcd projmat(_size,_size);
  //do {
  // _init_state();
  // for(size_t i=0; i<_size; i++) {
  //   size_t pos=_state.find(i+1);
  //   //if(_state[i]==1) pos=i;
  //   //else pos=i+_size;
  //   for(size_t j=0; j<_size; j++) {
  //     projmat(i, j) = redmat(pos, j);
  //   }
  // }
  //} while(fabs(projmat.determinant())<1e-12);

  //// build Green's function matrix
  //_gmat=redmat*(projmat.inverse());
  //_compute_jas_sum();
  //std::cout << _jas_sum << std::endl;
  //print_state();
  //(*_obsvec[0])(_state, _gmat, _par_lst_ptr->jas_vec()); 
}

// =============================================================================
// Private member functions
// =============================================================================

// output streams --------------------------------------------------------------
void SpinWavefunction::_setup_varfile() {
  _varfile_name=_varfile_name+".dat";
  _varfile_output.open(_varfile_name, std::ofstream::out | std::ios::trunc);  
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
                    << (*_par_lst_ptr)[i].val.real();
  }
  _varfile_output << '\n';
  _varfile_output.close();
}


void SpinWavefunction::_update_varfile() {
  _varfile_output.open(_varfile_name, std::ios::app);
  for(size_t i=0; i<_par_lst_ptr->npar(); i++) {
    _varfile_output << std::setw(FILE_COLUMN_WIDTH) 
                    << std::left
                    << std::setfill(' ')
                    << std::setprecision(FILE_PRECISION)
                    << (*_par_lst_ptr)[i].val.real();
  }
  _varfile_output << '\n';
  _varfile_output.close(); 
}


void SpinWavefunction::_update_obsfile(const size_t& i) {
  std::string curr_obs=_obsfile_name+std::to_string(i)+".dat";
  size_t nbins=(*_obsvec[0]).local_meas.nbins();
  _obsfile_output.open(curr_obs, std::ofstream::out | std::ios::trunc);
  _obsfile_output << i << '\n';
  for(auto it=_obsvec.begin(); it!=_obsvec.end(); it++) {
    _obsfile_output << std::setw(FILE_COLUMN_WIDTH)
                    << std::left 
                    << std::setfill(' ')
                    << (*it)->name;
  }
  for(size_t p=0; p<_par_lst_ptr->npar(); p++) {
    _obsfile_output << std::setw(FILE_COLUMN_WIDTH)
                    << std::left
                    << std::setfill(' ')
                    << (*_par_lst_ptr)[p].name;
  }
  _obsfile_output << '\n';
  for(size_t i=0; i<nbins; i++) {
    for(size_t b=0; b<_obsvec.size(); b++) {
      _obsfile_output << std::setw(FILE_COLUMN_WIDTH)
                      << std::left 
                      << std::setfill(' ')
                      << (*_obsvec[b]).local_meas(i); 
    }
    for(size_t p=0; p<_par_lst_ptr->npar(); p++) {
      _obsfile_output << std::setw(FILE_COLUMN_WIDTH)
                      << std::left 
                      << std::setfill(' ')
                      << (*_par_lst_ptr)[p].local_meas(i);  
    }
    _obsfile_output << '\n';
  }
  _obsfile_output.close();
  _obsfile_output.clear();
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
  for(size_t i=0; i<(*_par_lst_ptr).npar(); i++) {
    (*_par_lst_ptr)[i].local_meas.clear();
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
      if(j==lindex) d=1;
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
  _reinit_gmat();
  _compute_jas_sum(); 
}


void SpinWavefunction::_update_params(const double& delta) {
  
  size_t npar=_par_lst_ptr->npar();
  Eigen::MatrixXd s(npar,npar);
  Eigen::VectorXd f(npar);
  s=Eigen::MatrixXd::Zero(npar,npar);
  
  // computer S matrix and force vector
  for(size_t ka=0; ka<npar; ka++) {
    size_t na=(*_par_lst_ptr)[ka].local_meas.nvals();
    for(size_t kb=ka; kb<npar; kb++) {
      std::complex<double> avea=(*_par_lst_ptr)[ka].local_meas.ave();
      std::complex<double> aveb=(*_par_lst_ptr)[kb].local_meas.ave();
      std::complex<double> total=0.0;
      for(size_t i=0; i<na; i++) {
        std::complex<double> ai=(*_par_lst_ptr)[ka].local_meas[i];
        std::complex<double> bi=(*_par_lst_ptr)[kb].local_meas[i];
        total+=(ai-avea) * (bi-aveb);
      }
      s(ka,kb)=(total/(double)na).real();
      s(kb,ka)=s(ka,kb);
    } // kb loop
    std::complex<double> total=0.0;
    std::complex<double> avea=(*_par_lst_ptr)[ka].local_meas.ave();
    std::complex<double> avee=(*_obsvec[0]).local_meas.ave();
    for(size_t i=0; i<na; i++) {
      std::complex<double> ai=(*_par_lst_ptr)[ka].local_meas[i];
      // total += std::conj((*_obsvec[0]).local_meas[i]-avee)*(ai-avea);
      total += std::conj((*_obsvec[0]).local_meas[i])*(ai-avea);
    }
    f(ka)=-2.0*total.real()/(double)na;
  } // ka loop 

  // FULL SOLUTION
  Eigen::MatrixXd s_alt(npar, npar);
  Eigen::VectorXd f_alt(npar); 
  //for(size_t ka=0; ka<npar; ka++) {
  //  for(size_t kb=0; kb<npar; kb++) {
  //    s_alt(ka,kb)=s(ka,kb)/(std::sqrt(s(ka,ka)*s(kb,kb)));
  //  }
  //  s_alt(ka,ka)=1;
  //  f_alt(ka)=f(ka)/std::sqrt(s(ka,ka));
  //  if(f(ka)<1e-5) f_alt(ka)=0;
  //}
  for(size_t i=0; i<npar; i++) s(i,i)+=1e-3;
  Eigen::VectorXd da(npar);
  std::cout << "dE:" << '\t' 
            << -1.0*delta*f.transpose()*s.inverse()*f
            << std::endl;
  da=s.colPivHouseholderQr().solve(delta*f);
  for(size_t k=0; k<npar; k++) 
    (*_par_lst_ptr)[k].val+=da[k];
  std::cout << "dA:" << std::endl;
  std::cout << "---" << std::endl; 
  std::cout << da << std::endl;
  std::cout << "---" << std::endl; 
  // SUBSPACE SOLUTION
  // determine which indices to keep
  //Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solve_s;
  //solve_s.compute(s);
  //Eigen::VectorXd e_vals=solve_s.eigenvalues();
  //Eigen::MatrixXd e_vecs=solve_s.eigenvectors();
  //std::vector<size_t> keep_indices;
  //for(size_t i=0; i<npar; i++) if(fabs(e_vals(i))>1e-3) keep_indices.push_back(i);
  //size_t nnew=keep_indices.size();
  //if(nnew>0) {
  //  // modify parameters above threshold
  //  //s=(e_vecs.adjoint())*s*(e_vecs);
  //  s=Eigen::MatrixXd::Zero(npar,npar);
  //  for(size_t i=0; i<npar; i++) {
  //    s(i,i)=e_vals(i);
  //  }
  //  f=e_vecs*f;
  //  Eigen::MatrixXd s_new(nnew, nnew);
  //  Eigen::VectorXd f_new(nnew);
  //  Eigen::VectorXd dA_trans(nnew);
  //  Eigen::VectorXd dA(npar);
  //  for(size_t i=0; i<nnew; i++) {
  //    size_t i_index=keep_indices[i];
  //    for(size_t j=0; j<nnew; j++) {
  //      size_t j_index=keep_indices[j];
  //      s_new(i,j)=s(i_index,j_index);
  //    }
  //    f_new(i)=f(i_index);
  //  }
  //  // try preconditioning reduced matrix
  //  //Eigen::MatrixXd s_pc(nnew,nnew);
  //  //Eigen::VectorXd f_pc(nnew);
  //  //for(size_t ka=0; ka<nnew; ka++) {
  //  //  f_pc(ka)=f_new(ka)/std::sqrt(s_new(ka,ka));
  //  //  for(size_t kb=0; kb<nnew; kb++) {
  //  //    s_pc(ka,kb)=s_new(ka,kb)/(std::sqrt(s_new(ka,ka)*s_new(kb,kb)));
  //  //  }
  //  //}
  //  dA_trans=s_new.colPivHouseholderQr().solve(delta*f_new);
  //  for(size_t i=0; i<npar; i++) dA(i)=0;
  //  for(size_t i=0; i<nnew; i++) dA(keep_indices[i])=dA_trans(i);
  //  dA=e_vecs.adjoint()*dA;
  //  for(size_t k=0; k<npar; k++) (*_par_lst_ptr)[k].val+=dA[k];
  //} 
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
    // allocate reduced matrix and create initial state
    Eigen::MatrixXcd redmat(2*_size,_size);
    redmat=_aux_ham_ptr->get_reduced_matrix(_size);
    Eigen::MatrixXcd projmat(_size,_size);
    do {
     _init_state();
     for(size_t i=0; i<_size; i++) {
       size_t pos=_state.find(i+1);
       //if(_state[i]==1) pos=i;
       //else pos=i+_size;
       for(size_t j=0; j<_size; j++) {
         projmat(i, j) = redmat(pos, j);
       }
     }
    } while(fabs(projmat.determinant())<1e-12);
    // build Green's function matrix
    _gmat=redmat*(projmat.inverse());
    _compute_jas_sum();
    _clear_obs();
    for(size_t equil=0; equil<equil_steps; equil++) _sweep();  
    for(size_t simul=0; simul<simul_steps; simul++) {
      _sweep();
      //print_spins();
      // measure expectation values for variational parameters
      for(size_t i=0; i<_par_lst_ptr->npar(); i++) 
        (*_par_lst_ptr)[i](_state, _gmat);
      // measure local  expectation value of the energy
      (*_obsvec[0])(_state, _gmat, _par_lst_ptr->jas_vec());   
    } // simulation loop
    _update_params(delta);
    _aux_ham_ptr->init(_par_lst_ptr->aux_vec());
    _update_varfile();
    _update_obsfile(opt);
    std::cout << std::left 
              << "variational step: " 
              << opt 
              << " complete"
              << std::endl;
    std::cout << "energy: " 
              << (*_obsvec[0]).local_meas.ave() 
              << std::endl;
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


void SpinWavefunction::print_spins() {
  for(size_t i=0; i<_size; i++) 
   std::cout << std::setw(STATE_WIDTH) << std::left << std::setfill(' ')
             << _state[i];
  std::cout << std::endl;
}

} // namespace Wavefunction
} // namespace VMC
