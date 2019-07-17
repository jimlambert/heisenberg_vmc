#include <cstdlib>
#include <random>
#include <time.h>
#include <iostream>
#include <fstream>
#include <Eigen/Dense>
#include "IsingChainEnergy.h"
#include "HeisenbergChainEnergy.h"
#include "Wavefunction.h"
#include "LadderWavefunction.h"
#include "HoppingChainHamiltonian.h"
#include "PairingChainHamiltonian.h"
#include "PairingLadderHamiltonian.h"
#include "KjgLadderEnergy.h"
#include "BasisState.h"
#include "ParameterList.h"
#include "SpinSpin.h"
#include "LocalMeasurement.h"

using namespace std;

using VMC::Parameters::ONSITE;
using VMC::Parameters::HOPPING;
using VMC::Parameters::PAIRING;
using VMC::Parameters::SPIN;

int main(int argc, char* argv[]) {

  // Initial setup
  // ==========================================================================
  size_t L=2;
  size_t equil_steps=5000;
  size_t simul_steps=10000;
  size_t opt_steps=30;
  double df=0.1;
  VMC::ParamListSPtr par_lst_ptr=make_shared<VMC::Parameters::ParameterList>();
  std::random_device rd;
  std::mt19937 mteng(rd());
  std::uniform_real_distribution<double> rand_num(-1.0, 1.0);
  std::uniform_real_distribution<double> pos_rand_num(0.0, 1.0);
  std::uniform_real_distribution<double> neg_rand_num(-1.0, 0.0); 
  // ==========================================================================
  // Choose variational parameters
  // ==========================================================================
  
  // Heisenberg chain model variational parameters
  //par_lst_ptr->build_aux_param(HOPPING, "hopp1", rand_num(mteng), 0, 1, true, 100);
  //par_lst_ptr->build_aux_param(PAIRING, "pair1", rand_num(mteng), 0, 1, true, 100);
  //par_lst_ptr->build_aux_param(PAIRING, "pair2", rand_num(mteng), 0, 2, true, 100);
  //par_lst_ptr->build_aux_param(PAIRING, "pair3", rand_num(mteng), 0, 3, true, 100);

  // KJG ladder model variational parameters
  par_lst_ptr->build_aux_param(HOPPING, "hopp11", rand_num(mteng), 0, 0, 0, 1, true, 100);
  par_lst_ptr->build_aux_param(HOPPING, "hopp12", rand_num(mteng), 1, 0, 1, 1, true, 100);
  par_lst_ptr->build_aux_param(PAIRING, "pair11", rand_num(mteng), 0, 0, 0, 1, true, 100);
  //par_lst_ptr->build_aux_param(PAIRING, "pair12", rand_num(mteng), 0, 0, 0, 2, true, 100);
  //par_lst_ptr->build_aux_param(PAIRING, "pair13", rand_num(mteng), 0, 0, 0, 3, true, 100);
  par_lst_ptr->build_aux_param(PAIRING, "pair21", rand_num(mteng), 1, 0, 1, 1, true, 100);
  //par_lst_ptr->build_aux_param(PAIRING, "pair22", rand_num(mteng), 1, 0, 1, 2, true, 100);
  //par_lst_ptr->build_aux_param(PAIRING, "pair23", rand_num(mteng), 1, 0, 1, 3, true, 100);

  //par_lst_ptr->build_jas_param(SPIN, "spin1", 0.1, 0, 1, true, 100);
  //par_lst_ptr->build_jas_param(SPIN, "spin1", 0.1, 0, 2, true, 100);
  //par_lst_ptr->build_jas_param(SPIN, "spin2", pos_rand_num(mteng), 0, 2, true, 100);
  
  par_lst_ptr->report_aux_params();
  //par_lst_ptr->report_jas_params(); 
  // ==========================================================================
  // Choose energy to minimize
  // ==========================================================================
  //VMC::ObsUPtr enrg_ptr=make_unique<VMC::Observables::IsingChainEnergy>
  //                      ("Ising-Energy", 100, 1.0, 0.0); 
  //VMC::ObsUPtr enrg_ptr=make_unique<VMC::Observables::HeisenbergChainEnergy>
  //                      ("Heisenberg-Energy", 100, 1.0, 1.0, 1.0); 
  VMC::ObsUPtr enrg_ptr=make_unique<VMC::Observables::KjgLadderEnergy>
                        ("KJG-Energy", 100, 1.0, 0.0, 0.0, 0.0); 
  // ==========================================================================
  // Choose observables
  // ==========================================================================
  
  // ==========================================================================
  // Choose auxiliary Hamiltonian
  // ==========================================================================
  //VMC::AuxHamUPtr aux_ham_ptr=make_unique<VMC::HopChainHam>
  //                            (false, 2*L, par_lst_ptr->aux_vec());
  //VMC::AuxHamUPtr aux_ham_ptr=make_unique<VMC::BCSChainHam>
  //                            (true, 2*L, par_lst_ptr->aux_vec());
  VMC::AuxHamUPtr aux_ham_ptr=make_unique<VMC::BCSLadderHam>
                              (true, 2*L, par_lst_ptr->aux_vec()); 

  //cout << aux_ham_ptr->get_matrix() << endl;
  // ==========================================================================
  // Wavefunction
  // ========================================================================== 
 
  VMC::Wavefunctions::LadderWavefunction test_wave_func
  (
    L, 
    par_lst_ptr, 
    aux_ham_ptr, 
    enrg_ptr,
    "../vmc-dat/varparams",
    "../vmc-dat/n4obsvals"  
  ); 

  test_wave_func.optimize(equil_steps, simul_steps, opt_steps, df);
  //par_lst_ptr->report_aux_params();
  //par_lst_ptr->report_jas_params();
  
  
  return 0;
}
