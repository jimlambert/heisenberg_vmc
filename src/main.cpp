#include <cstdlib>
#include <random>
#include <time.h>
#include <iostream>
#include <fstream>
#include <Eigen/Dense>
#include "utils.h"
#include "SpinState.h"
#include "IsingChainEnergy.h"
#include "HeisenbergChainEnergy.h"
#include "Wavefunction.h"
#include "LadderWavefunction.h"
#include "HoppingChainHamiltonian.h"
#include "PairingChainHamiltonian.h"
#include "PairingLadderHamiltonian.h"
#include "KjgLadderEnergy.h"
#include "BasisState.h"
#include "SpinHilbert.h"
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
  size_t simul_steps=50000;
  size_t opt_steps=500;
  double df=0.01;
  VMC::ParamListSPtr par_lst_ptr=make_shared<VMC::Parameters::ParameterList>();
  std::random_device rd;
  std::mt19937 mteng(rd());
  std::uniform_real_distribution<double> rand_num(-1.0, 1.0);
  std::uniform_real_distribution<double> pos_rand_num(0.0, 1.0);
  std::uniform_real_distribution<double> neg_rand_num(-1.0, 0.0); 
  // ==========================================================================
  // Choose variational parameters
  // ==========================================================================
  
  // Ising chain model variational parameters
  par_lst_ptr->build_aux_param(ONSITE, "onsite", rand_num(mteng), 0, 1, true, 100);
  par_lst_ptr->build_aux_param(HOPPING, "hop1", rand_num(mteng), 0, 1, true, 100);
  par_lst_ptr->build_aux_param(HOPPING, "hop2", rand_num(mteng), 0, 2, true, 100);
  par_lst_ptr->build_aux_param(HOPPING, "hop3", rand_num(mteng), 0, 3, true, 100);
  par_lst_ptr->build_aux_param(HOPPING, "hop4", rand_num(mteng), 0, 4, true, 100);
  par_lst_ptr->build_aux_param(HOPPING, "hop5", rand_num(mteng), 0, 5, true, 100);
  par_lst_ptr->build_aux_param(PAIRING, "pair1", rand_num(mteng), 0, 1, true, 100);
  par_lst_ptr->build_aux_param(PAIRING, "pair2", rand_num(mteng), 0, 2, true, 100);
  par_lst_ptr->build_aux_param(PAIRING, "pair3", rand_num(mteng), 0, 3, true, 100);
  par_lst_ptr->build_aux_param(PAIRING, "pair4", rand_num(mteng), 0, 4, true, 100);
  par_lst_ptr->build_aux_param(PAIRING, "pair5", rand_num(mteng), 0, 5, true, 100);
  par_lst_ptr->build_aux_param(PAIRING, "pair6", rand_num(mteng), 0, 6, true, 100);

  // parameters where N=4 get's struck from TFIM
  //par_lst_ptr->build_aux_param(ONSITE, "onsite", -0.759878, 0, 1, true, 100);
  //par_lst_ptr->build_aux_param(HOPPING, "hop1", -0.0307682, 0, 1, true, 100);
  //par_lst_ptr->build_aux_param(PAIRING, "pair1", 0.448716, 0, 1, true, 100);
  //par_lst_ptr->build_aux_param(PAIRING, "pair2", -1.0158, 0, 2, true, 100);
  //par_lst_ptr->build_aux_param(PAIRING, "pair3", 0.447487, 0, 3, true, 100);
  
  // parameters start at E\approx -0.88 for TFIM
  //par_lst_ptr->build_aux_param(ONSITE, "onsite", 0.1, 0, 1, false, 100);
  //par_lst_ptr->build_aux_param(ONSITE, "onsite", 0.583116, 0, 1, true, 100);
  //par_lst_ptr->build_aux_param(HOPPING, "hop1", 0.589188, 0, 1, true, 100);
  //par_lst_ptr->build_aux_param(PAIRING, "pair3", -0.00321491, 0, 1, true, 100);
  
  // parameters start at E\approx 0.43 for TFIM
  //par_lst_ptr->build_aux_param(ONSITE, "onsite", 0.591611, 0, 1, true, 100);
  //par_lst_ptr->build_aux_param(HOPPING, "hop1", 0.580829, 0, 1, true, 100);
  //par_lst_ptr->build_aux_param(PAIRING, "pair3", 0.00584684, 0, 1, true, 100);
  
  // Heisenberg chain model variational parameters
  //par_lst_ptr->build_aux_param(HOPPING, "hopp1", rand_num(mteng), 0, 1, true, 100);
  //par_lst_ptr->build_aux_param(PAIRING, "pair1", rand_num(mteng), 0, 1, true, 100);
  //par_lst_ptr->build_aux_param(PAIRING, "pair2", rand_num(mteng), 0, 2, true, 100);
  //par_lst_ptr->build_aux_param(PAIRING, "pair3", rand_num(mteng), 0, 3, true, 100); 

  // KJG ladder model variational parameters
  //par_lst_ptr->build_aux_param(HOPPING, "hopp11", rand_num(mteng), 0, 0, 0, 1, true, 100);
  //par_lst_ptr->build_aux_param(HOPPING, "hopp12", rand_num(mteng), 1, 0, 1, 1, true, 100);
  //par_lst_ptr->build_aux_param(HOPPING, "hopp1r", rand_num(mteng), 0, 0, 1, 0, true, 100);
  //par_lst_ptr->build_aux_param(PAIRING, "pair11", rand_num(mteng), 0, 0, 0, 1, true, 100);
  //par_lst_ptr->build_aux_param(PAIRING, "pair12", rand_num(mteng), 0, 0, 0, 2, true, 100);
  //par_lst_ptr->build_aux_param(PAIRING, "pair13", rand_num(mteng), 0, 0, 0, 3, true, 100);
  //par_lst_ptr->build_aux_param(PAIRING, "pair21", rand_num(mteng), 1, 0, 1, 1, true, 100);
  //par_lst_ptr->build_aux_param(PAIRING, "pair22", rand_num(mteng), 1, 0, 1, 2, true, 100);
  //par_lst_ptr->build_aux_param(PAIRING, "pair23", rand_num(mteng), 1, 0, 1, 3, true, 100);
  //par_lst_ptr->build_aux_param(PAIRING, "pair1r", rand_num(mteng), 0, 0, 1, 0, true, 100);

  par_lst_ptr->build_jas_param(SPIN, "spin1", 0.1, 0, 1, true, 100);
  par_lst_ptr->build_jas_param(SPIN, "spin2", 0.1, 0, 2, true, 100);
  par_lst_ptr->build_jas_param(SPIN, "spin3", 0.1, 0, 3, true, 100);
  par_lst_ptr->build_jas_param(SPIN, "spin4", 0.1, 0, 4, true, 100);
  par_lst_ptr->build_jas_param(SPIN, "spin5", 0.1, 0, 5, true, 100);
  //par_lst_ptr->build_jas_param(SPIN, "spin2", pos_rand_num(mteng), 0, 2, true, 100);
  
  par_lst_ptr->report_aux_params();
  par_lst_ptr->report_jas_params(); 
  // ==========================================================================
  // Choose energy to minimize
  // ==========================================================================
  VMC::ObsUPtr enrg_ptr=make_unique<VMC::Observables::IsingChainEnergy>
                        ("Ising-Energy", 100, 1.0, 0.5); 
  //VMC::ObsUPtr enrg_ptr=make_unique<VMC::Observables::HeisenbergChainEnergy>
  //                      ("Heisenberg-Energy", 100, 0.5, 1.0, 1.0); 
  //VMC::ObsUPtr enrg_ptr=make_unique<VMC::Observables::KjgLadderEnergy>
  //                      ("KJG-Energy", 100, 1.0, 1.0, 0.0, 0.0); 
  // ==========================================================================
  // Choose observables
  // ==========================================================================
  
  // ==========================================================================
  // Choose auxiliary Hamiltonian
  // ==========================================================================
  VMC::AuxParamUVec test_noise;
  VMC::Utils::init_noise(2*L, test_noise);
  //VMC::AuxHamUPtr aux_ham_ptr=make_unique<VMC::HopChainHam>
  //                            (true, 2*L, par_lst_ptr->aux_vec());
  VMC::AuxHamUPtr aux_ham_ptr=make_unique<VMC::BCSChainHam>
                              (false, 2*L, par_lst_ptr->aux_vec(), 0.0, test_noise);
  //VMC::AuxHamUPtr aux_ham_ptr=make_unique<VMC::BCSLadderHam>
  //                            (false, 2*L, par_lst_ptr->aux_vec()); 

  //std::cout << aux_ham_ptr->get_matrix() << std::endl;
  //cout << "eigenvectors" << endl;
  //cout << aux_ham_ptr->get_eigenvectors() << endl;
  //cout << "eigenvalues" << endl;
  //cout << aux_ham_ptr->get_eigenvalues() << endl;
  
  // ==========================================================================
  // Wavefunction
  // ========================================================================== 

  //for(size_t i=0; i<par_lst_ptr->naux(); i++) {
  //  cout << (*par_lst_ptr).aux(i).name << endl;
  //  cout << (*par_lst_ptr).aux(i).mmat << endl;
  //}

  VMC::Wavefunctions::SpinWavefunction test_wave_func
  (
    L, 
    0,
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
