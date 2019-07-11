#include <cstdlib>
#include <random>
#include <time.h>
#include <iostream>
#include <fstream>
#include <Eigen/Dense>
#include "IsingChainEnergy.h"
#include "HeisenbergChainEnergy.h"
#include "Wavefunction.h"
#include "HoppingChainHamiltonian.h"
#include "PairingChainHamiltonian.h"
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

  size_t L=22;
  double df=0.1;
  VMC::ParamListSPtr par_lst_ptr=make_shared<VMC::Parameters::ParameterList>();

  std::random_device rd;
  std::mt19937 mteng(rd());
  std::uniform_real_distribution<double> rand_num(-1.0, 1.0);
  std::uniform_real_distribution<double> pos_rand_num(0.0, 1.0);
  std::uniform_real_distribution<double> neg_rand_num(-1.0, 0.0);
  
  //par_lst_ptr->build_aux_param(ONSITE, "onsite", rand_num(mteng), 0, L, true, 100);  
  //par_lst_ptr->build_aux_param(HOPPING, "hop1", rand_num(mteng), 0, 1, true, 100);  
  //par_lst_ptr->build_aux_param(HOPPING, "hop2", rand_num(mteng), 0, 2, true, 100);  
  //par_lst_ptr->build_aux_param(HOPPING, "hopL", rand_num(mteng), 0, L, true, 100);  
  //par_lst_ptr->build_aux_param(HOPPING, "hopL1", rand_num(mteng), 0, L+1, true, 100);  

  //for(size_t i=0; i<2*L; i++) {
  //  string name="onsite"+to_string(i);
  //  par_lst_ptr->build_aux_param(ONSITE, "onsite", rand_num(mteng), i, 0, false, 100);  
  //  
  //  string name_jas="spin"+to_string(i);
  //  if(i==0) continue;
  //  else if(i%2==0) 
  //    par_lst_ptr->build_jas_param(SPIN, name_jas, neg_rand_num(mteng), i, 1, false, 100);
  //  else  
  //    par_lst_ptr->build_jas_param(SPIN, name_jas, pos_rand_num(mteng), i, 1, false, 100);
  //}

  //par_lst_ptr->build_aux_param(ONSITE, "onsite", rand_num(mteng), 0, 0, false, 100);

  //for(size_t i=0; i<L; i++) {
  //  string name1="hop"+to_string(i);
  //  par_lst_ptr->build_aux_param(HOPPING, name1, rand_num(mteng), 0, i, true, 100);
  //
  //  string name2="pairing"+to_string(i);
  //  par_lst_ptr->build_aux_param(PAIRING, name2, rand_num(mteng), 0, i, true, 100);  
  //  //string name2="hopL"+to_string(i);
  //  //par_lst_ptr->build_aux_param(HOPPING, name2, rand_num(mteng), 0, L+i, true, 100);
  //}
  //par_lst_ptr->build_jas_param(SPIN, "spin1", pos_rand_num(mteng), 0, 1, true, 100);
 
  par_lst_ptr->build_aux_param(HOPPING, "hopp1", rand_num(mteng), 0, 1, true, 100);
  //par_lst_ptr->build_aux_param(HOPPING, "hopp2", rand_num(mteng), 0, 2, true, 100);
  //par_lst_ptr->build_aux_param(HOPPING, "hopp3", rand_num(mteng), 0, 3, true, 100);
  par_lst_ptr->build_aux_param(PAIRING, "pair1", rand_num(mteng), 0, 1, true, 100);
  par_lst_ptr->build_aux_param(PAIRING, "pair2", rand_num(mteng), 0, 2, true, 100);
  par_lst_ptr->build_aux_param(PAIRING, "pair3", rand_num(mteng), 0, 3, true, 100);


  //par_lst_ptr->build_jas_param(SPIN, "spin1", 0.1, 0, 1, true, 100);
  //par_lst_ptr->build_jas_param(SPIN, "spin1", 0.1, 0, 2, true, 100);
  //par_lst_ptr->build_jas_param(SPIN, "spin2", pos_rand_num(mteng), 0, 2, true, 100);
  
  par_lst_ptr->report_aux_params();
  //par_lst_ptr->report_jas_params(); 

  //VMC::ObsUPtr enrg_ptr=make_unique<VMC::Observables::IsingChainEnergy>
  //                      ("Ising-Energy", 100, 1.0, 0.0); 
  VMC::ObsUPtr enrg_ptr=make_unique<VMC::Observables::HeisenbergChainEnergy>
                        ("Heisenberg-Energy", 100, 1.0, 1.0, 1.0); 

  //VMC::AuxHamUPtr aux_ham_ptr=make_unique<VMC::HopChainHam>
  //                            (false, 2*L, par_lst_ptr->aux_vec());
  VMC::AuxHamUPtr aux_ham_ptr=make_unique<VMC::BCSChainHam>
                              (true, 2*L, par_lst_ptr->aux_vec());
  VMC::Wavefunctions::SpinWavefunction test_wave_func
  (
    L, 
    par_lst_ptr, 
    aux_ham_ptr, 
    enrg_ptr,
    "../vmc-dat/varparams",
    "../vmc-dat/n4obsvals"  
  ); 

  test_wave_func.optimize(1000, 5000, 30, df);
  //par_lst_ptr->report_aux_params();
  //par_lst_ptr->report_jas_params();
  
  
  return 0;
}
