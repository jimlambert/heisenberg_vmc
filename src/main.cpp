#include <cstdlib>
#include <random>
#include <time.h>
#include <iostream>
#include <fstream>
#include <Eigen/Dense>
#include "VarParam.h"
#include "IsingChainEnergy.h"
#include "HeisenbergChainEnergy.h"
#include "Wavefunction.h"
#include "HoppingChainHamiltonian.h"
#include "PairingChainHamiltonian.h"
#include "BasisState.h"
#include "ParameterList.h"
#include "SpinSpin.h"
#include "Simulation.h"
#include "LocalMeasurement.h"
#include "BCSHamiltonian.h"

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

  par_lst_ptr->build_aux_param(ONSITE, "onsite", rand_num(mteng), 0, 0, false, 100);

  for(size_t i=0; i<L; i++) {
    string name1="hop"+to_string(i);
    par_lst_ptr->build_aux_param(HOPPING, name1, rand_num(mteng), 0, i, true, 100);
  
    string name2="pairing"+to_string(i);
    par_lst_ptr->build_aux_param(PAIRING, name2, rand_num(mteng), 0, i, true, 100);  
    //string name2="hopL"+to_string(i);
    //par_lst_ptr->build_aux_param(HOPPING, name2, rand_num(mteng), 0, L+i, true, 100);
  }
  par_lst_ptr->build_jas_param(SPIN, "spin1", pos_rand_num(mteng), 0, 1, true, 100);
 
  //par_lst_ptr->build_aux_param(HOPPING, "hopp1", rand_num(mteng), 0, 1, true, 100);
  //par_lst_ptr->build_aux_param(PAIRING, "pair1", rand_num(mteng), 0, 1, true, 100);
  //par_lst_ptr->build_aux_param(PAIRING, "pair2", rand_num(mteng), 0, 2, true, 100);
  //par_lst_ptr->build_aux_param(PAIRING, "pair3", rand_num(mteng), 0, 3, true, 100);


  //par_lst_ptr->build_jas_param(SPIN, "spin1", 0.1, 0, 1, true, 100);
  //par_lst_ptr->build_jas_param(SPIN, "spin1", 0.1, 0, 2, true, 100);
  //par_lst_ptr->build_jas_param(SPIN, "spin2", pos_rand_num(mteng), 0, 2, true, 100);
  
  par_lst_ptr->report_aux_params();
  //par_lst_ptr->report_jas_params(); 

  VMC::ObsUPtr enrg_ptr=make_unique<VMC::Observables::IsingChainEnergy>
                        ("Ising-Energy", 100, 1.0, 0.0); 
  //VMC::ObsUPtr enrg_ptr=make_unique<VMC::Observables::HeisenbergChainEnergy>
  //                      ("Heisenberg-Energy", 100, 1.0, 1.0, 1.0); 

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
  
  
  //size_t L=6;
  //size_t equil=2000;
  //size_t simul=10000;
  //size_t binsize=100;
  //size_t vsteps = 50;
  //double df = 0.05;
  //VMC::ParamList_t params;
  //VMC::JspList_t jsparams;
  
  // BCS parameters ------------------------------------------------------------
  //VMC::BCSVarParam onsite(0.5, 0, VMC::Onsite, 2*L, binsize, "onsite");
  //VMC::BCSVarParam nnhop(0.5, 1, VMC::Hopping, 2*L, binsize, "nnhop");
  //VMC::BCSVarParam ospair(0.5, 0, VMC::Pairing, 2*L, binsize, "ospair");
  //VMC::BCSVarParam nnpair(0.5, 1, VMC::Pairing, 2*L, binsize, "nnpair");
  //VMC::BCSVarParam nnnpair(0.5, 2, VMC::Pairing, 2*L, binsize, "nnnpair");
  //VMC::BCSVarParam nnnnpair(0.5, 3, VMC::Pairing, 2*L, binsize, "nnnnpair");
  // ---------------------------------------------------------------------------
  //
  //// Jastrow parameters --------------------------------------------------------
  //VMC::JastrowParam nnjastrow(-0.1, 1, binsize, "nnJastrow");
  //VMC::JastrowParam nnnjastrow(0.1, 2, binsize, "nnnJastrow");
  //VMC::JastrowParam nnnnjastrow(-0.1, 3, binsize, "nnnnJastrow");
  //// ---------------------------------------------------------------------------
  //
  //params.push_back(onsite);
  //params.push_back(nnhop);
  //params.push_back(ospair);
  //params.push_back(nnpair);
  //params.push_back(nnnpair);
  //params.push_back(nnnnpair);
  //jsparams.push_back(nnjastrow);
  //jsparams.push_back(nnnjastrow);
  //jsparams.push_back(nnnnjastrow);
  //VMC::BCSChainHamiltonian testauxham(2*L, params);
  //cout << testauxham.get_hamiltonian() << endl;
  //cout << testauxham.get_eigenvals() << endl;
  //cout << testauxham.get_eigenvecs() << endl;
  //for(auto it=params.begin(); it!=params.end(); it++) {
  //  cout << it->name << endl;
  //  cout << "vmat: " << endl;
  //  cout << it->vmat << endl;
  //  cout << "mmat: " << endl;
  //  cout << it->mmat << endl;
  //  cout << "----" << endl; 
  //} 
  
  //VMC::HeisChainSim simulator(L, binsize, params);
  ////simulator.print_gmat();
  ////simulator.print_spinstate();
  ////simulator.print_operslist();
  ////std::cout << simulator._heisenergy() << std::endl;
  ////simulator._flipspin(0);
  ////simulator.print_gmat();
  ////simulator.print_spinstate();
  ////simulator.print_operslist();
  ////std::cout << simulator._heisenergy() << std::endl;
  //simulator.optimize(vsteps, equil, simul, df, "../dat/n4_optvals_ising"); 
  //simulator.print_params();  
  return 0;
}
