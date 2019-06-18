#include <cstdlib>
#include <time.h>
#include <iostream>
#include <fstream>
#include <Eigen/Dense>
#include "VarParam.h"
#include "HoppingChainHamiltonian.h"
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

  
  VMC::ParamListPtr par_lst_ptr=make_unique<VMC::Parameters::ParameterList>();

  par_lst_ptr->build_aux_param(ONSITE, "onsite", 0.5, 0, 0, true, 100);  
  par_lst_ptr->build_aux_param(HOPPING, "hop1", 0.1, 0, 1, true, 100);
  par_lst_ptr->build_aux_param(HOPPING, "hop2", 0.2, 0, 2, true, 100);
  //par_lst_ptr->build_aux_param(HOPPING, "hop3", 0.3, 0, 3, true, 100);
  par_lst_ptr->report_aux_params();  
  
  VMC::AuxHamUPtr aux_ham_ptr=make_unique<VMC::HopChainHam>
                         (true, 4, par_lst_ptr->aux_vec());
  cout << aux_ham_ptr->get_eigenvectors() << endl;
  cout << aux_ham_ptr->get_eigenvalues() << endl;
  for(size_t i=0; i<par_lst_ptr->naux(); i++) {
    cout << par_lst_ptr->aux(i).mmat << endl;
    cout << "----" << endl;
  } 
  //size_t L=10;
  //size_t equil=2000;
  //size_t simul=10000;
  //size_t binsize=100;
  //size_t vsteps = 50;
  //double df = 0.05;
  //double p1=0.5;
  //double p2=0.5;
  //double p3=0.5;
  //VMC::ParamList_t params;
  //VMC::JspList_t jsparams;
  //
  //// BCS parameters ------------------------------------------------------------
  //VMC::BCSVarParam onsite(1.0, 0, VMC::Onsite, 2*L, binsize, "onsite");
  //VMC::BCSVarParam nnhop(1.0, 1, VMC::Hopping, 2*L, binsize, "nnhop");
  //VMC::BCSVarParam ospair(0.5, 0, VMC::Pairing, 2*L, binsize, "ospair");
  //VMC::BCSVarParam nnpair(0.001, 1, VMC::Pairing, 2*L, binsize, "nnpair");
  //VMC::BCSVarParam nnnpair(0.027, 2, VMC::Pairing, 2*L, binsize, "nnnpair");
  //VMC::BCSVarParam nnnnpair(-0.02, 3, VMC::Pairing, 2*L, binsize, "nnnnpair");
  //// ---------------------------------------------------------------------------
  //
  //// Jastrow parameters --------------------------------------------------------
  //VMC::JastrowParam nnjastrow(-0.1, 1, binsize, "nnJastrow");
  //VMC::JastrowParam nnnjastrow(0.1, 2, binsize, "nnnJastrow");
  //VMC::JastrowParam nnnnjastrow(-0.1, 3, binsize, "nnnnJastrow");
  //// ---------------------------------------------------------------------------
  //
  ////params.push_back(onsite);
  //params.push_back(nnhop);
  ////params.push_back(ospair);
  //params.push_back(nnpair);
  //params.push_back(nnnpair);
  //params.push_back(nnnnpair);
  ////jsparams.push_back(nnjastrow);
  ////jsparams.push_back(nnnjastrow);
  ////jsparams.push_back(nnnnjastrow);
  ////VMC::BCSChainHamiltonian testauxham(2*L, params);
  ////cout << testauxham.get_hamiltonian() << endl;
  ////cout << testauxham.get_eigenvals() << endl;
  ////cout << testauxham.get_eigenvecs() << endl;
  ////for(auto it=params.begin(); it!=params.end(); it++) {
  ////  cout << it->name << endl;
  ////  cout << "vmat: " << endl;
  ////  cout << it->vmat << endl;
  ////  cout << "mmat: " << endl;
  ////  cout << it->mmat << endl;
  ////  cout << "----" << endl; 
  ////} 
  //
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
