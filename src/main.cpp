#include <cstdlib>
#include <time.h>
#include <iostream>
#include <fstream>
#include <Eigen/Dense>
#include "VarParam.h"
#include "Simulation.h"
#include "LocalMeasurement.h"
#include "BCSHamiltonian.h"

using namespace std;

int main(int argc, char* argv[]) {

  size_t L=8;
  size_t equil=10000;
  size_t simul=30000;
  size_t binsize=100;
  size_t vsteps = 100;
  double df = 0.01;
  double p1=0.5;
  double p2=0.5;
  double p3=0.5;
  VMC::ParamList_t params;
  VMC::JspList_t jsparams;
  
  // BCS parameters ------------------------------------------------------------
  VMC::BCSVarParam onsite(p1, 0, VMC::Onsite, 2*L, binsize, "onsite");
  VMC::BCSVarParam nnhop(p2, 1, VMC::Hopping, 2*L, binsize, "nnhop");
  VMC::BCSVarParam ospair(p3, 0, VMC::Pairing, 2*L, binsize, "ospair");
  // ---------------------------------------------------------------------------
  
  // Jastrow parameters --------------------------------------------------------
  VMC::JastrowParam nnjastrow(0.1, 1, binsize, "nnJastrow");
  VMC::JastrowParam nnnjastrow(-0.1, 2, binsize, "nnnJastrow");
  // ---------------------------------------------------------------------------
  
  params.push_back(onsite);
  params.push_back(nnhop);
  params.push_back(ospair);
  jsparams.push_back(nnjastrow);
  jsparams.push_back(nnnjastrow);
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
  
  VMC::HeisChainSim simulator(L, binsize, params, jsparams);
  //simulator.print_spinstate();
  //simulator.print_operslist();
  //std::cout << simulator.jsf() << std::endl;
  //cout << "----" << endl;
  //simulator._flipspin(0);
  //simulator.print_spinstate();
  //simulator.print_operslist();
  //std::cout << simulator.jsf() << std::endl;
  //cout << "----" << endl;
  //simulator.print_redmat();
  //cout << "----" << endl;
  //simulator.print_gmat();
  //simulator.print_spinstate();
  //simulator.print_operslist();
  //cout << "----" << endl;
  //simulator._exchange(0,1);
  //simulator.print_gmat();
  //cout << "---- reinit ----" << endl;
  //simulator._reinitgmat();
  //simulator.print_gmat();
  //simulator.print_spinstate();
  //simulator.print_operslist();
  simulator.optimize(vsteps, equil, simul, df, "./n4_optvals_ising"); 
  return 0;
}
