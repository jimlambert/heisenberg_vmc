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

  size_t L=2;
  size_t equil=10000;
  size_t simul=20000;
  size_t binsize=100;
  size_t vsteps = 200;
  double df = 0.01;
  double p1=1;
  double p2=1;
  double p3=3;
  double p4=1;
  double p5=1;
  VMC::ParamList_t params;
  VMC::BCSVarParam onsite(p1, 0, VMC::Onsite, 2*L, binsize, "onsite");
  VMC::BCSVarParam nnhop(p2, 1, VMC::Hopping, 2*L, binsize, "nnhop");
  VMC::BCSVarParam nnnhop(p3, 2, VMC::Hopping, 2*L, binsize, "nnnhop");
  VMC::BCSVarParam nnnnhop(p4, 3, VMC::Hopping, 2*L, binsize, "nnnnhop");
  VMC::BCSVarParam ospair(p3, 0, VMC::Pairing, 2*L, binsize, "ospair");
  VMC::BCSVarParam nnpair(p4, 1, VMC::Pairing, 2*L, binsize, "nnpair");
  VMC::BCSVarParam nnnpair(p5, 2, VMC::Pairing, 2*L, binsize, "nnnpair");
  params.push_back(onsite);
  params.push_back(nnhop);
  //params.push_back(nnnhop);
  //params.push_back(nnnnhop);
  params.push_back(ospair);
  params.push_back(nnpair);
  //params.push_back(nnnpair); 
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
  
  VMC::HeisChainSim simulator(L, binsize, params);
  //cout << "----" << endl;
  //simulator.print_redmat();
  //cout << "----" << endl;
  simulator.print_gmat();
  simulator.print_spinstate();
  simulator.print_operslist();
  std::cout << simulator._heisenergy() << std::endl;
  cout << "----" << endl;
  //simulator._flipspin(1);
  //simulator.print_gmat();
  //simulator.print_spinstate();
  //simulator.print_operslist();
  //std::cout << simulator._heisenergy() << std::endl;
  //simulator.optimize(vsteps, equil, simul, df, "./n4_optvals_ising"); 
  return 0;
}
