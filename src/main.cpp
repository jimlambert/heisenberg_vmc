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
  
  size_t L=10;
  size_t equil=10000;
  size_t simul=100000;
  size_t binsize=100;
  size_t vsteps = 300;
  double df = 0.001;
  double p1=0.5;
  double p2=0.5;
  double p3=0.5;
  VMC::ParamList_t params;
  VMC::BCSVarParam onsite(p1, 0, VMC::Onsite, 2*L, binsize, "onsite");
  VMC::BCSVarParam nnhop(p2, 1, VMC::Hopping, 2*L, binsize, "nnhop");
  VMC::BCSVarParam ospair(p3, 0, VMC::Pairing, 2*L, binsize, "ospair");
  VMC::BCSVarParam nnpair(p3, 1, VMC::Pairing, 2*L, binsize, "nnpair");
  params.push_back(onsite);
  params.push_back(nnhop);
  params.push_back(ospair);
  params.push_back(nnpair); 
  //VMC::BCSChainHamiltonian testauxham(2*L, params);
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
  
  VMC::HeisChainSim simulator(L, params);
  simulator.optimize(vsteps, equil, simul, df, "./n2_optvals"); 
  return 0;
}
