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
 
  //VMC::LocalMeasurement<std::complex<double> >  testmeas(100);

  //for(int i=0; i<1000; i++) {
  //  testmeas.push(i*0.5); 
  //} 
  //cout << testmeas.ave() << endl;

  size_t L=4;
  size_t equil=10000;
  size_t simul=100000;
  size_t binsize=100;
  size_t vsteps = 100;
  double df = 0.001;
  double p1=1;
  double p2=2;
  double p3=3;
  double p4=4;
  VMC::ParamList_t params;
  VMC::BCSVarParam onsite(p1, 0, VMC::Onsite, 2*L, binsize, "onsite");
  VMC::BCSVarParam nnhop(p2, 1, VMC::Hopping, 2*L, binsize, "nnhop");
  VMC::BCSVarParam ospair(p3, 0, VMC::Pairing, 2*L, binsize, "ospair");
  VMC::BCSVarParam nnpair(p4, 1, VMC::Pairing, 2*L, binsize, "nnpair");
  //VMC::BCSVarParam nnnpair(p3, 2, VMC::Pairing, 2*L, binsize, "nnnpair");
  params.push_back(onsite);
  params.push_back(nnhop);
  params.push_back(ospair);
  params.push_back(nnpair); 
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
  simulator.optimize(vsteps, equil, simul, df, "./n4_optvals_ising"); 
  return 0;
}
