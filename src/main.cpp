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
  if(argc==1) {
    cout << "missing input file" << endl;
    exit(1);
  }
  ifstream infile;
  infile.open(argv[1]); 
  size_t L=4;
  size_t equil=100;
  size_t simul=100;
  size_t vstep;
  size_t binsize=100;
  double df = 0.01;
  double p1, p2, p3;
  infile >> vstep;
  infile >> p1;
  infile >> p2;
  infile >> p3;
  VMC::ParamList_t params;
  VMC::BCSVarParam onsite(p1, 0, VMC::Onsite, 2*L, binsize, "onsite");
  VMC::BCSVarParam nnhop(p2, 1, VMC::Hopping, 2*L, binsize, "nnhop");
  VMC::BCSVarParam ospair(p3, 0, VMC::Pairing, 2*L, binsize, "ospair");
  VMC::BCSVarParam nnpair(p3, 1, VMC::Pairing, 2*L, binsize, "nnpair");
  params.push_back(onsite);
  params.push_back(nnhop);
  //params.push_back(ospair);
  params.push_back(nnpair);
  //VMC::BCSChainHamiltonian testaux(2*L, params);
  
  VMC::HeisChainSim simulator(L, params);
  //simulator.print_spinstate();
  //simulator.print_operslist();
  //simulator.print_gmat();
  //simulator._flipspin(0);
  //cout << "------" << endl;
  //simulator.print_spinstate();
  //simulator.print_operslist();
  //Vsimulator.print_gmat();
  simulator.optimize(vstep, equil, simul, df);
 
  return 0;
}
