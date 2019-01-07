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
  size_t equil=1000;
  size_t simul=100000;
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
  VMC::BCSVarParam nnhopp(p2, 1, VMC::Hopping, 2*L, binsize, "nnhop");
  VMC::BCSVarParam nnpair(p3, 1, VMC::Pairing, 2*L, binsize, "nnpair");
  params.push_back(onsite);
  params.push_back(nnhopp);
  params.push_back(nnpair);
  VMC::HeisChainSim simulator(L, params);
  simulator.print_spinstate();
  for(size_t i=0; i<vstep; i++) {
    simulator._sweep();
    simulator.print_spinstate();
  }
  //simulator.optimize(vstep, equil, simul, df);
 
  return 0;
}
