#include <cstdlib>
#include <time.h>
#include <iostream>
#include <Eigen/Dense>
#include "VarParam.h"
#include "Simulation.h"
#include "LocalMeasurement.h"
#include "BCSHamiltonian.h"

using namespace std;

int main() {
  size_t L=4;
  size_t equil=100000;
  size_t simul=100000;
  size_t vstep=1;
  size_t binsize=100;
  double df = 0.01;
  VMC::ParamList_t params;
  VMC::BCSVarParam onsite(0.5, 0, VMC::Onsite, 2*L, binsize, "onsite");
  VMC::BCSVarParam nnhopp(1.0, 1, VMC::Hopping, 2*L, binsize, "nnhop");
  VMC::BCSVarParam nnpair(1.0, 1, VMC::Pairing, 2*L, binsize, "nnpair");
  params.push_back(onsite);
  params.push_back(nnhopp);
  params.push_back(nnpair);
  VMC::HeisChainSim simulator(L, params);
  simulator.optimize(vstep, equil, simul, df);
 
  return 0;
}
