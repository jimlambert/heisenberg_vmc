#include <iostream>
#include <Eigen/Dense>
#include "VarParam.h"
#include "Simulation.h"
#include "BCSHamiltonian.h"

using namespace std;

int main() {
  int L=2;
  double val=2.0;
  int space=0;
  VMC::ParamList_t params;
  VMC::BCSVarParam onsite(val, space, VMC::Onsite, 2*L);
  VMC::BCSVarParam nnhopp(2.0, 1, VMC::Hopping, 2*L);
  VMC::BCSVarParam nnpair(3.0, 1, VMC::Pairing, 2*L);
  params.push_back(onsite);
  params.push_back(nnhopp);
  params.push_back(nnpair);
  VMC::HeisChainSim simulation(L, params);
  simulation.optimize();
  return 0;
}
