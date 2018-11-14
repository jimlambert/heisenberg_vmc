#include <iostream>
#include <Eigen/Dense>
#include "VarParam.h"
#include "BCSHamiltonian.h"

using namespace std;

int main(){
 
  int L = 8; 
  VMC::ParamList_t params;
  VMC::BCSVarParam onsite;
  VMC::BCSVarParam nnhop;
  VMC::BCSVarParam nnnhop;
  VMC::BCSVarParam pairing;
  onsite.val=2.0;
  onsite.space=0;
  onsite.type=VMC::Onsite; 
  nnhop.val=3.0;
  nnhop.space=1;
  nnhop.type=VMC::Hopping;
  nnnhop.val=4.0;
  nnnhop.space=2;
  nnnhop.type=VMC::Hopping;
  pairing.val=5.0;
  pairing.space=2;
  pairing.type=VMC::Pairing;
  params.push_back(onsite);
  //params.push_back(nnhop);
  //params.push_back(nnnhop);
  params.push_back(pairing);
  VMC::BCSChainHamiltonian bcsChain(L, 5, params);
  return 0;
}
