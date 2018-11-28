#include <iostream>
#include "VarParam.h"

namespace VMC {

BCSVarParam::BCSVarParam
(const double& v, const int& s, const ParamType& t, const size_t& L, 
 const size_t& bs) : localvals(bs), val(v), space(s), type(t) {
  mmat.resize(L,L);
  vmat.resize(L,L);
  vmat=Eigen::MatrixXd::Zero(L,L);
  switch(type) {
    case Onsite:
      for(size_t i=0; i<L; i++) {
        if(i<(L/2)) vmat(i, i)=-1.0;
        else vmat(i, i) = 1; 
      }
      break;
    case Hopping:
      for(size_t i=0; i<L; i++) {
        if(i<(L/2)) {
          vmat(i, (i+space)%(L/2))=1;
          vmat((i+space)%(L/2), i)=1;
        }
        else {
          vmat(i, (i+space)%(L/2)+(L/2))=-1;
          vmat((i+space)%(L/2)+(L/2), i)=-1;
        }
      } 
      break;
    case Pairing:
      for(size_t i=0; i<L/2; i++) {
        vmat(i, (i+space)%(L/2)+(L/2))=1;
        vmat((i+space)%(L/2)+(L/2), i)=1;
      }
      break;
  }
}
}
