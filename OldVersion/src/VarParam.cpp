#include <iostream>
#include <string>
#include <iostream>
#include "BasisState.h"
#include "VarParam.h"

namespace VMC {

//BCSVarParam::BCSVarParam
//(const double& v, const size_t& s, const ParamType& t, const size_t& L, 
// const size_t& bs, const std::string& n) 
//: val(v), space(s), type(t), lmeas(bs), name(n) {
//  mmat.resize(L,L);
//  vmat.resize(L,L);
//  vmat=Eigen::MatrixXd::Zero(L,L);
//  switch(type) {
//    case Onsite:
//      for(size_t i=0; i<L; i++) {
//        if(i<(L/2)) vmat(i, i)+=-1.0;
//        else vmat(i, i)+=1; 
//      }
//      break;
//    case Hopping:
//      for(size_t i=0; i<L; i++) {
//        if(i<(L/2)) {
//          vmat(i, (i+space)%(L/2))+=1;
//          vmat((i+space)%(L/2), i)+=1;
//        }
//        else {
//          vmat(i, (i+space)%(L/2)+(L/2))+=-1;
//          vmat((i+space)%(L/2)+(L/2), i)+=-1;
//        }
//      } 
//      break;
//    case Pairing:
//      //int ffac=std::pow(-1, s);
//      for(size_t i=0; i<L/2; i++) {
//        vmat(i, (i+space)%(L/2)+(L/2))+=1;
//        vmat((i+space)%(L/2)+(L/2), i)+=1;
//      }
//      break;
//  }
//}

BCSVarParam::BCSVarParam
(const double& v, const size_t& s, const ParamType& t, const size_t& L, 
 const size_t& bs, const std::string& n) 
: val(v), space(s), type(t), lmeas(bs), name(n) {
  mmat.resize(L,L);
  vmat.resize(L,L);
  int bcfac;
  if((L/2)%4==0) bcfac=-1;
  else           bcfac=1;
  vmat=Eigen::MatrixXd::Zero(L,L);
  switch(type) {
    case Onsite:
      for(size_t i=0; i<L; i++) {
        if(i<(L/2)) vmat(i, i)+=-1.0;
        else vmat(i, i)+=1; 
      }
      break;
    case Hopping:
      for(size_t i=0; i<L; i++) {
        if(i<(L/2)) {
          int cfac=1;
          size_t j=(i+space)%(L/2);
          if((i+space)>=L/2) cfac=bcfac;
          vmat(i,j)+=cfac*1;
          vmat(j,i)+=cfac*1;
        }
        else {
          int cfac=1;
          size_t j=(i+space)%(L/2)+(L/2);
          if((i+space)>=L) cfac=bcfac;
          vmat(i,j)+=-1*cfac;
          vmat(j,i)+=-1*cfac;
        }
      } 
      break;
    case Pairing:
      //int ffac=std::pow(-1, s);
      for(size_t i=0; i<L/2; i++) {
        int cfac=1;
        size_t j=(i+space)%(L/2)+(L/2);
        if((i+space)>=(L/2)) cfac=bcfac;
        vmat(i, (i+space)%(L/2)+(L/2))+=1*cfac;
        vmat((i+space)%(L/2)+(L/2), i)+=1*cfac;
      }
      break;
  }
}

}
