// -----------------------------------------------------------------------------
// Structure for variational parameters. Each paramter can be a hopping
// parameter or a pairing parameter. The ParamList is updated during optimizing
// then passed to the auxiliary Hamiltonian where is it used to udpate the
// Hamiltonian and construct a new groundstate.
// -----------------------------------------------------------------------------

#ifndef VARPARAM_H
#define VARPARAM_H

#include <Eigen/Dense>
#include <vector>

namespace VMC{

enum ParamType {Onsite, Hopping, Pairing};

struct BCSVarParam{
  double val;
  int space;
  ParamType type;
  Eigen::MatrixXd _vmat;
  Eigen::MatrixXd _mmat;
  BCSVarParam(const double& v, const int& s, type& t, const size_t L) 
  : val(v), space(s) type(t) {
    _vmat.resize(L,L)=MatrixXd::Zero(L,L);
    switch(type) {
      case Onsite:
        for(size_t i=0; i<2*L; i++) {
          if(i<(_size/2)) _vmat(i, i)=-1.0;
          else _vmat(i, i) = 1; 
        }
        break
      case Hopping:
        for(size_t i=0; i<_size; i++) {
          if(i<(_size/2)) {
            _vmat(i, (i+space)%(_size/2))=1;
          }
          else {
          
          }
        } 
        break
      case Pairing:

        break
    }  
  }
};

typedef std::vector<BCSVarParam> ParamList_t;

}

#endif // VARPARAM_H
