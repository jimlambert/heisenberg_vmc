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
  Eigen::MatrixXd _mmat;
};

typedef std::vector<BCSVarParam> ParamList_t;

}

#endif // VARPARAM_H
