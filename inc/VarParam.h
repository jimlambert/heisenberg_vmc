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
#include <complex>
#include "LocalMeasurement.h"

namespace VMC{

enum ParamType {Onsite, Hopping, Pairing};

struct BCSVarParam{
  LocalMeasurement<std::complex<double> > lmeas; // local measurements of associated operator 
  std::complex<double> val; // current value of variational parameter
  int space; // spacing for operator
  ParamType type; // type of parameter chosen from ParamType
  std::string name;
  Eigen::MatrixXd vmat; // matrix defining structure of associated operator
  Eigen::MatrixXcd mmat; // matrix used to measure local effect of parameter
  BCSVarParam(const double&, const int&, const ParamType&, const size_t&, 
              const size_t&, const std::string&);
};

typedef std::vector<BCSVarParam> ParamList_t;

}

#endif // VARPARAM_H
