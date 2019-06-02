// -----------------------------------------------------------------------------
// Structure for variational parameters. Each paramter can be a hopping
// parameter or a pairing parameter. The ParamList is updated during optimizing
// then passed to the auxiliary Hamiltonian where is it used to udpate the
// Hamiltonian and construct a new groundstate.
// -----------------------------------------------------------------------------

#ifndef VARPARAM_H
#define VARPARAM_H

#include <Eigen/Dense>
#include <string>
#include <vector>
#include <complex>
#include <memory>
#include "BasisState.h"
#include "LocalMeasurement.h"

namespace VMC{

enum ParamType {Onsite, Hopping, Pairing};

struct BCSVarParam {
  double val; // current value of variational parameter
  size_t space; // spacing for operator
  LocalMeasurement<std::complex<double> > lmeas; // local measurements of associated operator 
  ParamType type; // type of parameter chosen from ParamType
  std::string name;
  Eigen::MatrixXd vmat; // matrix defining structure of associated operator
  Eigen::MatrixXcd mmat; // matrix used to measure local effect of parameter
  BCSVarParam(const double&, const size_t&, const ParamType&, const size_t&, 
              const size_t&, const std::string&);
};

typedef std::vector<BCSVarParam> ParamList_t;

struct JastrowParam {
  double val; // value of variational parameter
  size_t space; // spacing between S_i^z and S_j^z
  LocalMeasurement<std::complex<double> > lmeas; // local measurement values
  std::string name;
  JastrowParam
  (const double& v, const size_t& s, const size_t& bs, const std::string& n) 
    : val(v), space(s), lmeas(bs), name(n) {}
};

typedef std::vector<JastrowParam> JspList_t;

} // VMC namespace

#endif // VARPARAM_H
