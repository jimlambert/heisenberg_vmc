#ifndef ABSTRACT_PARAMETER_H
#define ABSTRACT_PARAMETER_H


#include <Eigen/Dense>
#include <vector>
#include <complex>
#include <memory>
#include <string>
#include "BasisState.h"
#include "LocalMeasurement.h"

// VMC::Parameters structs
// =============================================================================
//  Variational parameters are constructed through an inheritance structure that
// allows us to take advantage of polymorphism when accessing them. The
// structure is as follows:
//                              Parameter
//                                  |
//                 -----------------------------------
//                 |                                 |
//             Auxiliary                          Jastrow
//                                                   |
//                                             -------------    
//                                             |           |
//                                          SpinSpin      ...
//
// More parameters may be added under the appropriate branch of this structure,
// such as density-density Jastrow factors.
//
//  Parameter
//  ---------
//  
//      The paremeter class contains information that is universal to all
//    variational parameters such as the name, value, and sites on which that
//    parameter acts. Moreover, it contains a local_meas container to recored 
//    the values of this variational parameter during the simulation.
//      Finally, the Parameter class contains a pure virtual methods,
//    
//      virtual void operator() (const BasisState&)=0
//
//    which determine the value of this variational parameter in the current basis
//    state. 
//
//  AuxiliaryParameter
//  ------------------
//
//      The auxiliary parameters are those which act as terms in the auxiliary
//    Hamiltonian. Determining how these terms need to be modified in order to
//    minimize the energy involves a low order perturbation theory calculation. 
//    In particular, the matrix vmat encodes the general structure of the
//    variational parameter and its action on the raising and lowering operators
//    in the auxiliary Hamiltonian. 
//      The matrix mmat is the matrix which will actually come into play when
//    determining the expected value of the operator corresponding to a shift in
//    the variational parameter. This matrix is recalculated for each new
//    variational step.
//      Both mmat and vmat are defined when the auxiliary Hamiltonian is
//    constructed. 
//      In the case of the Onsite parameter there is a redundant member variable,
//    since such a parameter only acts on one particular position. In such 
//    cases it is convention to use site_i in the code as the variable defining 
//    the site on which this parameter acts.
//  
//  JastrowParameter
//  ----------------
//
//      These parameters enter the variational computation as prefactors on the
//    variational wavefunction. They are generally diagonal in the sampling
//    basis to avoid computational inefficiency. The Jastrow parameters also
//    contain a boolean, trans_inv which determines whether or not these
//    variables are translation invariant. This is useful when meaeasuring their
//    value on a particular basis state
//
// =============================================================================

namespace VMC {
namespace Parameters {

// Parameter types are denoted by enums
// =============================================================================

enum ParameterType {AUXILIARY, JASTROW};
enum ParameterSubtype {ONSITE, HOPPING, PAIRING, SPIN};

// =============================================================================
// Parent struct forming the basis of all variational parameters
// =============================================================================

struct Parameter {
  ParameterType           type;
  std::string             name;
  std::complex<double>    val;
  size_t                  site_i;
  size_t                  site_j;
  bool                    trans_inv; 
  MeasCd                  local_meas;
  Parameter(
    ParameterType      pt, // set parameter type
    const std::string& n,  // parameter name 
    const double&      v,  // value of parameter
    const size_t&      si, // parameter site i
    const size_t&      sj, // parameter site j
    const bool&        ti, // translation invariant (true)
    const size_t&      bs  // binsize for local_meas 
  ) : type(pt), 
      name(n), 
      val(v), 
      site_i(si), 
      site_j(sj), 
      trans_inv(ti),
      local_meas(bs) 
  {}
  ParameterType get_type() {return type;}
  virtual ParameterSubtype get_subtype()=0;
  virtual void operator() (BasisState&)=0;
};

// =============================================================================

} // Parameters namespace

// Typedefs for Parameters namespace
// =============================================================================

typedef std::unique_ptr<Parameters::Parameter> ParamUPtr;
typedef std::shared_ptr<Parameters::Parameter> ParamSPtr;
typedef std::vector<ParamUPtr>                 ParamsUVec;
typedef std::vector<ParamSPtr>                 ParamsSVec;

} // VMC namespace

#endif // ABSTRACT_PARAMETER_H
