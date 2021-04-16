#ifndef AUXILIARY_PARAMETER_H
#define AUXILIARY_PARAMETER_H

#include "AbstractParameter.h"

namespace VMC {
namespace Parameters {

// Auxiliary parameter struct 
// =============================================================================

struct AuxiliaryParameter : Parameter {
  ParameterSubtype        subtype;       
  Eigen::MatrixXd         vmat;
  Eigen::MatrixXcd        mmat;
  bool                    vinit; // set true when vmat initialized
  bool                    minit; // set true when mmat initialized
  
  AuxiliaryParameter(
    ParameterSubtype   apt, // auxiliary parameter type
    const std::string& n,   // parameter name 
    const double&      v,   // value of parameter
    const size_t&      s,   // parameter site 
    const size_t&      d,   // change in position
    const bool&        ti,  // translation invariant (true)
    const size_t&      bs   // binsize for local_meas 
  ) : Parameter(AUXILIARY, n, v, s, d, ti, bs), 
      subtype(apt), 
      vinit(false), 
      minit(false)
  {}
  
  AuxiliaryParameter(
    ParameterSubtype   apt, // auxiliary parameter type
    const std::string& n,   // parameter name 
    const double&      v,   // value of parameter
    const size_t&      x1,  // parameter site 
    const size_t&      y1,  // parameter site 
    const size_t&      x2,  // parameter site 
    const size_t&      y2,  // parameter site 
    const bool&        ti,  // translation invariant (true)
    const size_t&      bs   // binsize for local_meas 
  ) : Parameter(AUXILIARY, n, v, x1, y1, x2, y2, ti, bs), 
      subtype(apt), 
      vinit(false), 
      minit(false)
  {}
  ParameterSubtype get_subtype() {return subtype;}
  void operator() (BasisState&, Eigen::MatrixXcd&);
};

// =============================================================================

} // namespace Parameter 

// Typedefs for Parameters namespace
// =============================================================================

typedef std::unique_ptr<Parameters::AuxiliaryParameter> AuxParamUPtr;
typedef std::shared_ptr<Parameters::AuxiliaryParameter> AuxParamSPtr;
typedef std::vector<AuxParamUPtr>                       AuxParamUVec;
typedef std::vector<AuxParamSPtr>                       AuxParamSVec;
typedef std::vector<Parameters::AuxiliaryParameter>     AuxParamVec;

} // namespace VMC 

#endif // AUXILIARY_PARAMETER_H
