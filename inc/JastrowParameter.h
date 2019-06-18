#ifndef JASTROW_PARAMETER_H
#define JASTROW_PARAMETER_H

#include "AbstractParameter.h"

namespace VMC {
namespace Parameters {

// Jastrow parameter parent struct
// =============================================================================

struct JastrowParameter : Parameter {
  ParameterSubtype     subtype;
  JastrowParameter(
    ParameterSubtype   jpt, // subtype for Jastrow parameter
    const std::string& n,   // parameter name 
    const double&      v,   // value of parameter
    const size_t&      s,   // parameter site 
    const size_t&      d,   // change in position
    const bool&        ti,  // translation invariant (true)
    const size_t&      bs   // binsize for local_meas 
  ) : Parameter(JASTROW, n, v, s, d, ti, bs), subtype(jpt) {}
  ParameterSubtype get_subtype() {return subtype;}
};

// =============================================================================

} // namespace Parameters 

// Typedefs for Parameters namespace
// =============================================================================

typedef std::unique_ptr<Parameters::JastrowParameter> JasParamUPtr;
typedef std::shared_ptr<Parameters::JastrowParameter> JasParamSPtr;
typedef std::vector<JasParamUPtr>                     JasParamUVec;
typedef std::vector<JasParamSPtr>                     JasParamSVec;

// =============================================================================

} // namespace VMC 

#endif // JASTROW_PARAMETER_H
