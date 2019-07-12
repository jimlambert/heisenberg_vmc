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
  JastrowParameter(
    ParameterSubtype   jpt, // subtype for Jastrow parameter
    const std::string& n,   // parameter name 
    const double&      v,   // value of parameter
    const size_t&      x1,  // parameter site 
    const size_t&      y1,  // parameter site 
    const size_t&      x2,  // parameter site 
    const size_t&      y2,  // parameter site 
    const bool&        ti,  // translation invariant (true)
    const size_t&      bs   // binsize for local_meas 
  ) : Parameter(JASTROW, n, v, x1, y1, x2, y2, ti, bs), subtype(jpt) {}
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
