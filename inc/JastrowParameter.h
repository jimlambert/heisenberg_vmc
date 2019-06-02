#ifndef JASTROW_PARAMERTER_H
#define JASTROW_PARAMERTER_H

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
    const size_t&      si,  // parameter site i
    const size_t&      sj,  // parameter site j
    const bool&        ti,  // translation invariant (true)
    const size_t&      bs   // binsize for local_meas 
  ) : Parameter(JASTROW, n, v, si, sj, ti, bs), subtype(jpt) {}
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

} // namespace VMC 

#endif // JASTROW_PARAMERTER_H
