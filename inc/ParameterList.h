#ifndef PARAMETER_LIST_H
#define PARAMETER_LIST_H

#include <vector>
#include <memory>
#include "AuxiliaryParameter.h"
#include "JastrowParameter.h"

namespace VMC {
namespace Parameters {

// Parameter list class use by the wavefunction
// =============================================================================

class ParameterList {
  private:
    size_t       _naux;
    size_t       _njas;
    AuxParamUVec _aux_params;
    JasParamUVec _jas_params;
  public:
    ParameterList() : _naux(0), _njas(0) {}
    
    // push externally created parameters
    void push_aux_param(AuxParamUPtr);
    void push_jas_param(JasParamUPtr);

    // build parameters directly
    void build_aux_param(
      ParameterSubtype   apt, // auxiliary parameter type
      const std::string& n,   // parameter name 
      const double&      v,   // value of parameter
      const size_t&      s,   // parameter site 
      const size_t&      d,   // change in position
      const bool&        ti,  // translation invariant (true)
      const size_t&      bs   // binsize for local_meas 
    );
    //void build_jas_param(
    //  ParameterSubtype   jpt, // subtype for Jastrow parameter
    //  const std::string& n,   // parameter name 
    //  const double&      v,   // value of parameter
    //  const size_t&      s,   // parameter site i
    //  const size_t&      d,   // parameter site j
    //  const bool&        ti,  // translation invariant (true)
    //  const size_t&      bs   // binsize for local_meas 
    //);

    // access functions
    size_t njas() {return _njas;}
    size_t naux() {return _naux;}
    size_t npar() {return _njas+_naux;}
    AuxParamUVec& aux_vec();
    JasParamUVec& jas_vec();
    AuxParamUPtr& aux(const size_t&);
    JasParamUPtr& jas(const size_t&);
    
    // output functions
    void report_params() {
      report_aux_params();
      report_jas_params();
    }
    void report_aux_params();
    void report_jas_params();
};

// =============================================================================

} // namespace Parameters

typedef std::unique_ptr<Parameters::ParameterList> ParamListPtr; 

} // namespace VMC

#endif // PARAMETER_LIST_H
