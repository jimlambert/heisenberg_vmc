#ifndef PARAMETER_LIST_H
#define PARAMETER_LIST_H

#include <vector>
#include <memory>
#include "AuxiliaryParameter.h"
#include "JastrowParameter.h"

namespace VMC {
namespace Parameters {

// Parameter list class used by the wavefunction
// =============================================================================
// This class is used to store the variational parameters used in the VMC.
// All parameters are stored at the end of a unique_ptr. The parameters may be
// constructed through a pointer that is outside the ParameterList and then
// added using the push_... functions, or the parameters may be built directly
// inside the the parameter list sing the build_... functions (recommended).
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
    void build_jas_param(
      ParameterSubtype   jpt, // subtype for Jastrow parameter
      const std::string& n,   // parameter name 
      const double&      v,   // value of parameter
      const size_t&      s,   // parameter site i
      const size_t&      d,   // parameter site j
      const bool&        ti,  // translation invariant (true)
      const size_t&      bs   // binsize for local_meas 
    );

    // access functions
    size_t njas() {return _njas;}
    size_t naux() {return _naux;}
    size_t npar() {return _njas+_naux;}
    AuxParamUVec& aux_vec(){return _aux_params;}
    JasParamUVec& jas_vec(){return _jas_params;}
    AuxiliaryParameter& aux(const size_t& i){return *_aux_params[i];}
    JastrowParameter&   jas(const size_t& i){return *_jas_params[i];}
   
    // useful when constructing S_kk'
    Parameter& operator[](const size_t&);
   
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

typedef std::unique_ptr<Parameters::ParameterList> ParamListUPtr; 
typedef std::shared_ptr<Parameters::ParameterList> ParamListSPtr; 

} // namespace VMC

#endif // PARAMETER_LIST_H
