#ifndef PARAMETER_LIST_H
#define PARAMETER_LIST_H

#include <vector>
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
    AuxParamSVec _aux_params;
    JasParamSVec _jas_params;
  public:
    ParameterList() : _naux(0), _njas(0) {}
    void push_aux_param(JasParamUPtr);
    void push_jas_param(AuxParamUPtr);
    size_t njas() {return _njas;}
    size_t naux() {return _naux;}
    size_t npar() {return _njas+_naux;}
    AuxParamSVec aux_vec();
    JasParamSVec jas_vec();
    AuxParamSPtr aux(const size_t&);
    JasParamSPtr jas(const size_t&);
    void report_params();
    void report_jas_params();
    void report_aux_params();
};

// =============================================================================

} // namespace Parameters

typedef std::unique_ptr<Parameters::ParameterList> ParamListPtr; 

} // namespace VMC

#endif // PARAMETER_LIST_H
