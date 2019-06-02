#ifndef WAVEFUNCTION_H
#define WAVEFUNCTION_H

// VMC::Wavefunctions
// =============================================================================
//  This namespace includes the wavefuntion class for the VMC. In order to allow
// transparent and customizable construction, some of the wavefunctions
// constituent pieces are passed as unique_ptrs which are then moved into the
// ownership of the wavefunction. The purpose of this is to allow the
// wavefunction class to be as versatile as possible. 
// =============================================================================

#include <memory>
#include "ParameterList.h"
#include "AuxiliaryHamiltonian.h"
#include "BasisState.h"
#include "VarParam.h"

namespace VMC {
namespace Wavefunctions {

class Wavefunction {
  private:
    ParameterList _param_list;
    AuxHamPtr     _aux_ham;
    BasisState    _state;              
  public:
};


}
}
#endif // WAVEFUNCTION_H
