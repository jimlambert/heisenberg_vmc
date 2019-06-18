#ifndef WAVEFUNCTION_H
#define WAVEFUNCTION_H

#include <iostream>
#include <string>
#include <memory>
#include <fstream>
#include <iomanip>
#include <random>
#include "iomacros.h"
#include "ParameterList.h"
#include "AbstractHamiltonian.h"
#include "AbstractObservable.h"
#include "BasisState.h"

namespace VMC {
namespace Wavefunctions {

// VMC::Wavefunctions
// =============================================================================
//  This namespace includes the wavefuntion class for the VMC. In order to allow
// transparent and customizable construction, some of the wavefunctions
// constituent pieces are passed as unique_ptrs which are then moved into the
// ownership of the wavefunction. The purpose of this is to allow the
// wavefunction class to be as versatile as possible. 
//  The auxiliary Hamiltonian must already have been built with the variational
// parameters passed to the wavefunction (this is likely not an ideal setup). 
// =============================================================================

class SpinWavefunction {
  private: 
    // random number stuff -----------------------------------------------------
    std::random_device _rd;
    std::mt19937      _mteng{_rd()};
    std::uniform_real_distribution<double> _rnum{0.0, 1.0};
    std::uniform_int_distribution<int>*    _rpos;

    // member variables --------------------------------------------------------
    size_t            _size; // size of the spin model
    double            _delta; // change in variational parameters
    double            _jas_sum; // running total of Jastrow factors 
    Eigen::MatrixXcd  _gmat; // Green's function matrix 
    ParamListUPtr     _par_lst_ptr; // pointer to parameter list
    AuxHamSPtr        _aux_ham_ptr; // pointer to auxiliary Hamiltonian
    BasisState        _state; // spin state with order of operators    
    ObsUVec           _obsvec; // functors for observables, _obsvec[0]=energy
    
    // output streams ----------------------------------------------------------
    std::ofstream _var_param_output;
    std::ofstream _obs_output;
    std::string   _var_param_file_name;
    std::string   _obs_file_name;
    
    // setup functions ---------------------------------------------------------
    void _init_state();
    void _reinit_gmat();
    void _compute_jas_sum();

    // update functions --------------------------------------------------------
    void _flipspin(const size_t&);
    void _exchange(const size_t&, const size_t&);
    void _sweep();
    void _update_params(const double&);
  public:
    SpinWavefunction(
      const size_t&,      // system size
      const double&,      // delta in the variational update
      ParamListUPtr&,     // unique pointer to parameter list
      AuxHamUPtr&,        // unique pointer to auxiliary Hamiltonian
      ObsUPtr&,           // unique pointer to energy functor
      const std::string&, // name of file for variational parameters
      const std::string&  // name of file for observables
    );

    void add_obs(ObsUPtr& ptr){_obsvec.push_back(std::move(ptr));}

    void optimize(const size_t&, const size_t&, const size_t&);
    void sample(const size_t&, const size_t&, const size_t&);

    // output functions -----------------------------------------------------------
    void print_gmat() {std::cout << _gmat << std::endl;}
    void print_state(); 
};

// =============================================================================

} // namespace Wavefunctions
} // namespace VMC

#endif // WAVEFUNCTION_H
