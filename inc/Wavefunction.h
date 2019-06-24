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
    std::random_device                     _rd;
    std::mt19937                           _mteng{_rd()};
    std::uniform_real_distribution<double> _rand_num{0.0, 1.0};
    std::uniform_int_distribution<int>*    _rand_pos;

    // member variables --------------------------------------------------------
    size_t            _size; // size of the spin model
    double            _jas_sum; // running total of Jastrow factors 
    Eigen::MatrixXcd  _gmat; // Green's function matrix 
    ParamListSPtr     _par_lst_ptr; // pointer to parameter list
    AuxHamSPtr        _aux_ham_ptr; // pointer to auxiliary Hamiltonian
    BasisState        _state; // spin state with order of operators    
    ObsUVec           _obsvec; // functors for observables, _obsvec[0]=energy
    
    // output streams ----------------------------------------------------------
    std::ofstream _varfile_output;
    std::ofstream _obsfile_output;
    std::string   _varfile_name;
    std::string   _obsfile_name;
    void          _setup_varfile();    
    void          _update_varfile();
    void          _update_obsfile(const size_t&);

    // setup functions ---------------------------------------------------------
    void   _init_state();
    void   _reinit_gmat();
    void   _clear_obs();
    void   _compute_jas_sum();
    double _dj_sum(const size_t&, const int&);

    // update functions --------------------------------------------------------
    void _flipspin(const int&);
    void _exchange(const int&, const int&);
    void _sweep();
    void _update_params(const double&);
  public:
    SpinWavefunction(
      const size_t&,      // system size
      ParamListSPtr&,     // unique pointer to parameter list
      AuxHamUPtr&,        // unique pointer to auxiliary Hamiltonian
      ObsUPtr&,           // unique pointer to energy functor
      const std::string&, // name of file for variational parameters
      const std::string&  // name of file for observables
    );

    void add_obs(ObsUPtr& ptr){_obsvec.push_back(std::move(ptr));}

    void optimize(const size_t&, const size_t&, const size_t&, const double&);
    void sample(const size_t&, const size_t&, const size_t&);

    // output functions -----------------------------------------------------------
    void print_gmat() {std::cout << _gmat << std::endl;}
    void print_state(); 
    void print_spins();
};

// =============================================================================

} // namespace Wavefunctions
} // namespace VMC

#endif // WAVEFUNCTION_H
