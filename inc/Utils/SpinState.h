#ifndef VMC_UTILS_SPIN_STATE_H
#define VMC_UTILS_SPIN_STATE_H

#include <random>
#include <vector>
#include "typedefs.h"
#include "SpinHilbert.h"

namespace VMC {
namespace Utils {

class SpinState {


  private:
     
    u_int _phys_size;
    u_int _extd_size;
    Hilbert::SpinHilbert* _hilbert_space_ptr; 
    std::vector<double> _spins_list;
    std::vector<u_int> _opers_list;

  public:

    SpinState
    (
      const u_int& size, 
      Hilbert::SpinHilbert& hilbert_space
    ) : _phys_size(size), _hilbert_space_ptr(&hilbert_space) {

      // for each spin state we need another copy of the lattice      
      _extd_size=_phys_size * (_hilbert_space_ptr->nstates());
     
      // generate a random spin state
      random_state(); 
    }  

    // Generate a random state
    void random_state();
   
    // Find the extended position of the lth electron 
    u_int find(const u_int&) const;

    // This function retrieves the position of the spinless fermion representing
    // the spin at physical position i.
    u_int get_extd_position(const u_int&);
    u_int get_extd_position(const u_int&, const double&);

    // Return state properties
    double max_spin() const { return _hilbert_space_ptr->max_spin(); }
    u_int phys_size() const { return _phys_size; }
    u_int extd_size() const { return _extd_size; }
    
    // Access function for _spins_list
    const double& spin_list(const u_int i) const { return _spins_list[i]; }
    double& spin_list(const u_int i) { return _spins_list[i]; }

    // Access function for _opers_list
    const u_int& oper_list(const u_int i) const { return _opers_list[i]; }
    u_int& oper_list(const u_int i) { return _opers_list[i]; } 

    // Functions to display internal state
    void report_spin_state() const;
    void report_oper_state() const;
};

} // namespace Utils
} // namespace VMC

#endif // VMC_UTILS_SPIN_STATE_H
