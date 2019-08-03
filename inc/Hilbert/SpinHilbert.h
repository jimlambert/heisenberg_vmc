#ifndef VMC_HILBERT_SPIN_HILBERT_H
#define VMC_HILBERT_SPIN_HILBERT_H

#include <random>
#include <vector>
#include "typedefs.h"

namespace VMC {
namespace Hilbert {

class SpinHilbert {

  
  private:
    
    std::random_device _rd;
    std::mt19937 _mteng{_rd()}; 
    random_double_t _random_double{0.0, 1.0};
    
    double _max_spin;
    size_t _nstates;
    std::vector<double> _states;
  
  public:
    
    SpinHilbert(const double& max_spin);
    
    // returns a random allowed spin
    double random_spin();

    double state(const size_t& i) { return _states[i]; }
      
    // returns the maximum spin;
    double max_spin() const { return _max_spin; }
    size_t nstates() const { return _nstates; }
};

} // namespace Utils
} // namespace VMC

#endif //  VMC_HILBERT_SPIN_HILBERT_H

