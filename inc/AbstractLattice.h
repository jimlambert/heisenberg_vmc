#ifndef ABSTRACT_LATTICE_H
#define ABSTRACT_LATTICE_H

#include "BasisState.h"
#include <vector>

namespace VMC {
namespace Lattices {

class AbstractLattice {
  private:
    size_t              _dim;    // dimension of lattice
    size_t              _nsites; // number of sites in lattice
    size_t              _nbasis; // number of sites in basis
    std::vector<size_t> _ndim;   // number of sites along each dimension
    std::vector<bool>   _pbc;    // boundary conditions along each dimension 
    BasisState          _state;  // current basis state

    // these functions take a site in the basis and the coordinates of a unit
    // cell to a single site in the state vector and the inverse 
    virtual size_t _lattice2coordinate(const size_t&, const std::vector<int>&)=0;
    virtual size_t _coordinate2lattice(const size_t&)=0;
  public:
    


};

} // namespace VMC
} // namespace Lattices

#endif // ABSTRACT_LATTICE_H
