#ifndef COORDINATE_PAIR_H
#define COORDINATE_PAIR_H

// Coordinate
// =============================================================================
// class that is used to define the actiion of variational parameters etc...
// =============================================================================

namespace VMC {
namespace Utils {

struct CoordinatePair {
  std::vector<int> coord1;
  std::vector<int> coord2;
};


} // namespace Utils 
} // namespace VMC

#endif // COORDINATE_PAIR_H
