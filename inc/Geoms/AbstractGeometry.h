#ifndef ABSTRACT_GEOMETRY_H
#define ABSTRACT_GEOMETRY_H

namespace VMC {
namespace Geometries {

struct AbstractGeometry {
  std::vector<size_t> ndims;
  AbstractGeometry(const std::vector<size_t>& dim_vec) {
    for(auto it= dim_vec.begin(); it!=dim_vec.end(); it++) 
      ndims.push_back(*it);
  }
  size_t operator(std::vector<size_t>&)()=0;
};

}
}

#endif // ABSTRACT_GEOMETRY_H
