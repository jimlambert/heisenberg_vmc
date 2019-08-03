#ifndef VMC_UTILS_UPDATE_INFO_H
#define VMC_UTILS_UPDATE_INFO_H

#include <vector>

namespace VMC {
namespace Utils {

struct UpdateInfo {
  std::vector<int> spin_change;
  std::vector<int> position_indices;
  std::vector<int> electron_indices;
};

} // namespace Utils
} // namespace VMC

#endif // VMC_UTILS_UPDATE_INFO_H
