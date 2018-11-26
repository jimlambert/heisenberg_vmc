#ifndef LOCAL_MEASUREMENT_H
#define LOCAL_MEASUREMENT_H

#include <vector>

namespace VMC {

class LocalMeasurement {
  private:
    std::vector<double> _localvals;
    size_t _nmeas;
    double _total;
  public:
    LocalMeasurement() {}
    void push(double val) {
      _localvals.push_back(val);
      _total += val;
      _nmeas += 1;
    }
    double ave() {return _total/(double)_nmeas;}
    double operator [](const size_t& i) const {return _localvals[i];}
};  
}

#endif // LOCAL_MEASUREMENT_H
