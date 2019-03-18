#ifndef LOCAL_MEASUREMENT_H
#define LOCAL_MEASUREMENT_H

#include <vector>

namespace VMC {

template <class T>
class LocalMeasurement {
  private:
    //std::vector<double> _binvals; // binned values
    std::vector<T> _vals; // directly measured values
    size_t _binsize;
    size_t _nmeas;
    size_t _nbins;
    T _total;
    T _bintotal;
  public:
    LocalMeasurement(const size_t& b) : _binsize(b), _nmeas(0), _nbins(0), 
                                        _total(0.0) {}
    void push(T val) {
      _vals.push_back(val);
      _total += val;
      _nmeas += 1;
    }
    size_t nbins() {return _nbins;}
    size_t nmeas() {return _nmeas;}
    //double ave() {return (double)_bintotal/(double)_nbins;}
    T ave() {return (T)_total/(T)_nmeas;}
    //double operator [](const size_t& i) const {return _binvals[i];}
    T operator [](const size_t& i) const {return _vals[i];}
};  
}

#endif // LOCAL_MEASUREMENT_H
