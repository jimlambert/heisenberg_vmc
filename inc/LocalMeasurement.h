#ifndef LOCAL_MEASUREMENT_H
#define LOCAL_MEASUREMENT_H

#include <vector>

namespace VMC {

class LocalMeasurement {
  private:
    //std::vector<double> _binvals; // binned values
    std::vector<double> _vals; // directly measured values
    size_t _binsize;
    size_t _nmeas;
    size_t _nbins;
    double _total;
    double _bintotal;
  public:
    LocalMeasurement(const size_t& b) : _binsize(b), _nmeas(0), _nbins(0), 
                                        _total(0.0) {}
    void push(double val) {
      //if(_nmeas<_binsize-1) {
      //  _total+=val;
      //  _nmeas+=1;
      //} 
      //else {
      //  _binvals.push_back((double)_total/(double)_nmeas);
      //  _bintotal+=(double)_total/(double)_nmeas;
      //  _nmeas=0;
      //  _total=0;
      //  _nbins+=1;
      //}
      _vals.push_back(val);
      _total += val;
      _nmeas += 1;
    }
    size_t nbins() {return _nbins;}
    size_t nmeas() {return _nmeas;}
    //double ave() {return (double)_bintotal/(double)_nbins;}
    double ave() {return (double)_total/(double)_nmeas;}
    //double operator [](const size_t& i) const {return _binvals[i];}
    double operator [](const size_t& i) const {return _vals[i];}
};  
}

#endif // LOCAL_MEASUREMENT_H
