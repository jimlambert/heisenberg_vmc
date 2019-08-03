#ifndef LOCAL_MEASUREMENT_H
#define LOCAL_MEASUREMENT_H

#include <vector>
#include <iostream>

namespace VMC {

template <class T>
class LocalMeasurement {
  private:
    using vec_t = std::vector<T>;
    vec_t _binaves; // average bin values
    vec_t _vals; // average bin values
    size_t _bsize;
    size_t _nmeas;
    size_t _nbins;
    size_t _nvals;
    T _total;
  public:
    typedef typename vec_t::iterator vecIter;
    typedef typename vec_t::const_iterator constVecIter;
    LocalMeasurement(const size_t& b) : _bsize(b), _nmeas(0), _nbins(0), 
                                        _nvals(0), _total(0.0) {}
    void push(T val) {
      if(_nmeas<(_bsize-1)) {
        _vals.push_back(val);
        _total += val;
        _nmeas += 1;
        _nvals += 1;
      }
      else {
        _binaves.push_back((_total + val)/(double)_bsize);
        _vals.push_back(val);
        _total = 0.0;
        _nmeas = 0;
        _nvals += 1;
        _nbins += 1;
      }
    }
    
    void clear() { 
      _total=0.0;
      _nmeas=0;
      _nbins=0;
      _nvals=0;
      _binaves.clear();
      _vals.clear();
    }
   
    // access members 
    size_t bsize() {return _bsize;}
    size_t nbins() {return _nbins;}
    size_t nmeas() {return _nmeas;}
    size_t nvals() {return _nvals;}
    T ave() {
      T total = 0.0;
      for(auto it=_binaves.begin(); it!=_binaves.end(); it++) total += *it;
      return (T)total/(T)_nbins;
    }

    T var() {
      T total=0.0;
      T average=ave();
      for(auto it=_binaves.begin(); it!=_binaves.end(); it++) {
        total+=(*it-average)*(*it-average);
      }
      return std::sqrt((T)(total)/(T)_nbins);
    }
    
    // define access operators including iterators to underlying container. 
    T operator [](const size_t& i) const {
      if(i>=_vals.size()) {
        std::cout << "ERROR: overflow on measurement vals" << std::endl;
        exit(1);
      }
      else if(i<0) {
        std::cout << "ERROR: underflow on measurement vals" << std::endl;
        exit(1);
      }
      return _vals[i];
    }
    T operator ()(const size_t& i) const {
      if(i>=_binaves.size()) {
        std::cout << "ERROR: overflow on measurement bins" << std::endl;
        exit(1);
      } 
      else if(i<0) {
        std::cout << "ERROR: underflow on measurement bins" << std::endl;
        exit(1);
      }
      return _binaves[i];
    }
    vecIter begin(){return _binaves.begin();}
    vecIter end(){return _binaves.end();}
    constVecIter begin() const {return _binaves.begin();}
    constVecIter end() const {return _binaves.end();}
};  

typedef LocalMeasurement<std::complex<double> > MeasCd;
typedef LocalMeasurement<double>                MeasD;

} // namespace VMC

#endif // LOCAL_MEASUREMENT_H
