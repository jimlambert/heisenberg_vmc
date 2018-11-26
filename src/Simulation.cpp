#include <algorithm>
#include <complex>
#include <random>
#include <vector>
#include <fstream>
#include <iostream>
#include <Eigen/Dense>
#include "LocalMeasurement.h"
#include "Simulation.h"
#include "VarParam.h"

namespace VMC {
namespace SIMULATION {

HeisenbergChainSimulator::HeisenbergChainSimulator
(const size_t& N, ParamList_t& params) : _size(N), _auxham(2*N, params) {
  _rpos = new std::uniform_int_distribution<>(0, _size-1);
  _auxham.solve();
  _spinstate.resize(_size);
  _operslist.resize(2*_size);
  Eigen::MatrixXd projmat(_size, _size);
  Eigen::MatrixXd redmat(2*_size, _size);
  redmat=_auxham.get_reduced_matrix(_size);
  do {
    _genstate();
    for(size_t i=0; i<_size; i++) {
      int pos;
      if(_spinstate[i]==1) pos=i;
      else pos=i+_size;
      for(size_t j=0; j<_size; j++) {
        projmat(i, j) = redmat(pos, j);
      }
    }
  } while(projmat.determinant()==0);
  _gmat = redmat*projmat.inverse();
}

void HeisenbergChainSimulator::_genstate() {
  std::fill(_spinstate.begin(), _spinstate.end(), 0);
  std::fill(_operslist.begin(), _operslist.end(), 0);
  for(size_t i=0; i<_size; i++) {
    double r = _rnum(_mteng);
    if(r < 0.5) { 
      _spinstate[i] = 1;
      _operslist[i] = i+1;
    } 
    else {
      _spinstate[i] = -1;
      _operslist[i+_size] = i+1;
    }  
  }
}

void HeisenbergChainSimulator::_reinitgmat() {
  Eigen::MatrixXd projmat(_size, _size);
  Eigen::MatrixXd redmat(2*_size, _size);
  redmat=_auxham.get_reduced_matrix(_size);
  for(size_t i=0; i<_size; i++) {
    int pos;
    if(_spinstate[i]==1) pos=i;
    else pos=i+_size;
    for(size_t j=0; j<_size; j++) {
      projmat(i, j) = redmat(pos, j);
    }
  }
  _gmat = redmat*projmat.inverse();
}

size_t HeisenbergChainSimulator::_flipspin() {
  double rnum=_rnum(_mteng);
  int rpos=(*_rpos)(_mteng);
  size_t exipos;
  size_t exnpos;
  size_t lindex; // index of position in W array
  int ds;
  if(_spinstate[rpos]==1){
    exipos=rpos;
    exnpos=exipos+_size;
    ds=-2;
  }
  else{
    exipos=rpos+_size;
    exnpos=exipos-_size;
    ds=2;
  }
  lindex=_operslist[exipos]-1;
  double amp=abs(_gmat(exnpos, lindex));
  // accept flip according to green function
  if(rnum<amp*amp){
    _operslist[exnpos]=_operslist[exipos];
    _operslist[exipos]=0;
    _spinstate[rpos]+=ds;
    Eigen::VectorXd a(2*_size);
    Eigen::VectorXd b(_size);
    for(size_t i=0; i<2*_size; i++) a(i)=_gmat(i, lindex);
    for(size_t i=0; i<_size; i++) {
      double d=0;
      if(i==lindex) d=1;
      b(i)= -1.0*(_gmat(exnpos, i) - d)/_gmat(exnpos, lindex);
    }
    _gmat=_gmat+a*b.transpose();
    return 1;
  }
  return 0;
}

void HeisenbergChainSimulator::_sweep() {
  size_t counter=0;
  // flip until N electron moves have occured.
  // while(counter<_size) counter+=_flipspin();
  for(size_t i=0; i<_size; i++) counter += _flipspin();
  _reinitgmat();
}

double HeisenbergChainSimulator::_isingenergy() {
  double total=0;
  for(size_t i=0; i<_size; i++) {
    if(i<(_size-1)) total+=_spinstate[i]*_spinstate[i+1];
    else total+=_spinstate[0]*_spinstate[i];
  }
  return total;
}

void HeisenbergChainSimulator::optimize() {
  _equilenergy.open("equilenergy.dat");
  for(size_t i=0; i<10000; i++) {
    std::cout << i << std::endl;
    _sweep();
    _equilenergy << _isingenergy() << '\n';
  }  
  _equilenergy.close();
}

void HeisenbergChainSimulator::print_spinstate() {
  for(size_t i=0; i<_size; i++) std::cout << _spinstate[i] << '\t';
  std::cout << std::endl; 
}

void HeisenbergChainSimulator::print_operslist() {
  for(size_t i=0; i<2*_size; i++) std::cout << _operslist[i] << '\t';
  std::cout << std::endl; 
}

}
}
