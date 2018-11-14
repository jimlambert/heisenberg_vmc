#include <random>
#include <vector>
#include <Eigen/Dense>
#include "Simulation.h"
#include "VarParam.h"

namespace VMC {
namespace SIMULATION {

HeisenbergChainSimulator::HeisenbergChainSimulator
(const size_t& N, const ParamList_t& params) : _size(N), _auxham(2*N, params) {
  _rpos = new std::uniform_int_distribution<>(0, _size-1);
  _spinstate.resize(_size);
  _operslist.resize(2*_size);
  Eigen::MatrixXd _projmat(_size, _size);
  Eigen::MatrixXd _redcmat(2*_size, _size);
  do {
    _genstate();
  } while(_projmat.determinant()==0);
}

void HeisenbergChainSimulator::_genstate() {
  for(size_t i=0; i<_size; i++) {
    double r = _rnum(_mteng);
    if(r < 0.5) { 
      _spinstate[i] = 1;
      _operslist[i] = i;
    } 
    else {
      _spinstate[i] = -1;
      _operslist[i+_size] = i;
    }  
  }
}

}
}
