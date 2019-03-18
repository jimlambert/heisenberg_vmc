#ifndef SIMULATION_H
#define SIMULATION_H

#include <vector>
#include <string>
#include <Eigen/Dense>
#include <random>
#include <fstream>
#include "LocalMeasurement.h"
#include "BCSHamiltonian.h"
#include "VarParam.h"

namespace VMC {
namespace SIMULATION {

class HeisenbergChainSimulator {
  private:
    //setup random number generation
    std::random_device _rd;
    std::mt19937 _mteng{_rd()};
    std::uniform_real_distribution<double> _rnum{0.0, 1.0};
    std::uniform_int_distribution<int>* _rpos; // random position
    
    // simulation parameters
    size_t _size;
    ParamList_t _params; // variational parameter list
    std::vector<int> _spinstate; // spin state
    std::vector<size_t> _operslist; // positions of creation operators
    BCSChainHamiltonian _auxham; // auxiliary Hamiltonian
    Eigen::MatrixXcd _gmat; // Green's function matrix
    LocalMeasurement<double> _el{1000}; // local energy

    // help functions
    void _genstate();   // generate a random state
    void _reinitgmat(); // reinitialize _gmat every sweep
    void _updateparams(const double&); 
    double _isingenergy();
  public:
    HeisenbergChainSimulator(const size_t&, ParamList_t&);
    // optimize function accepts the number of variatonal steps, the number of
    // equilibrations per step, the number of configurations to sample per step,
    // and the size of each variational step.
    size_t _flipspin(); // single spin flip operation
    size_t _flipspin(const size_t&); // single spin flip operation
    void optimize(const size_t&, const size_t&, const size_t&, const double&);
    void _sweep();
    void print_spinstate();
    void print_operslist();
    void print_gmat(){std::cout << _gmat << std::endl;}
};

}
typedef SIMULATION::HeisenbergChainSimulator HeisChainSim;
}

#endif // SIMULATION_H
