#ifndef SIMULATION_H
#define SIMULATION_H

#include <vector>
#include <Eigen/Dense>
#include <random>
#include "BCSHamiltonian.h"
#include "VarParam.h"

namespace VMC {
namespace SIMULATION {

class HeisenbergChainSimulator {
  private:
    std::random_device _rd;
    std::mt19937 _mteng{_rd()};
    std::uniform_real_distribution<double> _rnum{0.0, 1.0};
    std::uniform_int_distribution<int>* _rpos; // random position
    size_t _size;
    ParamList_t params;
    std::vector<int> _spinstate; // spin state
    std::vector<size_t> _operslist; // positions of creation operators
    BCSChainHamiltonian _auxham; // auxiliary Hamiltonian
    Eigen::MatrixXd _gmat; // Green's function matrix
    void _genstate();   // generate a random state
    void _reinitgmat(); // reinitialize _gmat every few steps
  public:
    HeisenbergChainSimulator(const size_t&, const ParamList_t&);
    void _flipspin();   // single spin flip operation
    void print_spinstate();
    void print_operslist();
};

}
typedef SIMULATION::HeisenbergChainSimulator HeisChainSim;
}

#endif // SIMULATION_H
