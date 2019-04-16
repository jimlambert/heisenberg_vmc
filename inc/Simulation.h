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
    ParamList_t _params; // variational parameter list from BCS wavefunction
    JspList_t _jsparams; // list of variational parameters from Jastrow factors
    double _jsf; // Jastrow spin factor which needs to be tracked
    std::vector<int> _spinstate; // spin state in the S_z basis
    std::vector<size_t> _operslist; // positions of creation operators
    BCSChainHamiltonian _auxham; // auxiliary Hamiltonian
    Eigen::MatrixXcd _gmat; // Green's function matrix
    LocalMeasurement<std::complex<double> > _el; // local energy

    // help functions
    void _genstate();   // generate a random state
    bool _updateparams(const double&); 
  public:
    
    size_t _exchange(const size_t&, const size_t&); // flip two spins together
    size_t _flipspin(const size_t&); // single spin flip operation
    void _reinitgmat(); // reinitialize _gmat every sweep
    void _sweep();
    
    HeisenbergChainSimulator(const size_t&, const size_t&, ParamList_t&);
    HeisenbergChainSimulator(const size_t&, const size_t&, ParamList_t&, JspList_t&);
    
    void optimize(const size_t&, const size_t&, const size_t&, 
                  const double&, const std::string);
    std::complex<double> _isingenergy();
    std::complex<double> _heisenergy();
    
    void print_auxham(){_auxham.print_matrix();}
    void print_eigvecs(){std::cout << _auxham.get_eigenvecs() << std::endl;}
    void print_eigvals(){std::cout << _auxham.get_eigenvals() << std::endl;}
    void print_redmat(){std::cout << _auxham.get_reduced_matrix(_size) << std::endl;}
    void print_spinstate();
    void print_operslist();
    void print_gmat(){std::cout << _gmat << std::endl;}
    double jsf(){return _jsf;}
};

}
typedef SIMULATION::HeisenbergChainSimulator HeisChainSim;
}

#endif // SIMULATION_H
