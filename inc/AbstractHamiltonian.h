#ifndef ABSTRACT_HAMILTONIAN_H
#define ABSTRACT_HAMILTONIAN_H

#include <memory>
#include "AuxiliaryParameter.h"

namespace VMC {
namespace AuxiliaryHamiltonians {

// Parent class for auxiliary Hamiltonians
// ============================================================================= 

class AuxiliaryHamiltonian {
  public:
    virtual void solve()=0;
    virtual void set_vmats(AuxParamUVec&)=0;
    virtual void set_mmats(AuxParamUVec&)=0;
    virtual void init(AuxParamUVec&)=0;
    virtual Eigen::MatrixXcd get_reduced_matrix(const size_t&)=0;
    virtual Eigen::MatrixXcd get_matrix()=0;
    virtual Eigen::MatrixXcd get_eigenvectors()=0;
    virtual Eigen::VectorXd  get_eigenvalues()=0;
};

// ============================================================================= 

} // namespace AuxiliaryHamiltonians

typedef AuxiliaryHamiltonians::AuxiliaryHamiltonian AuxHam;
typedef std::unique_ptr<AuxHam>                     AuxHamUPtr;
typedef std::shared_ptr<AuxHam>                     AuxHamSPtr;

} // namespace VMC

#endif // ABSTRACT_HAMILTONIAN_H
