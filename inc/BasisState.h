// -----------------------------------------------------------------------------
// CLASS: BASISSTATE
// -----------------------------------------------------------------------------
// BasisState class for VMC algorithm. This class wraps to basis vectors,
// _spinstate, and _operslist. _spinstate is accessed by [], and _operslist is
// accessed by ().
// -----------------------------------------------------------------------------

#ifndef BASISSTATE_H
#define BASISSTATE_H

#include<vector>
#include<iostream>
#include<Eigen/Dense>

namespace VMC {
  
class BasisState {
  typedef std::vector<int> state_t;
  typedef state_t::reference stateref_t;
  typedef state_t::const_reference const_stateref_t;
  typedef state_t::size_type sindex_t;
  typedef std::vector<size_t> list_t;
  typedef list_t::reference listref_t;
  typedef list_t::const_reference const_listref_t;
  typedef list_t::size_type lindex_t;
  private:
    size_t            _size;      // system size in reduced bases
    state_t           _spinstate; // spin state in the S_z basis
    list_t            _operslist; // positions of creation operators
    Eigen::MatrixXcd* _gmat_ptr;  // pointer to current Greens function
  public:
    
    BasisState(const size_t&);
    size_t find(const size_t&); // method returns position of lindex
 
    // access operators for spin state
    // -----------------------------------------------------------------------
    stateref_t       operator [](const sindex_t& i) {return _spinstate[i];}
    const_stateref_t operator [](const sindex_t& i) const 
      {return _spinstate[i];} 
    // -----------------------------------------------------------------------   
    // access operators for operator list
    // -----------------------------------------------------------------------
    listref_t        operator ()(const lindex_t& i) {return _operslist[i];}
    const_listref_t  operator ()(const lindex_t& i) const 
      {return _operslist[i];}
    // -----------------------------------------------------------------------
    size_t size() {return _size;} 
    Eigen::MatrixXcd* gmat_ptr() {return _gmat_ptr;}
};

} // VMC namespace

#endif // BASSISTATE_H
