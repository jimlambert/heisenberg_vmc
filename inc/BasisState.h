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
      size_t _L;
      state_t _spinstate; // spin state in the S_z basis
      std::vector<size_t> _operslist; // positions of creation operators
    public:
      BasisState(const size_t&);
      size_t find(const size_t& lindex); // method returns position of lindex
      // access operators ------------------------------------------------------
      stateref_t operator [](const sindex_t& i) {return _spinstate[i];}
      const_stateref_t operator [](const sindex_t& i) const 
      {return _spinstate[i];} 
      listref_t operator ()(const lindex_t& i) {return _operslist[i];}
      const_listref_t operator()(const lindex_t& i) const {return _operslist[i];}
      // -----------------------------------------------------------------------
  };
}

#endif // BASSISTATE_H
