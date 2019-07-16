#ifndef BASISSTATE_H
#define BASISSTATE_H

#include<vector>
#include<iostream>
#include<Eigen/Dense>

namespace VMC {

// BasisState class
// =============================================================================
// The only real point of having a separate class for this is that it keeps the
// operlist and spinlist in the same place. Both of these things are initialized
// by the wavefunction and whatever degrees of freedom you want are stored here.  
// =============================================================================

class BasisState {
  
  typedef std::vector<int>         state_t;
  typedef state_t::reference       stateref_t;
  typedef state_t::const_reference const_stateref_t;
  typedef state_t::size_type       sindex_t;
  typedef std::vector<size_t>      list_t;
  typedef list_t::reference        listref_t;
  typedef list_t::const_reference  const_listref_t;
  typedef list_t::size_type        lindex_t;
  
  private:
    size_t              _size;      // system size in reduced bases
    state_t             _spinstate; // spin state in the S_z basis
    list_t              _operslist; // positions of creation operators
  public:
    
    BasisState(const size_t&);
    size_t find (const size_t&) const; // method returns position of lindex
    void refresh(); // resest lists to have all zeros

    // access for spin state ---------------------------------------------------
    stateref_t       operator [](const sindex_t& i) {return _spinstate[i];}
    const_stateref_t operator [](const sindex_t& i) const 
      {return _spinstate[i];} 
    stateref_t site(const sindex_t& i) {return _spinstate[i];}
    const_stateref_t site(const sindex_t& i) const {return _spinstate[i];}
    // access for operator list ------------------------------------------------
    listref_t        operator ()(const lindex_t& i) {return _operslist[i];}
    const_listref_t  operator ()(const lindex_t& i) const 
      {return _operslist[i];}
    listref_t list(const lindex_t& i) {return _operslist[i];}
    const_listref_t list(const lindex_t& i) const {return _operslist[i];}
    // access for members ------------------------------------------------------
    size_t size() const {return _size;}
};

// =============================================================================

} // VMC namespace

#endif // BASSISTATE_H
