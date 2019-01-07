#include <algorithm>
#include <complex>
#include <random>
#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <Eigen/Dense>
#include <cmath>
#include "LocalMeasurement.h"
#include "Simulation.h"
#include "VarParam.h"

namespace VMC {
namespace SIMULATION {

HeisenbergChainSimulator::HeisenbergChainSimulator
(const size_t& N, ParamList_t& params) 
: _size(N), _auxham(2*N, params) {
  // THIS SOLUTION SUCKS
  for(size_t i=0; i<params.size(); i++) _params.push_back(params[i]);
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
  _auxham.setopers(_params);
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

void HeisenbergChainSimulator::_updateparams(const double& df) {
  size_t N=_params.size();
  Eigen::MatrixXd S(N,N);
  Eigen::VectorXd F(N); 
  Eigen::VectorXd dA(N);
  for(size_t ka=0; ka<N; ka++) {  
    size_t na=_params[ka].lmeas.nmeas();
    for(size_t kb=ka; kb<N; kb++) {
      size_t nb=_params[kb].lmeas.nmeas();
      if(nb==0) {
        std::cout << "ERROR: local measurements missing for parameter, "
                  << _params[kb].name << std::endl;
        exit(1);
      }
      if(na!=nb) {
        std::cout << "ERROR: local measurements, " << _params[kb].name << " and "
                  << _params[ka].name << " are different sizes." << std::endl;
        exit(1);
      }
      double avea=_params[ka].lmeas.ave();
      double aveb=_params[kb].lmeas.ave();
      for(size_t i=0; i<na; i++) {
        double ai=_params[ka].lmeas[i];
        double bi=_params[ka].lmeas[i];
        S(ka, kb)+=(ai-avea) * (bi-aveb);
      }
      S(ka,kb)=S(ka,kb)/(double)(_params[ka].lmeas.nmeas());
      kb++;
    }
    for(size_t i=0; i<na; i++) F(ka)+=_el[i]*(_params[ka].lmeas[i] - _params[ka].lmeas.ave());
    F(ka)=-2*F(ka)/(double)na;
    ka++;
  }
  // preconditioning
  Eigen::MatrixXd S_pc(N,N);
  Eigen::VectorXd F_pc(N); 
  for(size_t ka=0; ka<N; ka++) {
    F_pc(ka)=F(ka)/std::sqrt(S(ka,ka));
    for(size_t kb=0; kb<N; kb++) {
      S_pc(ka,kb)=S(ka,kb)/(std::sqrt(S(ka,ka)*S(kb,kb)));
    }
  }
  F=df*F;
  dA = S_pc.inverse()*F_pc;
  for(size_t k=0; k<N; k++) _params[k].val+=dA[k]/S(k,k); 
}

size_t HeisenbergChainSimulator::_flipspin(const size_t& rpos) {
  size_t iexpos; // initial position in extended space
  size_t fexpos; // final position in extended space
  size_t lindex; // index of position in W array
  int ds; // change in spin direction
  // determine initial and final positions on extended lattice
  if(_spinstate[rpos]==1){
    iexpos=rpos;
    fexpos=iexpos+_size;
    ds=-2;
  }
  else{
    iexpos=rpos+_size;
    fexpos=iexpos-_size;
    ds=2;
  }
  // determine position of associated creation operator
  lindex=_operslist[iexpos]-1;
  double amp=abs(_gmat(fexpos, lindex)); // extract Green's function
  // accept flip according to green function
  double rnum=_rnum(_mteng);
  if(rnum<amp*amp){
    // update spin state and extended state
    _operslist[fexpos]=_operslist[iexpos];
    _operslist[iexpos]=0;
    _spinstate[rpos]+=ds;
    // update green's function matrix
    Eigen::VectorXd a(2*_size);
    Eigen::VectorXd b(_size);
    for(size_t i=0; i<2*_size; i++) a(i)=_gmat(i, lindex);
    for(size_t i=0; i<_size; i++) {
      double d=0;
      if(i==lindex) d=1;
      b(i)= -1.0*(_gmat(fexpos, i) - d)/_gmat(fexpos, lindex);
    }
    _gmat=_gmat+a*b.transpose();
    return 1;
  }
  return 0;
}

void HeisenbergChainSimulator::_sweep() {
  // flip until N electron moves have occured.
  // while(counter<_size) counter+=_flipspin();
  for(size_t i=0; i<_size; i++) {
    int rpos=(*_rpos)(_mteng);
    size_t n=_flipspin(rpos);
  } 
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

void HeisenbergChainSimulator::optimize
(const size_t& vsteps, const size_t& equil, const size_t& simul, 
 const double& df) {

  // open file for variational parameters
  std::ofstream var_params;
  std::string fvar_params="var_params.dat";
  var_params.open(fvar_params);
  for(auto it=_params.begin(); it!=_params.end(); it++) 
    var_params << std::setw(10) << std::left << it->name;
  var_params << '\n';
  for(auto it=_params.begin(); it!=_params.end(); it++) 
    var_params << std::setw(10) << std::left << it->val;
  var_params << '\n';

  for(size_t vstep=0; vstep<vsteps; vstep++) {
    // open file for local observables
    std::ofstream measvals;
    std::string fmeasvals="measvals"+std::to_string(vstep)+".dat";
    measvals.open(fmeasvals);  
    // first equilibrate this particular set
    for(size_t i=0; i<equil; i++) _sweep();
    std::cout << "equilibration " << vstep << " done." << std::endl;
    // sample current wavefunction
    for(size_t i=0; i<simul; i++) { 
      _sweep();
      double e=_isingenergy();
      _el.push(e);
      measvals << e << '\n';
      // Calculate O_k(x) for each variational parameter
      for(auto it=_params.begin(); it!=_params.end(); it++) {
        Eigen::MatrixXd redMat(_size, 2*_size); // N_e X 2L matrix
        Eigen::MatrixXd okmat; // operator to measure
        for(size_t l=0; l<_size; l++)
        for(size_t r=0; r<2*_size; r++) {
          auto itstart=_operslist.begin();
          auto itend=_operslist.end();
          size_t lpos=std::distance(itstart, std::find(itstart, itend, l+1));
          redMat(l,r) = it->mmat(lpos, r);  
        }
        // compute actual local value
        okmat=redMat*_gmat;
        // add measurement to local values for this operator
        it->lmeas.push(okmat.trace());
      }
    }
    measvals.close();
    measvals.clear();
    // need to actually adjust variational energies now. First construct matrix
    // S_kk' and force vector f_k 
    _updateparams(df);
    // update variational parameters
    std::cout << "outputing new variational parameters" << std::endl;
    for(auto it=_params.begin(); it!=_params.end(); it++) 
      var_params << std::setw(10) << std::left << it->val; 
    var_params << '\n';
    std::cout << "reinitialzing auxiliary Hamiltonian" << std::endl;
    _auxham.init(_params);
  }
  var_params.close();
  var_params.clear();
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
