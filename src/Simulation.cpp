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

#define COLUMN_WITH 20

HeisenbergChainSimulator::HeisenbergChainSimulator
(const size_t& N, const size_t& bs, ParamList_t& params) 
: _size(N), _el(bs), _auxham(2*N, params) {
  // THIS SOLUTION SUCKS
  for(size_t i=0; i<params.size(); i++) _params.push_back(params[i]);
  _rpos = new std::uniform_int_distribution<>(0, _size-1);
  _auxham.solve();
  _spinstate.resize(_size);
  _operslist.resize(2*_size);
  Eigen::MatrixXcd projmat(_size, _size);
  Eigen::MatrixXcd redmat(2*_size, _size);
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
  } while(abs(projmat.determinant())==0);
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
  Eigen::MatrixXcd projmat(_size, _size);
  Eigen::MatrixXcd redmat(2*_size, _size);
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
  Eigen::VectorXcd dA(N);
  S = Eigen::MatrixXd::Zero(N, N);
  for(size_t ka=0; ka<N; ka++) {  
    size_t na=_params[ka].lmeas.nvals();
    for(size_t kb=ka; kb<N; kb++) {
      size_t nb=_params[kb].lmeas.nvals();
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
      std::complex<double> avea=_params[ka].lmeas.ave();
      std::complex<double> aveb=_params[kb].lmeas.ave();
      std::complex<double> total=0.0;
      for(size_t i=0; i<na; i++) {
        std::complex<double> ai=_params[ka].lmeas[i];
        std::complex<double> bi=_params[kb].lmeas[i];
        total+=(ai-avea) * (bi-aveb);
      }
      S(ka,kb)=(total/(double)na).real();
      S(kb,ka)=S(ka,kb);
    }
    std::complex<double> total=0.0;
    for(size_t i=0; i<na; i++) 
      total+=std::conj(_el[i])*(_params[ka].lmeas[i] - _params[ka].lmeas.ave());
    F(ka)=-2.0*total.real()/(double)na;
  } 
  // ----------------------
  // preconditioning
  // ----------------------
  Eigen::MatrixXd S_pc(N,N);
  Eigen::VectorXd F_pc(N); 
  for(size_t ka=0; ka<N; ka++) {
    F_pc(ka)=F(ka)/std::sqrt(S(ka,ka));
    for(size_t kb=0; kb<N; kb++) {
      S_pc(ka,kb)=S(ka,kb)/(std::sqrt(S(ka,ka)*S(kb,kb)));
    }
  }
  F_pc=df*F_pc;
  dA = S_pc.inverse()*F_pc;
  //std::cout << F << std::endl;
  //dA = df*F;
  for(size_t k=0; k<N; k++) _params[k].val+=dA[k]/S(k,k);
  //dA = S.inverse()*F*df; 
  //for(size_t k=0; k<N; k++) _params[k].val+=dA[k]/S(k,k);
}

size_t HeisenbergChainSimulator::_flipspin() {
  int rpos=(*_rpos)(_mteng);
  size_t iexpos;
  size_t nexpos;
  size_t lindex; // index of position in W array
  int ds; // change in spin direction
  // determine initial and final positions on extended lattice
  if(_spinstate[rpos]==1){
    iexpos=rpos;
    nexpos=iexpos+_size;
    ds=-2;
  }
  else{
    iexpos=rpos+_size;
    nexpos=iexpos-_size;
    ds=2;
  }
  // determine position of associated creation operator
  lindex=_operslist[iexpos]-1;
  double amp=fabs(_gmat(nexpos, lindex)); // extract Green's function
  lindex=_operslist[iexpos]-1;
  // accept flip according to green function
  double rnum=_rnum(_mteng);
  if(rnum<amp*amp){
    // update spin state and extended state
    _operslist[nexpos]=_operslist[iexpos];
    _operslist[iexpos]=0;
    _spinstate[rpos]+=ds;
    // update green's function matrix
    //Eigen::VectorXcd a(2*_size);
    //Eigen::VectorXcd b(_size);
    //for(size_t i=0; i<2*_size; i++) a(i)=_gmat(i, lindex);
    //for(size_t i=0; i<_size; i++) {
    //  double d=0;
    //  if(i==lindex) d=1;
    //  b(i)= -1.0*(_gmat(nexpos, i) - d)/_gmat(nexpos, lindex);
    //}
    //_gmat=_gmat+a*b.transpose();
    Eigen::MatrixXcd newgmat(2*_size, _size);
    for(size_t i=0; i<2*_size; i++)
    for(size_t j=0; j<_size; j++) {
      double d=0;
      if(i==lindex) d=1;
      std::complex<double> den = _gmat(nexpos,lindex);
      std::complex<double> num = _gmat(i,lindex);
      newgmat(i,j)=_gmat(i,j)-(num/den)*(_gmat(nexpos,j)-d); 
    }
    _gmat=newgmat;
    return 1;
  }
  return 0;
}

size_t HeisenbergChainSimulator::_flipspin(const size_t& rpos) {
  size_t iexpos; // initial position in extended basis
  size_t nexpos; // new position in extended basis
  size_t lindex; // index of position in W array
  int ds; // change in spin
  // -----------------------------------
  // determine new location of spin in extended state
  // -----------------------------------
  if(_spinstate[rpos]==1){
    iexpos=rpos;
    nexpos=iexpos+_size;
    ds=-2;
  }
  else{
    iexpos=rpos+_size;
    nexpos=iexpos-_size;
    ds=2;
  }
  // ----------------------------------
  // accept flip according to Green's function
  // ----------------------------------
  lindex=_operslist[iexpos]-1;
  double amp=abs(_gmat(nexpos, lindex));
  double rnum = _rnum(_mteng);
  if(rnum<amp*amp){
    // -----------------------
    // swap operator positions and update spin state
    // -----------------------
    _operslist[nexpos]=_operslist[iexpos]; 
    _operslist[iexpos]=0;
    _spinstate[rpos]+=ds;
    // -----------------------
    // prepare vectors to update Green's function matrix
    // -----------------------
    //Eigen::VectorXcd a(2*_size);
    //Eigen::VectorXcd b(_size);
    //for(size_t i=0; i<2*_size; i++) a(i)=_gmat(i, lindex);
    //for(size_t i=0; i<_size; i++) {
    //  double d=0;
    //  if(i==lindex) d=1; // implementing the delta function
    //  b(i)= -1.0*(_gmat(nexpos, i) - d)/_gmat(nexpos, lindex);
    //}
    //_gmat=_gmat+a*b.transpose(); // update Green's function here
    Eigen::MatrixXcd newgmat(2*_size, _size);
    for(size_t i=0; i<2*_size; i++)
    for(size_t j=0; j<_size; j++) {
      double d=0;
      if(i==lindex) d=1;
      std::complex<double> den = _gmat(nexpos,lindex);
      std::complex<double> num = _gmat(i,lindex);
      newgmat(i,j)=_gmat(i,j)-(num/den)*(_gmat(nexpos,j)-d); 
    }
    _gmat=newgmat;
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

std::complex<double> HeisenbergChainSimulator::_isingenergy() {
  std::complex<double> total=0.0;
  for(size_t i=0; i<_size; i++) {
    if(i<(_size-1)) total+=_spinstate[i]*_spinstate[i+1];
    else total+=_spinstate[0]*_spinstate[i];
  }
  return total;
}

std::complex<double> HeisenbergChainSimulator::_heisenergy() {
  std::complex<double> total=0.0;
  for(size_t i=0; i<_size; i++) {
    size_t j;
    if(i<(_size-1)) j = i+1;
    else j=0;
    total += 0.25*_spinstate[i]*_spinstate[j];
    if(_spinstate[i]==_spinstate[j]) continue;
    else if(_spinstate[i]==1) {
      size_t iexpos1=i;   
      size_t iexpos2=j+_size;   
      size_t nexpos1=i+_size;   
      size_t nexpos2=j;   
      size_t lindex1=_operslist[iexpos1]-1;
      size_t lindex2=_operslist[iexpos2]-1;
      std::complex<double> det=_gmat(nexpos1, lindex1)*_gmat(nexpos2, lindex2)
        -_gmat(nexpos2, lindex1)*_gmat(nexpos1, lindex2);
      total += 0.5 * det;
    }
    else {
      size_t iexpos1=i+_size;   
      size_t iexpos2=j;   
      size_t nexpos1=i;   
      size_t nexpos2=j+_size;   
      size_t lindex1=_operslist[iexpos1]-1;
      size_t lindex2=_operslist[iexpos2]-1;
      std::complex<double> det=_gmat(nexpos1, lindex1)*_gmat(nexpos2, lindex2)
        -_gmat(nexpos2, lindex1)*_gmat(nexpos1, lindex2);
      total += 0.5 * det;
    }
  }
  return total;
}

void HeisenbergChainSimulator::optimize
(const size_t& vsteps, const size_t& equil, const size_t& simul, 
 const double& df, const std::string ofname) {

  // ------------------------------
  // open file for variational parameters, prepare header
  // ------------------------------
  std::ofstream var_params;
  std::string fvar_params="var_params.dat";
  var_params.open(fvar_params);
  for(auto it=_params.begin(); it!=_params.end(); it++) 
    var_params << std::setw(COLUMN_WITH) << std::left << it->name;
  var_params << '\n';
  for(auto it=_params.begin(); it!=_params.end(); it++) 
    var_params << std::setw(COLUMN_WITH) << std::left << it->val;
  var_params << '\n';

  // -----------------------------
  // Main loop for variational steps
  // -----------------------------
  for(size_t vstep=0; vstep<vsteps; vstep++) {
    // open file for local observables
    std::ofstream measfile;
    std::string measfileName=ofname+std::to_string(vstep)+".dat";
    measfile.open(measfileName);  
    measfile << vstep << '\n';
    measfile << std::setw(COLUMN_WITH) << std::left << "e_L";
    for(auto it=_params.begin(); it!=_params.end(); it++) {
      measfile << std::setw(COLUMN_WITH) << std::left << it->name; 
    }
    measfile << '\n'; 
    // first equilibrate this particular set
    for(size_t i=0; i<equil; i++){
      _sweep();
    }
    std::cout << "equilibration " << vstep << " done." << std::endl;    
    // clear any old measurements
    for(auto it=_params.begin(); it!=_params.end(); it++) {
      it->lmeas.clear();
    }
    _el.clear();
    // -------------------
    // Loop for sampling wavefunction
    // -------------------
    for(size_t i=0; i<simul; i++) { 
      // start with a sweep
      _sweep();
      std::complex<double> e=_isingenergy();
      //std::complex<double> e=_heisenergy();
      _el.push(e); // record energy
      // -------------------------
      // Loop through O_k(x) for each variational parameter
      // -------------------------
      for(auto it=_params.begin(); it!=_params.end(); it++) {
        Eigen::MatrixXcd redMat(_size, 2*_size); // N_e X 2L matrix
        Eigen::MatrixXcd okmat; // operator to measure
        for(size_t l=0; l<_size; l++)
        for(size_t r=0; r<2*_size; r++) {
          auto itstart=_operslist.begin();
          auto itend=_operslist.end();
          size_t lpos=std::distance(itstart, std::find(itstart, itend, l+1));
          redMat(l,r) = it->mmat(lpos, r);  
        }
        // compute actual local value
        okmat=redMat*_gmat;
        std::complex<double> okval = okmat.trace();
        it->lmeas.push(okval);
      }
    }
    for(size_t i=0; i<_el.nbins(); i++) {
      measfile << std::setw(COLUMN_WITH) << std::left << _el(i);
      for(auto pit=_params.begin(); pit!=_params.end(); pit++) 
        measfile << std::setw(COLUMN_WITH) << std::left << pit->lmeas(i);
      measfile << '\n';
    }
    measfile.close();
    measfile.clear();
    // update variational parameters
    _updateparams(df);
    std::cout << "outputing new variational parameters" << std::endl;
    for(auto it=_params.begin(); it!=_params.end(); it++) 
      var_params << std::setw(COLUMN_WITH) << std::left << it->val; 
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
