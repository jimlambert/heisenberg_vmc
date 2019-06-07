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
#include <limits>
#include "LocalMeasurement.h"
#include "Simulation.h"
#include "VarParam.h"

namespace VMC {
namespace SIMULATION {

#define COLUMN_WIDTH 20
#define OUTPUT_WIDTH 10

HeisenbergChainSimulator::HeisenbergChainSimulator
(const size_t& N, const size_t& bs, ParamList_t& params) 
: _size(N), _auxham(2*N, params), _el(bs) {
  // THIS SOLUTION SUCKS -------------------------------------------------------
  for(size_t i=0; i<params.size(); i++) _params.push_back(params[i]);
  // ---------------------------------------------------------------------------
  _jsf=1.0;
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
  } while(fabs(projmat.determinant())<1e-12);
  //print_spinstate();
  //print_operslist();
  //std::cout << redmat << std::endl;
  //std::cout << "----" << std::endl;
  //std::cout << projmat << std::endl;
  _gmat = redmat*(projmat.inverse());
}

HeisenbergChainSimulator::HeisenbergChainSimulator
(const size_t& N, const size_t& bs, ParamList_t& params, JspList_t& jsparams) 
: _size(N), _auxham(2*N, params), _el(bs) {
  // THIS SOLUTION SUCKS -------------------------------------------------------
  for(size_t i=0; i<params.size(); i++) _params.push_back(params[i]);
  for(size_t i=0; i<jsparams.size(); i++) _jsparams.push_back(jsparams[i]);
  // ---------------------------------------------------------------------------
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
  } while(fabs(projmat.determinant())<1e-12);
  _gmat = redmat*projmat.inverse();
  // calculate initial value of Jastrow factor ---------------------------------
  double total=0.0;
  for(size_t i=0; i<_size; i++) 
  for(size_t j=0; j<_jsparams.size(); j++) {
    double v=_jsparams[j].val;
    size_t dx=_jsparams[j].space;
    total += 0.25*v*(double)_spinstate[i]*(double)_spinstate[(i+dx)%_size];
  }
  _jsf = std::exp(total);
  // ---------------------------------------------------------------------------
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

bool HeisenbergChainSimulator::_updateparams(const double& df) {
  
  size_t N=_params.size()+_jsparams.size();
  size_t Np = _params.size();
  size_t Nj = _jsparams.size();
  Eigen::MatrixXd S(N,N);
  Eigen::VectorXd F(N); 
  Eigen::VectorXd dA(N);
  S = Eigen::MatrixXd::Zero(N, N);
  
  for(size_t ka=0; ka<Np; ka++) {  
    size_t na=_params[ka].lmeas.nvals();
    // Loop for BCS params -----------------------------------------------------
    for(size_t kb=ka; kb<Np; kb++) {
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
    // Loop for Jastrow factors ------------------------------------------------
    for(size_t j=0; j<_jsparams.size(); j++) {
      std::complex<double> avea=_params[ka].lmeas.ave();
      std::complex<double> aveb=_jsparams[j].lmeas.ave();
      std::complex<double> total=0.0;
      for(size_t i=0; i<na; i++) {
        std::complex<double> ai=_params[ka].lmeas[i];
        std::complex<double> bi=_jsparams[j].lmeas[i];
        total+=(ai-avea) * (bi-aveb);
      }
      S(ka, Np+j)+=(total/(double)na).real();
      S(Np+j, ka)+=S(ka, Np+j);
    }
    // Loop for BCS forces -----------------------------------------------------
    std::complex<double> total=0.0;
    std::complex<double> avea=_params[ka].lmeas.ave();
    for(size_t i=0; i<na; i++) 
      total+=std::conj(_el[i])*(_params[ka].lmeas[i] - avea);
    F(ka)= -2.0*total.real()/(double)na;
  }
  
  for(size_t ja=0; ja<Nj; ja++) {
    size_t na = _jsparams[ja].lmeas.nvals();
    for(size_t jb=ja; jb<Nj; jb++) {
      std::complex<double> avea=_jsparams[ja].lmeas.ave();
      std::complex<double> aveb=_jsparams[jb].lmeas.ave();
      std::complex<double> total=0.0;
      for(size_t i=0; i<na; i++) {
        std::complex<double> ai=_jsparams[ja].lmeas[i];
        std::complex<double> bi=_jsparams[jb].lmeas[i];
        total+=(ai-avea) * (bi-aveb);
      }
      S(Np+ja,Np+jb)=(total/(double)na).real();
      S(Np+jb,Np+ja)=S(Np+ja,Np+jb);
    }
    std::complex<double> total=0.0;
    std::complex<double> avea=_jsparams[ja].lmeas.ave();
    for(size_t i=0; i<na; i++) 
      total+=std::conj(_el[i])*(_jsparams[ja].lmeas[i] - avea);
    F(Np+ja)= -2.0*total.real()/(double)na;
  }   

  // preconditioning -----------------------------------------------------------
  //Eigen::MatrixXd S_pc(N,N);
  //Eigen::VectorXd F_pc(N); 
  //for(size_t ka=0; ka<N; ka++) {
  //  F_pc(ka)=F(ka)/std::sqrt(S(ka,ka));
  //  if(fabs(F(ka)) < 1e-6) F_pc(ka) = 0;
  //  for(size_t kb=0; kb<N; kb++) {
  //    S_pc(ka,kb)=S(ka,kb)/(std::sqrt(S(ka,ka)*S(kb,kb)));
  //    //if(ka==kb) S_pc(ka,kb)=1+1e-3;
  //  }
  //}
  // ---------------------------------------------------------------------------
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solveS;
  solveS.compute(S);
  Eigen::VectorXd e=solveS.eigenvalues();
  Eigen::MatrixXd U=solveS.eigenvectors();
  std::vector<size_t> keep_indices;
  for(size_t i=0; i<N; i++) {
    if(e(i)>1e-3) keep_indices.push_back(i);
  }
  size_t Nnew=keep_indices.size();
  if(Nnew==0) return true; 
  Eigen::MatrixXd Snew(Nnew,Nnew);
  Eigen::VectorXd Fnew(Nnew);
  for(size_t i=0; i<Nnew; i++) {
    size_t i_index=keep_indices[i];
    for(size_t j=0; j<Nnew; j++) {
      size_t j_index=keep_indices[j];
      Snew(i,j)=S(i_index, j_index);
    }
    Fnew(i)=F(keep_indices[i]);
  }
  //dA = df*S_pc.inverse()*F_pc;
  dA = df*Snew.inverse()*Fnew;
  for(size_t k=0; k<Nnew; k++) {
    //if(fabs(S(k,k))<1e-3) continue;
    //_params[k].val+=dA[k]/std::sqrt(S(k,k));
    _params[keep_indices[k]].val+=dA[k];
  }  
  //for(size_t k=0; k<Nj; k++) {
  //  if(fabs(F_pc[k+Np])<1e-5) continue;
  //  _jsparams[k].val+=dA[k+Np]/std::sqrt(S(k+Np,k+Np));
  //}
  return false;
}

size_t HeisenbergChainSimulator::_flipspin(const size_t& rpos) {
  size_t iexpos; // initial position in extended basis
  size_t nexpos; // new position in extended basis
  size_t lindex; // index of position in W array
  size_t nferms=0; // number of fermion swaps
  int ds; // change in spin
  // determine new location of spin in extended state --------------------------
  if(_spinstate[rpos]==1){
    iexpos=rpos;
    nexpos=iexpos+_size;
    ds=-2;
    for(size_t i=iexpos+1; i<nexpos; i++)
      if(_operslist[i]!=0) nferms += 1;
  }
  else{
    iexpos=rpos+_size;
    nexpos=iexpos-_size;
    ds=2;
    for(size_t i=nexpos+1; i<iexpos; i++)
      if(_operslist[i]!=0) nferms += 1;
  }
  // ---------------------------------------------------------------------------
  double dj=0.0;
  //double newjsf;
  //for(auto it=_jsparams.begin(); it!=_jsparams.end(); it++) {
  //  size_t dr=it->space;
  //  size_t rl=(rpos+(_size-dr))%_size;
  //  size_t rr=(rpos+(_size+dr))%_size;
  //  double v=it->val;
  //  dj+=0.25*(double)(_spinstate[rl]+_spinstate[rr])*v*(double)ds;
  //}
  //newjsf=_jsf*std::exp(dj);
  // ---------------------------------------------------------------------------
  
  // accept flip according to Green's function ---------------------------------
  lindex=_operslist[iexpos]-1;
  //double amp=(std::exp(dj))*fabs(_gmat(nexpos, lindex));
  double amp=fabs(_gmat(nexpos, lindex));
  double rnum=_rnum(_mteng);
   
  if(rnum<amp*amp){
    // swap operator positions and update spin state
    _operslist[nexpos]=_operslist[iexpos]; 
    _operslist[iexpos]=0;
    _spinstate[rpos]+=ds;

    // update gmat -------------------------------------------------------------
    Eigen::MatrixXcd newgmat(2*_size, _size);
    for(size_t i=0; i<2*_size; i++)
    for(size_t j=0; j<_size; j++) {
      double d=0;
      if(j==lindex) d=1.0;
      std::complex<double> den = _gmat(nexpos,lindex);
      std::complex<double> num = _gmat(i,lindex);
      newgmat(i,j)=_gmat(i,j)-(num/den)*(_gmat(nexpos,j)-d); 
    }
    _gmat=newgmat;
    // -------------------------------------------------------------------------    
    // update Jastrow factor
    //_jsf=newjsf;
    _ffac = std::pow(-1.0, nferms)*_ffac;
    return 1;
  }
  return 0;
}

void HeisenbergChainSimulator::_sweep() {
  // flip until N electron moves have occured.
  // while(counter<_size) counter+=_flipspin();
  for(size_t i=0; i<2*_size; i++) {
    int rpos=(*_rpos)(_mteng);
    size_t n=0;
    n = _flipspin(rpos);
  } 
  _reinitgmat();
}

std::complex<double> HeisenbergChainSimulator::_isingenergy() {
  std::complex<double> total=0.0;
  for(size_t i=0; i<_size; i++) {
    if(i<(_size-1)) total+=0.25*(double)_spinstate[i]*(double)_spinstate[i+1];
    else total+=0.25*(double)_spinstate[0]*(double)_spinstate[i];
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
      //size_t iexpos1=i;   
      //size_t iexpos2=j+_size;   
      //size_t nexpos1=j;   
      //size_t nexpos2=i+_size;   
      size_t lindex1=_operslist[iexpos1]-1;
      size_t lindex2=_operslist[iexpos2]-1;
      std::complex<double> det=_gmat(nexpos1, lindex1)*_gmat(nexpos2, lindex2)
        -_gmat(nexpos2, lindex1)*_gmat(nexpos1, lindex2);
      double jsum=0.0;
      for(auto it=_jsparams.begin(); it!=_jsparams.end(); it++) {
        size_t dr=it->space;
        double v=it->val;
        size_t ri=(i+(_size-dr))%_size;
        size_t li=(i+(_size+dr))%_size;
        size_t rj=(j+(_size-dr))%_size;
        size_t lj=(j+(_size+dr))%_size;
        if(li!=j) jsum+=0.25*v*(_spinstate[ri]+_spinstate[li])*(-2.0);
        else jsum+=0.25*v*_spinstate[ri]*(-2.0);
        if(rj!=i) jsum+=0.25*v*(_spinstate[rj]+_spinstate[lj])*(-2.0);
        else jsum+=0.25*v*_spinstate[lj]*(-2.0);
      }
      //std::cout << "up-down to down-up" << std::endl;
      //std::cout << det << std::endl;
      //det = std::exp(jsum)*det;
      total -= 0.5 * det; 
    }
    else {
      size_t iexpos1=i+_size;   
      size_t iexpos2=j;   
      size_t nexpos1=i;   
      size_t nexpos2=j+_size;   
      //size_t iexpos1=i+_size;   
      //size_t iexpos2=j;   
      //size_t nexpos1=j+_size;   
      //size_t nexpos2=i;   
      size_t lindex1=_operslist[iexpos1]-1;
      size_t lindex2=_operslist[iexpos2]-1;
      std::complex<double> det=_gmat(nexpos1, lindex1)*_gmat(nexpos2, lindex2)
        -_gmat(nexpos2, lindex1)*_gmat(nexpos1, lindex2);
      double jsum=0.0;
      for(auto it=_jsparams.begin(); it!=_jsparams.end(); it++) {
        size_t dr=it->space;
        double v=it->val;
        size_t ri=(i+(_size-dr))%_size;
        size_t li=(i+(_size+dr))%_size;
        size_t rj=(j+(_size-dr))%_size;
        size_t lj=(j+(_size+dr))%_size;
        if(li!=j) jsum+=0.25*v*(_spinstate[ri]+_spinstate[li])*(-2.0);
        else jsum+=0.25*v*_spinstate[ri]*(-2.0);
        if(rj!=i) jsum+=0.25*v*(_spinstate[rj]+_spinstate[lj])*(-2.0);
        else jsum+=0.25*v*_spinstate[lj]*(-2.0);
      }
      //det = std::exp(jsum)*det;
      //std::cout << "down-up to up-down" << std::endl;
      //std::cout << det << std::endl;
      total -= 0.5 * det;
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
    var_params << std::setw(COLUMN_WIDTH) << std::left << it->name;
  for(auto it=_jsparams.begin(); it!=_jsparams.end(); it++) 
    var_params << std::setw(COLUMN_WIDTH) << std::left << it->name;
  var_params << '\n';
  for(auto it=_params.begin(); it!=_params.end(); it++) 
    var_params << std::setw(COLUMN_WIDTH) << std::left << it->val;
  for(auto it=_jsparams.begin(); it!=_jsparams.end(); it++) 
    var_params << std::setw(COLUMN_WIDTH) << std::left << it->val;
  var_params << '\n';

  // Loop for variational steps ------------------------------------------------
  for(size_t vstep=0; vstep<vsteps; vstep++) {
    // open file for local observables
    std::ofstream measfile;
    std::string measfileName=ofname+std::to_string(vstep)+".dat";
    measfile.open(measfileName);  
    measfile << vstep << '\n';
    measfile << std::setw(COLUMN_WIDTH) << std::left << "e_L";
    for(auto it=_params.begin(); it!=_params.end(); it++) {
      measfile << std::setw(COLUMN_WIDTH) << std::left << it->name; 
    }
    for(auto it=_jsparams.begin(); it!=_jsparams.end(); it++) {
      measfile << std::setw(COLUMN_WIDTH) << std::left << it->name; 
    }
    measfile << '\n'; 
    _ffac = 1.0;
    // equilibrate -------------------------------------------------------------
    for(size_t i=0; i<equil; i++) _sweep();
    std::cout << "equilibration " << vstep << " done." << std::endl;   
    // clear any old measurements ----------------------------------------------
    for(auto it=_params.begin(); it!=_params.end(); it++) it->lmeas.clear();  
    for(auto it=_jsparams.begin(); it!=_jsparams.end(); it++) it->lmeas.clear();
    _el.clear();
    // -------------------------------------------------------------------------
    
    // Sample wavefunction -----------------------------------------------------
    for(size_t i=0; i<simul; i++) { 
      // start with a sweep
      _sweep();
      //std::complex<double> e=_isingenergy();
      std::complex<double> e=_heisenergy();
      _el.push(e); // record energy

      // Measure variational parameters ----------------------------------------
      for(auto it=_params.begin(); it!=_params.end(); it++) {
        Eigen::MatrixXcd redMat(_size, 2*_size); // N_e X 2L matrix
        Eigen::MatrixXcd okmat; // operator to measure
        for(size_t l=0; l<_size; l++){   
          auto itstart=_operslist.begin();
          auto itend=_operslist.end();
          size_t lpos=std::distance(itstart, std::find(itstart, itend, l+1));
          for(size_t r=0; r<2*_size; r++) redMat(l,r) = it->mmat(lpos, r);  
        }
        // compute actual local value
        okmat=redMat*_gmat;
        std::complex<double> okval = okmat.trace();
        it->lmeas.push(okval);
      }
      for(auto it=_jsparams.begin(); it!=_jsparams.end(); it++) {
        size_t dx=it->space;
        double total=0.0;
        for(size_t r=0; r<_size; r++) 
          total+=0.25*_spinstate[r]*_spinstate[(r+dx)%_size];    
        it->lmeas.push(total/(double)_size);
      }
      // -----------------------------------------------------------------------
    }
    // -------------------------------------------------------------------------
    for(size_t i=0; i<_el.nbins(); i++) {
      measfile << std::setw(COLUMN_WIDTH) << std::left << _el(i);
      for(auto pit=_params.begin(); pit!=_params.end(); pit++) 
        measfile << std::setw(COLUMN_WIDTH) << std::left << pit->lmeas(i);
      for(auto jit=_jsparams.begin(); jit!=_jsparams.end(); jit++) {
        measfile << std::setw(COLUMN_WIDTH) << std::left << jit->lmeas(i);
      }
      measfile << '\n';
    }
    measfile.close();
    measfile.clear();
    // update variational parameters
    _updateparams(df);
    std::cout << "outputing new variational parameters" << std::endl;
    for(auto it=_params.begin(); it!=_params.end(); it++) 
      var_params << std::setw(COLUMN_WIDTH) << std::left << it->val; 
    for(auto it=_jsparams.begin(); it!=_jsparams.end(); it++) 
      var_params << std::setw(COLUMN_WIDTH) << std::left << it->val; 
    var_params << '\n';
    std::cout << "reinitialzing auxiliary Hamiltonian" << std::endl;
    _auxham.init(_params);
    _reinitgmat();
  }
  // ---------------------------------------------------------------------------
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

void HeisenbergChainSimulator::print_params() {
  for(auto it=_params.begin(); it!=_params.end(); it++)
    std::cout << std::setw(OUTPUT_WIDTH) << std::left 
              << it->name << std::setw(3) << " : " 
              << std::setw(OUTPUT_WIDTH) << std::left << it->val
              << std:: endl;
  for(auto it=_jsparams.begin(); it!=_jsparams.end(); it++)
    std::cout << std::setw(OUTPUT_WIDTH) << std::left 
              << it->name << std::setw(3) << " : " 
              << std::setw(OUTPUT_WIDTH) << std::left << it->val
              << std::endl;
}

}
}
