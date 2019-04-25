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

#define COLUMN_WITH 20

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
  _gmat = redmat*projmat.inverse();
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
    total += 0.25*v*_spinstate[i]*_spinstate[(i+dx)%_size];
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
      S(ka, Np+j)=(total/(double)na).real();
      S(Np+j, ka)=S(ka, Np+j);
    }
    // Loop for BCS forces -----------------------------------------------------
    std::complex<double> total=0.0;
    std::complex<double> avea=_params[ka].lmeas.ave();
    for(size_t i=0; i<na; i++) 
      total+=std::conj(_el[i])*(_params[ka].lmeas[i] - avea);
    F(ka)=-2.0*total.real()/(double)na;
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
    F(Np+ja)=-2.0*total.real()/(double)na;
  }   

  // preconditioning -----------------------------------------------------------
  Eigen::MatrixXd S_pc(N,N);
  Eigen::VectorXd F_pc(N); 
  for(size_t ka=0; ka<N; ka++) {
    F_pc(ka)=F(ka)/std::sqrt(S(ka,ka));
    if(fabs(F(ka)) < 1e-10) F_pc(ka) = 0;
    for(size_t kb=0; kb<N; kb++) {
      S_pc(ka,kb)=S(ka,kb)/(std::sqrt(S(ka,ka)*S(kb,kb)));
      if(fabs(S(ka,kb)) < 1e-14) S_pc(ka, kb) = 0;
      if(ka==kb) S_pc(ka, kb) += 1e-3;
    }
  }
  // ---------------------------------------------------------------------------
  //Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solveS;
  //solveS.compute(S_pc);
  //std::cout << solveS.eigenvalues() << std::endl;
  //std::cout << S << std::endl;
  //std::cout << "----" << std::endl;
  //std::cout << S_pc << std::endl;
  //std::cout << "----" << std::endl;
  //std::cout << S_pc.inverse() << std::endl;
  //std::cout << "----" << std::endl;
  //std::cout << F << std::endl;
  //std::cout << "----" << std::endl;
  //std::cout << F_pc << std::endl;

  if(fabs(S_pc.determinant())<1e-12) {
    std::cout << "CONVERGED: S IS SINGULAR" << std::endl;
    return true;
  }
  else {
    dA = df*S_pc.inverse()*F_pc;
    //std::cout << F << std::endl;
    //dA = df*F;
    for(size_t k=0; k<Np; k++) {
      if(fabs(F_pc[k])<1e-12) continue;
      _params[k].val+=dA[k]/std::sqrt(S(k,k));
    }  
    for(size_t k=0; k<Nj; k++) {
      if(fabs(F_pc[k])<1e-12) continue;
      _jsparams[k].val+=dA[k]/std::sqrt(S(k+Np,k+Np));
    }
    return false;
  }
  return false;
  //dA = df*F;
  //dA = S.inverse()*F*df; 
  //for(size_t k=0; k<N; k++) _params[k].val+=dA[k];
}

size_t HeisenbergChainSimulator::_exchange(const size_t& i, const size_t& j) { 
  int rpos1=i;
  int rpos2=j;
  int ds1;
  int ds2;
  size_t iexpos1, iexpos2, nexpos1, nexpos2, lindex1, lindex2;
  // generate a random second position different from the first
  do rpos2=(*_rpos)(_mteng); while(rpos2==rpos1);
  if(_spinstate[rpos1]==1) {
    iexpos1=rpos1;
    nexpos1=rpos1+_size;
    ds1=-2;
  }
  else {
    iexpos1=rpos1+_size;
    nexpos1=rpos1;
    ds1=2;
  }
  if(_spinstate[rpos2]==1) {
    iexpos2=rpos2;
    nexpos2=rpos2+_size;
    ds2=-2;
  }
  else {
    iexpos2=rpos2+_size;
    nexpos2=rpos2;
    ds2=2;
  }
  lindex1=_operslist[iexpos1]-1;
  lindex2=_operslist[iexpos2]-1;
  double amp=fabs(_gmat(nexpos1, lindex1)*_gmat(nexpos2, lindex2)
        -_gmat(nexpos2, lindex1)*_gmat(nexpos1, lindex2));
  double rnum=_rnum(_mteng);
  if(rnum<amp*amp) {
    // update stored state -----------------------------------------------------
    _operslist[nexpos1]=_operslist[iexpos1];
    _operslist[nexpos2]=_operslist[iexpos2];
    _operslist[iexpos1]=0;
    _operslist[iexpos2]=0;
    _spinstate[rpos1]+=ds1;
    _spinstate[rpos2]+=ds2;
    // update Green's function matrix ------------------------------------------
    Eigen::MatrixXcd newgmat(2*_size, _size);
    std::complex<double> d11=_gmat(nexpos1, lindex1);
    std::complex<double> d12=_gmat(nexpos1, lindex2);
    std::complex<double> d21=_gmat(nexpos2, lindex1);
    std::complex<double> d22=_gmat(nexpos2, lindex2);
    //for(size_t i=0; i<2*_size; i++) 
    //for(size_t j=0; j<_size; j++) {
    //  std::complex<double> n1=_gmat(i, lindex1);
    //  std::complex<double> n2=_gmat(i, lindex2);
    //  std::complex<double> p1=_gmat(nexpos1, j);
    //  std::complex<double> p2=_gmat(nexpos2, j);
    //  std::complex<double> id1=0, id2=0;
    //  if(lindex1==j) id1=1.0; 
    //  if(lindex2==j) id2=1.0;
    //  if(fabs(n1)<1e-12) {
    //    d11=1;
    //    d12=1;
    //  }
    //  if(fabs(n2)<1e-12) {
    //    d21=1;
    //    d22=1;
    //  }
    //  newgmat(i,j) = _gmat(i,j) -(n1/d11)*(p1-id1)-(n1/d12)*(p2-id2)
    //                            -(n2/d21)*(p1-id1)-(n2/d22)*(p2-id2); 
    //}
    //_gmat=newgmat;
    _reinitgmat();
    return 1;
  }
  return 0;
}

size_t HeisenbergChainSimulator::_flipspin(const size_t& rpos) {
  size_t iexpos; // initial position in extended basis
  size_t nexpos; // new position in extended basis
  size_t lindex; // index of position in W array
  int ds; // change in spin
  // determine new location of spin in extended state --------------------------
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
  // ---------------------------------------------------------------------------
  double dj=0.0;
  double newjsf;
  for(auto it=_jsparams.begin(); it!=_jsparams.end(); it++) {
    size_t dr=it->space;
    size_t rl=(rpos+(_size-dr))%_size;
    size_t rr=(rpos+(_size+dr))%_size;
    double v=it->val;
    dj+=0.25*(_spinstate[rl]+_spinstate[rr])*v*ds;
  }
  newjsf=_jsf*std::exp(dj);
  // ---------------------------------------------------------------------------
  
  // accept flip according to Green's function ---------------------------------
  lindex=_operslist[iexpos]-1;
  double amp=(newjsf/_jsf)*fabs(_gmat(nexpos, lindex));
  double rnum = _rnum(_mteng);
   
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
      if(j==lindex) d=1;
      std::complex<double> den = _gmat(nexpos,lindex);
      std::complex<double> num = _gmat(i,lindex);
      newgmat(i,j)=_gmat(i,j)-(num/den)*(_gmat(nexpos,j)-d); 
    }
    _gmat=newgmat;
    // -------------------------------------------------------------------------    

    // update Jastrow factor
    _jsf=newjsf;
    _ffac = -1.0*_ffac;
    return 1;
  }
  return 0;
}

void HeisenbergChainSimulator::_sweep() {
  // flip until N electron moves have occured.
  // while(counter<_size) counter+=_flipspin();
  for(size_t i=0; i<_size; i++) {
    int rpos=(*_rpos)(_mteng);
    size_t n=0;
    n = _flipspin(rpos);
    //int rpos1=(*_rpos)(_mteng);
    //int rpos2;
    //do rpos2=(*_rpos)(_mteng); while(rpos2==rpos1);
    //n = _exchange(rpos1, rpos2);
    //if(rnum<0.5) n=_flipspin(rpos);
    //else n=_exchange();
  } 
  _reinitgmat();
}

std::complex<double> HeisenbergChainSimulator::_isingenergy() {
  std::complex<double> total=0.0;
  for(size_t i=0; i<_size; i++) {
    if(i<(_size-1)) total+=0.25*_spinstate[i]*_spinstate[i+1];
    else total+=0.25*_spinstate[0]*_spinstate[i];
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
      // We'll just use the bad solution for the time being
      std::vector<int> newstate;
      double newjsf;
      double newjsum=0.0;
      for(size_t k=0; k<_spinstate.size(); k++) 
        newstate.push_back(_spinstate[k]);
      newstate[i]=newstate[i]-2;
      newstate[j]=newstate[j]+2;
      for(size_t k=0; k<_size; k++) 
      for(size_t l=0; l<_jsparams.size(); l++) { 
        double v=_jsparams[l].val;
        size_t dx=_jsparams[l].space;
        newjsum += 0.25*v*newstate[k]*newstate[(k+dx)%_size];
      }
      newjsf=std::exp(newjsum);
      det = (newjsf/_jsf)*det;
      total += 0.5 * det*_ffac*(-1.0);
      std::cout << "det: " << det << std::endl;
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
      std::vector<int> newstate;
      double newjsf;
      double newjsum=0.0;
      for(size_t k=0; k<_spinstate.size(); k++) 
        newstate.push_back(_spinstate[k]);
      newstate[i]=newstate[i]+2;
      newstate[j]=newstate[j]-2;
      for(size_t k=0; k<_size; k++) 
      for(size_t l=0; l<_jsparams.size(); l++) { 
        double v=_jsparams[l].val;
        size_t dx=_jsparams[l].space;
        newjsum += 0.25*v*newstate[k]*newstate[(k+dx)%_size];
      }
      newjsf=std::exp(newjsum);
      std::cout << "det: " << det << std::endl;
      det = (newjsf/_jsf)*det;
      total += 0.5 * det*_ffac*(-1.0);
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
  for(auto it=_jsparams.begin(); it!=_jsparams.end(); it++) 
    var_params << std::setw(COLUMN_WITH) << std::left << it->name;
  var_params << '\n';
  for(auto it=_params.begin(); it!=_params.end(); it++) 
    var_params << std::setw(COLUMN_WITH) << std::left << it->val;
  for(auto it=_jsparams.begin(); it!=_jsparams.end(); it++) 
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
    for(auto it=_jsparams.begin(); it!=_jsparams.end(); it++) {
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
    for(auto it=_jsparams.begin(); it!=_jsparams.end(); it++) {
      it->lmeas.clear();
    }
    _el.clear();
    // -------------------
    // Loop for sampling wavefunction
    // -------------------
    for(size_t i=0; i<simul; i++) { 
      // start with a sweep
      _sweep();
      //std::complex<double> e=_isingenergy();
      std::complex<double> e=_heisenergy();
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
      for(auto it=_jsparams.begin(); it!=_jsparams.end(); it++) {
        size_t dx=it->space;
        double total=0.0;
        for(size_t r=0; r<_size; r++) 
          total+=_spinstate[r]*_spinstate[(r+dx)%_size];    
        it->lmeas.push(total/(double)_size);
      }
    }
    for(size_t i=0; i<_el.nbins(); i++) {
      measfile << std::setw(COLUMN_WITH) << std::left << _el(i);
      for(auto pit=_params.begin(); pit!=_params.end(); pit++) 
        measfile << std::setw(COLUMN_WITH) << std::left << pit->lmeas(i);
      for(auto jit=_jsparams.begin(); jit!=_jsparams.end(); jit++) {
        measfile << std::setw(COLUMN_WITH) << std::left << jit->lmeas(i);
      }
      measfile << '\n';
    }
    measfile.close();
    measfile.clear();
    // update variational parameters
    _updateparams(df);
    std::cout << "outputing new variational parameters" << std::endl;
    for(auto it=_params.begin(); it!=_params.end(); it++) 
      var_params << std::setw(COLUMN_WITH) << std::left << it->val; 
    for(auto it=_jsparams.begin(); it!=_jsparams.end(); it++) 
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
