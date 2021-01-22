#ifndef __CHEBY_HPP__
#define __CHEBY_HPP__

#include <Eigen/Dense>

namespace Cheby{


template <int _n_>
class Cheby{
private:
  Eigen::Matrix<double,_n_,1> c_,abci_;
  double a_,b_;
public:
  Cheby(){}
  int getN(){return _n_;}
  void setRange(const double a, const double b){
    a_ = a; 
    b_ = b;
    if(a_>=b_){
      std::cout << "a>=b in Cheby::setRange" << std::endl;
    }
  }
  Eigen::Matrix<double,_n_,1> getAbcissa(){
    double bma = 0.5*(b_-a_);
    double bpa = 0.5*(b_+a_);
    double y;
    double dn=static_cast<double>(_n_);
    for(int i=0;i<_n_;i++){
      y = cos(3.141592653589793116*(static_cast<double>(i)+0.5)/dn);
      abci_(i)=y*bma+bpa;
    }
    return abci_;
  }
  Eigen::Matrix<double,_n_,1> getAbcissa(const double a, const double b){
    setRange(a,b);
    return getAbcissa();
  }
  void setFunVal(const Eigen::Matrix<double,Eigen::Dynamic,1> &fn){
    double dn=static_cast<double>(_n_);
    double fac = 2.0/dn;
    double su;
    double dj,dk;
    for(int j=0;j<_n_;j++){
      su = 0.0;
      dj = static_cast<double>(j);
      for(int k=0;k<_n_;k++){
        dk = static_cast<double>(k);
        su +=  fn[k]*cos(3.141592653589793116*dj*(dk+0.5)/dn);
      }
      c_[j] = su*fac;
    }
  }
  double eval(const double x) const {
    double y = (2.0*x-a_-b_)/(b_-a_);
    double y2 = 2.0*y;
    double d = 0.0;
    double dd = 0.0;
    double sv;
    for(int j=_n_-1;j>0;j--){
      sv = d;
      d = y2*d - dd + c_[j];
      dd = sv;
    }
    return(y*d - dd + 0.5*c_[0]);
  }
  
   void deriv(Cheby<_n_> &ret) const {
    
    ret.a_ = a_;
    ret.b_ = b_;
    ret.abci_ = abci_;
    double dn = static_cast<double>(_n_);
    ret.c_[_n_-1] = 0.0;
    ret.c_[_n_-2] = 2.0*(dn-1.0)*c_[_n_-1];
    for(int j=_n_-2;j>0;j--){
      ret.c_[j-1]=ret.c_[j+1]+2.0*c_[j]*static_cast<double>(j);
    }
    ret.c_ *= 2.0/(b_-a_);
  }
  Cheby<_n_> deriv() const {
    Cheby<_n_> ret;
    deriv(ret);
    return(ret);
  }
  
  
  void integral(Cheby<_n_> &ret) const {
    ret.a_ = a_;
    ret.b_ = b_;
    ret.abci_ = abci_;
    double su = 0.0;
    double fac = 1.0;
    double con = 0.25*(b_-a_);
    for(int j=1; j<_n_-1;j++){
      ret.c_[j] = con*(c_[j-1]-c_[j+1])/static_cast<double>(j);
      su = su + fac*ret.c_[j];
      fac = -fac;
    }
    ret.c_[_n_-1] = con*ret.c_[_n_-2]/static_cast<double>(_n_-1);
    su = su + fac*ret.c_[_n_-1];
    ret.c_[0] = 2.0*su;
  }
  Cheby<_n_> integral() const{
    Cheby<_n_> ret;
    integral(ret);
    return(ret);
  }
  // solves f(x)=y for x\in [a,b], when f is the function represented by this 
  // derivative is a Cheby-object of f'(x)
  double solve(const double y, const Cheby<_n_> &derivative) const {
#ifdef _PDP_DEBUG_
    std::cout << "Cheby::solve : a = " << a_ << " b = " << b_ << " y = " << y << std::endl;
#endif
    double ffl = eval(a_)-y;
    double ffr = eval(b_)-y;
    if(ffl*ffr>0.0){
      std::cout << "no unique root in Cheby::solve" << std::endl;
      return(0.5*(a_+b_));
    }
    double lb=a_;
    double rb=b_;
    double fft,tb;
    // bisection search
    while(rb-lb>1.0e-7*(b_-a_)){
      tb = 0.5*(rb+lb);
      fft = eval(tb)-y;
      if(fft*ffl>=0.0){
        lb = tb;
        ffl = fft;
      } else {
        rb = tb;
        ffr = fft;
      }
    }
    
    // refine solution Newton iterations
    tb = 0.5*(lb+rb);
    double res;
    for(int i=0;i<4;i++){
      res = eval(tb) - y;
      tb = tb - res/derivative.eval(tb);
      tb = fmin(rb,fmax(lb,tb));
    }
    
    if(fabs(eval(tb)-y)>1.0e-6*fmax(1.0,fabs(y))){
      std::cout << "Cheby::solve : residual at returned solution = " << res << std::endl; 
    }
#ifdef _PDP_DEBUG_
    std::cout << "Cheby::solve : last residual = " << eval(tb)-y << std::endl;
#endif
    return(tb);
  }
  
  double solve(const double y) const {
    return solve(y,deriv());
  }
  
  
};

} // namespace Cheby


#endif


