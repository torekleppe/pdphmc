#ifndef __TRANSFORMED_PRIOR_HPP__
#define __TRANSFORMED_PRIOR_HPP__
namespace transformed_prior{



template <class var>
class uniform_logit{
  double a_;
  double b_;
  var lpar_;
public:
  uniform_logit(const double a, const double b,const var lpar): a_(a),b_(b),lpar_(lpar) {
    if(a>=b){
      _ERROR_LOG_.push("bad parameters in uniform_logit");
    }
  }
  var par() const {
    if(doubleValue(lpar_)>0.0){
      return(a_+(b_-a_)/(1.0+exp(-lpar_)));
    } else {
      var tmp = exp(lpar_);
      return(a_+(b_-a_)*tmp/(1.0+tmp));
    }
  }
  var implied_prior_lpdf() const {
    var tmp = (doubleValue(lpar_)>0.0) ? -lpar_ : lpar_;
    return(tmp - 2.0*log(1.0+exp(tmp)));
  }
  
};


template <class var>
class uniform_probit{
  double a_;
  double b_;
  var lpar_;
public:
  uniform_probit(const double a, const double b,const var lpar): a_(a),b_(b),lpar_(lpar) {
    if(a>=b){
      _ERROR_LOG_.push("bad parameters in uniform_logit");
    }
  }
  var par() const { return a_ + (b_-a_)*fast_spec_funs::fast_std_pnorm(lpar_);}
  var implied_prior_lpdf() const { return -0.5*lpar_*lpar_ -0.9189385332046729; }
  
};
}

#endif






