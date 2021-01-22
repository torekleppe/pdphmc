#ifndef __ADTARGET_HPP__
#define __ADTARGET_HPP__

/*
 * This class provides AD gradient of the original target distribution
 * 
 * 
 */

template <class _target_type_>
class ADtarget{
private:
  _target_type_ t_;
  int d_;
  int dg_;
  
  VectorXad dpar_;
  
  double ret_;
  stan::math::var dret_; 
  
public:
  inline int copies() const {return 1;}
  // run the setup method of the target distribution and allocates some memory
  inline void setup(){
    t_.setup();
    d_ = t_.dim();
    dg_ = t_.dimGenerated();
    dpar_.resize(d_);
  }
  
  inline int dim() const {return d_;}
  inline int dimGenerated() const {return dg_;}
  
  // run the eval function with double variables, i.e. no gradient
  inline double feval(const Eigen::VectorXd &par, Eigen::VectorXd &gen){
    return(t_.eval(par,gen));
  }
  
  inline double feval(const Eigen::VectorXd &par, Eigen::Ref<Eigen::VectorXd> gen){
    return(t_.eval(par,gen));
  }
  
  // gradient eval
  double geval(const Eigen::VectorXd &par, 
               Eigen::Ref< Eigen::VectorXd > gen, 
               Eigen::Ref< Eigen::VectorXd > grad){
      dpar_ = par;
      
      try{
        dret_ = t_.eval(dpar_,gen);
      } 
      catch(...){
        stan::math::recover_memory();
        std::cout << "Bad function eval" << std::endl;
        return std::numeric_limits<double>::quiet_NaN();
      }

      ret_ = dret_.val();
      dret_.grad();
      for(int i=0;i<d_;i++) grad.coeffRef(i)=dpar_.coeff(i).adj();
      stan::math::recover_memory();
      return(ret_);
    
  }
};

/*
 * 
 * Two independent copies of the original target
 */
template <class _target_type_>
class ADtarget_cp2{
  ADtarget<_target_type_> adt_;
  int d_;
  int dg_;
  int ed_;
  int edg_;
public:
  inline int copies() const {return 2;}
  inline void setup(){
    adt_.setup();
    d_ = adt_.dim();
    ed_ = 2*d_;
    dg_ = adt_.dimGenerated();
    edg_ = 2*dg_;
  }
  inline int dim() const {return ed_;}
  inline int dimGenerated() const {return edg_;}
  
  inline double feval(const Eigen::VectorXd &par, Eigen::VectorXd &gen){
    double ret = adt_.feval(par.head(d_),gen.head(dg_));
    ret += adt_.feval(par.tail(d_),gen.tail(dg_));
    return(ret);
  }
  
  inline double feval(const Eigen::VectorXd &par, Eigen::Ref<Eigen::VectorXd> gen){
    double ret = adt_.feval(par.head(d_),gen.head(dg_));
    ret += adt_.feval(par.tail(d_),gen.tail(dg_));
    return(ret);
  }
  
  double geval(const Eigen::VectorXd &par, 
               Eigen::Ref< Eigen::VectorXd > gen, 
               Eigen::Ref< Eigen::VectorXd > grad){
    double ret = adt_.geval(par.head(d_),gen.head(dg_),grad.head(d_));
    ret += adt_.geval(par.tail(d_),gen.tail(dg_),grad.tail(d_));
    return(ret);
  }

};



#endif



