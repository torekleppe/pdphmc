#ifndef __CIP_Wishart_HPP__
#define __CIP_Wishart_HPP__
#include <Eigen/Dense>


namespace PDPHMC_CIP_Wishart{

int CIP_Wishart_Vdim(const int n){
  return(std::round(0.5*static_cast<double>(n*(n-1))));
}


template <class df_type__, class scale_type__,class lambda_type__, class V_type__> 
class CIP_Wishart{
  bool dr_;
  int n_;
  df_type__ df_;
  lambda_type__ logDetJac_;
  
  
  Eigen::Matrix<scale_type__,Eigen::Dynamic,Eigen::Dynamic> scale_,scaleChol_;
  Eigen::Matrix<lambda_type__,Eigen::Dynamic,1> lambda_;
  Eigen::Matrix<typename boost::math::tools::promote_args<df_type__,
                                                          scale_type__,
                                                          lambda_type__, 
                                                          V_type__>::type,
                                                          Eigen::Dynamic,
                                                          Eigen::Dynamic> LP_,F_;
public:
  CIP_Wishart(const double df,
              const Eigen::Matrix<double,Eigen::Dynamic,1> scale,
              const Eigen::Matrix<lambda_type__,Eigen::Dynamic,1> lambda,
              const Eigen::Matrix<V_type__,Eigen::Dynamic,1> V,
              const bool dynamicRescale=true) : dr_(dynamicRescale),df_(df),lambda_(lambda){
    
    n_ = scale.size();
    scale_ = scale.asDiagonal();
    scaleChol_ = scale.cwiseSqrt().asDiagonal();
    LP_.resize(n_,n_);
    F_.resize(n_,n_);
    LP_.setIdentity(n_,n_);
    if(lambda.size() != n_){
      _ERROR_LOG_.push("size mismatch in argument lambda in CIP_Wishart_diagScale");
    } else if(V.size() != CIP_Wishart_Vdim(n_)){
      _ERROR_LOG_.push("size mismatch in argument V in CIP_Wishart_diagScale");
    } else {
      // build the 
      int l = 1;
      int f = 0;
      logDetJac_ = 0.0;
      LP_.coeffRef(n_-1,n_-1) = exp(0.5*lambda.coeff(n_-1));
      for(int j=n_-2;j>=0;j--){
        LP_.coeffRef(j,j) = exp(0.5*lambda.coeff(j));
        if(dr_){
          LP_.col(j).tail(l) = V.segment(f,l); 
          logDetJac_ -= 0.5*static_cast<double>(l)*lambda.coeff(j);
        } else {
          LP_.col(j).tail(l) = LP_.coeff(j,j)*V.segment(f,l);  
        }
        f+=l;
        l++;
      }
    }
    F_ = scaleChol_.diagonal().asDiagonal().inverse()*LP_;
  }
  
  
  // gives logdensity so that P will have a wishart(df,scale)-distrbution 
  // Note, normalization factor depending on df and scale missing
  typename boost::math::tools::promote_args<df_type__,
                                            scale_type__,
                                            lambda_type__, 
                                            V_type__>::type
  lpdf() const {
    typename boost::math::tools::promote_args<df_type__,scale_type__,lambda_type__, V_type__>::type ret = 0.0;
    df_type__ tmp = 0.5*(df_+static_cast<double>(n_)+1.0);
    // prior on lambda
    for(int i=0;i<n_;i++){
      ret += (tmp-static_cast<double>(i+1))*lambda_.coeff(i) - 0.5*exp(lambda_.coeff(i))/pow(scaleChol_.coeff(i,i),2);
    }
    // prior on V
    for(int i=0;i<n_-1;i++){
      ret += -0.5*F_.col(i).tail(n_-i-1).squaredNorm();
    }
    // rescaling Jacobian
    if(dr_) ret += logDetJac_;
    return(ret);
  }
  
  Eigen::Matrix<typename boost::math::tools::promote_args<df_type__,
                                                          scale_type__,
                                                          lambda_type__, 
                                                          V_type__>::type,
                                                          Eigen::Dynamic,
                                                          Eigen::Dynamic> 
  L() const {return(LP_);}
  
  template <class T__>
  void LTsolveInPlace(Eigen::Matrix<T__,Eigen::Dynamic,Eigen::Dynamic> &arg){
    LP_.template triangularView<Eigen::Lower>().transpose().solveInPlace(arg);
  }
  
  template <class T__>
  void LTsolveInPlace(Eigen::Matrix<T__,Eigen::Dynamic,1> &arg){
    LP_.template triangularView<Eigen::Lower>().transpose().solveInPlace(arg);
  }
  
  template <class T__>
  Eigen::Matrix<T__,Eigen::Dynamic,Eigen::Dynamic> 
  LTsolve(const Eigen::Matrix<T__,Eigen::Dynamic,Eigen::Dynamic> &lhs){
    Eigen::Matrix<T__,Eigen::Dynamic,Eigen::Dynamic> ret = lhs;
    LTsolveInPlace(ret);
    return(ret);
  }
  
  template <class T__>
  Eigen::Matrix<T__,Eigen::Dynamic,1> 
  LTsolve(const Eigen::Matrix<T__,Eigen::Dynamic,1> &lhs){
    Eigen::Matrix<T__,Eigen::Dynamic,1> ret = lhs;
    LTsolveInPlace(ret);
    return(ret);
  }
  
  
  // this function should be avoided as far as possible, i.e. use P_double for generated quantities
  Eigen::Matrix<typename boost::math::tools::promote_args<df_type__,
                                                          scale_type__,
                                                          lambda_type__, 
                                                          V_type__>::type,
                                                          Eigen::Dynamic,
                                                          Eigen::Dynamic> 
  P() const {
    return(LP_.template triangularView<Eigen::Lower>()*LP_.transpose());
  }
  
  Eigen::MatrixXd P_double() const{
    Eigen::MatrixXd ret = doubleValue(LP_);
    ret *= ret.template triangularView<Eigen::Lower>().transpose();
    return(ret);
  }
  
  
};


/*
 * to do: DR, AD-types for df and scale (requires normalization factor of logProb)
 * 
 */

} // namespace
#endif


