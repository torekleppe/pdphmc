

template <class _massmatrix_type_>
class lambda_constant{
private:
    double beta_;
    double fac_; 
public:
  lambda_constant() : beta_(1.0), fac_(2.0) {};
  
  inline double lambda(const Eigen::VectorXd &q, 
                       const Eigen::VectorXd &p, 
                       const Eigen::VectorXd &grad,
                       const _massmatrix_type_ &M){
    return 1.0/beta_;
  }
  // return draw \propto lambda(p)*N(p|0,M)
  inline void momentumUpdate(rng &r,
                             const Eigen::VectorXd &q, 
                             Eigen::VectorXd &p, 
                             const Eigen::VectorXd &grad,
                             const _massmatrix_type_ &M){
    r.rnorm(p);
    M.inPlaceApplySqrtM(p);
  }
  // return parameters, mainly for diagnostics purposes
  inline double getPar(const int par) const {
    if(par==0){
      return beta_;
    } else if(par==1) {
      return fac_; 
    } else {
      return 0.0;
    }
  }
  inline void setPar(const int parNum, const double parVal){
    if(parNum==0){
      beta_ = parVal;
    } else if(parNum==1){
      fac_ = parVal;
    }
  }
  inline void adapt(const Eigen::VectorXd &diagRow){
    // adaptation is performed using an exponential moving average in log(time_opt)
    double propnext = _LAMBDA_EMA_ALPHA_*log(diagRow.coeff(4)) + (1.0-_LAMBDA_EMA_ALPHA_)*log(beta_/fac_);
    beta_ = fmin(2.0*beta_,fmax(0.1*beta_,fac_*exp(propnext)));
  }
  inline void reset(){
    beta_ = 1.0;
  }
};

template <class _massmatrix_type_>
class lambda_constant_partial{
private:
  double beta_;
  double phi_;
  Eigen::VectorXd innov_;
public:
  lambda_constant_partial() : beta_(1.0), phi_(0.9) {};
  
  inline double lambda(const Eigen::VectorXd &q, 
                       const Eigen::VectorXd &p, 
                       const Eigen::VectorXd &grad,
                       const _massmatrix_type_ &M){
    return 1.0/beta_;
  }
  // return draw \propto lambda(p)*N(p|0,M)
  inline void momentumUpdate(rng &r,
                             const Eigen::VectorXd &q, 
                             Eigen::VectorXd &p, 
                             const Eigen::VectorXd &grad,
                             const _massmatrix_type_ &M){
    innov_ = r.rnorm(p.size());
    M.inPlaceApplySqrtM(innov_);
    p = phi_*p + sqrt(1.0-pow(phi_,2))*innov_;
  }
  // return parameters, mainly for diagnostics purposes
  inline double getPar(const int par) const {
    if(par==0){
      return beta_;
    } else if(par==1){
      return phi_;
    } else {
      return 0.0; 
    }
  }
  inline void setPar(const int parNum, const double parVal){
    if(parNum==0) {
      beta_ = parVal; 
    } else if(parNum==1) {
      phi_ = parVal;
    }
  }
  inline void adapt(const Eigen::VectorXd &diagRow){
    // adaptation is performed using an exponential moving average in log(time_opt)
    double propnext = _LAMBDA_EMA_ALPHA_*log(diagRow.coeff(4)) + (1.0-_LAMBDA_EMA_ALPHA_)*log(beta_/1.442695);
    beta_ = fmin(2.0*beta_,fmax(0.1*beta_,exp(propnext)*1.442695));
  }
  
  inline void reset(){
    beta_ = 1.0;
  }
  
};

template <class _massmatrix_type_>
class lambda_arclength{
private:
  double beta_;
  double fac_;
public:
  lambda_arclength() : beta_(1.0), fac_(2.0) {};
  
  inline double lambda(const Eigen::VectorXd &q, 
                       const Eigen::VectorXd &p, 
                       const Eigen::VectorXd &grad,
                       const _massmatrix_type_ &M){
    return sqrt(p.dot(M.applyMinv(p)))/beta_;
  }
  /* p = sqrt(M)*z*(r/norm(z))
   *  where z \sim N(0,I_d), r^2 \sim \chi^2_{d+1}
   */
  inline void momentumUpdate(rng &r,
                             const Eigen::VectorXd &q, 
                             Eigen::VectorXd &p, 
                             const Eigen::VectorXd &grad,
                             const _massmatrix_type_ &M){
    double scal = sqrt(r.rgamma(0.5*static_cast<double>(p.size()+1),2.0)); // \chi^2_{d+1}
    r.rnorm(p);
    scal/= p.norm();
    p *= scal;
    M.inPlaceApplySqrtM(p);
  }
  
  // return parameters, mainly for diagnostics purposes
  inline double getPar(const int par) const {
    if(par==0){
      return beta_;
    } else if(par==1) {
      return fac_; 
    } else {
      return 0.0;
    }
  }
  inline void setPar(const int parNum, const double parVal){
    if(parNum==0){
      beta_ = parVal;
    } else if(parNum==1){
      fac_ = parVal;
    }
  }
  
  inline void adapt(const Eigen::VectorXd &diagRow){
    // 
    double propnext = _LAMBDA_EMA_ALPHA_*log(diagRow.coeff(3)*diagRow.coeff(5)) + (1.0-_LAMBDA_EMA_ALPHA_)*log(beta_/fac_);
    beta_ = fmin(2.0*beta_,fmax(0.1*beta_,exp(propnext)*fac_));
  }
  inline void reset(){
    beta_ = 1.0;
  }
  
};


