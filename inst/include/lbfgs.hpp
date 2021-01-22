
#define __LBFGS_DEBUG__ 0


template <class _ADtarget_type_>
class lbfgs{
  _ADtarget_type_* t_;
  bool maximize_;
  int d_;
  int m_;
  int eflag_;
  Eigen::MatrixXd x_curr_,p_curr_,gens_curr_,grad_curr_;
  Eigen::MatrixXd x_try_,p_try_,gens_try_,grad_try_;
  double f_curr_,f_try_;
  Eigen::MatrixXd sk_;
  Eigen::MatrixXd yk_;
  Eigen::VectorXd rhok_;
  int iterSinceRestart_;
  
  
  double lineSearch(const double alpha0){
    const double c1 = 1.0e-4;
    double alpha = alpha0;
    double alphaprop;
    Eigen::MatrixXd dirder = grad_curr_.transpose()*p_curr_;
    // locate upper bound on search line
    for(int i=0;i<30;i++){
      x_try_ = x_curr_ + alpha*p_curr_;
      f_try_ = (*t_).feval(x_try_.col(0),gens_try_.col(0));

      if(maximize_){
        f_try_ = -f_try_;
        grad_try_ = -grad_try_;
      }
      
#ifdef __LBFGS_DEBUG__    
      std::cout << "linsesearch with alpha = " << alpha << " fval = " << f_try_ << std::endl;
#endif
      
      if(!isfinite(f_try_) || !grad_curr_.col(0).array().isFinite().all()){
        alpha *= 0.5;
      } else if(f_try_<= f_curr_ + c1*alpha*dirder(0,0)){
#ifdef __LBFGS_DEBUG__  
        std::cout << "linsesearch exited OK" << std::endl;
#endif
        return alpha;
        
      } else {
        // do backtracking based on quadratic polynomial
        alphaprop = 0.5*dirder(0,0)*pow(alpha,2)/(f_curr_ + alpha*dirder(0,0) - f_try_);
#ifdef __LBFGS_DEBUG__  
        std::cout << "alphaprop = " << alphaprop << std::endl;
#endif       
        alpha = fmin(0.9*alpha,fmax(0.1*alpha,alphaprop));
      }
    }
    
    return -1.0;
    
  }
  
  
  Eigen::VectorXd parabolic_roots(const double a, const double fa,
                                    const double b, const double fb,
                                    const double c, const double fc){
    Eigen::VectorXd pcoeffs(3),ret(3);
    double t1 = b * b;
    double t2 = c * c;
    double t3 = b - c;
    double t4 = a - b;
    double t5 = a - c;
    double t6 = -fb + fc;
    double t7 = -fa + fc;
    double t8 = fa - fb;
    t4 = 0.1e1 / t4;
    t5 = 0.1e1 / t5;
    double t9 = 0.1e1 / t3;
    pcoeffs[0] = (((b * fc - c * fb) * a + fb * t2 - fc * t1) * a + b * c * fa * t3) * t4 * t9 * t5;
    pcoeffs[1] = -(a * a * t6 - t1 * t7 - t2 * t8) * t4 * t9 * t5;
    pcoeffs[2] = (a * t6 - b * t7 - c * t8) * t4 * t9 * t5;
    
    double dd = -4.0*pcoeffs(0)*pcoeffs(2)+ pow(pcoeffs(1),2);
    if(dd>=0 && fabs(pcoeffs(2))>1.0e-8){
      dd = sqrt(dd);
      ret(0) = 0.0;
      ret(1) = 0.5*(-pcoeffs(1)+dd)/pcoeffs(2);
      ret(2) = 0.5*(-pcoeffs(1)-dd)/pcoeffs(2);
    } else {
      ret(0) = 1.0;
      ret(1) = 0.0;
      ret(2) = 0.0;
    }
    return(ret);
  }
  
  double lineSearchExact(const double alpha0){
    const double c1 = 1.0e-4;
    const double c2 = 0.9;
    double dd0 = grad_curr_.col(0).dot(p_curr_.col(0));
    double dd1;
    double alpha = alpha0;
    bool evalGood;
    if(dd0>0.0){
      std::cout << "current direction has positive directional derivative!!!" << std::endl;
      return(-1.0);
    }
    
    
    x_try_ = x_curr_ + alpha*p_curr_;
    f_try_ = (*t_).geval(x_try_.col(0),gens_curr_.col(0),grad_try_.col(0));
    if(maximize_){
      f_try_ = -f_try_;
      grad_try_ = -grad_try_;
    }
    evalGood = isfinite(f_try_) && grad_try_.col(0).array().isFinite().all();
    dd1 = grad_try_.col(0).dot(p_curr_.col(0));
    
    bool wolfe1 = f_try_<= f_curr_ + c1*alpha*dd0;
    bool wolfe2 = dd1 >= c2*dd0;
    
    if(wolfe1 && wolfe2 && evalGood){
#ifdef __LBFGS_DEBUG__       
      std::cout << "both wolfe conditions OK for intial alpha = " << alpha << std::endl;
#endif
      return(alpha);
    } 
#ifdef __LBFGS_DEBUG__     
    std::cout << "wolfe conditions not OK for initial alpha" << std::endl;
#endif
    
    // shorten initial bracket if first eval failed
    if(!evalGood) alpha *= 0.05;
    
    // first try to bracket minimizer
    double a,c,old;
    double fa,fc,fold;
    
    a = 0.0;
    fa = f_curr_;
    if(dd1<0.0 || !evalGood){
      // search further
      for(int biter=0;biter<20;biter++){
        alpha *= 2.0;
        x_try_ = x_curr_ + alpha*p_curr_;
        f_try_ = (*t_).geval(x_try_.col(0),gens_curr_.col(0),grad_try_.col(0));
        if(maximize_){
          f_try_ = -f_try_;
          grad_try_ = -grad_try_;
        }
        evalGood = isfinite(f_try_) && grad_try_.col(0).array().isFinite().all();
        dd1 = grad_try_.col(0).dot(p_curr_.col(0));
        if(dd1>=0.0 && evalGood){
#ifdef __LBFGS_DEBUG__ 
          std::cout << "found initial bracket with alpha = " << alpha << std::endl;
#endif
          break;
        } else if(dd1<0.0 && evalGood) {
          a = alpha;
          fa = dd1;
        } else {
          alpha *= 0.1;
        }
      }
    } else {
#ifdef __LBFGS_DEBUG__ 
      std::cout << "minimizer already bracketed" << std::endl;
#endif
      a = 0.0;
      fa = dd0;
    }
    
    // gets called if above loop is not exited by the break
    if(dd1<0.0 || !evalGood){
      std::cout << "line search failed, could not bracket minimizer" << std::endl;
      return(-1.0);
    }
    
    // do bisection search
    
    c = alpha;
    fc = dd1;
    double aa,bb,cc,pcand;
    double ddThresh = fmin(1.0,1.0e-1*c2*fabs(dd0));
    Eigen::VectorXd pret(3);
    for(int kk=0;kk<30;kk++){

      if(kk>0 && kk%2==0){
        // linear interpolation root if parabolic fails
        alpha = (a*fc-fa*c)/(fc-fa);
        // try parabolic
        pret = parabolic_roots(a,fa,c,fc,old,fold);
        if(pret(0)<0.1){
          // parabolic has real roots, now check if any of these are in current interval
          pcand = -1.0;
          if(pret(2)<c && pret(2)>a) pcand = pret(2);
          if(pret(1)<c && pret(1)>a) pcand = pret(1); // pick the largest root if possible
          if(pcand>0.0) alpha = pcand;
        } 
        
      } else {
        alpha = 0.5*(a+c);
      }
      x_try_ = x_curr_ + alpha*p_curr_;
      f_try_ = (*t_).geval(x_try_.col(0),gens_curr_.col(0),grad_try_.col(0));
      if(maximize_){
        f_try_ = -f_try_;
        grad_try_ = -grad_try_;
      }
      evalGood = isfinite(f_try_) && grad_try_.col(0).array().isFinite().all();
      if(!evalGood){
#ifdef __LBFGS_DEBUG__ 
        std::cout << "bad eval within bracket, exiting" << std::endl;
        return(-1.0);
#endif
      }
      dd1 = grad_try_.col(0).dot(p_curr_.col(0));
#ifdef __LBFGS_DEBUG__ 
      std::cout << "alpha = " << alpha << "\t dd1 = " << dd1 << std::endl;
#endif
      if(fabs(dd1)<ddThresh || c-a<1.0e-4+1.0e-3*(c+a)){
#ifdef __LBFGS_DEBUG__ 
        std::cout << "linesearch successful with alpha = " << alpha << std::endl;
#endif
        return(alpha);
      }
      if(dd1>0.0){
        old = c;
        fold = fc;
        c = alpha;
        fc = dd1;
      } else {
        old = a;
        fold = fa;
        a = alpha;
        fa = dd1;
      }
    }
    return -1.0;
  }
  
  
  
  
public:
  lbfgs(_ADtarget_type_& t,
          const int m, 
          const bool maximize, 
          const Eigen::VectorXd q0,
          const int print,
          const double gradTol = 1.0e-4,
          const double xTol = 1.0e-4,
          const int maxIter = 1000) : m_(m), maximize_(maximize) {
    eflag_ = 0;
    t_ = &t;
    // allocate storage
    d_ = q0.size();
    if(d_ != (*t_).dim()){
      std::cout << "bad dimension of initial point in lbfgs" << std::endl;
      _ERROR_LOG_.push("bad dimension of initial point in lbfgs");
      return;
    }
    Eigen::VectorXd qq(d_),rr(d_);
    x_curr_.resize(d_,1);
    grad_curr_.resize(d_,1);
    grad_try_.resize(d_,1);
    x_curr_.col(0) = q0;
    gens_curr_.resize((*t_).dimGenerated(),1);
    gens_try_.resize((*t_).dimGenerated(),1);
    sk_.resize(d_,m_);
    yk_.resize(d_,m_);
    rhok_.resize(m_);
    sk_.setZero();
    yk_.setZero();
    rhok_.setZero();
    // initial evaluation
    f_curr_ = (*t_).geval(x_curr_.col(0),gens_curr_.col(0),grad_curr_.col(0));
    if(!isfinite(f_curr_) || ! grad_curr_.col(0).array().isFinite().all()){
      _ERROR_LOG_.push("bad initial point in lbfgs"); 
      eflag_ = -1;
      return;
    }
    if(maximize_){
      f_curr_ = -f_curr_;
      grad_curr_ = -grad_curr_;
    }
    double maxgrad = grad_curr_.col(0).cwiseAbs().maxCoeff();
    if(print>0) std::cout << "initial point OK, fval = " << f_curr_ << " , max grad =" << maxgrad << std::endl;
    p_curr_ = -grad_curr_;
    
    double alpha,alpha0;
    double maxP,maxXdiff;
    double gamk,bfgs_beta;
    int writec,readc;
    iterSinceRestart_=0;
    Eigen::VectorXd bfgs_alp(m_);
    for(int iter = 1; iter<maxIter;iter++){
      iterSinceRestart_++;

      
      alpha = lineSearchExact(1.0);  
      if(alpha<=0.0){
        std::cout << "optimizer not making progress, exiting" << std::endl;
        eflag_ = -2;
        return;
      }
      
      maxgrad = grad_try_.col(0).cwiseAbs().maxCoeff();
      if(iter % print == 0){
        std::cout << "iteration # " << iter << " \t fval =  " << f_try_ 
                  << " \t max grad = " << maxgrad
                  << std::endl;
      }
      
     
      
      // BFGS stuff
      
      writec = (iterSinceRestart_-1) % m_;
      sk_.col(writec) = x_try_ - x_curr_;
      maxXdiff = sk_.col(writec).cwiseAbs().maxCoeff();
      yk_.col(writec) = grad_try_ - grad_curr_;
      gamk = sk_.col(writec).dot(yk_.col(writec));
      rhok_(writec) = 1.0/gamk;
      gamk = gamk/(yk_.col(writec).squaredNorm());
      x_curr_ = x_try_;
      grad_curr_ = grad_try_;
      f_curr_ = f_try_;
      
      // check termination criterion
      if(maxgrad < gradTol){
        eflag_ = 0;
        if(print>0) std::cout << "optimization exited: max grad element < gradTol" << std::endl;
        break;
      }
      
      std::cout << "maxXdiff : " << maxXdiff << std::endl;
      if(maxXdiff<xTol){
        eflag_ = 1;
        if(print>0) std::cout << "optimization exited: max difference in x < xTol" << std::endl;
        break;
      }
      
      
      // L-BFGS recurisons
      qq = grad_curr_;
      
      for(int k = 0;k<std::min(iterSinceRestart_,m_);k++){
        
        readc = (m_ + writec - k) % m_ ;
        bfgs_alp(readc) = rhok_(readc)*sk_.col(readc).dot(qq);
        qq -= bfgs_alp(readc)*yk_.col(readc);
      }
      rr = gamk*qq;
      for(int k=std::min(iterSinceRestart_,m_)-1;k>=0;k--){
        readc = (m_ + writec - k) % m_ ;
        bfgs_beta = rhok_(readc)*yk_.col(readc).dot(rr);
        rr += (bfgs_alp(readc) - bfgs_beta)*sk_.col(readc);
      }
      p_curr_ = -rr;
      
      // check 
      
      maxP = p_curr_.cwiseAbs().maxCoeff();
      if(maxP<xTol && iter>10){
        eflag_ = 2;
        std::cout << "optimization exited: max search direction < xTol" <<  std::endl;
        break;
      }
/*
      if(iter%(5*m_)==0 || (maxP<xTol && iter<=10)){
        std::cout << "restarting" << std::endl;
        iterSinceRestart_ = 0;
        p_curr_ = -grad_try_;
      } */
    } // main bfgs iterations loop
  } // constructor
  
  int exitflag() const {return eflag_;}
  Eigen::VectorXd optX() const {
    if(eflag_<0){
      std::cout << "warning, attempting to access failed lbfgs object" << std::endl;
    }
    return x_curr_.col(0);
  }
  Eigen::VectorXd optGrad() const {
    if(eflag_<0){
      std::cout << "warning, attempting to access failed lbfgs object" << std::endl;
    }
    return grad_curr_.col(0);
  }
  
  
  /*
   * the diagonal elements of the BFGS updates
   * 
   */ 
/*  Eigen::VectorXd bfgsDiag(){
    Eigen::VectorXd D(d_);
    int writec = (iterSinceRestart_-1) % m_;
    double gamk = 1.0/rhok_(writec);
    gamk = gamk/(yk_.col(writec).squaredNorm());
    D.setConstant(d_,gamk);
    int readc;
    double rho;
    double ydy;
    for(int k=std::min(iterSinceRestart_,m_)-1;k>=0;k--){
      readc = (m_ + writec - k) % m_;
      rho = rhok_(readc);
      ydy = yk_.col(readc).array().square().matrix().dot(D);
      D += rho*(1.0+rho*ydy)*sk_.col(readc).array().square().matrix()
        - 2.0*rho*sk_.col(readc).cwiseProduct(D.cwiseProduct(yk_.col(readc)));
    }
    return(D);
    
  }
  */
  
}; // class lbfgs




