
#ifndef __EXTRKN_HPP__
#define __EXTRKN_HPP__


/*
 * The following class implements the order 6,4,6 Runge Kutta Nystöm triple
 * RKN6(4)6FD taken from Dormand, Prince, Runge-Kutta-Nystrom Triples,
 * Comput. Math. .Applic.Vol. 13, No. 12, pp. 937-949, 1987.
 * 
 * It provides an embedded order 6 and 4 pair, along with dense
 * formula that allow evaluation at each time-point within the
 * current step, with error also of order 6
 * 
 */
template < class _ADtarget_type_, 
           class _lambda_type_,
           class _massmatrix_type_>
class extRKN64{
private:
  _ADtarget_type_* t_;
  _lambda_type_* lam_;
  _massmatrix_type_* M_;
  rng* r_;
  
  double absTol_;
  double relTol_;
  double energyTol_;
  int d_;
  int dg_;
  
  double eps_;
  
  Eigen::MatrixXd fs_;
  Eigen::MatrixXd grads_;
  Eigen::MatrixXd gradsSq_;
  Eigen::MatrixXd gens_;
  Eigen::MatrixXd qs_;
  
  Eigen::VectorXd p0_,p1_,ptmp_;
  Eigen::VectorXd qd0_,qd1_;
  Eigen::VectorXd r0_,r1_,r1low_,rtmp_;
  Eigen::VectorXd sqGradInt1_,sqGradInt1Low_;
  Eigen::VectorXd M1_,M2_;
  Eigen::VectorXd q1low_,qd1low_;
  Eigen::VectorXd absDev_,targErr_,adaptRet_;
  
  Eigen::Matrix<double,1,6> funvals_;
  Eigen::Matrix<double,1,6> lams_;
  
  double H0_,H1_;
  double Ekin0_,Ekin1_;
  double Lam0_,Lam1_,Lam1low_;
  
  double err_;
  double Herr_;
  double errOld_;
  double PI_beta_;
  
  /*
   * Polynomials used for interpolation purposes
   */
  // for position
  Eigen::Matrix<double,6,1> bqstar(const double s) const {
    Eigen::Matrix<double,6,1> ret;
    ret << 
      0.5+(-1.2450142450142450142+(1.5161443494776828110+(-.90669515669515669516+.21367521367521367521*s)*s)*s)*s,
      0.0,
      (1.9859371988705031985+(-3.5776528722949209438+(2.5384011164109503556-.65839370179676583254*s)*s)*s)*s,
      (-1.8559358276926614325+(5.5475867510523291130+(-5.1641335539080937350+1.5949081677229964647*s)*s)*s)*s,
      (2.1786492374727668845+(-6.9269873191441818893+(7.2733366851013909837-2.5138260432378079437*s)*s)*s)*s,
      (-1.0636363636363636364+(3.4409090909090909091+(-3.7409090909090909091+1.3636363636363636364*s)*s)*s)*s;
    
    return(ret);
  }
  // for position first derivative
  Eigen::Matrix<double,6,1> bqdstar(const double s) const {
    Eigen::Matrix<double,6,1> ret;
    ret << 
      1.0+(-3.7350427350427350427+(6.0645773979107312441+(-4.5334757834757834758+1.2820512820512820513*s)*s)*s)*s,
      0.0,
      (5.9578115966115095956+(-14.310611489179683776+(12.692005582054751778-3.9503622107805949953*s)*s)*s)*s,
      (-5.5678074830779842974+(22.190347004209316452+(-25.820667769540468676+9.5694490063379787879*s)*s)*s)*s,
      (6.5359477124183006536+(-27.707949276576727557+(36.366683425506954919-15.082956259426847662*s)*s)*s)*s,
      (-3.1909090909090909091+(13.763636363636363636+(-18.704545454545454545+8.1818181818181818182*s)*s)*s)*s;
    
    return(ret);
  }
  
  // derivative wrt s of bqdstar
  Eigen::Matrix<double,6,1> dbqdstar(const double s) const {
    Eigen::Matrix<double,6,1> ret;
    ret << -3.7350427350427350427+(12.129154795821462488+(-13.600427350427350427+5.1282051282051282052*s)*s)*s,
           0.0,
           5.9578115966115095956+(-28.621222978359367552+(38.076016746164255334-15.801448843122379981*s)*s)*s,
           -5.5678074830779842974+(44.380694008418632904+(-77.462003308621406028+38.277796025351915152*s)*s)*s,
           6.5359477124183006536+(-55.415898553153455114+(109.10005027652086476-60.331825037707390648*s)*s)*s,
             -3.1909090909090909091+(27.527272727272727272+(-56.113636363636363635+32.727272727272727272*s)*s)*s;
      
    return(ret);
  }
  
public:
  extRKN64() : absTol_(1.0e-3), relTol_(1.0e-3), energyTol_(1.0e-3), eps_(0.1), d_(0), dg_(0),
   PI_beta_(0.04) {}
  void setAbsTol(const double absTol){absTol_=absTol;}
  double getAbsTol() const {return absTol_;}
  void setRelTol(const double relTol){relTol_=relTol;}
  double getRelTol() const {return relTol_;}
  void setup(_ADtarget_type_& t,
             _lambda_type_& lam,
             _massmatrix_type_& M,
             rng& r){
    t_ = &t;
    lam_ = &lam;
    M_ = &M;
    r_ = &r;
    d_ = (*t_).dim();
    dg_ = (*t_).dimGenerated();
    fs_.resize(d_,6);
    grads_.resize(d_,6);
    
    gens_.resize(dg_,6);
    qs_.resize(d_,6);
    
    p0_.resize(d_);
    p1_.resize(d_);
    ptmp_.resize(d_);
    qd0_.resize(d_);
    qd1_.resize(d_);
    
    r0_.resize(dg_);
    r1_.resize(dg_);
    rtmp_.resize(dg_);
    
    r1low_.resize(dg_);
    q1low_.resize(d_);
    qd1low_.resize(d_);
    
    
  }
  
  inline double denseLam(const double t) const {
    return(Lam0_ + t*(lams_.dot(bqdstar(t/eps_))));
  }
  // derivative of denseLam wrt t, which gives interpolated values of \lambda
  inline double denseDLam(const double t) const {
    double s = t/eps_;
    return(lams_.dot(bqdstar(s)) + s*lams_.dot(dbqdstar(s)));
  }
  
  inline Eigen::VectorXd denseP(const double t) const {
    return(p0_ + t*(grads_*bqdstar(t/eps_)));
  }
  inline Eigen::VectorXd denseR(const double t) const {
    return(r0_ + t*(gens_*bqdstar(t/eps_)));
  }
  inline Eigen::VectorXd denseSqGradInt(const double t) const {
    return(t*(gradsSq_*bqdstar(t/eps_)));
  }
  inline Eigen::VectorXd denseDR(const double t) const {
    double s = t/eps_;
    return(gens_*bqdstar(s) + gens_*(s*dbqdstar(s)));
  }
 
  inline Eigen::VectorXd denseQ(const double t) const {
    return(qs_.col(0) + t*(qd0_+fs_*(t*bqstar(t/eps_))));
  }
  
  // solves \Lambda = u for time t
  double solveLam(const double u) const {
    if(u<Lam0_ || u>Lam1_ || eps_<1.0e-14){
      std::cout << "bad problem in extRKN64::solveLam" << std::endl;
      return(0.5*eps_);
    }
    // binary search
    double ttry;
    double tleft = 0.0;
    double tright = eps_;
    while(tright-tleft>eps_*1.0e-4){
      ttry = 0.5*(tleft+tright);
      if(denseLam(ttry)-u>0.0){
        tright = ttry;
      } else {
        tleft = ttry; 
      }
     /* std::cout << tleft << "  " << tright << std::endl;
      std::cout << denseLam(0.5*(tleft+tright))-u << std::endl;*/
    }
    // do some additional Newton iterations
    ttry = 0.5*(tleft+tright);
    for(int i=0;i<4;i++){
      ttry -= (denseLam(ttry)-u)/denseDLam(ttry);
      ttry = fmin(tright,fmax(ttry,tleft));
      //std::cout << "residual " << denseLam(ttry)-u << std::endl;
    }
    double res = denseLam(ttry)-u;
    if(fabs(res)>1.0e-8*u){
      std::cout << "residual in extRKN64::solveLam = " << res << std::endl;
    }
    
    return(ttry);
  }
  
  inline void setEps(const double eps){eps_ = fabs(eps);}
  inline double getEps() const {return eps_;}
  inline double getH0() const {return H0_;}
  inline double getH1() const {return H1_;}
  inline double getEkin0() const {return Ekin0_;}
  inline double getEkin1() const {return Ekin1_;}
  inline double getLp0() const {return funvals_(0);}
  inline double getLp1() const {return funvals_(5);}
  inline double getEnergyErr() const {return Herr_;} 
  inline double getLam1() const {return Lam1_;}
  inline Eigen::VectorXd getP1() const {return p1_;}
  inline Eigen::VectorXd getQ1() const {return qs_.col(5);}
  inline Eigen::VectorXd getQ0() const {return qs_.col(0);}
  inline Eigen::VectorXd getR1() const {return r1_;}
  inline Eigen::VectorXd getSqGradInt1() const {return sqGradInt1_;} 
  inline Eigen::VectorXd getM1() const {return M1_;} 
  inline Eigen::VectorXd getM2() const {return M2_;} 
  inline Eigen::VectorXd getGen0() const {return gens_.col(0);}
  inline void getISGsample(Eigen::VectorXd &sample) const {sample = (1.0/eps_)*sqGradInt1_;}
  inline void getM1sample(Eigen::VectorXd &sample) const {sample = (1.0/eps_)*M1_;}
  inline void getM2sample(Eigen::VectorXd &sample) const {sample = (1.0/eps_)*M2_;}
  
  
  
  inline void adaptEps(){
    if(err_>0.0){
      eps_ *= fmin(4.0,fmax(0.1,0.9*pow(err_,-(0.2-0.75*PI_beta_))*pow(errOld_,PI_beta_)));
      errOld_ = err_;
    } else {
      eps_ *= 0.1;
      errOld_ = 1.0; // switch to P-controller
    }
  }
  
  inline void energyAdaptEps(){
    if(isfinite(Herr_) && err_> 0.0 && Herr_<100.0 ){
      eps_ *= fmin(1.5,fmax(0.1,0.9*pow(Herr_,-0.1428571)));  // -1/7
    } else {
      eps_ *= 0.1;
    }
  }
  
  int firstEval(const Eigen::VectorXd& q0, 
                Eigen::VectorXd& p0,
                const bool updateMomentum,
                const bool computeGrad=true){
    errOld_ = 1.0;
    qs_.col(0) = q0;
    
    
    r0_.setZero();
    
    
    Lam0_ = 0.0;
    if(computeGrad){
      funvals_(0) = (*t_).geval(q0,gens_.col(0),grads_.col(0));
    }
    
    if(updateMomentum){
      (*lam_).momentumUpdate(*r_,q0,p0,grads_.col(0),(*M_));
    }
    
    p0_ = p0;
   
    (*M_).applyMinv(p0_,qd0_); // qd0 = Minv*p0
    Ekin0_ = 0.5*(p0_.dot(qd0_));
    H0_ = -funvals_(0) + Ekin0_;
    
    (*M_).applyMinv(grads_.col(0),fs_.col(0)); 
    
    lams_(0) = (*lam_).lambda(q0,p0,grads_.col(0),(*M_));
    
    return(!(isfinite(funvals_(0)) && isfinite(lams_(0))));
  }
  
  
  
  /*
   * The RKN prepNext() method assumes that this call has been preceeded by
   * a call to step()
   * 
   */
  void prepNext(const bool updateStepSize){
    // update time step size
    
    if(updateStepSize) adaptEps(); 
    
    // maintain FSAL
    qs_.col(0) = qs_.col(5);
    fs_.col(0) = fs_.col(5);
    grads_.col(0) = grads_.col(5);
    gens_.col(0) = gens_.col(5);
    funvals_(0) = funvals_(5);
    lams_(0) = lams_(5);
    
    // book keeping
    H0_ = H1_;
    Ekin0_ = Ekin1_;
    p0_ = p1_;
    qd0_ = qd1_;
    r0_ = r1_;
    
    Lam0_ = Lam1_;
   
  }
  
  void energyPrepNext(const bool updateStepSize){
    if(updateStepSize) energyAdaptEps();
    prepNext(false);
  }
  
  
  /*
   * The RKN step() method assumes that this call is preceeded by either firstEval()
   * or step() and prepNext()
   */
  
  
  Eigen::VectorXd massAdjust(){
    absDev_ = (qs_.col(5)-q1low_).cwiseAbs();
    targErr_ = ((relTol_*(qs_.col(5).cwiseAbs().cwiseMax(qs_.col(0).cwiseAbs()))).array()+absTol_).matrix();
    adaptRet_ = absDev_.cwiseQuotient(targErr_);
 //   std::cout << adaptRet_ << std::endl; 
    absDev_ = (qd1_-qd1low_).cwiseAbs();
    targErr_ = ((relTol_*(qd1_.cwiseAbs().cwiseMax(qd0_.cwiseAbs()))).array()+absTol_).matrix();
    adaptRet_ = (absDev_.cwiseQuotient(targErr_)).cwiseMax(adaptRet_);
   // std::cout << adaptRet_ << std::endl; 
   /* absDev_ = (sqGradInt1_-sqGradInt1Low_).cwiseAbs();
    targErr_ = ((relTol_*(sqGradInt1_.cwiseAbs())).array()+absTol_).matrix();
    adaptRet_ = (absDev_.cwiseQuotient(targErr_)).cwiseMax(adaptRet_);
    std::cout << adaptRet_ << std::endl; 
    std::cout << sqGradInt1_  << "\n ld \n" << sqGradInt1Low_ << std::endl;*/
    adaptRet_ = adaptRet_.array().pow(-0.4).matrix();
    
    //ptmp_ = ((1.0/eps_)*sqGradInt1_).cwiseInverse().array().sqrt().matrix();
    //ptmp_ = ((relTol_*(qs_.col(5).cwiseAbs().cwiseMax(qs_.col(0).cwiseAbs()))).array()+absTol_).matrix();
    //ptmp_ = ((qs_.col(5)-q1low_).cwiseAbs().cwiseQuotient(ptmp_)).array().pow(-0.4).matrix();
    std::cout << "correction max = " << adaptRet_.maxCoeff() << "\t  min = " << adaptRet_.minCoeff() << std::endl;
    return(adaptRet_);
  }
  
  
  double step(bool intSqGrad, bool intM2){
    double epssq = pow(eps_,2);
    
    // first eval
    qs_.col(1) = qs_.col(0) + (0.12929590313670441529*eps_)*qd0_+ 
      epssq*0.008358715283968025328*fs_.col(0);
    funvals_(1) = (*t_).geval(qs_.col(1),gens_.col(1),grads_.col(1));
    (*M_).applyMinv(grads_.col(1),fs_.col(1));
    if(!isfinite(funvals_(1)) || ! fs_.col(1).array().isFinite().all()){
      err_ = -1.0;
      return(err_);
    }
    // second eval
    qs_.col(2) = qs_.col(0) + (0.25859180627340883058*eps_)*qd0_ + 
      (epssq*0.011144953711957367104)*fs_.col(0) + 
      (epssq*0.022289907423914734209)*fs_.col(1);
    funvals_(2) = (*t_).geval(qs_.col(2),gens_.col(2),grads_.col(2));
    (*M_).applyMinv(grads_.col(2),fs_.col(2));
    if(!isfinite(funvals_(2)) || ! fs_.col(2).array().isFinite().all()){
      err_ = -2.0;
      return(err_);
    }
    // third eval
    qs_.col(3) = qs_.col(0) + (0.67029708261548005830*eps_)*qd0_ +
      (epssq*0.14547474280109178590)*fs_.col(0) + 
      (epssq*(-0.22986064052264747311))*fs_.col(1) +
      (epssq*0.30903498720296753653)*fs_.col(2);
    funvals_(3) = (*t_).geval(qs_.col(3),gens_.col(3),grads_.col(3));
    (*M_).applyMinv(grads_.col(3),fs_.col(3));
    if(!isfinite(funvals_(3)) || ! fs_.col(3).array().isFinite().all()){
      err_ = -3.0;
      return(err_);
    } 
    // fourth eval
    qs_.col(4) = qs_.col(0) + (0.9*eps_)*qd0_ + 
      (epssq*(-0.20766826295078995434))*fs_.col(0) + 
      (epssq*0.68636678429251431227)*fs_.col(1) + 
      (epssq*(-0.19954927787234925220))*fs_.col(2) + 
      (epssq*0.12585075653062489426)*fs_.col(3);
    funvals_(4) = (*t_).geval(qs_.col(4),gens_.col(4),grads_.col(4));
    (*M_).applyMinv(grads_.col(4),fs_.col(4));
    if(!isfinite(funvals_(4)) || ! fs_.col(4).array().isFinite().all()){
      err_ = -4.0;
      return(err_);
    } 
    // fifth eval 
    qs_.col(5) = qs_.col(0) + eps_*qd0_ + 
      (epssq*0.078110161443494776828)*fs_.col(0) + 
      (epssq*0.28829174118976677768)*fs_.col(2) +
      (epssq*0.12242553717457041018)*fs_.col(3) + 
      (epssq*0.011172560192168035305)*fs_.col(4);
    funvals_(5) = (*t_).geval(qs_.col(5),gens_.col(5),grads_.col(5));
    (*M_).applyMinv(grads_.col(5),fs_.col(5));
    if(!isfinite(funvals_(5)) || ! fs_.col(5).array().isFinite().all()){
      err_ = -5.0;
      return(err_);
    } 
    /* Done evals; qs_.col(5) is the high order position approximation
     * 
     */
   if(intSqGrad) gradsSq_ = grads_.array().square().matrix();
   
    
    // high order q-dot (used for advancing the integrator)
    qd1_ = qd0_ + (eps_*0.078110161443494776828)*fs_.col(0) + 
      (eps_*0.38884347870598260272)*fs_.col(2) +
      (eps_*0.37132075792884226740)*fs_.col(3) + 
      (eps_*0.11172560192168035305)*fs_.col(4) + 
      (eps_*0.05)*fs_.col(5);
    // high order momentum
    
    (*M_).applyM(qd1_,p1_);
    Ekin1_ = 0.5*qd1_.dot(p1_);
    H1_ = -funvals_(5) + Ekin1_;
    
    // high order integrated generated quantities
    r1_ = r0_ + (eps_*0.078110161443494776828)*gens_.col(0) + 
      (eps_*0.38884347870598260272)*gens_.col(2) +
      (eps_*0.37132075792884226740)*gens_.col(3) + 
      (eps_*0.11172560192168035305)*gens_.col(4) + 
      (eps_*0.05)*gens_.col(5);
    
    if(intSqGrad){
      
      sqGradInt1_ = (eps_*0.078110161443494776828)*gradsSq_.col(0) + 
        (eps_*0.38884347870598260272)*gradsSq_.col(2) +
        (eps_*0.37132075792884226740)*gradsSq_.col(3) + 
        (eps_*0.11172560192168035305)*gradsSq_.col(4) + 
        (eps_*0.05)*gradsSq_.col(5);
      /*
      sqGradInt1Low_ = (eps_*0.054605887939221272555)*gradsSq_.col(0) +
        (eps_*0.46126678590362684429)*gradsSq_.col(2) + 
        (eps_*0.19588085947931265629)*gradsSq_.col(3) +
        (eps_*0.38824646667783922686)*gradsSq_.col(4) + 
        (eps_*(-0.1))*gradsSq_.col(5); */
    }
    if(intM2){ // first and second order integrated moments 
      M1_ = (eps_*0.078110161443494776828)*qs_.col(0) + 
        (eps_*0.38884347870598260272)*qs_.col(2) +
        (eps_*0.37132075792884226740)*qs_.col(3) + 
        (eps_*0.11172560192168035305)*qs_.col(4) + 
        (eps_*0.05)*qs_.col(5);
      
      M2_ = (eps_*0.078110161443494776828)*qs_.col(0).array().square().matrix() + 
        (eps_*0.38884347870598260272)*qs_.col(2).array().square().matrix() +
        (eps_*0.37132075792884226740)*qs_.col(3).array().square().matrix() + 
        (eps_*0.11172560192168035305)*qs_.col(4).array().square().matrix() + 
        (eps_*0.05)*qs_.col(5).array().square().matrix();
    }
    
    // low order approximations for error estimates
    q1low_ = qs_.col(0) + eps_*qd0_ +
      (epssq*1.0588592603704182782)*fs_.col(0) + 
      (epssq*(-2.4067513719244520532))*fs_.col(1) +
      (epssq*1.8478921115540337750)*fs_.col(2);
    
    qd1low_ = qd0_ + (eps_*0.054605887939221272555)*fs_.col(0) +
      (eps_*0.46126678590362684429)*fs_.col(2) + 
      (eps_*0.19588085947931265629)*fs_.col(3) +
      (eps_*0.38824646667783922686)*fs_.col(4) + 
      (eps_*(-0.1))*fs_.col(5);  
    
    r1low_ = r0_ + (eps_*0.054605887939221272555)*gens_.col(0) +
      (eps_*0.46126678590362684429)*gens_.col(2) + 
      (eps_*0.19588085947931265629)*gens_.col(3) +
      (eps_*0.38824646667783922686)*gens_.col(4) + 
      (eps_*(-0.1))*gens_.col(5);  
    
    
    
    // compute integrated lambda:
    (*M_).applyM(qd0_+(0.12929590313670441529*eps_)*fs_*bqdstar(0.12929590313670441529),ptmp_); //interpolated p
    
    lams_(1) = (*lam_).lambda(qs_.col(1),ptmp_,grads_.col(1),*M_);
    
    (*M_).applyM(qd0_+(0.25859180627340883058*eps_)*fs_*bqdstar(0.25859180627340883058),ptmp_);//interpolated p
    lams_(2) = (*lam_).lambda(qs_.col(2),ptmp_,grads_.col(2),*M_);
    
    (*M_).applyM(qd0_+(0.67029708261548005830*eps_)*fs_*bqdstar(0.67029708261548005830),ptmp_);//interpolated p
    lams_(3) = (*lam_).lambda(qs_.col(3),ptmp_,grads_.col(3),*M_);
    
    (*M_).applyM(qd0_+(0.9*eps_)*fs_*bqdstar(0.9),ptmp_);//interpolated p
    lams_(4) = (*lam_).lambda(qs_.col(4),ptmp_,grads_.col(4),*M_);
    
    lams_(5) = (*lam_).lambda(qs_.col(5),p1_,grads_.col(5),*M_);
    
    Lam1_ = Lam0_ + (eps_*0.078110161443494776828)*lams_(0) + 
      (eps_*0.38884347870598260272)*lams_(2) +
      (eps_*0.37132075792884226740)*lams_(3) + 
      (eps_*0.11172560192168035305)*lams_(4) + 
      (eps_*0.05)*lams_(5);
    
    Lam1low_ = Lam0_ + (eps_*0.054605887939221272555)*lams_(0) +
      (eps_*0.46126678590362684429)*lams_(2) + 
      (eps_*0.19588085947931265629)*lams_(3) +
      (eps_*0.38824646667783922686)*lams_(4) + 
      (eps_*(-0.1))*lams_(5);

    // error estimation
    // from integrated lambda
    //double sqrtd = sqrt(static_cast<double>(d_));
    //double sqrtdg = sqrt(static_cast<double>(dg_));
    err_ = fabs(Lam1_-Lam1low_)/(absTol_+relTol_*fmax(fabs(Lam1_),fabs(Lam0_)));   
    // from position
    ptmp_ = ((relTol_*(qs_.col(5).cwiseAbs().cwiseMax(qs_.col(0).cwiseAbs()))).array()+absTol_).matrix();
    err_ = fmax(err_,(qs_.col(5)-q1low_).cwiseAbs().cwiseQuotient(ptmp_).maxCoeff());
    //err_ = fmax(err_,(qs_.col(5)-q1low_).cwiseQuotient(ptmp_).norm()/sqrtd);
    // from position dot
    ptmp_ = ((relTol_*(qd1_.cwiseAbs().cwiseMax(qd0_.cwiseAbs()))).array()+absTol_).matrix();
    err_ = fmax(err_,(qd1_-qd1low_).cwiseAbs().cwiseQuotient(ptmp_).maxCoeff());
    //err_ = fmax(err_,(qd1_-qd1low_).cwiseQuotient(ptmp_).norm()/sqrtd);
    
    // from integrated generated quantities
    rtmp_ = ((relTol_*(r1_.cwiseAbs().cwiseMax(r0_.cwiseAbs()))).array()+absTol_).matrix();
    err_ = fmax(err_,(r1_-r1low_).cwiseAbs().cwiseQuotient(rtmp_).maxCoeff());
    //err_ = fmax(err_,(r1_-r1low_).cwiseQuotient(rtmp_).norm()/sqrtdg);
    // energy-based error estimate
    Herr_ = fabs(H1_-H0_)/energyTol_;
    
    return(err_);  
  }
  
  void dump() const {
    std::cout << "step size eps = " << eps_ << std::endl;
     
    std::cout << "q1 = \n" << qs_.col(5) << std::endl;
    std::cout << "p1 = \n" << p1_ << std::endl;
    std::cout << "r1 = \n" << r1_ << std::endl;
    std::cout << "sqGradInt1 = \n" << sqGradInt1_ << std::endl;
    
    std::cout << "Lam1 = \n" << Lam1_ << std::endl;
    
    std::cout << "q1low = \n" << q1low_ << std::endl;
    std::cout << "qd1low = \n" << qd1low_ << std::endl;
    std::cout << "r1low = \n" << r1low_ << std::endl;
    
    std::cout << "Lam1low = \n" << Lam1low_ << std::endl;
    
    std::cout << "err = \n" << err_ << std::endl;
    std::cout << "H1 = \n " << H1_ << std::endl;
    std::cout << "Herr = \n" << H1_-H0_ << std::endl;
    
  }
  
  
};


#endif




