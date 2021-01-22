#ifndef __PDPSAMPLER_HPP__
#define __PDPSAMPLER_HPP__

#include <boost/math/special_functions/gamma.hpp>

template <class _model_type_, 
          template <typename> class _ADwrapper_,
          template <typename,typename,typename> class _integrator_template_,
          template <typename> class _lambda_template_,
          class _massmatrix_type_>
class PDPsampler{
private:
  _ADwrapper_<_model_type_> adt_;
  _integrator_template_< _ADwrapper_< _model_type_ >, _lambda_template_<_massmatrix_type_>, _massmatrix_type_ > rkn_;
  //ADtarget<_model_type_> adt_;
  //_integrator_template_< ADtarget< _model_type_ >, _lambda_template_<_massmatrix_type_>, _massmatrix_type_ > rkn_;
  _lambda_template_<_massmatrix_type_> lambda_;
  _massmatrix_type_ M_;
  rng r_;
  
  Eigen::VectorXd qcurr_,qnew_,qsample_;
  Eigen::VectorXd pcurr_,pnew_;
  
  // dimensions
  int d_;
  int dg_;
  
  // MCMC stuff
  double time_; // process time
  double eventTime_;
  double u_; // exponential(1) rw
  double Tmax_; // length of simulation time
  
  Eigen::VectorXi storePars_;
  Eigen::VectorXd r_last_sampling_,r_last_trajectory_;
  
  Eigen::VectorXd rtmp_,ISG_sample_,M1_sample_,M2_sample_;
  int sampleCount_;
  int nSamples_;
  double samplesSpacing_;
  double invSamplesSpacing_;
  
  // diagnostics-related storage
  int maxTrajSteps_;
  int nstep_;
  int nacc_;
  int exitType_;
  int nevents_;
  double energyErr_;
  Eigen::MatrixXd intDiag_;
  
  
  
  // initial trajectories storage
  int initMassAdapt_;
  int initMomentumRenormalization_;
  int maxInitTraj_;
  double Lpdiff_;
  
  // lambda-adaptation related storage
  bool lambdaAdapt_;
  double time_opt_;
  double Lam_opt_;
  
  
  // print options
  std::string printPrefix_;
  
  double nextSampleTime() const {return static_cast<double>(sampleCount_)*samplesSpacing_;} 
  
  
  
  
  
  /*
   * trajectories intended for transient first part of simulation
   * -Uses less stringent error control 
   * -Does not record adaptation information
   */
/*  int initTrajectory(){
    int eflag = rkn_.firstEval(qcurr_,pcurr_,true);
    if(eflag != 0){
      std::cout << "bad first eval in PDPsampler::initTrajectory" << std::endl;
      return(1);  
    }
    
    double absTol = rkn_.getAbsTol();
    double relTol = rkn_.getRelTol();
    rkn_.setAbsTol(1.0e-3);
    rkn_.setRelTol(1.0e-3);
    
    // initial mass matrix adaptation
    double err;

    if(M_.allowsDiagAdaptation() && initMassAdapt_>0){
      Eigen::VectorXd p_std = pcurr_;
      Eigen::VectorXd MinvAdjust(d_);
      Eigen::VectorXd newMinv(d_);
      M_.inPlaceApplySqrtInvM(p_std);
      bool massAdaptOK = false;
      double epsRef = 1.0;
      if(d_>1) epsRef = 3.0*pow(1.0e-3,0.2)*pow(log(static_cast<double>(d_)),-0.1); // se maple work sheet ref_eps
      epsRef *= 0.5;
      rkn_.setEps(epsRef);
      std::cout << "Locating good mass matrix diagonal with epsRef = " << epsRef << std::endl;
      for(int iter=0;iter < initMassAdapt_;iter++){
        err = rkn_.step(true);
        std::cout << "step err: " << err << std::endl;
        
        if(err>0.0){
          MinvAdjust = rkn_.massAdjust();
          M_.applyMinv(MinvAdjust,newMinv);
          M_.acceptMinvDiag(newMinv,10.0,0.1);
        } else { 
          newMinv.setZero(); // effectively causes the Minv *= 0.1
          M_.acceptMinvDiag(newMinv,1.0,0.1);
        }
        
        std::cout << "Minv max : " << M_.MinvDiagMax() << "\t Minv min : " << M_.MinvDiagMin() << std::endl;
        pcurr_ = p_std;
        M_.inPlaceApplySqrtM(pcurr_);
        rkn_.firstEval(qcurr_,pcurr_,false,false);
        if(err>0.0 && MinvAdjust.minCoeff()>0.95 && MinvAdjust.maxCoeff()<1.05){
          std::cout << "Mass matrix determination successfull" << std::endl;
          massAdaptOK = true;
          break;
        }
      }
      if(! massAdaptOK){
        M_.reset();
      } else {
        rkn_.setEps(epsRef*0.9); // ensure first step gets accepted with high prob
      }
    }
    
    
    
    bool trajDone = false;
    bool stepGood;
    
    nstep_ = 0;
    nacc_ = 0;
    double Lpcurr;
    double Lpmax = rkn_.getLp0();
    double NUT; 
    
    double Lp0 = rkn_.getLp0();
    double Lp1 = Lp0;
    
    
    while(nstep_< maxTrajSteps_){
      stepGood = false;
      while(rkn_.getEps()>1.0e-14 && nstep_< maxTrajSteps_){
        err = rkn_.step(false);
        intDiag_.coeffRef(nstep_,0) = rkn_.getEps();
        intDiag_.coeffRef(nstep_,1) = err;
        intDiag_.coeffRef(nstep_,2) = rkn_.getLp1();
        nstep_++;
        if(err<1.0 && err>0.0){
          stepGood = true;
          nacc_++;
          break;
        }
        rkn_.adaptEps();
      }
      if(! stepGood){
        std::cout << "integrator is not making progress, rejecting trajectory" << std::endl;
        std::cout << "nstep = " << nstep_ << std::endl;
        std::cout << "nacc = " << nacc_ << std::endl;
        rkn_.dump();
        rkn_.setAbsTol(absTol);
        rkn_.setRelTol(relTol);
        exitType_ = -1;
        return(2);
      }
      
      // first stop criterion based on potential energy decrease
      
      Lpcurr = rkn_.getLp1();
      Lpmax = fmax(Lpmax,Lpcurr);
      if(Lpcurr < Lpmax - static_cast<double>(d_)){
        qcurr_ = rkn_.getQ0();
        Lp1 = rkn_.getLp0();
        r_.rnorm(pcurr_); // resample momentum in case partial refreshments are used
        break;
      }
      
      // second stop criteron based on NUT
      
      NUT = rkn_.getP1().dot(rkn_.getQ1()-qcurr_); // NUT criterion
      if(NUT<=0.0 && nacc_>10){
        qcurr_ = rkn_.getQ1();
        Lp1 = rkn_.getLp1();
        r_.rnorm(pcurr_); // resample momentum in case partial refreshments are used
        break;
      }
      rkn_.prepNext(true);
    } // main time stepping loop
    
    // put back error tolerances
    rkn_.setAbsTol(absTol);
    rkn_.setRelTol(relTol);
    Lpdiff_ = Lp1-Lp0;
    return(0);
  }
  
 */


  int initTrajectoryAdapt(){
    int eflag = rkn_.firstEval(qcurr_,pcurr_,true);
    if(eflag != 0){
      std::cout << "bad first eval in PDPsampler::initTrajectory" << std::endl;
      return(1);  
    }
    
    double absTol = rkn_.getAbsTol();
    double relTol = rkn_.getRelTol();
    rkn_.setAbsTol(1.0e-3);
    rkn_.setRelTol(1.0e-3);
    
    double err;
    
    bool stepGood;
    bool momentumModified;
    
    nstep_ = 0;
    nacc_ = 0;
    double Lpcurr;
    double Lpmax = rkn_.getLp0();
    double NUT; 
    
    
    
    double Lp0 = rkn_.getLp0();
    double Lp1 = Lp0;
    
    Eigen::VectorXd p_std(d_);
    
    // thresholds for renormalization of momentum (possibly move to avoid recomputing)
    double halfChiSqThresh;
    double halfChiSqMean = 0.5*static_cast<double>(d_);
    if(initMomentumRenormalization_>0){
      halfChiSqThresh = boost::math::gamma_p_inv(halfChiSqMean, 0.999);
    }
    
    while(nstep_< maxTrajSteps_){
      stepGood = false;
      while(rkn_.getEps()>1.0e-14 && nstep_< maxTrajSteps_){
        err = rkn_.step(M_.acceptsISG() && initMassAdapt_>0,M_.acceptsM2() && initMassAdapt_>0);
        intDiag_.coeffRef(nstep_,0) = rkn_.getEps();
        intDiag_.coeffRef(nstep_,1) = err;
        intDiag_.coeffRef(nstep_,2) = rkn_.getLp1();
        nstep_++;
        if(err<1.0 && err>0.0){
          stepGood = true;
          nacc_++;
          break;
        }
        rkn_.adaptEps();
      }
      if(! stepGood){
        std::cout << "integrator is not making progress, rejecting trajectory" << std::endl;
        std::cout << "nstep = " << nstep_ << std::endl;
        std::cout << "nacc = " << nacc_ << std::endl;
        rkn_.dump();
        rkn_.setAbsTol(absTol);
        rkn_.setRelTol(relTol);
        exitType_ = -1;
        return(2);
      }
      
      
      
      if(M_.acceptsISG() && initMassAdapt_>0 ){
        rkn_.getISGsample(ISG_sample_);
        M_.adaptPushISG(ISG_sample_,rkn_.getEps());
      }
      
      if(M_.acceptsM2() && initMassAdapt_>0 ){
        rkn_.getM1sample(M1_sample_);
        rkn_.getM2sample(M2_sample_);
        M_.adaptPushM2(M1_sample_,M2_sample_,rkn_.getEps());
      }
      
      
      // first stop criterion based on potential energy decrease
      
      Lpcurr = rkn_.getLp1();
      Lpmax = fmax(Lpmax,Lpcurr);
      if(Lpcurr < Lpmax - static_cast<double>(d_)){
        qcurr_ = rkn_.getQ0();
        Lp1 = rkn_.getLp0();
        r_.rnorm(pcurr_); // resample momentum in case partial refreshments are used
        break;
      }
      
      // second stop criteron based on NUT
      
      NUT = rkn_.getP1().dot(rkn_.getQ1()-qcurr_); // NUT criterion
      if(NUT<=0.0 && nacc_>10){
        qcurr_ = rkn_.getQ1();
        Lp1 = rkn_.getLp1();
        r_.rnorm(pcurr_); // resample momentum in case partial refreshments are used
        break;
      }
      
      if(nacc_>300){
        qcurr_ = rkn_.getQ1();
        Lp1 = rkn_.getLp1();
        r_.rnorm(pcurr_); // resample momentum in case partial refreshments are used
        break;
      }
      
      
      momentumModified = false;
      // mass adaptation
      if((M_.acceptsISG() || M_.acceptsM2()) && initMassAdapt_>0 && nacc_ % 30==0){
        p_std = rkn_.getP1();
        M_.inPlaceApplySqrtInvM(p_std);
        M_.adaptUpdate();
        M_.inPlaceApplySqrtM(p_std);
        momentumModified = true;
      }
      
      // reduce kinetic energy by scaling if the momentum is much larger than
      // what is expected under the Boltzmann distribution
      if(initMomentumRenormalization_ > 0 && rkn_.getEkin1() > halfChiSqThresh ){
        if(! momentumModified) p_std = rkn_.getP1();
        p_std *= sqrt(halfChiSqMean/rkn_.getEkin1());
        std::cout << "momentum renormalized with factor " << sqrt(halfChiSqMean/rkn_.getEkin1()) << std::endl;
      }
      
      
      if(momentumModified){
        eflag = rkn_.firstEval(rkn_.getQ1(),p_std,false);
      } else {
        rkn_.prepNext(true);
      }
      
      
      
    } // main time stepping loop
    
    // put back error tolerances
    rkn_.setAbsTol(absTol);
    rkn_.setRelTol(relTol);
    Lpdiff_ = Lp1-Lp0;
    return(0);
  }
  
  
  int trajectory(const bool warmup){
#ifdef _USE_R_RNGS_
    r_.seed(nevents_);
#endif
    double ltime = 0.0; //local time
    // first evaluation, including momentum resampling
    int eflag = rkn_.firstEval(qcurr_,pcurr_,true);
    if(eflag != 0){
      std::cout << "bad first eval in PDPsampler::trajectory" << std::endl;
      return(1);  
    }
    
    // store initial position configurations
    if(sampleCount_==0){
      for(int i=0;i<storePars_.size();i++) pointSamples_.coeffRef(i,0) = qcurr_.coeff(storePars_.coeff(i));
      pointSamples_.col(0).tail(dg_) = rkn_.getGen0();
      sampleCount_++;
    }
    
    // intra trajectory dependence of integrals (firstEval() above resets r and rsq to zero at ltime=0)
    r_last_sampling_.setZero();
    
    
    
    // sample exponential(1) rw
    u_ = -log(r_.runif());
    
    // stop conditions
    bool Lam_done = false;
    bool Adapt_done = ! (warmup && lambdaAdapt_) ;
    bool eventThisIter = false;
    
    // storage related to propose/accept mechanism
    double err;
    bool stepGood;
    
    // diagnostics info
    nstep_ = 0;
    nacc_ = 0;
    exitType_ = -1;
    double H0 = rkn_.getH0();
    energyErr_ = -999.0;
    
    // assorted storage
    double tinc;
    double time1;
    double tinterp;
    eventTime_ = 0.0;
    
    // adaptation_related
    double adapt_criterion;
    time_opt_ = 0.0;
    Lam_opt_ = 0.0;
    
    // main time integration loop
    while((! Lam_done || ! Adapt_done) ){
      // propose/accept mechanism
      stepGood = false;
      while(rkn_.getEps()>1.0e-14 && nstep_< maxTrajSteps_){
        err = rkn_.step(warmup && M_.acceptsISG(),warmup && M_.acceptsM2());
        
        intDiag_.coeffRef(nstep_,0) = rkn_.getEps();
        intDiag_.coeffRef(nstep_,1) = err;
        intDiag_.coeffRef(nstep_,2) = rkn_.getLp1();
        nstep_++;
        
        if(err>0.0 && err<1.0){
          stepGood = true;
          nacc_++;
          break;
        }
        rkn_.adaptEps();
      }
      if(! stepGood){
        std::cout << "integrator is not making progress, rejecting trajectory" << std::endl;
        std::cout << "nstep = " << nstep_ << std::endl;
        rkn_.dump();
        if(nstep_>=maxTrajSteps_ && warmup && rkn_.getEps()>1.0e-14){
          qnew_ = rkn_.getQ1();
          pnew_ = rkn_.getP1();
          exitType_ = -2;
          return(0);
        } else {
          exitType_ = -1;
          return(2);
        }
      }
      
      // push ISG sample
      if(warmup && M_.acceptsISG()){
        rkn_.getISGsample(ISG_sample_);
        M_.adaptPushISG(ISG_sample_,rkn_.getEps());
      }
      // push M2 sample
      if(M_.acceptsM2() && initMassAdapt_>0 ){
        rkn_.getM1sample(M1_sample_);
        rkn_.getM2sample(M2_sample_);
        M_.adaptPushM2(M1_sample_,M2_sample_,rkn_.getEps());
      }
      
      tinc = rkn_.getEps();
      
      // check sampler stop criterion
      if(! Lam_done){
        if(rkn_.getLam1()>=u_){ // if true, \Lambda>u
          // work out exact time of event where \Lambda = u
          tinc = rkn_.solveLam(u_); // event happened at ltime+tinc
          eventTime_ = ltime + tinc;
          qnew_ = rkn_.denseQ(tinc);
          pnew_ = rkn_.denseP(tinc);
          Lam_done = true;
          exitType_ = 0;
          eventThisIter = true;
        }
      }
      
      // check adaptation criterion
      if(! Adapt_done){
        adapt_criterion = rkn_.getP1().dot(rkn_.getQ1()-qcurr_); // NUT criterion
        
        if(adapt_criterion<=0.0){
          time_opt_ = ltime+rkn_.getEps();
          Lam_opt_ = rkn_.getLam1();
          Adapt_done = true;
        }
      }
      
      /*
       *  stop simulation if Lam criterion reqires too much time relative to the adaptation criterion
       *  mainly relevant in the first few iterations of the warmup phase
       */
      if( warmup && Adapt_done && (!Lam_done) && ltime + rkn_.getEps() > 4.0*time_opt_ && time_opt_>0.0){
        qnew_ = rkn_.getQ1();
        pnew_ = rkn_.getP1();
        
        eventTime_ = ltime+tinc;  
        Lam_done = true;
        exitType_ = 1;
        eventThisIter = true;
      }
      
      /*
       * Store samples of the process
       * 
       */
      
      // furthest time to be sampled
      if(Lam_done){
        time1 = time_ + eventTime_;
      } else {
        time1 = time_ + ltime + rkn_.getEps();  
      }
      
      //std::cout << nextSampleTime() << "  " << time1 << "  " << Lam_done << std::endl;
      
      while(nextSampleTime()<=time1 && sampleCount_<=nSamples_ ){
        
        tinterp = nextSampleTime() - (time_+ltime); // time_ + ltime is left mesh point
        
        // pointwise samples
        qsample_ = rkn_.denseQ(tinterp);
        for(int i=0;i<storePars_.size();i++) pointSamples_.coeffRef(i,sampleCount_) = qsample_.coeff(storePars_.coeff(i));
        pointSamples_.col(sampleCount_).tail(dg_) = rkn_.denseDR(tinterp);
        // exclude first few trajectories for mass matrix adaptation
        if(warmup && M_.acceptsSample()) M_.adaptPushSample(qsample_);
        
        // integrated samples
        rtmp_ = rkn_.denseR(tinterp);
        intSamples_.col(sampleCount_) = invSamplesSpacing_*(rtmp_-r_last_sampling_+r_last_trajectory_);
        r_last_sampling_ = rtmp_;
        r_last_trajectory_.setZero();
        
        sampleCount_++;
      }
      
      // prepare for next trajectory (runs only in the iteration where the event happened)
      if(eventThisIter){
        if(exitType_==0){
          // regular exit
          r_last_trajectory_ += (rkn_.denseR(tinc)-r_last_sampling_);
          
        } else {
          // trajectory cut short during warmup
          r_last_trajectory_ += (rkn_.getR1()-r_last_sampling_);
          
        }
        energyErr_ = rkn_.getH1()-H0;
      }
      
      // update local time
      
      ltime += rkn_.getEps();
      
      
      
      // check if we are done with the complete simulation
      if(fmin(ltime+time_,time1)>Tmax_){
        std::cout << printPrefix_ << "PDP simulation done" << std::endl;
        eventTime_ = ltime;
        Lam_done = true;
        Adapt_done = true;
        exitType_ = 2;
        energyErr_ = rkn_.getH1()-H0;
      }
      if(!Lam_done || !Adapt_done) rkn_.prepNext(true);
      eventThisIter = false;
      
    } // main time stepping loop
    
    return(0);
  }
  
  
  std::vector<std::string> propListName_;
  std::vector<int> propListDim_;
  std::vector<bool> propListWrite_;
  std::vector<bool> propListRead_;
  
  void fillPropertyList() {
    //list.push_back("1 # W # A # seed");  //0
    propListName_.push_back("seed");
    propListDim_.push_back(1);
    propListWrite_.push_back(true);
    propListRead_.push_back(false);
    //list.push_back("1 # R # A # dim");  //1
    propListName_.push_back("dim");
    propListDim_.push_back(1);
    propListWrite_.push_back(false);
    propListRead_.push_back(true);
    //list.push_back("1 # R # A # targetDim"); //2
    propListName_.push_back("targetDim");
    propListDim_.push_back(1);
    propListWrite_.push_back(false);
    propListRead_.push_back(true);
    //list.push_back("1 # R # A # targetDimGenerated"); //3
    propListName_.push_back("targetDimGenerated");
    propListDim_.push_back(1);
    propListWrite_.push_back(false);
    propListRead_.push_back(true);
    //list.push_back("1 # RW # P # absTol"); //4
    propListName_.push_back("absTol");
    propListDim_.push_back(1);
    propListWrite_.push_back(true);
    propListRead_.push_back(true);
    //list.push_back("1 # RW # P # relTol"); //5
    propListName_.push_back("relTol");
    propListDim_.push_back(1);
    propListWrite_.push_back(true);
    propListRead_.push_back(true);
    //list.push_back("3 # RW # APP # lambdaAdapt"); //6
    propListName_.push_back("lambdaAdapt");
    propListDim_.push_back(3);
    propListWrite_.push_back(true);
    propListRead_.push_back(true);
    //list.push_back("1 # RW # A # maxInitTraj"); // 7
    propListName_.push_back("maxInitTraj");
    propListDim_.push_back(1);
    propListWrite_.push_back(true);
    propListRead_.push_back(true);
    //list.push_back("3 # R # AAA # CPUtime"); // 8
    propListName_.push_back("CPUtime");
    propListDim_.push_back(3);
    propListWrite_.push_back(false);
    propListRead_.push_back(true);
    //list.push_back("1 # RW # P # maxTrajSteps"); // 9
    propListName_.push_back("maxTrajSteps");
    propListDim_.push_back(1);
    propListWrite_.push_back(true);
    propListRead_.push_back(true);
    //list.push_back("1 # R # A # massallowsFixedSubvector"); // 10
    propListName_.push_back("massallowsFixedSubvector");
    propListDim_.push_back(1);
    propListWrite_.push_back(false);
    propListRead_.push_back(true);
    //list.push_back("D # W # A # massFixedSubvector"); // 11
    propListName_.push_back("massFixedSubvector");
    propListDim_.push_back(d_);
    propListWrite_.push_back(true);
    propListRead_.push_back(false);
    //list.push_back("D # R # A # lastMiDiag"); // 12
    propListName_.push_back("lastMiDiag");
    propListDim_.push_back(d_);
    propListWrite_.push_back(false);
    propListRead_.push_back(true);
    //list.push_back("1 # RW # A # initMomentumRenormalization"); // 13
    propListName_.push_back("initMomentumRenormalization");
    propListDim_.push_back(1);
    propListWrite_.push_back(true);
    propListRead_.push_back(true);
    //list.push_back("1 # RW # A # initMassAdapt"); // 14
    propListName_.push_back("initMassAdapt");
    propListDim_.push_back(1);
    propListWrite_.push_back(true);
    propListRead_.push_back(true);
  }
  
  
  
public:
  
  // storage for samples
  Eigen::MatrixXd pointSamples_;
  Eigen::MatrixXd intSamples_;
  Eigen::MatrixXd eventSamples_;
  
  // storage for diagnostics
  Eigen::MatrixXd diagnostics_;
  // storage for CPU time
  Eigen::Matrix<double,3,1> CPUtime_;
  
  PDPsampler() : lambdaAdapt_(true),  initMassAdapt_(1), initMomentumRenormalization_(0), maxInitTraj_(30), maxTrajSteps_(20000), printPrefix_("") {}
  inline void seed(const int seed){r_.seed(seed);}
  inline int getDim() const {return d_;}
  inline int getDimGenerated() const {return dg_;}
  inline int getTargetCopies() const {return adt_.copies();}
  inline int getTargetDim() const {return adt_.dim();}
  inline int getTargetDimGenerated() const {return adt_.dimGenerated();}
  inline void setAbsTol(const double absTol){rkn_.setAbsTol(absTol);}
  inline void setRelTol(const double relTol){rkn_.setRelTol(relTol);}
  inline void setPrintPrefix(const std::string& prefix){printPrefix_="[ " + prefix + " ]\t";}
  inline void lambdaAdaptOff(const double lambdaPar0, const double lambdaPar1){
    lambdaAdapt_ = false;
    lambda_.setPar(0,lambdaPar0);
    lambda_.setPar(1,lambdaPar1);
  }
  inline void lambdaAdaptOn(){ lambdaAdapt_ = true;}
  
  
  
  
  
  void setProperty(const int whichProp, const Eigen::VectorXd theProp){
    switch(whichProp){
    case 0:
      r_.seed(static_cast<int>(theProp[0]));
      break;
    case 4:
      rkn_.setAbsTol(theProp[0]);
      break;
    case 5:
      rkn_.setRelTol(theProp[0]);
      break;
    case 6:
      lambdaAdapt_ = fabs(theProp[0])>1.0e-12;
      lambda_.setPar(0,theProp[1]);
      lambda_.setPar(1,theProp[2]);
      break;
    case 7:
      maxInitTraj_ = static_cast<int>(theProp[0]);
      break;
    case 9:
      maxTrajSteps_ = static_cast<int>(theProp[0]);
      break;
    case 11:
      M_.fixedSubvector(theProp);
      break;
    case 13:
      initMomentumRenormalization_ = static_cast<int>(theProp[0]);
      break;
    case 14:
      initMassAdapt_ = static_cast<int>(theProp[0]);
      break;
    default:
      std::cout << "Warning: property # " << whichProp <<" cannot be set " << std::endl;
    break;
    }
  }
  Eigen::VectorXd getProperty(const int whichProp) const {
    Eigen::VectorXd ret;
    switch(whichProp){
    case 1:
      ret.resize(1);
      ret(0) = static_cast<double>(d_);
      break;
    case 2:
      ret.resize(1);
      ret(0) = static_cast<double>(adt_.dim());
      break;
    case 3:
      ret.resize(1);
      ret(0) = static_cast<double>(adt_.dimGenerated());
      break;
    case 4:
      ret.resize(1);
      ret(0) = rkn_.getAbsTol();
      break;
    case 5:
      ret.resize(1);
      ret(0) = rkn_.getRelTol();
      break;
    case 6:
      ret.resize(3);
      ret(0) = static_cast<double>(lambdaAdapt_);
      ret(1) = lambda_.getPar(0);
      ret(2) = lambda_.getPar(1);
      break;
    case 7:
      ret.resize(1);
      ret(0) = static_cast<double>(maxInitTraj_);
      break;
    case 8:
      ret = CPUtime_;
      break;
    case 9:
      ret.resize(1);
      ret(0) = static_cast<double>(maxTrajSteps_);
      break;
    case 10:
      ret.resize(1);
      ret(0) = static_cast<double>(M_.allowsFixedSubvector());
      break;
    case 12:
      ret = M_.MinvDiag();
      break;
    case 13:
      ret.resize(1);
      ret(0) = static_cast<double>(initMomentumRenormalization_);
      break;
    case 14:
      ret.resize(1);
      ret(0) = static_cast<double>(initMassAdapt_);
      break;
    default:
      std::cout << "Warning: property # " << whichProp << " cannot be read! " << std::endl;
    break;
    }
    return(ret);
  }
  
  
  void setProperty(const std::string &whichProp, const Eigen::VectorXd theProp){
    int propNum;
    auto its = std::find(propListName_.begin(),propListName_.end(),whichProp);
    if(its != propListName_.end()){
      //std::cout << "found property " << whichProp << std::endl;
      propNum = std::distance(propListName_.begin(), its);
      if(propListWrite_[propNum]){
        if(theProp.size() == propListDim_[propNum]){
          setProperty(propNum,theProp);
        } else {
          std::cout << "WARNING : setProperty : wrong dimension of property " << whichProp << std::endl;
        }
      } else {
        std::cout << "WARNING : setProperty : attempting to read from write only property " << whichProp << std::endl;
      }
    } else {
      std::cout << "WARNING : setProperty : no property named " << whichProp << ", ignored" << std::endl;
    }
  }
  
  void setProperty(const std::string &whichProp, const double theProp){
    Eigen::VectorXd thePropVec(1);
    thePropVec(0) = theProp;
    setProperty(whichProp,thePropVec);
  }
  
  void setProperty(const std::string &whichProp, const int theProp){
    setProperty(whichProp,static_cast<double>(theProp));
  }
  
  
  Eigen::VectorXd getProperty(const std::string &whichProp) const {
    int propNum;
    auto its = std::find(propListName_.begin(),propListName_.end(),whichProp);
    if(its != propListName_.end()){
      //std::cout << "found property " << whichProp << std::endl;
      propNum = std::distance(propListName_.begin(), its);
      if(propListRead_[propNum]){
        return(getProperty(propNum));
      } else {
        std::cout << "WARNING : attempting to write to read only property " << whichProp << std::endl;
      }
    } else {
      std::cout << "WARNING : getProperty no property named " << whichProp << ", ignored" << std::endl;
    }
    Eigen::VectorXd ret;
    return(ret); // only returned upon error
  }
  
  
  
  void setup(){
    adt_.setup();
    d_ = adt_.dim();
    dg_ = adt_.dimGenerated();
    M_.setup(d_);
    rkn_.setup(adt_,lambda_,M_,r_);
    r_.seed(0);
    qcurr_.resize(d_);
    qnew_.resize(d_);
    pcurr_.resize(d_);
    pnew_.resize(d_);
    rtmp_.resize(dg_);
    r_last_sampling_.resize(dg_);
    r_last_trajectory_.resize(dg_);
    if(M_.acceptsISG()) ISG_sample_.resize(d_);
    if(M_.acceptsM2()){
      M1_sample_.resize(d_);
      M2_sample_.resize(d_);
    }
    fillPropertyList();
  }
  
  double target(const Eigen::VectorXd &q, Eigen::VectorXd &gen){
    if(q.size()!= d_ || gen.size() != dg_){
      std::cout << " bad argument to PDPsampler::target"  << std::endl;
      return(0.0);
    } else {
      return(adt_.feval(q,gen.col(0)));
    }
  }
  
  double gradient(const Eigen::VectorXd &q,
                  Eigen::VectorXd &gen,
                  Eigen::VectorXd &grad){
    if(q.size()!= d_ || gen.size() != dg_ || grad.size() != d_){
      std::cout << " bad argument to PDPsampler::target"  << std::endl;
      return(0.0);
    } else {
      return(adt_.geval(q,gen.col(0),grad.col(0)));
    }
  }
  
  
  int run(const int nSamples, 
          const double Tmax,
          const double warmupFrac,
          const Eigen::Matrix<int,Eigen::Dynamic,1> &storePars,
          const Eigen::Matrix<double,Eigen::Dynamic,1 > &q0){

    
    // in case multiple runs are done on the same instance 
    M_.reset();
    if(lambdaAdapt_) lambda_.reset();
    
    // sample initial config in all copies
    qcurr_ = q0;
   
    r_.rnorm(pcurr_); // in case only partial refreshes are done
    M_.inPlaceApplySqrtM(pcurr_);
    time_ = 0.0;
    
    Tmax_ = Tmax;
    
    // storePars
    storePars_ = storePars;
    
    
    sampleCount_ = 0;
    nSamples_ = nSamples;
    pointSamples_.resize(storePars_.size()+dg_,nSamples+1);
    intSamples_.resize(dg_,nSamples+1);
    intSamples_.col(0).setZero();
    samplesSpacing_ = Tmax/static_cast<double>(nSamples);
    invSamplesSpacing_ = 1.0/samplesSpacing_;
    
    diagnostics_.resize(_DIAG_BLOCK_SIZE_,_DIAG_WIDTH_);
    diagnostics_.setZero();
    
#ifdef __RCPP_DEFINED__
    Progress _progress_object(0, false); 
#endif
    
    
    intDiag_.resize(maxTrajSteps_,4);
    intDiag_.setZero();
    r_last_trajectory_.setZero();
    
    
    
    int eflag;
    nevents_ = 0;
    int samples0 = 0;
    int sampleCountBackup;
    bool warmup; 
    
    // progress
    double printCounter = 1.0;
    int oldEventsCount = 0; 
    int oldSampleCount = 0;
    
    // timing
    auto start_time = std::chrono::high_resolution_clock::now();
    bool warmupTimed = false;
    
    
    // do some intial low-fi trajectories
    if(maxInitTraj_>0) std::cout << printPrefix_ << "locating good initial point" << std::endl;
    for(int i=0;i<maxInitTraj_;i++) {
      eflag = initTrajectoryAdapt();
      std::cout << printPrefix_ << "lp-diff " << Lpdiff_ << std::endl;
      if(eflag >0 && i==0){
        std::cout << printPrefix_ << "giving up, please provide another initial q" << std::endl;
        return(1);
      }
      if(Lpdiff_<0.0 && i>4){
        std::cout << printPrefix_ << "done locating inital point" << std::endl;
        break;
      }
#ifdef __RCPP_DEFINED__
      // stop simulation if control-c is given (R-only)
      if (Progress::check_abort() ){
        std::cout << printPrefix_ << "Warning computation exited early, results may be incomplete!!!" << std::endl;
        break;
      }
#endif
    }
    
    
    
#ifdef __STORE_EVENT_Q__
    int eventSamplesCount = 1;
    eventSamples_.resize(d_,_DIAG_BLOCK_SIZE_);
    eventSamples_.col(0) = qcurr_;
#endif
    
    auto last_print = std::chrono::high_resolution_clock::now();
    
    while(time_ < Tmax){
      
      while((0.1*printCounter)*Tmax_ <= time_){
        auto tmp_time = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> since_last_print = tmp_time-last_print;
        if(since_last_print.count()>0.5){
          std::cout << printPrefix_ << "time = " << time_ << ", done with " << static_cast<int>(100.0*time_/Tmax_) << "% of run" << std::endl;
          std::cout << printPrefix_ << "events per time unit: " << static_cast<double>(nevents_-oldEventsCount)/(0.1*Tmax_) << std::endl;
          std::cout << printPrefix_ << "samples per trajectory: " << static_cast<double>(sampleCount_-oldSampleCount)/
            static_cast<double>(nevents_-oldEventsCount) << std::endl;
          last_print = tmp_time;
        }
        oldEventsCount = nevents_;
        oldSampleCount = sampleCount_;
        printCounter += 1.0;
      }
      
      // is current trajectory a warmup trajectory
      warmup = time_<warmupFrac*Tmax_;
      
      if(! warmup && !warmupTimed){
        auto sampling_start = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> warmup_elapsed = sampling_start-start_time;
        CPUtime_(0) = warmup_elapsed.count();
        if(warmupFrac>0.0) std::cout << printPrefix_ << "warmup done, required " << CPUtime_(0) << " seconds" << std::endl;
        warmupTimed = true;
      }
      
      
      // trajectory
      sampleCountBackup = sampleCount_;
      eflag = trajectory(warmup);
      
      if(eflag >0 && nevents_==0){
        std::cout << printPrefix_ << "giving up, please provide another initial q" << std::endl;
        return(1);
      }
      
      
      // store diagnostics info
      diagnostics_.coeffRef(nevents_,0) = time_; // start time
      diagnostics_.coeffRef(nevents_,1) = eventTime_; // interevent duration
      diagnostics_.coeffRef(nevents_,2) = static_cast<double>(exitType_); // exit type, positive indicates successful
      diagnostics_.coeffRef(nevents_,3) = Lam_opt_;
      diagnostics_.coeffRef(nevents_,4) = time_opt_;
      diagnostics_.coeffRef(nevents_,5) = lambda_.getPar(0);
      diagnostics_.coeffRef(nevents_,6) = lambda_.getPar(1);
      diagnostics_.coeffRef(nevents_,7) = static_cast<double>(sampleCount_-samples0); //# samples during trajectory
      diagnostics_.coeffRef(nevents_,8) = static_cast<double>(nstep_); // total number of steps
      diagnostics_.coeffRef(nevents_,9) = static_cast<double>(nacc_); // number of accepted steps
      diagnostics_.coeffRef(nevents_,10) = energyErr_;
      diagnostics_.coeffRef(nevents_,11) = static_cast<double>(warmup);
      diagnostics_.coeffRef(nevents_,12) = M_.MinvDiagMin();
      diagnostics_.coeffRef(nevents_,13) = M_.MinvDiagMax();
      
      
      
      
      
      // advance chain if the trajectory was good
      if(eflag==0){
        
        // adaptation of the lambda class
        if(warmup && lambdaAdapt_) lambda_.adapt(diagnostics_.row(nevents_));
        // adaptation of the mass matrix class
        if(warmup) {
          M_.inPlaceApplySqrtInvM(pnew_);
          M_.adaptUpdate();
          M_.inPlaceApplySqrtM(pnew_);
        }
        
        qcurr_ = qnew_;
        pcurr_ = pnew_;
        time_ += eventTime_;
        
#ifdef __STORE_EVENT_Q__
          if(eventSamplesCount>eventSamples_.cols()) eventSamples_.conservativeResize(d_,eventSamples_.cols()+_DIAG_BLOCK_SIZE_);
          eventSamples_.col(eventSamplesCount) = qcurr_;
          eventSamplesCount++;
#endif
        
      } else {
        // throw out recorded samples for failed trajectories. Note slight
        // inaccuracy in the integrated samples as these are not backed up
        sampleCount_ = sampleCountBackup;
      }
      
      
      
      
      // bookkeeping stuff
      nevents_++;
      samples0 = sampleCount_;
      
      // resize diagnostics storage if neccesary
      if(nevents_==diagnostics_.rows()){
        diagnostics_.conservativeResize(diagnostics_.rows()+_DIAG_BLOCK_SIZE_,_DIAG_WIDTH_);
      }
      
#ifdef __RCPP_DEFINED__
      // stop simulation if control-c is given (R-only)
      if (Progress::check_abort() ){
        std::cout << printPrefix_ << "Warning computation exited early, results may be incomplete!!!" << std::endl;
        break;
      }
#endif
    }
    
    diagnostics_.conservativeResize(nevents_,_DIAG_WIDTH_);
#ifdef __STORE_EVENT_Q__
    eventSamples_.conservativeResize(d_,eventSamplesCount);
#endif
    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_total = end_time-start_time;
    CPUtime_(2) = elapsed_total.count();
    CPUtime_(1) = CPUtime_(2) - CPUtime_(0);
    std::cout << printPrefix_ << "total time = " << CPUtime_(2) << " seconds " << std::endl;
    
    //std::cout << intSamples_.transpose() << std::endl;
    
    return eflag;
  }
  
  void diagnosticsToFile(const std::string filename){
    Eigen::IOFormat CSV(8, Eigen::DontAlignCols, ", ", "\n");
    std::ofstream file;
    file.open(filename);
    file << "time , eventTime , exitType , Lam_opt , time_opt , lambdaPar0 , lambdaPar1, nsamples, nSteps , nAccSteps , HErr , warmup , Minv_min, Minv_max " << std::endl;
    file << diagnostics_.topRows(nevents_).format(CSV) << std::endl;
    file.close();
  }
  
  void integratorDiagnosticsToFile(const std::string filename){
    Eigen::IOFormat CSV(8, Eigen::DontAlignCols, ", ", "\n");
    std::ofstream file;
    file.open(filename);
    file << " eps , err , empty1, empty2 " << std::endl;
    file << intDiag_.topRows(nstep_).format(CSV) << std::endl;
    file.close();
  }
  
  void pointSamplesToFile(const std::string filename){
    Eigen::IOFormat CSV(8, Eigen::DontAlignCols, ", ", "\n");
    std::ofstream file;
    file.open(filename);
    file << pointSamples_.transpose().format(CSV) << std::endl;
    file.close();
  }
  
  void samplesToFile(const int csvPrec, 
                     const bool point, //otherwise integrated samples
                     const std::vector<std::string> colhead,
                     const std::string filename){
    Eigen::IOFormat CSV(csvPrec, Eigen::DontAlignCols, ", ", "\n");
    std::ofstream file;
    file.open(filename);
    file << "\"";
    for(int c=0;c<colhead.size()-1;c++) file << colhead[c] << "\" , \"" ;
    file << colhead[colhead.size()-1] << "\"" << std::endl;
    if(point){
      file << pointSamples_.transpose().format(CSV) << std::endl;
    } else {
      file << intSamples_.transpose().format(CSV) << std::endl;
    }
    file.close();
  }
  
  /* 
   void isgSamplesToFile(const std::string filename){
   Eigen::IOFormat CSV(4, Eigen::DontAlignCols, ", ", "\n");
   std::ofstream file;
   file.open(filename);
   file << isgSamples_.transpose().format(CSV) << std::endl;
   file.close();
   }
   */
  void intSamplesToFile(const std::string filename){
    Eigen::IOFormat CSV(4, Eigen::DontAlignCols, ", ", "\n");
    std::ofstream file;
    file.open(filename);
    file << intSamples_.transpose().format(CSV) << std::endl;
    file.close();
  }
#ifdef __STORE_EVENT_Q__  
  void eventSamplesToFile(const std::string filename){
    
      Eigen::IOFormat CSV(8, Eigen::DontAlignCols, ", ", "\n");
      std::ofstream file;
      file.open(filename);
      file << eventSamples_.transpose().format(CSV) << std::endl;
      file.close();
  }
#endif
  /* LEGACY R interface, not used in current R-package
   * 
   * 
   * 
   * List of all possible fields that can be read and written to
   * using the {set,get}Property functions.
   * Each field needs two 4 pieces of information divided by a "#"
   * 1: dimension
   * 2: R=read, W=write, RW=both
   * 3: A=all values possible, P=positive, N=negative
   * 4: Name of field (available in R)
   
   std::vector<std::string> propertyList() const {
   std::vector<std::string> list;
   list.push_back("1 # W # A # seed");  //0
   list.push_back("1 # R # A # dim");  //1
   list.push_back("1 # R # A # targetDim"); //2
   list.push_back("1 # R # A # targetDimGenerated"); //3
   list.push_back("1 # RW # P # absTol"); //4
   list.push_back("1 # RW # P # relTol"); //5
   list.push_back("3 # RW # APP # lambdaAdapt"); //6
   list.push_back("1 # RW # A # maxInitTraj"); // 7
   list.push_back("3 # R # AAA # CPUtime"); // 8
   list.push_back("1 # RW # P # maxTrajSteps"); // 9
   list.push_back("1 # R # A # massallowsFixedSubvector"); // 10
   list.push_back("D # W # A # massFixedSubvector"); // 11
   list.push_back("D # R # A # lastMiDiag"); // 12
   list.push_back("1 # RW # A # initMomentumRenormalization"); // 13
   list.push_back("1 # RW # A # initMassAdapt"); // 14
   return(list);
   }
   
   */
  
  
}; // end class 


#endif



