#ifndef _MASSMATRIX_HPP_
#define _MASSMATRIX_HPP_

#define _MASS_MIN_SAMPLES_ 25

class identityMass{
private: 
  int d_;
public:
  identityMass() {}
  void setup(const int d){
    d_=d;
  }
  void reset(){}
  inline double MinvDiagMax() const { return 1.0;}
  inline double MinvDiagMin() const { return 1.0;}
  inline Eigen::VectorXd MinvDiag() const { 
    Eigen::VectorXd ret(d_);
    ret.setOnes();
    return(ret);
  }
  inline bool acceptsSample() const {return false;}
  inline bool acceptsISG() const {return false;}
  inline bool acceptsM2() const {return false;}
  inline bool allowsDiagAdaptation() const {return false;}
  inline bool allowsFixedSubvector() const {return false;}
  inline void fixedSubvector(const Eigen::VectorXd &massFixed){
    std::cout << "WARNING, function fixedSubvector should not be called for this class" << std::endl;
  }
  void acceptMinvDiag(const Eigen::VectorXd &newDiag, const double maxChange, const double minChange){
    std::cout << "WARNING, acceptMinvDiag: this function should not be called !!!" << std::endl; 
  }
  void adaptUpdate(){}
  void adaptPushSample(const Eigen::VectorXd &sample){}
  void adaptPushM2(const Eigen::VectorXd &m1, const Eigen::VectorXd &m2, const double eps){
    std::cout << "mass::adaptPushM2 : this function should not be called" << std::endl;
  }
  
  void adaptPushISG(const Eigen::VectorXd &isg, const double eps){}
  
  inline void applyMinv(const Eigen::VectorXd &input, Eigen::Ref< Eigen::VectorXd > output) const {
    output = input;
  }
  inline Eigen::VectorXd applyMinv(const Eigen::VectorXd &arg) const {
    return(arg);
  }
  inline void applyM(const Eigen::VectorXd &input, Eigen::Ref< Eigen::VectorXd > output) const {
    output = input;
  }
  inline Eigen::VectorXd applyM(const Eigen::VectorXd &arg) const {
    return(arg);
  }
  inline void inPlaceApplySqrtM(Eigen::VectorXd &inout) const {}
  inline void inPlaceApplySqrtInvM(Eigen::VectorXd &inout) const {}
  
  inline Eigen::VectorXd applySqrtM(const Eigen::VectorXd &arg) const {
    return(arg);
  }
};


class diagMassFixed{
private: 
  // diagonal of current mass matrix inverse (\approx cov(q))
  Eigen::VectorXd Minv_;
public:
  diagMassFixed() {}
  void setup(const int d){
    Minv_.resize(d);
    Minv_.setOnes();
  }
  void reset(){  }
  inline double MinvDiagMax() const { return Minv_.maxCoeff();}
  inline double MinvDiagMin() const { return Minv_.minCoeff();}
  inline Eigen::VectorXd MinvDiag() const {return Minv_;}
  inline bool acceptsSample() const {return false;}
  inline bool acceptsISG() const {return false;}
  inline bool acceptsM2() const {return false;}
  inline bool allowsDiagAdaptation() const {return false;}
  inline bool allowsFixedSubvector() const {return true;}
  inline void fixedSubvector(const Eigen::VectorXd &massFixed){
    for(int i=0;i<Minv_.size();i++){
      if(massFixed.coeff(i)>0.0) Minv_.coeffRef(i) = massFixed.coeff(i);
    }  
  }
  void acceptMinvDiag(const Eigen::VectorXd &newDiag, const double maxChange, const double minChange){
   std::cout << "WARNING, acceptMinvDiag: this function should not be called !!!" << std::endl; 
  }
  
  void adaptUpdate(){}
  void adaptPushSample(const Eigen::VectorXd &sample){}
  
  void adaptPushISG(const Eigen::VectorXd &isg, const double eps){}
  
  void adaptPushM2(const Eigen::VectorXd &m1, const Eigen::VectorXd &m2, const double eps){
    std::cout << "mass::adaptPushM2 : this function should not be called" << std::endl;
  }
  
  
  inline void applyMinv(const Eigen::VectorXd &input, Eigen::Ref< Eigen::VectorXd > output) const {
    output = input;
    output.array() *= Minv_.array();
  }
  inline Eigen::VectorXd applyMinv(const Eigen::VectorXd &arg) const {
    Eigen::VectorXd ret(arg.size());
    applyMinv(arg,ret);
    return(ret);
  }
  inline void applyM(const Eigen::VectorXd &input, Eigen::Ref< Eigen::VectorXd > output) const {
    output = input;
    output.array() /= Minv_.array();
  }
  inline Eigen::VectorXd applyM(const Eigen::VectorXd &arg) const {
    Eigen::VectorXd ret(arg.size());
    applyM(arg,ret);
    return(ret);
  }
  
  
  inline void inPlaceApplySqrtM(Eigen::VectorXd &inout) const {
    inout.array() /= Minv_.array().sqrt();
  }
  
  inline void inPlaceApplySqrtInvM(Eigen::VectorXd &inout) const {
    inout.array() *= Minv_.array().sqrt();
  }
  
  inline Eigen::VectorXd applySqrtM(const Eigen::VectorXd &arg) const {
    Eigen::VectorXd ret = arg;
    inPlaceApplySqrtM(ret);
    return(ret);
  }
};


/*
 *  Fixed mass matrix intended for two independent copies of the target
 * 
 */
class ncorrMassFixed_cp2{
private: 
  int origDim_;
  double rho_;
  Eigen::VectorXd Minv_;
  Eigen::VectorXd Miwt1_;
  Eigen::VectorXd Miwt2_;
  Eigen::VectorXd sqrtMi_;
  Eigen::VectorXd wtSqrtMi_;
  
  void updateWeights(){
    Miwt1_ = (1.0/(1.0-std::pow(rho_,2)))*Minv_;
    Miwt2_ = -(rho_/(1.0-std::pow(rho_,2)))*Minv_;
    sqrtMi_.array() = Minv_.array().sqrt();
    wtSqrtMi_ = (1.0/sqrt(1.0-pow(rho_,2)))*sqrtMi_;
  }
  
public:
  ncorrMassFixed_cp2() : rho_(-0.99) {}
  void setup(const int d){
    origDim_ = d/2;
    Minv_.resize(origDim_);
    Minv_.setOnes();
    updateWeights();
  }
  void reset(){  }
  inline double MinvDiagMax() const { return Minv_.maxCoeff();}
  inline double MinvDiagMin() const { return Minv_.minCoeff();}
  inline Eigen::VectorXd MinvDiag() const {return Minv_;}
  inline bool acceptsSample() const {return false;}
  inline bool acceptsISG() const {return false;}
  inline bool acceptsM2() const {return false;}
  inline bool allowsDiagAdaptation() const {return false;}
  inline bool allowsFixedSubvector() const {return true;}
  inline void fixedSubvector(const Eigen::VectorXd &massFixed){
    for(int i=0;i<Minv_.size();i++){
      if(massFixed.coeff(i)>0.0) Minv_.coeffRef(i) = massFixed.coeff(i);
    }
    updateWeights();
  }
  void acceptMinvDiag(const Eigen::VectorXd &newDiag, const double maxChange, const double minChange){
    std::cout << "WARNING, acceptMinvDiag: this function should not be called !!!" << std::endl; 
  }
  void adaptUpdate(){
  }
  
  void adaptPushSample(const Eigen::VectorXd &sample){
    std::cout << "WARNING, adaptPushSample : should not be called" << std::endl;
  }
  
  void adaptPushISG(const Eigen::VectorXd &isg, const double eps){
    std::cout << "WARNING, adaptPushISG : should not be called" << std::endl;
  }
  
  void adaptPushM2(const Eigen::VectorXd &m1, const Eigen::VectorXd &m2, const double eps){
    std::cout << "mass::adaptPushM2 : this function should not be called" << std::endl;
  }
  
  
  inline void applyMinv(const Eigen::VectorXd &input, Eigen::Ref< Eigen::VectorXd > output) const {
    output.head(origDim_).array() = 
      Miwt1_.array()*input.head(origDim_).array() + Miwt2_.array()*input.tail(origDim_).array();
    output.tail(origDim_).array() = 
      Miwt1_.array()*input.tail(origDim_).array() + Miwt2_.array()*input.head(origDim_).array();
  }
  inline Eigen::VectorXd applyMinv(const Eigen::VectorXd &arg) const {
    Eigen::VectorXd ret(arg.size());
    applyMinv(arg,ret);
    return(ret);
  }
  inline void applyM(const Eigen::VectorXd &input, Eigen::Ref< Eigen::VectorXd > output) const {
    output.head(origDim_).array() = 
      input.head(origDim_).array()/Minv_.array() + 
      rho_*input.tail(origDim_).array()/Minv_.array();
    output.tail(origDim_).array() =
      input.tail(origDim_).array()/Minv_.array() +
      rho_*input.head(origDim_).array()/Minv_.array();
      
  }
  inline Eigen::VectorXd applyM(const Eigen::VectorXd &arg) const {
    Eigen::VectorXd ret(arg.size());
    applyM(arg,ret);
    return(ret);
  }
  
  
  inline void inPlaceApplySqrtM(Eigen::VectorXd &inout) const {
    inout.head(origDim_).array() /= sqrtMi_.array();
    inout.tail(origDim_).array() /= wtSqrtMi_.array();
    inout.tail(origDim_) += rho_*inout.head(origDim_);
  }
  
  inline void inPlaceApplySqrtInvM(Eigen::VectorXd &inout) const {
    inout.head(origDim_).array() *= wtSqrtMi_.array();
    inout.tail(origDim_).array() *= sqrtMi_.array();
    inout.tail(origDim_) -= rho_*inout.head(origDim_);
  }
  
  inline Eigen::VectorXd applySqrtM(const Eigen::VectorXd &arg) const {
    Eigen::VectorXd ret = arg;
    inPlaceApplySqrtM(ret);
    return(ret);
  }
};




class diagMassVarS{
private: 
  // diagonal of current mass matrix inverse (\approx cov(q))
  Eigen::VectorXd Minv_;
  
  // storage related to adaptation
  Eigen::VectorXd mean_;
  Eigen::VectorXd var_;
  Eigen::VectorXd tmpVec_;
  int numPushed_;
public:
  diagMassVarS() {}
  void setup(const int d){
    Minv_.resize(d);
    Minv_.setOnes();
    mean_.resize(d);
    var_.resize(d);
    tmpVec_.resize(d);
    numPushed_ = 0;
  }
  void reset(){
    Minv_.setOnes();
    numPushed_ = 0;
  }
  inline double MinvDiagMax() const { return Minv_.maxCoeff();}
  inline double MinvDiagMin() const { return Minv_.minCoeff();}
  inline Eigen::VectorXd MinvDiag() const {return Minv_;}
  inline bool acceptsSample() const {return true;}
  inline bool acceptsISG() const {return false;}
  inline bool acceptsM2() const {return false;}
  inline bool allowsDiagAdaptation() const {return true;}
  inline bool allowsFixedSubvector() const {return false;}
  inline void fixedSubvector(const Eigen::VectorXd &massFixed){
    std::cout << "WARNING, function fixedSubvector should not be called for this class" << std::endl;
  }
  void acceptMinvDiag(const Eigen::VectorXd &newDiag, const double maxChange, const double minChange){
    for(int i=0;i<Minv_.size();i++){
      Minv_(i) = 0.05*Minv_(i)+0.95*fmin(maxChange*Minv_(i),fmax(minChange*Minv_(i),newDiag(i)));
    }
  }
  
  void adaptUpdate(){
    if(numPushed_>_MASS_MIN_SAMPLES_){
    Minv_ = var_.cwiseMax(0.01*Minv_).cwiseMin(100.0*Minv_);
    } 
  }
  void adaptPushSample(const Eigen::VectorXd &sample){
    if(numPushed_==0){
      mean_ = sample;
      var_.setZero();
    } else {
      tmpVec_ = sample-mean_;
      mean_ += _MASS_EMA_ALPHA_*tmpVec_;
      var_ = (1.0-_MASS_EMA_ALPHA_)*(var_ + _MASS_EMA_ALPHA_*(tmpVec_.array().square().matrix()));
    }
    numPushed_++;
  }
  
  void adaptPushISG(const Eigen::VectorXd &isg, const double eps){
    std::cout << "diagMassVar::adaptPushISG : this function should not be called" << std::endl;
  }
  
  void adaptPushM2(const Eigen::VectorXd &m1, const Eigen::VectorXd &m2, const double eps){
    std::cout << "mass::adaptPushM2 : this function should not be called" << std::endl;
  }
  
  
  inline void applyMinv(const Eigen::VectorXd &input, Eigen::Ref< Eigen::VectorXd > output) const {
    output = input;
    output.array() *= Minv_.array();
  }
  inline Eigen::VectorXd applyMinv(const Eigen::VectorXd &arg) const {
    Eigen::VectorXd ret(arg.size());
    applyMinv(arg,ret);
    return(ret);
  }
  inline void applyM(const Eigen::VectorXd &input, Eigen::Ref< Eigen::VectorXd > output) const {
    output = input;
    output.array() /= Minv_.array();
  }
  inline Eigen::VectorXd applyM(const Eigen::VectorXd &arg) const {
    Eigen::VectorXd ret(arg.size());
    applyM(arg,ret);
    return(ret);
  }
  
  
  inline void inPlaceApplySqrtM(Eigen::VectorXd &inout) const {
    inout.array() /= Minv_.array().sqrt();
  }
  
  inline void inPlaceApplySqrtInvM(Eigen::VectorXd &inout) const {
    inout.array() *= Minv_.array().sqrt();
  }
  
  inline Eigen::VectorXd applySqrtM(const Eigen::VectorXd &arg) const {
    Eigen::VectorXd ret = arg;
    inPlaceApplySqrtM(ret);
    return(ret);
  }
};


/*
class diagMassISG{
private: 
  // diagonal of current mass matrix inverse (\approx cov(q))
  Eigen::VectorXd Minv_;
  
  // storage related to adaptation
  Eigen::VectorXd mean_;
  Eigen::VectorXd var_;
  Eigen::VectorXd logISG_;
  Eigen::VectorXd tmpVec_;
  int numPushed_;
public:
  diagMassISG() {}
  void setup(const int d){
    Minv_.resize(d);
    Minv_.setOnes();
    mean_.resize(d);
    var_.resize(d);
    logISG_.resize(d);
    tmpVec_.resize(d);
    numPushed_ = 0;
  }
  void reset(){
    Minv_.setOnes();
    numPushed_ = 0;
  }
  inline double MinvDiagMax() const { return Minv_.maxCoeff();}
  inline double MinvDiagMin() const { return Minv_.minCoeff();}
  inline bool acceptsSample() const {return false;}
  inline bool acceptsISG() const {return true;}
  inline bool allowsDiagAdaptation() const {return true;}
  inline bool allowsFixedSubvector() const {return false;}
  inline void fixedSubvector(const Eigen::VectorXd &massFixed){
    std::cout << "WARNING, function fixedSubvector should not be called for this class" << std::endl;
  }
  void acceptMinvDiag(const Eigen::VectorXd &newDiag, const double maxChange, const double minChange){
    for(int i=0;i<Minv_.size();i++){
      Minv_(i) = 0.5*Minv_(i)+0.5*fmin(maxChange*Minv_(i),fmax(minChange*Minv_(i),newDiag(i)));
    }
  }
  
  
  void adaptUpdate(){
    if(numPushed_>_MASS_MIN_SAMPLES_){
      tmpVec_ = (-mean_ - 0.5*var_).array().exp().matrix(); //- (0.5/static_cast<double>(numPushed_-1))*var_
      Minv_ = tmpVec_.cwiseMax(0.5*Minv_).cwiseMin(2.0*Minv_);
    } 
  }
  void adaptPushSample(const Eigen::VectorXd &sample){
    std::cout << "diagMassISG::adaptPushSample : this function should not be called" << std::endl;
  }
  
  void adaptPushISG(const Eigen::VectorXd &isg, const double eps){
    if(numPushed_==0){
      mean_ = isg.array().log().matrix();
      var_.setZero();
    } else {
      logISG_ = isg.array().log().matrix();
      tmpVec_ = logISG_-mean_;
      mean_ += _MASS_EMA_ALPHA_*tmpVec_;
      var_ = (1.0-_MASS_EMA_ALPHA_)*(var_ + _MASS_EMA_ALPHA_*(tmpVec_.array().square().matrix()));
    }
    numPushed_++;
    
  }
  
  inline void applyMinv(const Eigen::VectorXd &input, Eigen::Ref< Eigen::VectorXd > output) const {
    output = input;
    output.array() *= Minv_.array();
  }
  inline Eigen::VectorXd applyMinv(const Eigen::VectorXd &arg) const {
    Eigen::VectorXd ret(arg.size());
    applyMinv(arg,ret);
    return(ret);
  }
  inline void applyM(const Eigen::VectorXd &input, Eigen::Ref< Eigen::VectorXd > output) const {
    output = input;
    output.array() /= Minv_.array();
  }
  inline Eigen::VectorXd applyM(const Eigen::VectorXd &arg) const {
    Eigen::VectorXd ret(arg.size());
    applyM(arg,ret);
    return(ret);
  }
  
  
  inline void inPlaceApplySqrtM(Eigen::VectorXd &inout) const {
    inout.array() /= Minv_.array().sqrt();
  }
  
  inline void inPlaceApplySqrtInvM(Eigen::VectorXd &inout) const {
    inout.array() *= Minv_.array().sqrt();
  }
  
  inline Eigen::VectorXd applySqrtM(const Eigen::VectorXd &arg) const {
    Eigen::VectorXd ret = arg;
    inPlaceApplySqrtM(ret);
    return(ret);
  }
};
*/

class diagMassISG{
private: 
  // diagonal of current mass matrix inverse (\approx cov(q))
  Eigen::VectorXd Minv_;
  
  // storage related to adaptation
  Eigen::VectorXd mean_;
  Eigen::VectorXd var_;
  Eigen::VectorXd logISG_;
  Eigen::VectorXd tmpVec_;
  Eigen::VectorXd fixedMass_;
  Eigen::VectorXi fixedInds_;
  int numPushed_;
public:
  diagMassISG() {}
  void setup(const int d){
    Minv_.resize(d);
    Minv_.setOnes();
    mean_.resize(d);
    var_.resize(d);
    logISG_.resize(d);
    tmpVec_.resize(d);
    numPushed_ = 0;
  }
  void reset(){
    Minv_.setOnes();
    numPushed_ = 0;
  }
  inline double MinvDiagMax() const { return Minv_.maxCoeff();}
  inline double MinvDiagMin() const { return Minv_.minCoeff();}
  inline Eigen::VectorXd MinvDiag() const {return Minv_;}
  inline bool acceptsSample() const {return false;}
  inline bool acceptsISG() const {return true;}
  inline bool acceptsM2() const {return false;}
  inline bool allowsDiagAdaptation() const {return true;}
  inline bool allowsFixedSubvector() const {return true;}
  
  inline void fixedSubvector(const Eigen::VectorXd &fixedMass){
    fixedMass_ = fixedMass;
    // count number of fixed elements
    int numFixed = (fixedMass_.array()>0.0).count();
    fixedInds_.resize(numFixed);
    
    int k = 0;
    for(int j = 0;j<fixedMass_.size();j++){
      if(fixedMass_.coeff(j)>0.0){
        fixedInds_.coeffRef(k) = j;
        k++;
      }
    }
    for(int i=0; i<fixedInds_.size();i++){
      Minv_.coeffRef(fixedInds_.coeff(i)) = fixedMass_.coeff(fixedInds_.coeff(i));
    }
  }
  
  void acceptMinvDiag(const Eigen::VectorXd &newDiag, const double maxChange, const double minChange){
    for(int i=0;i<Minv_.size();i++){
      Minv_(i) = 0.5*Minv_(i)+0.5*fmin(maxChange*Minv_(i),fmax(minChange*Minv_(i),newDiag(i)));
    }
    for(int i=0; i<fixedInds_.size();i++){
      Minv_.coeffRef(fixedInds_.coeff(i)) = fixedMass_.coeff(fixedInds_.coeff(i));
    }
  }
  
  
  void adaptUpdate(){
    if(numPushed_>_MASS_MIN_SAMPLES_){
      tmpVec_ = (-mean_ - 0.5*var_).array().exp().matrix(); //- (0.5/static_cast<double>(numPushed_-1))*var_
      Minv_ = tmpVec_.cwiseMax(0.5*Minv_).cwiseMin(2.0*Minv_);
      for(int i=0; i<fixedInds_.size();i++){
        Minv_.coeffRef(fixedInds_.coeff(i)) = fixedMass_.coeff(fixedInds_.coeff(i));
      }
    } 
  }
  
  void adaptPushSample(const Eigen::VectorXd &sample){
    std::cout << "diagMassISG::adaptPushSample : this function should not be called" << std::endl;
  }
  
  
  void adaptPushM2(const Eigen::VectorXd &m1, const Eigen::VectorXd &m2, const double eps){
    std::cout << "mass::adaptPushM2 : this function should not be called" << std::endl;
  }
  
  
  void adaptPushISG(const Eigen::VectorXd &isg, const double eps){
    if(numPushed_==0){
      mean_ = isg.array().log().matrix();
      var_.setZero();
    } else {
      logISG_ = isg.array().log().matrix();
      tmpVec_ = logISG_-mean_;
      mean_ += _MASS_EMA_ALPHA_*tmpVec_;
      var_ = (1.0-_MASS_EMA_ALPHA_)*(var_ + _MASS_EMA_ALPHA_*(tmpVec_.array().square().matrix()));
    }
    numPushed_++;
    
  }
  
  inline void applyMinv(const Eigen::VectorXd &input, Eigen::Ref< Eigen::VectorXd > output) const {
    output = input;
    output.array() *= Minv_.array();
  }
  inline Eigen::VectorXd applyMinv(const Eigen::VectorXd &arg) const {
    Eigen::VectorXd ret(arg.size());
    applyMinv(arg,ret);
    return(ret);
  }
  inline void applyM(const Eigen::VectorXd &input, Eigen::Ref< Eigen::VectorXd > output) const {
    output = input;
    output.array() /= Minv_.array();
  }
  inline Eigen::VectorXd applyM(const Eigen::VectorXd &arg) const {
    Eigen::VectorXd ret(arg.size());
    applyM(arg,ret);
    return(ret);
  }
  
  
  inline void inPlaceApplySqrtM(Eigen::VectorXd &inout) const {
    inout.array() /= Minv_.array().sqrt();
  }
  
  inline void inPlaceApplySqrtInvM(Eigen::VectorXd &inout) const {
    inout.array() *= Minv_.array().sqrt();
  }
  
  inline Eigen::VectorXd applySqrtM(const Eigen::VectorXd &arg) const {
    Eigen::VectorXd ret = arg;
    inPlaceApplySqrtM(ret);
    return(ret);
  }
};

/*
 * Mass matrix estimating precision matrix diagonal based on integrated squared gradients
 * 
 */
class diagMassISG_P{
private: 
  // diagonal of current mass matrix inverse (\approx cov(q))
  Eigen::VectorXd Minv_;
  
  // storage related to adaptation
  Eigen::VectorXd mean_;
  double T_;
//  Eigen::VectorXd var_;
//  Eigen::VectorXd logISG_;
  Eigen::VectorXd tmpVec_;
  Eigen::VectorXd fixedMass_;
  Eigen::VectorXi fixedInds_;
  int numPushed_;
public:
  diagMassISG_P() {}
  void setup(const int d){
    Minv_.resize(d);
    Minv_.setOnes();
    mean_.resize(d);
    
    numPushed_ = 0;
  }
  void reset(){
    Minv_.setOnes();
    numPushed_ = 0;
  }
  inline double MinvDiagMax() const { return Minv_.maxCoeff();}
  inline double MinvDiagMin() const { return Minv_.minCoeff();}
  inline Eigen::VectorXd MinvDiag() const {return Minv_;}
  inline bool acceptsSample() const {return false;}
  inline bool acceptsISG() const {return true;}
  inline bool acceptsM2() const {return false;}
  inline bool allowsDiagAdaptation() const {return true;}
  inline bool allowsFixedSubvector() const {return true;}
  
  inline void fixedSubvector(const Eigen::VectorXd &fixedMass){
    fixedMass_ = fixedMass;
    // count number of fixed elements
    int numFixed = (fixedMass_.array()>0.0).count();
    fixedInds_.resize(numFixed);
    
    int k = 0;
    for(int j = 0;j<fixedMass_.size();j++){
      if(fixedMass_.coeff(j)>0.0){
        fixedInds_.coeffRef(k) = j;
        k++;
      }
    }
    for(int i=0; i<fixedInds_.size();i++){
      Minv_.coeffRef(fixedInds_.coeff(i)) = fixedMass_.coeff(fixedInds_.coeff(i));
    }
  }
  
  void acceptMinvDiag(const Eigen::VectorXd &newDiag, const double maxChange, const double minChange){
    for(int i=0;i<Minv_.size();i++){
      Minv_(i) = 0.5*Minv_(i)+0.5*fmin(maxChange*Minv_(i),fmax(minChange*Minv_(i),newDiag(i)));
    }
    for(int i=0; i<fixedInds_.size();i++){
      Minv_.coeffRef(fixedInds_.coeff(i)) = fixedMass_.coeff(fixedInds_.coeff(i));
    }
  }
  
  
  void adaptUpdate(){
    if(numPushed_>_MASS_MIN_SAMPLES_){
      tmpVec_ =  mean_.cwiseInverse();
      Minv_ = tmpVec_.cwiseMax(0.5*Minv_).cwiseMin(2.0*Minv_);
      for(int i=0; i<fixedInds_.size();i++){
        Minv_.coeffRef(fixedInds_.coeff(i)) = fixedMass_.coeff(fixedInds_.coeff(i));
      }
    } 
  }
  
  void adaptPushSample(const Eigen::VectorXd &sample){
    std::cout << "diagMassISG_G::adaptPushSample : this function should not be called" << std::endl;
  }
  
  void adaptPushISG(const Eigen::VectorXd &isg, const double eps){
    if(numPushed_==0){
      mean_ = isg;
      T_ = eps;
    } else {
      T_ += eps;
      mean_ = (1.0-eps/T_)*mean_ + (eps/T_)*isg;
    }
    numPushed_++;
  }
  
  void adaptPushM2(const Eigen::VectorXd &m1, const Eigen::VectorXd &m2, const double eps){
    std::cout << "mass::adaptPushM2 : this function should not be called" << std::endl;
  }
  
  
  inline void applyMinv(const Eigen::VectorXd &input, Eigen::Ref< Eigen::VectorXd > output) const {
    output = input;
    output.array() *= Minv_.array();
  }
  inline Eigen::VectorXd applyMinv(const Eigen::VectorXd &arg) const {
    Eigen::VectorXd ret(arg.size());
    applyMinv(arg,ret);
    return(ret);
  }
  inline void applyM(const Eigen::VectorXd &input, Eigen::Ref< Eigen::VectorXd > output) const {
    output = input;
    output.array() /= Minv_.array();
  }
  inline Eigen::VectorXd applyM(const Eigen::VectorXd &arg) const {
    Eigen::VectorXd ret(arg.size());
    applyM(arg,ret);
    return(ret);
  }
  
  
  inline void inPlaceApplySqrtM(Eigen::VectorXd &inout) const {
    inout.array() /= Minv_.array().sqrt();
  }
  
  inline void inPlaceApplySqrtInvM(Eigen::VectorXd &inout) const {
    inout.array() *= Minv_.array().sqrt();
  }
  
  inline Eigen::VectorXd applySqrtM(const Eigen::VectorXd &arg) const {
    Eigen::VectorXd ret = arg;
    inPlaceApplySqrtM(ret);
    return(ret);
  }
};


class diagMassVARI{
private: 
  // diagonal of current mass matrix inverse (\approx cov(q))
  Eigen::VectorXd Minv_;
  
  // storage related to adaptation
  Eigen::VectorXd M1_;
  Eigen::VectorXd M2_;
  double T_;
  //  Eigen::VectorXd var_;
  //  Eigen::VectorXd logISG_;
  Eigen::VectorXd tmpVec_;
  Eigen::VectorXd fixedMass_;
  Eigen::VectorXi fixedInds_;
  int numPushed_;
public:
  diagMassVARI() {}
  void setup(const int d){
    Minv_.resize(d);
    Minv_.setOnes();
    M1_.resize(d);
    M2_.resize(d);
    
    numPushed_ = 0;
  }
  void reset(){
    Minv_.setOnes();
    numPushed_ = 0;
  }
  inline double MinvDiagMax() const { return Minv_.maxCoeff();}
  inline double MinvDiagMin() const { return Minv_.minCoeff();}
  inline Eigen::VectorXd MinvDiag() const {return Minv_;}
  inline bool acceptsSample() const {return false;}
  inline bool acceptsISG() const {return false;}
  inline bool acceptsM2() const {return true;}
  inline bool allowsDiagAdaptation() const {return true;}
  inline bool allowsFixedSubvector() const {return true;}
  
  inline void fixedSubvector(const Eigen::VectorXd &fixedMass){
    fixedMass_ = fixedMass;
    // count number of fixed elements
    int numFixed = (fixedMass_.array()>0.0).count();
    fixedInds_.resize(numFixed);
    
    int k = 0;
    for(int j = 0;j<fixedMass_.size();j++){
      if(fixedMass_.coeff(j)>0.0){
        fixedInds_.coeffRef(k) = j;
        k++;
      }
    }
    for(int i=0; i<fixedInds_.size();i++){
      Minv_.coeffRef(fixedInds_.coeff(i)) = fixedMass_.coeff(fixedInds_.coeff(i));
    }
  }
  
  void acceptMinvDiag(const Eigen::VectorXd &newDiag, const double maxChange, const double minChange){
    for(int i=0;i<Minv_.size();i++){
      Minv_(i) = 0.5*Minv_(i)+0.5*fmin(maxChange*Minv_(i),fmax(minChange*Minv_(i),newDiag(i)));
    }
    for(int i=0; i<fixedInds_.size();i++){
      Minv_.coeffRef(fixedInds_.coeff(i)) = fixedMass_.coeff(fixedInds_.coeff(i));
    }
  }
  
  
  void adaptUpdate(){
    if(numPushed_>_MASS_MIN_SAMPLES_){
      tmpVec_ =  M2_ - M1_.array().square().matrix();
      Minv_ = tmpVec_.cwiseMax(0.5*Minv_).cwiseMin(2.0*Minv_);
      for(int i=0; i<fixedInds_.size();i++){
        Minv_.coeffRef(fixedInds_.coeff(i)) = fixedMass_.coeff(fixedInds_.coeff(i));
      }
    } 
  }
  
  void adaptPushSample(const Eigen::VectorXd &sample){
    std::cout << "diagMassVARI::adaptPushSample : this function should not be called" << std::endl;
  }
  
  void adaptPushISG(const Eigen::VectorXd &isg, const double eps){
    std::cout << "mass::adaptPushISG : this function should not be called" << std::endl;
   /* if(numPushed_==0){
      mean_ = eps*isg;
      T_ = eps;
    } else {
      T_ += eps;
      mean_ = (1.0-eps/T_)*mean_ + (eps/T_)*isg;
    }
    numPushed_++;*/
  }
  
  void adaptPushM2(const Eigen::VectorXd &m1, const Eigen::VectorXd &m2, const double eps){
    
    double w1,w2;
    if(numPushed_==0){
      M1_ = m1;
      M2_ = m2;
      T_ = eps;
    } else {
      T_ += eps;
      w2 = eps/T_;
      w1 = 1.0-w2;
      M1_ = w1*M1_ + w2*m1;
      M2_ = w1*M2_ + w2*m2;
    }
    numPushed_++;
  }
  
  
  inline void applyMinv(const Eigen::VectorXd &input, Eigen::Ref< Eigen::VectorXd > output) const {
    output = input;
    output.array() *= Minv_.array();
  }
  inline Eigen::VectorXd applyMinv(const Eigen::VectorXd &arg) const {
    Eigen::VectorXd ret(arg.size());
    applyMinv(arg,ret);
    return(ret);
  }
  inline void applyM(const Eigen::VectorXd &input, Eigen::Ref< Eigen::VectorXd > output) const {
    output = input;
    output.array() /= Minv_.array();
  }
  inline Eigen::VectorXd applyM(const Eigen::VectorXd &arg) const {
    Eigen::VectorXd ret(arg.size());
    applyM(arg,ret);
    return(ret);
  }
  
  
  inline void inPlaceApplySqrtM(Eigen::VectorXd &inout) const {
    inout.array() /= Minv_.array().sqrt();
  }
  
  inline void inPlaceApplySqrtInvM(Eigen::VectorXd &inout) const {
    inout.array() *= Minv_.array().sqrt();
  }
  
  inline Eigen::VectorXd applySqrtM(const Eigen::VectorXd &arg) const {
    Eigen::VectorXd ret = arg;
    inPlaceApplySqrtM(ret);
    return(ret);
  }
};









/*
class diagMassISGW{
private: 
  double stepWeight_;
  double maxMi_;
  double minMi_;
  // diagonal of current mass matrix inverse (\approx cov(q))
  Eigen::VectorXd Minv_;
  
  // storage related to adaptation
  Eigen::VectorXd mean_;
  //Eigen::VectorXd var_;
  Eigen::VectorXd ISG_;
  Eigen::VectorXd tmpVec_;
  
  Eigen::VectorXd ISGsum_;
  Eigen::VectorXd ISGWsum_;
  double epsSum_;
  
  int numPushed_;
  int numUpdated_;
  
public:
  diagMassISGW() : stepWeight_(1.0), maxMi_(1000.0), minMi_(1.0/1000.0) {}
  void setup(const int d){
    Minv_.resize(d);
    Minv_.setOnes();
    mean_.resize(d);
    //var_.resize(d);
    ISG_.resize(d);
    tmpVec_.resize(d);
    ISGsum_.resize(d);
    ISGWsum_.resize(d);
    ISGsum_.setZero();
    ISGWsum_.setZero();
    epsSum_ = 0.0;
    
    numPushed_ = 0;
    numUpdated_ = 0;
  }
  void reset(){
    Minv_.setOnes();
    numPushed_ = 0;
    numUpdated_ = 0;
  }
  inline double MinvDiagMax() const { return Minv_.maxCoeff();}
  inline double MinvDiagMin() const { return Minv_.minCoeff();}
  inline bool acceptsSample() const {return false;}
  inline bool acceptsISG() const {return true;}
  inline bool allowsDiagAdaptation() const {return true;}
  inline bool allowsFixedSubvector() const {return false;}
  inline void fixedSubvector(const Eigen::VectorXd &massFixed){
    std::cout << "WARNING, function fixedSubvector should not be called for this class" << std::endl;
  }
  void acceptMinvDiag(const Eigen::VectorXd &newDiag, const double maxChange, const double minChange){
    for(int i=0;i<Minv_.size();i++){
      Minv_(i) = 0.5*Minv_(i)+0.5*fmin(maxChange*Minv_(i),fmax(minChange*Minv_(i),newDiag(i)));
    }
  }
  
  
  void adaptUpdate(){
    std::cout << "numPushed " << numPushed_ << std::endl; 
    if(numPushed_>_MASS_MIN_SAMPLES_){
      
      if(numUpdated_==0){
        mean_ = ((((1.0/epsSum_)*ISGsum_).array().pow(1.0-stepWeight_))*
          ((1.0/static_cast<double>(numPushed_))*ISGWsum_.array().pow(stepWeight_))).matrix();//.array().log().matrix();
        mean_ = mean_.cwiseMax(minMi_).cwiseMin(maxMi_);
        //var_.setZero();
      } else {
        ISG_ = ((((1.0/epsSum_)*ISGsum_).array().pow(1.0-stepWeight_))*
          ((1.0/static_cast<double>(numPushed_))*ISGWsum_.array().pow(stepWeight_))).matrix();
        ISG_ = ISG_.cwiseMax(minMi_).cwiseMin(maxMi_);//.array().log().matrix();
        tmpVec_ = ISG_-mean_;
        mean_ += _MASS_EMA_ALPHA_*tmpVec_;
        //var_ = (1.0-_MASS_EMA_ALPHA_)*(var_ + _MASS_EMA_ALPHA_*(tmpVec_.array().square().matrix()));
        tmpVec_ = mean_.cwiseInverse(); // (-mean_ - 0.5*var_).array().exp().matrix(); //- (0.5/static_cast<double>(numPushed_-1))*var_
        Minv_ = tmpVec_.cwiseMax(0.5*Minv_).cwiseMin(2.0*Minv_);
      }
      numUpdated_++;
      
      epsSum_ = 0.0;
      ISGsum_.setZero();
      ISGWsum_.setZero();
      numPushed_ = 0;
    }
    
    std::cout << "Minv:\n" << Minv_.maxCoeff() << "  " << Minv_.minCoeff() << std::endl;
    std::cout << "ISG:\n" << ISG_.maxCoeff() << "  " << ISG_.minCoeff() << std::endl;
  }
  
  void adaptPushSample(const Eigen::VectorXd &sample){
    std::cout << "diagMassISG::adaptPushSample : this function should not be called" << std::endl;
  }
  
  void adaptPushISG(const Eigen::VectorXd &isg, const double eps){
    ISGWsum_ += isg;
    ISGsum_ += eps*isg;
    epsSum_ += eps;
    numPushed_++;
  }
  
  inline void applyMinv(const Eigen::VectorXd &input, Eigen::Ref< Eigen::VectorXd > output) const {
    output = input;
    output.array() *= Minv_.array();
  }
  inline Eigen::VectorXd applyMinv(const Eigen::VectorXd &arg) const {
    Eigen::VectorXd ret(arg.size());
    applyMinv(arg,ret);
    return(ret);
  }
  inline void applyM(const Eigen::VectorXd &input, Eigen::Ref< Eigen::VectorXd > output) const {
    output = input;
    output.array() /= Minv_.array();
  }
  inline Eigen::VectorXd applyM(const Eigen::VectorXd &arg) const {
    Eigen::VectorXd ret(arg.size());
    applyM(arg,ret);
    return(ret);
  }
  
  
  inline void inPlaceApplySqrtM(Eigen::VectorXd &inout) const {
    inout.array() /= Minv_.array().sqrt();
  }
  
  inline void inPlaceApplySqrtInvM(Eigen::VectorXd &inout) const {
    inout.array() *= Minv_.array().sqrt();
  }
  
  inline Eigen::VectorXd applySqrtM(const Eigen::VectorXd &arg) const {
    Eigen::VectorXd ret = arg;
    inPlaceApplySqrtM(ret);
    return(ret);
  }
};


*/


#endif



