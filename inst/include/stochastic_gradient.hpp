#ifndef _stochastic_gradient_hpp_
#define _stochastic_gradient_hpp_

#include <cmath>
#include <iostream>
#include <Eigen/Dense>
#include "rng.hpp"


/*
 * Stochastic gradients wrapper class for the normal linear model 
 * y = x*beta + epsilon, where epsilon \sim iid N(0,sigma^2)
 * Note, the class only have pointers to y and x, which must be stored
 * elsewhere.
 */
class sg_lm{
  
  rng r_;
  
  const Eigen::VectorXd *y_;
  const Eigen::MatrixXd *x_;
  int n_; // total number of observations
  int p_; // number of covariate columns
  int m_; // minibatch size
  Eigen::Matrix<std::size_t,Eigen::Dynamic,1> inds_;
  
  double scaleFac_;
  double dm_;
  
  void sampleInds(){
    if(inds_.size()!=m_) inds_.resize(m_);
    for(int i=0;i<m_;i++){
      // improve this for large n (not storeable exactly as double)
      inds_.coeffRef(i) = static_cast<std::size_t>(std::floor(r_.runif()*static_cast<double>(n_)));
    }
  }
  
  void getDims(){
    n_ = (*y_).size();
    p_ = (*x_).cols();
    if(n_ != (*x_).rows()){
      std::cout << "sg_lm : Warning, dimension of y is not equal to number of rows in x" << std::endl;
    }
    dm_ = static_cast<double>(m_);
    scaleFac_ = static_cast<double>(n_)/dm_;
  }
  
  void checkDims(){
    if(n_!=(*y_).size() || p_ != (*x_).cols()) getDims();
  }
  
public:
  sg_lm(const Eigen::VectorXd &y,
        const Eigen::MatrixXd &x) : y_(&y),x_(&x),m_(30) {
    
    checkDims();
    r_.seed(1);
    sampleInds();
  }
  
  void resampleMiniBatch(){
    checkDims();
    sampleInds();
  }
  void miniBatchSize(const int m){
    m_ = m;
    getDims();
    sampleInds();
  }
  void seed(const int seed){
    r_.seed(seed);
    checkDims();
    sampleInds();
  }
  
  template <class T>
  T lpdf(const Eigen::Matrix<T,Eigen::Dynamic,1> beta,
         const T sigma) {
    checkDims(); // be careful here as data, and hence dimensions of y and x is added quite late
    T eta;
    T devs = 0.0;
    T sigmaSq = square(sigma);
    for(int i=0;i<m_;i++){
      eta = (*x_).row(inds_.coeff(i)).dot(beta);
      devs += square((*y_).coeff(inds_.coeff(i))-eta);
    }
    return(scaleFac_*(-0.5*devs/sigmaSq - 
           0.5*dm_*log(6.2831853071796*sigmaSq)));
  }
}; // end class sg_lm

#endif