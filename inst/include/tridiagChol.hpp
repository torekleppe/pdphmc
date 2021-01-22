#include <iostream>
#include <Eigen/Dense>
/*
 * This header provides Cholesky algorithms, including triangular solvers for symmetric
 *  *TRIDIAGONAL* matrices
 */

namespace tridiag_Chol{
using namespace stan::math;
template <class var>
class tridiagChol{
private:
  Eigen::Matrix<var,Eigen::Dynamic,1> L_;
  int n_;
public:
  /*
   *  tridiagonal n x n matrix A with diag(A) = [d1n, d2nm1, d2nm1, ... , d2nm1, d1n]
   *  and each of the first off-diagonals equal to od
   */
  tridiagChol(const int n, const var d1n, const var d2nm1, const var od) : n_(n) {
    L_.resize(2*n-1);
    var tmp;
    if(doubleValue(d1n)<0.0){
      std::cout << "Matrix in symTridiagChol not SPD, iteration 1, exiting" << std::endl;
      throw 0;
    }
    L_.coeffRef(0) = sqrt(d1n);
    for(int t=2; t<n;t++){
      L_.coeffRef(n+t-2) = od/L_.coeff(t-2);
      tmp = d2nm1-square(L_.coeff(n+t-2));
      if(doubleValue(tmp)<0.0){
        std::cout << "Matrix in symTridiagChol not SPD, iteration " << t << ", exiting" << std::endl;
        throw 0;
        
      }
      L_.coeffRef(t-1) = sqrt(tmp);
    }
    L_.coeffRef(2*(n-1)) = od/L_.coeff(n-2);
    tmp = d1n - square(L_.coeff(2*(n-1)));
    if(doubleValue(tmp)<0.0){
      std::cout << "Matrix in symTridiagChol not SPD, iteration n, exiting" << std::endl;
      throw 0;
    }
    L_.coeffRef(n-1) = sqrt(tmp);
  }
  
  /*
   * tridiagonal n x n matrix A with diag(A) = diag
   *  and each of the first off-diagonals equal to od
   */
  tridiagChol(const Eigen::Matrix<var,Eigen::Dynamic,1> &diag, const var od){
    n_ = diag.size();
    L_.resize(2*n_-1);
    var tmp;
    if(doubleValue(diag.coeff(0))<0.0){
      std::cout << "Matrix in symTridiagChol not SPD, iteration 1, exiting" << std::endl;
      return;
    }
    L_.coeffRef(0) = sqrt(diag.coeff(0));
    for(int t=2;t<=n_;t++){
      L_.coeffRef(n_+t-2) = od/L_.coeff(t-2);
      tmp = diag.coeff(t-1) - square(L_.coeff(n_+t-2));
      if(doubleValue(tmp)<0.0){
        std::cout << "Matrix in symTridiagChol not SPD, iteration " << t << ", exiting" << std::endl;
        return;
      }
      L_.coeffRef(t-1) = sqrt(tmp);
    }
  }
  
  /*
   * tridiagonal n x n matrix A with diag(A) = diag (length = n)
   *  and the first off-diagonals equal to od (length = n-1)
   * 
   */
  tridiagChol(const Eigen::Matrix<var,Eigen::Dynamic,1> &diag,
                 const Eigen::Matrix<var,Eigen::Dynamic,1> &od){
    n_ = diag.size();
    L_.resize(2*n_-1);
    var tmp;
    if(od.size() != n_-1){
      std::cout << "Dimension mismatch in symTridiagChol, exiting" << std::endl;
      return;
    }
    if(doubleValue(diag.coeff(0))<0.0){
      std::cout << "Matrix in symTridiagChol not SPD, iteration 1, exiting" << std::endl;
      return;
    }
    L_.coeffRef(0) = sqrt(diag.coeff(0));
    for(int t=2;t<=n_;t++){
      L_.coeffRef(n_+t-2) = od.coeff(t-2)/L_.coeff(t-2);
      tmp = diag.coeff(t-1) - square(L_.coeff(n_+t-2));
      if(doubleValue(tmp)<0.0){
        std::cout << "Matrix in symTridiagChol not SPD, iteration " << t << ", exiting" << std::endl;
        return;
      }
      L_.coeffRef(t-1) = sqrt(tmp);
    }
  }
  
  
  
  /*
   * Solves L x = b for x
   * 
   */
  template <class T__>
  void L_solve(const Eigen::Matrix<T__,Eigen::Dynamic,1> &b,
               Eigen::Matrix<typename boost::math::tools::promote_args<var,T__>::type,Eigen::Dynamic,1> &x) const {
    if(x.size() != n_ || b.size() != n_){
      std::cout << "Dimension mismatch in L_solve" << std::endl;
      return;
    }
    x.coeffRef(0) = b.coeff(0)/L_.coeff(0);
    for(int i=1;i<n_;i++) x.coeffRef(i) = (b.coeff(i) - x.coeff(i-1)*L_.coeff(n_+i-1))/L_.coeff(i);
  }
  template <class T__>
  Eigen::Matrix<typename boost::math::tools::promote_args<var,T__>::type,Eigen::Dynamic,1>
  L_solve(const Eigen::Matrix<T__,Eigen::Dynamic,1> &b) const {
    Eigen::Matrix<typename boost::math::tools::promote_args<var,T__>::type,Eigen::Dynamic,1> ret(n_);
    L_solve(b,ret);
    return(ret);
  }
  
  /*
   * Solves L^T x = b for x 
   * 
   */
  template <class T__>
  void LT_solve(const Eigen::Matrix<T__,Eigen::Dynamic,1> &b,
                Eigen::Matrix<typename boost::math::tools::promote_args<var,T__>::type,Eigen::Dynamic,1> &x) const {
    if(x.size() != n_ || b.size() != n_){
      std::cout << "Dimension mismatch in LT_solve" << std::endl;
      return;
    }  
    x.coeffRef(n_-1) = b.coeff(n_-1)/L_.coeff(n_-1);
    for(int i=n_-2; i>=0;i--) x.coeffRef(i) = (b.coeff(i) - x.coeff(i+1)*L_.coeff(n_+i))/L_.coeff(i);
  }
  
  template <class T__>
  void LT_solve(const Eigen::Matrix<T__,Eigen::Dynamic,1> &b,
                Eigen::Ref<Eigen::Matrix<typename boost::math::tools::promote_args<var,T__>::type,Eigen::Dynamic,1> > x) const {
    if(x.size() != n_ || b.size() != n_){
      std::cout << "Dimension mismatch in LT_solve" << std::endl;
      return;
    }  
    x.coeffRef(n_-1) = b.coeff(n_-1)/L_.coeff(n_-1);
    for(int i=n_-2; i>=0;i--) x.coeffRef(i) = (b.coeff(i) - x.coeff(i+1)*L_.coeff(n_+i))/L_.coeff(i);
  }
  
  
  template <class T__>
  Eigen::Matrix<typename boost::math::tools::promote_args<var,T__>::type,Eigen::Dynamic,1>
  LT_solve(const Eigen::Matrix<T__,Eigen::Dynamic,1> &b) const {
    Eigen::Matrix<typename boost::math::tools::promote_args<var,T__>::type,Eigen::Dynamic,1> ret(n_);
    LT_solve(b,ret);
    return(ret);
  }
  
  template <class T__>
  Eigen::Matrix<typename boost::math::tools::promote_args<var,T__>::type,Eigen::Dynamic,1>
  LT_solve(const Eigen::Ref< Eigen::Matrix<T__,Eigen::Dynamic,1> > b) const {
    Eigen::Matrix<typename boost::math::tools::promote_args<var,T__>::type,Eigen::Dynamic,1> ret(n_);
    LT_solve(b,ret);
    return(ret);
  }
  
  
  /*
   * Solves A x = b for x where A = LL^T
   * 
   */
  template <class T__>
  void solve(const Eigen::Matrix<T__,Eigen::Dynamic,1> &b,
                     Eigen::Matrix<typename boost::math::tools::promote_args<var,T__>::type,Eigen::Dynamic,1> &x) const {
    if(x.size() != n_ || b.size() != n_){
      std::cout << "Dimension mismatch in symTridiagChol::solve" << std::endl;
      return;
    }  
    LT_solve(L_solve(b),x);
  }
  template <class T__>
  Eigen::Matrix<typename boost::math::tools::promote_args<var,T__>::type,Eigen::Dynamic,1>
  solve(const Eigen::Matrix<T__,Eigen::Dynamic,1> &b) const {
    Eigen::Matrix<typename boost::math::tools::promote_args<var,T__>::type,Eigen::Dynamic,1> ret(n_);
    solve(b,ret);
    return(ret);
  }
  
  /*
   * Log-determinant of Cholesky factor L
   */
  var logDetL() const {
    var ret = 0.0;
    for(int i=0;i<n_;i++) ret+=log(L_.coeff(i));
    return(ret);
  }
  
  /*
   * Log-determinat of orignal matrix A=LL^T
   */
  var logDetA() const {return 2.0*logDetL();}
  
  /*
   * Prints Cholesky factor to std::cout, 
   * For debugging only, do not run for large matrices!
   */
  void dump(){
    Eigen::MatrixXd dense(n_,n_);
    dense.setZero();
    for(int i=0;i<n_;i++) dense(i,i) = doubleValue(L_(i));
    for(int i=0;i<n_-1;i++) dense(i+1,i) = doubleValue(L_(i+n_));
    std::cout << "Cholesky factor of " << n_ << " x " << n_ << " matrix: " << std::endl;
    std::cout << dense << std::endl;
  }
  
  
  
  
};


} // namespace



