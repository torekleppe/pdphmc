#ifndef __identityPlusRank1Chol_HPP__
#define __identityPlusRank1Chol_HPP__
#include <Eigen/Dense>

/*
 * This header provides O(d) routines for Cholesky factorizations L*L of matrices
 * on the form
 * 
 *  L*L^T = I + w*w^T = A
 * 
 * where w is a (column) vector
 * Notice that the factorization is not calculated explicitly, which gives rise to fast algorithms
 * 
 * In the algorithms below (using a matlab-like notation)
 *  S(0) = 1.0
 *  S(i) = 1.0 + sum(w(1:i).^2) = S(i-1)+w(i)^2, i=1,...,n
 *  
 *  D(i) = sqrt(S(i)/S(i-1)), i=1,...,n
 *  
 *  Then the cholesky factor has the form
 *  
 *  L(i,i) = D(i), i=1,...,n
 *  L(i,j) = w(i)*w(j)/(S(j-1)*D(j))
 *  
 *  This notation is used throughout the algorithms below
 */




namespace identityPlusRank1Chol{

// solve L^T*x = y for x
void LTsolve(const Eigen::Matrix<double,Eigen::Dynamic,1> &w,
             const Eigen::Matrix<double,Eigen::Dynamic,1> &y,
             Eigen::Matrix<double,Eigen::Dynamic,1> &x){
  int d = w.size();
  Eigen::Matrix<double,Eigen::Dynamic,1> S(d+1);
  Eigen::Matrix<double,Eigen::Dynamic,1> D(d);
  
  S.coeffRef(0) = 1.0;
  for(int i=0;i<d;i++){
    S.coeffRef(i+1) = S.coeff(i) + pow(w.coeff(i),2);
    D.coeffRef(i) = sqrt(S.coeff(i+1)/S.coeff(i));
  }
  
  double R = 0.0; // R_i = \sum_\{j=i+1}^d r_i, r_i =w_i*y_i/(d_i*s_{i-1})
  double tmp;
  
  for(int i=d-1;i>=0;i--){
    tmp = y.coeff(i)/D.coeff(i);
    x.coeffRef(i) = tmp - w.coeff(i)*R;
    R += tmp*w.coeff(i)/S.coeff(i);
  }
}

Eigen::Matrix<double,Eigen::Dynamic,1> LTsolve(const Eigen::Matrix<double,Eigen::Dynamic,1> &w,
                                               const Eigen::Matrix<double,Eigen::Dynamic,1> &y){
  Eigen::Matrix<double,Eigen::Dynamic,1> x(w.size());
  LTsolve(w,y,x);
  return(x);
}

// solves L*x = y for x
void Lsolve(const Eigen::Matrix<double,Eigen::Dynamic,1> &w,
            const Eigen::Matrix<double,Eigen::Dynamic,1> &y,
            Eigen::Matrix<double,Eigen::Dynamic,1> &x){
  
  int d = w.size();
  double S;
  double Sold = 1.0;
  double D;
  double R = 0.0;
  
  for(int i=0;i<d;i++){
    S = Sold + pow(w.coeff(i),2);
    D = sqrt(S/Sold);
    x.coeffRef(i) = (y.coeff(i)-w.coeff(i)*R/Sold)/D;
    R += w.coeff(i)*y.coeff(i);
    Sold = S;
  }
}
Eigen::Matrix<double,Eigen::Dynamic,1> Lsolve(const Eigen::Matrix<double,Eigen::Dynamic,1> &w,
                                              const Eigen::Matrix<double,Eigen::Dynamic,1> &y){
  Eigen::Matrix<double,Eigen::Dynamic,1> x(w.size());
  Lsolve(w,y,x);
  return(x);
}
// solves A*x = L*L^T*x = y for x using the Sherman-Morrison formula
void solve(const Eigen::Matrix<double,Eigen::Dynamic,1> &w,
           const Eigen::Matrix<double,Eigen::Dynamic,1> &y,
           Eigen::Matrix<double,Eigen::Dynamic,1> &x){
  x = y - (w.dot(y)/(1.0+w.squaredNorm()))*w;
}
Eigen::Matrix<double,Eigen::Dynamic,1> solve(const Eigen::Matrix<double,Eigen::Dynamic,1> &w,
                                             const Eigen::Matrix<double,Eigen::Dynamic,1> &y){
  Eigen::Matrix<double,Eigen::Dynamic,1> x(w.size());
  solve(w,y,x);
  return(x);
}

} // namespace



#endif

