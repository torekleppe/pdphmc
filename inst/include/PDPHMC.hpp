#ifndef __PDPHMC_HPP__
#define __PDPHMC_HPP__

#include <cmath>
#include <iostream>
#include <fstream>
#include <string>

#include <stan/math.hpp>
#include <Eigen/Dense>

#include <algorithm>
#include <vector>
#include <random>
#include <chrono>
#include <iterator>


/*
 * Used with either var=double or var=stan::math::var
 */
#define VectorXv Eigen::Matrix<var,Eigen::Dynamic,1> 
#define MatrixXv Eigen::Matrix<var,Eigen::Dynamic,Eigen::Dynamic>

/*
 * Vector of stan::math::var
 */

#define VectorXad Eigen::Matrix<stan::math::var,Eigen::Dynamic,1> 

/*
 * overloads to get value from ADtypes
 */
inline double doubleValue(const double var){return var;}
inline double doubleValue(const stan::math::var var){return var.val();}
inline Eigen::VectorXd doubleValue(const Eigen::VectorXd &var){return var;}
inline Eigen::VectorXd doubleValue(const Eigen::Matrix<stan::math::var,Eigen::Dynamic,1> &var){
  Eigen::VectorXd ret(var.size());
  for(int i=0;i<var.size();i++) ret.coeffRef(i) = var.coeff(i).val();
  return ret;
}
inline Eigen::MatrixXd doubleValue(const Eigen::MatrixXd &var){return var;}
inline Eigen::MatrixXd doubleValue(const Eigen::Matrix<stan::math::var,Eigen::Dynamic,Eigen::Dynamic> &var){
  Eigen::MatrixXd ret(var.rows(),var.cols());
  for(int j=0;j<var.cols();j++) for(int i=0;i<var.rows();i++) ret.coeffRef(i,j) = var.coeff(i,j).val();
  return ret;
}
  

/*
 * 
 * Macros used for passing data from R
 * 
 */

#include "codeTransform.hpp"

/*#define VectorXv Eigen::Matrix<var,Eigen::Dynamic,1> 
#define MatrixXv Eigen::Matrix<var,Eigen::Dynamic,Eigen::Dynamic>

#define DATA_SCALAR(x)  double x 
#define DATA_VECTOR(x)  Eigen::VectorXd x 
#define DATA_MATRIX(x)  Eigen::MatrixXd x 
#define DATA_ISCALAR(x)  static int x
#define DATA_IVECTOR(x)  Eigen::VectorXi x 
#define DATA_IMATRIX(x)  Eigen::MatrixXi x
*/
/*
#define PARAMETER_SCALAR(x) var x
#define PARAMETER_VECTOR(x,d1) VectorXv x(d1)
#define PARAMETER_MATRIX(x,d1,d2) MatrixXv x(d1,d2)

#define GENERATED_SCALAR(x) var x
#define GENERATED_VECTOR(x,d1) VectorXv x(d1)
#define GENERATED_MATRIX(x,d1,d2) MatrixXv x(d1,d2)
*/


/*
 *  constants
 */
//#define _MAX_STEPS_PER_TRAJECTORY_ 20000
#define _DIAG_BLOCK_SIZE_ 10000
#define _DIAG_WIDTH_ 14
#define _LAMBDA_EMA_ALPHA_ 0.02
#define _MASS_EMA_ALPHA_ 0.01
#define _STORE_INTEGRATED_SQUARED_GRADIENT_ 1 // set to zero except when debugging isg mechanism

/*
 * Error and warning handling
 */
#include "errorLog.hpp"
PDPHMC_errorLog::errorLog _PDPHMC_Err_("error");
PDPHMC_errorLog::errorLog _PDPHMC_Warn_("warning");
#define _ERROR_LOG_ _PDPHMC_Err_
#define _WARNING_LOG_ _PDPHMC_Warn_

#include "rng.hpp"


#ifdef __RCPP_DEFINED__ 
#include <progress.hpp>  // ensure that computations can be stopped using command-c if running under Rcpp
#endif

#include "ADtarget.hpp"
//#include "lbfgs.hpp"
#include "extRKN.hpp"
#include "lambda.hpp"
#include "massmatrix.hpp"
#include "PDPsampler.hpp"
#include "fast_spec_funs.hpp"
#include "transformed_prior.hpp"
#include "CIP_Wishart.hpp"
#include "tridiagChol.hpp"
#include "stochastic_gradient.hpp"


#include "json_wrap.hpp"


#endif


