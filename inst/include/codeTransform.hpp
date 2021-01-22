#ifndef __CODETRANSFORM_HPP__
#define __CODETRANSFORM_HPP__



#include <string>
#include <vector>
#include <iostream>

#define DATA_SCALAR(varname) double varname                                            
#define DATA_VECTOR(varname) Eigen::VectorXd varname
#define DATA_MATRIX(varname) Eigen::MatrixXd varname

#define DATA_ISCALAR(varname) int varname                                            
#define DATA_IVECTOR(varname) Eigen::VectorXi varname
#define DATA_IMATRIX(varname) Eigen::MatrixXi varname

#define PARAMETER_SCALAR(x) var x
#define PARAMETER_VECTOR(x,d1) VectorXv x(d1)
#define PARAMETER_MATRIX(x,d1,d2) MatrixXv x(d1,d2)
    
#define GENERATED_SCALAR(x) var x
#define GENERATED_VECTOR(x,d1) VectorXv x(d1)
#define GENERATED_MATRIX(x,d1,d2) MatrixXv x(d1,d2)

    
#endif
