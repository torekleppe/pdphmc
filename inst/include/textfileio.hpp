

#ifndef __TEXTFILEIO_HPP__
#define __TEXTFILEIO_HPP__

#include <iostream>
#include <fstream>
#include <string>
#include <Eigen/Dense>

/*
 * write corresponding file in R using the function
 * 
 */

Eigen::MatrixXd readVecMatrix(const std::string filename){
  int rows = 0;
  int cols = 0;
  Eigen::MatrixXd ret(0,0);
  std::ifstream file;
  std::string line;
  file.open(filename);
  if(file.is_open()){
    getline(file,line);
    rows = std::round(std::stod(line));
    getline(file,line);
    cols = std::round(std::stod(line));
    ret.resize(rows,cols);
    
    for(int j=0;j<cols;j++){
      for(int i=0;i<rows;i++){
        if(! file.eof()){
          getline(file,line);
          ret.coeffRef(i,j) = std::stod(line);
        } else{
          std::cout << "file ended before matrix was completed" << std::endl;
          ret.resize(0,0);
          return(ret);
        }  
      }
    }
    return(ret);
  } else {
    std::cout << "failed to open file" << std::endl;
  }
  return(ret);
}





#endif

