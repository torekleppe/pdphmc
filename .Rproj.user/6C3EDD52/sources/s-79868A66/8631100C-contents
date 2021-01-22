/* 
 * example model specification, the target is independent Gaussian
 */
 

INCLUDE_BLOCK{
  
}

DATA_BLOCK{
  DATA_VECTOR(y);
}
SETUP_BLOCK{
  
}
VARIABLES_BLOCK{
  
  PARAMETER_SCALAR(mu);
  PARAMETER_SCALAR(logTau); // log of precsion
  
  GENERATED_SCALAR(sigma);
  GENERATED_SCALAR(sigmaSquared);

}
MODEL_BLOCK{
  
  // transformation of the basic parameters
  var sd = exp(-0.5*logTau); // standard deviation
  
  // notice, stan math library is available using namespace stan::math
  target = stan::math::normal_lpdf(mu,0,1000); // mu prior, (flat prior on logTau)
  target += stan::math::normal_lpdf(y,mu,sd); // likelihood function
  
  // generated quantities
  sigma = sd;
  sigmaSquared = sigma*sigma;
}

