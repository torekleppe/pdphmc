/* 
 * example model specification, the target is independent Gaussian
 */
 

INCLUDE_BLOCK{
  
}

DATA_BLOCK{
  
}
SETUP_BLOCK{
  
}
VARIABLES_BLOCK{
  
  PARAMETER_VECTOR(x,4);
  PARAMETER_SCALAR(y);
  PARAMETER_MATRIX(z,2,3);
  
  GENERATED_SCALAR(x12);
  GENERATED_SCALAR(chi2);
  GENERATED_VECTOR(gx,4);
  
  
}
MODEL_BLOCK{
  
  target = -0.5*x.dot(x) - 0.5*(y-10.0)*(y-10.0)*100.0;
  target -= 0.5*z.row(0).dot(z.row(0)) + 0.5*z.row(1).dot(z.row(1));
  
  gx = doubleValue(x);
  x12 = doubleValue(x(0));
  chi2 = doubleValue(target);
  
}

