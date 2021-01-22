INCLUDE_BLOCK{
  using namespace transformed_prior;
  using namespace tridiag_Chol;
  using namespace stan::math;
}

DATA_BLOCK{
    
  DATA_MATRIX(yy);
  DATA_MATRIX(yinv);
  DATA_SCALAR(ldets);
    
    
  double s0;
  double k0;
  double dk;
  double dT;
  int k;
  int hsize;
  int T;
}
SETUP_BLOCK{
    s0 = 0.25;
    k0 = 4.0;
    k = yy.cols();
    dk = static_cast<double>(k);
    hsize = std::round(0.5*static_cast<double>(k*(k-1)));
    dT = static_cast<double>(yy.rows())/dk;
    T = std::round(dT);
}
VARIABLES_BLOCK{
	
	PARAMETER_VECTOR(lphi,k);
	PARAMETER_VECTOR(mu,k);
	PARAMETER_VECTOR(lsigmaSq,k);
  PARAMETER_VECTOR(h,hsize);
  PARAMETER_SCALAR(lnu);
  PARAMETER_VECTOR(zz,T*k);
	
	
  GENERATED_VECTOR(phi_g,k);
  GENERATED_VECTOR(mu_g,k);
  GENERATED_VECTOR(sigma_g,k);
  GENERATED_VECTOR(h_g,hsize);
  GENERATED_SCALAR(nu_g);
  GENERATED_VECTOR(firstx,k);
  GENERATED_VECTOR(firstz,k);
  GENERATED_VECTOR(lastz,k);
  
}
MODEL_BLOCK{
  //
  // transformation of basic parameters and their priors
  //
  VectorXv phi(k);
  VectorXv sigma(k);
  var nu;
  // the phis
  var tmp0;
  for(int i=0;i<k;i++){
    //tmp0 = exp(lphi.coeff(i));
    //phi.coeffRef(i) = -1.0+2.0*tmp0/(1.0+tmp0); // scaled logit transformation
    uniform_logit<var> tphi(-1.0,1.0,lphi.coeff(i));
    phi.coeffRef(i) = tphi.par();
    target += tphi.implied_prior_lpdf();
   // target += lphi.coeff(i) - 2.0*log(1.0+tmp0); // uniform(0,1) prior on phi
  }
  phi_g = doubleValue(phi);
  
  // the mus
  target += normal_lpdf(mu,0.0,5.0);
  mu_g = doubleValue(mu);
  
  // the sigmas
  for(int i=0;i<k;i++){
    sigma.coeffRef(i) = exp(0.5*lsigmaSq.coeff(i));
    target += -0.5*k0*lsigmaSq.coeff(i) - 0.5*s0*k0/pow(sigma.coeff(i),2);
  }
  sigma_g = doubleValue(sigma);
  
  // the h-tildes
  target += normal_lpdf(h,0.0,10.0);
  h_g = doubleValue(h);
  MatrixXv H(k,k);
  H.setZero();
  int f=0;
  int l=k-1;
  for(int i=0;i<k-1;i++){
    H.coeffRef(i,i) = 1.0;
    H.col(i).tail(l) = h.segment(f,l);
    f += l;
    l -= 1;
  }
  H.coeffRef(k-1,k-1) = 1.0;
  
  
  // nu
  uniform_probit<var> tnu(10.0,60.0,lnu);
  nu = tnu.par();
  target += tnu.implied_prior_lpdf();
  nu_g = doubleValue(nu);
  
  //std::cout << "nu = " << doubleValue(nu) << std::endl;
  
  
  // generate x-processes form zz
  MatrixXv xx(T,k);
  VectorXv ztmp(T);
  VectorXv meantmp(T);
  var ssq;
  var d1n,dn,od,ytilde,pm1n,pm;
  var xxsum = 0.0;
  for(int p=0;p<k;p++){
    // build matrix
    ssq = square(sigma.coeff(p));
    d1n = 1.0/ssq + 0.5*nu;
    dn = (1.0+square(phi.coeff(p)))/ssq + 0.5*nu;
    od = -phi.coeff(p)/ssq;
    // lhs vector for initial iterate
    pm1n = mu.coeff(p)*(1.0-phi.coeff(p))/pow(sigma.coeff(p),2);
    pm = mu.coeff(p)*pow((phi.coeff(p)-1.0)/sigma.coeff(p),2);
    // initial time
    ytilde = (H.col(p).transpose()*yinv.block(0,0,k,k)).dot(H.col(p));
    meantmp.coeffRef(0) = pm1n + 0.5*nu*log(nu/ytilde);
    for(int t=1;t<T-1;t++){
      ytilde = (H.col(p).transpose()*yinv.block(t*k,0,k,k)).dot(H.col(p));
      meantmp.coeffRef(t) = pm + 0.5*nu*log(nu/ytilde);
    }  
    ytilde = (H.col(p).transpose()*yinv.block((T-1)*k,0,k,k)).dot(H.col(p));
    meantmp.coeffRef(T-1) = pm1n + 0.5*nu*log(nu/ytilde); 
    
    
    tridiagChol<var> L(T,d1n,dn,od); // cholesky factorization
    ztmp = zz.segment(T*p,T);
    xx.col(p) = L.solve(meantmp) + L.LT_solve(ztmp);
    target -= L.logDetL();
    xxsum += xx.col(p).sum();
    firstz(p) = doubleValue(ztmp(0));
    lastz(p) = doubleValue(ztmp(T-1));
    firstx(p) = doubleValue(xx(0,p));
  }
  
  
 // std::cout << "done choleskys, target = " << doubleValue(target) << std::endl;
  
  
  // prior on latent processes
  
  for(int p=0;p<k;p++){
    // intial time has marginal distribution
    target += normal_lpdf(xx.coeff(0,p),mu.coeff(p),sigma.coeff(p)/sqrt(1.0-pow(phi.coeff(p),2)));
 //   std::cout << "xloop" << std::endl;
    for(int t=1;t<T;t++){ 
      target += normal_lpdf(xx.coeff(t,p),mu.coeff(p)+phi.coeff(p)*(xx.coeff(t-1,p)-mu.coeff(p)),
                            sigma.coeff(p));
    }
    
  }
  
  //std::cout << "done xx" << std::endl;
  
  // common factor in measurement equation
  var mvgamma = 0.0;
  for(int s=1;s<=k;s++) mvgamma += stan::math::lgamma(0.5*(nu+1.0-static_cast<double>(s)));
  target += -0.5*(nu+dk+1.0)*ldets - dT*(0.5*nu*dk*log(2.0) + mvgamma) + 0.5*nu*xxsum;
  
  
  // remaining parts of measurement equation
  for(int t=0;t<T;t++){
    target += -0.5*((H*(xx.row(t).array().exp().matrix().asDiagonal())*H.transpose())*yinv.block(t*k,0,k,k)).trace();
    
  }
  
}

