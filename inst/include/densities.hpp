#include <stan/math.hpp>


/*
 * density of x, so that logit^{-1}(x) = \exp(x)/(1+exp(x)) \sim beta(alpha,beta) (pdf \propo x^{a-1}(1-x)^{b-1})
 */

template <class Tx_, class Ta_, class Tb_>
inline typename boost::math::tools::promote_args<Tx_,Ta_,Tb_>::type beta_logit_lpdf(const Tx_ x, const Ta_ a, const Tb_ b){
  Tx_ y = stan::math::inv_logit(x);
  return a*log(y) + b*log(1.0-y) - stan::math::lbeta(a,b);
}

template <class Tx_, class Ta_, class Tb_>
inline typename boost::math::tools::promote_args<Tx_,Ta_,Tb_>::type beta_logit_lpdf_kernel(const Tx_ x, const Ta_ a, const Tb_ b){
  Tx_ y = stan::math::inv_logit(x);
  return a*log(y) + b*log(1.0-y);
}

/*
 *  density of x, so that \exp(x) \sim \Gamma(a,b) (note scale parameterization: pdf \propto x^{a-1}\exp(-x/b))
 */

template <class Tx_, class Ta_, class Tb_>
inline typename boost::math::tools::promote_args<Tx_,Ta_,Tb_>::type gamma_log_lpdf(const Tx_ x, const Ta_ a, const Tb_ b){
  return a*x - exp(x)/b - a*log(b) - stan::math::lgamma(a);
}

template <class Tx_, class Ta_, class Tb_>
inline typename boost::math::tools::promote_args<Tx_,Ta_,Tb_>::type gamma_log_lpdf_kernel(const Tx_ x, const Ta_ a, const Tb_ b){
  return a*x - exp(x)/b;
}


