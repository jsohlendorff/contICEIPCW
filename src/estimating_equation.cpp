// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

inline vec expit_vec(const vec& x) {
  return 1.0 / (1.0 + exp(-x));
}

// weighted crossproduct: X' diag(w) X
inline mat crossprod_weighted(const mat& X, const vec& w) {
  return X.t() * (X.each_col() % w);
}

// weighted score: X' v
inline vec crossprod_vec(const mat& X, const vec& v) {
  return X.t() * v;
}

// [[Rcpp::export]]
arma::vec estimating_equation_cpp(
    const arma::mat& X,
    const arma::vec& Y,
    std::string model_type,
    int maxit = 100,
    double tol = 1e-8
) {
  
  const int p = X.n_cols;
  vec beta(p, fill::zeros);
  
  bool is_oipcw     = model_type.find("oipcw") != std::string::npos;
  bool is_nls_expit = model_type == "nls_expit";
  bool is_nls_probit= model_type == "nls_probit";
  bool is_expit     = model_type.find("expit") != std::string::npos;
  
  for(int iter = 0; iter < maxit; iter++) {
    
    vec eta = X * beta;
    vec mu, pdf, F, w;
    mat J;
    
    // ==============================
    // OIPCW score: X'(Y - mu)
    // ==============================
    if(is_oipcw) {
      
      if(is_expit) {
        mu = expit_vec(eta);
        w  = mu % (1.0 - mu);
      } else {
        mu  = normcdf(eta);
        w   = normpdf(eta);
      }
      
      F = crossprod_vec(X, Y - mu);
      J = -crossprod_weighted(X, w);
    }
    
    // ==============================
    // Nonlinear LS – expit
    // ==============================
    else if(is_nls_expit) {
      
      mu  = expit_vec(eta);
      pdf = mu % (1.0 - mu);
      
      vec r = Y - mu;
      
      F = crossprod_vec(X, r % pdf);
      
      vec W = -pdf % pdf - r % pdf % (1.0 - 2.0*mu);
      J = crossprod_weighted(X, W);
    }
    
    // ==============================
    // Nonlinear LS – probit
    // ==============================
    else if(is_nls_probit) {
      
      mu  = normcdf(eta);
      pdf = normpdf(eta);
      
      vec r = Y - mu;
      
      F = crossprod_vec(X, r % pdf);
      
      vec W = -pdf % pdf - r % eta % pdf;
      J = crossprod_weighted(X, W);
    }
    else {
      stop("Unknown model_type");
    }
    
    // Newton step
    vec step = solve(J, F, solve_opts::fast);
    
    // Step halving for stability
    double step_factor = 1.0;
    vec beta_new = beta - step;
    
    while(!beta_new.is_finite() && step_factor > 1e-6) {
      step_factor *= 0.5;
      beta_new = beta - step_factor * step;
    }
    
    if(max(abs(beta_new - beta)) < tol)
      return beta_new;
    
    beta = beta_new;
  }
  
  warning("Did not converge");
  return beta;
}
