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
    arma::vec beta,
    std::string solve_opts = "fast",
    int maxit = 100,
    double tol = 1e-8,
    bool verbose = false) {
  
  const int p = X.n_cols;
  // check that beta has the right length
  if(beta.n_elem != p) {
    stop("Length of beta must match number of columns in X");
  }
  // vec beta(p, fill::zeros);
  
  bool is_oipcw_expit = model_type == "oipcw_expit";
  bool is_oipcw_probit = model_type == "oipcw_probit";
  bool is_nls_expit = model_type == "nls_expit";
  bool is_nls_probit = model_type == "nls_probit";
  
  for(int iter = 0; iter < maxit; iter++) {
    
    vec eta = X * beta;
    vec mu, pdf, F, w;
    mat J;
    
    // ==============================
    // OIPCW score: expit
    // ==============================
    if (is_oipcw_expit) {
      mu = expit_vec(eta);
      w  = mu % (1.0 - mu);
      
      F = crossprod_vec(X, Y - mu);
      J = -crossprod_weighted(X, w);
    }

    // ==============================
    // OIPCW score: probit
    // ==============================
    else if (is_oipcw_probit) {
      mu  = normcdf(eta);
      w  = mu % (1.0 - mu);
      pdf = normpdf(eta) / w;
      
      vec r = Y - mu;

      F = crossprod_vec(X, r % pdf);
      
      vec d_pdf = -eta % normpdf(eta) / w 
	- pdf % ( (1.0 - 2.0*mu) % normpdf(eta) / w );
  
      vec weight = pdf % normpdf(eta) + r % d_pdf;
  
      J = -crossprod_weighted(X, weight);

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
    // Fast but less stable: solve J step = F
    vec step;
    if (solve_opts == "fast") {
      step = solve(J, F, solve_opts::fast);
    } else if (solve_opts == "refine") {
      step = solve(J, F, solve_opts::refine);
    } else if (solve_opts == "equilibrate") {
      step = solve(J, F, solve_opts::equilibrate);
    } else if (solve_opts == "allow_ugly") {
      step = solve(J, F, solve_opts::allow_ugly);
    } else if (solve_opts == "no_approx") {
      step = solve(J, F, solve_opts::no_approx);
    } else if (solve_opts == "force_sym") {
      step = solve(J, F, solve_opts::force_sym);
    } else if (solve_opts == "force_approx") {
      step = solve(J, F, solve_opts::force_approx);
    } else {
      stop("Unknown solve_opts");
    }

    if (verbose) {
      Rcpp::Rcout << "step is " << step.t() << std::endl;
    }
    
    // Step halving for stability
    double step_factor = 1.0;
    vec beta_new = beta - step;
    
    while(!beta_new.is_finite() && step_factor > 1e-6) {
      step_factor *= 0.5;
      beta_new = beta - step_factor * step;
      if (verbose) {
	Rcpp::Rcout << "beta_new is " << beta_new.t() << std::endl;
      }
    }
    
    if(max(abs(beta_new - beta)) < tol)
      return beta_new;
    
    beta = beta_new;
  }
  
  warning("Did not converge");
  return beta;
}
