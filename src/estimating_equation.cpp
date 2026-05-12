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
    arma::vec beta,
    arma::vec offset,
    Nullable<arma::vec> weights_ = R_NilValue,
    int maxit = 100,
    double tol = 1e-8,
    bool verbose = false) {
  
  const uword p = X.n_cols;
  const uword n = X.n_rows;
  
  // check that beta has the right length
  if(beta.n_elem != p) {
    stop("Length of beta must match number of columns in X");
  }
  if(Y.n_elem != n) {
    stop("Length of Y must match number of rows in X");
  }
  if(offset.n_elem != n) {
    stop("Length of offset must match number of rows in X");
  }

  vec weights;
  bool has_weights = weights_.isNotNull();
  if (has_weights) {
    weights = as<arma::vec>(weights_);
    if(weights.n_elem != n) {
      stop("Length of weights_ must match number of rows in X");
    }
  }
    
  for(int iter = 0; iter < maxit; iter++) {
    
    vec eta = X * beta + offset;
    vec mu, pdf, F, w;
    mat J;
    
    // ==============================
    // OIPCW score: expit
    // ==============================
    mu = expit_vec(eta);
    w  = mu % (1.0 - mu);

    if (has_weights) {
      F = crossprod_vec(X, weights % (Y - mu));
      J = -crossprod_weighted(X, weights % w);
    } else {
      F = crossprod_vec(X, Y - mu);
      J = -crossprod_weighted(X, w);
    }
    
    // Newton step
    // Fast but less stable: solve J step = F
    vec step;
    bool solvable = false;
    solvable = solve(step, J, F, solve_opts::no_approx);
    
    if (!solvable) {
      warning("Linear system could not be solved, stopping iteration");
      return beta;
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
