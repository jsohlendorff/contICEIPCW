### compare_to_reference.R --- 
#----------------------------------------------------------------------
## Author: Johan Sebastian Ohlendorff
## Created: Mar 24 2026 (17:21) 
## Version: 
## Last-Updated: Mar 24 2026 (17:31) 
##           By: Johan Sebastian Ohlendorff
##     Update #: 9
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

#' Compare treatment groups to a reference group using influence curve-based SEs
#' @description
#' This function takes a reference group and one or more treatment groups, and computes the difference in estimates between each treatment group and the reference group, along with standard errors and confidence intervals based on the influence curve (IC) differences.
#' @param reference_group The name of the reference group (must match `treatment_name
#' in the input groups).
#' @param ... One or more group objects (output from `estimate_effects`) to
#' compare against the reference group. Each group must have a `treatment_name`, `ic` with `ic` and `time_horizons`, and `results` with `time_horizon` and `estimate`.
#' @return A data.table with columns: `treatment`, `reference_group`, `
#' time_horizon`, `estimate` (difference in estimates), `se` (standard error), `lower` and `upper` (95% confidence interval bounds).
#' @details
#' The function computes the standard error of the difference in estimates between each treatment group and the
#' reference group by calculating the standard deviation of the difference in influence curves (ICs) divided by the square root of the sample size. It assumes that the ICs for the treatment and reference groups are aligned by time horizon and have the same length.
#' @export
#' 
compare_to_reference <- function(reference_group, ...) {
  groups <- list(...)
  
  # Convert list of groups into named list by treatment_name
  group_list <- stats::setNames(groups, sapply(groups, function(x) x$treatment_name))
  
  # Extract reference group components
  ref <- group_list[[reference_group]]
  ref_ic <- ref$ic$ic
  ref_results <- ref$results
  
  # Prepare output collector
  out_list <- list()
  
  # Loop over all non-reference groups
  for (tr in names(group_list)) {
    if (tr == reference_group) next
    
    g <- group_list[[tr]]
    
    # merge results by time horizon
    merged <- merge(
      g$results,
      ref_results,
      by = "time_horizon",
      suffixes = c("_tr", "_ref")
    )
    
    # Compute SE from IC difference for each time horizon
    # Need to filter influence curves by matching time horizon
    dt_list <- lapply(merged$time_horizon, function(th) {
      
      # indices for correct time horizon in the stacked IC vector
      idx_tr <- which(g$ic$time_horizons == th)
      idx_ref <- which(ref$ic$time_horizons == th)
      
      ic_tr_h <- g$ic$ic[[idx_tr]]
      ic_ref_h <- ref_ic[[idx_ref]]
      
      if (length(ic_tr_h) != length(ic_ref_h)) {
        stop("IC vectors for treatment and reference must have same length per horizon.")
      }
      
      diff_ic <- ic_tr_h - ic_ref_h
      N <- length(diff_ic)
      se <- sd(diff_ic) / sqrt(N)
      
      estimate_diff <- merged$estimate_tr[merged$time_horizon == th] -
                       merged$estimate_ref[merged$time_horizon == th]
      
      data.table(
        treatment = tr,
        reference_group = reference_group,
        time_horizon = th,
        estimate = estimate_diff,
        se = se,
        lower = estimate_diff - 1.96 * se,
        upper = estimate_diff + 1.96 * se
      )
    })
    
    out_list[[tr]] <- rbindlist(dt_list)
  }
  
  rbindlist(out_list)
}

######################################################################
### compare_to_reference.R ends here
