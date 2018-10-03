## This R package implements the methods proposed in
## Dzemski, Andreas: An empirical model of dyadic link formation in
## a network with unobserved heterogeneity, Review of Economics and Statistics, forthcoming

## Copyright (C) 2018  Andreas Dzemski

## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.

## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.

## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <https://www.gnu.org/licenses/>.

#' Compute specification test based on transitive triangles
#'
#' @param fit_ML returned object from ML estimation (\code{MLE_stage1})
#' @param fit_ttest_theta fitted \code{ttest_theta} object
#' @param gammaS vector of sender effects  (default: estimated value in \code{fit_ML})
#' @param gammaR vetor of receiver effects (default: estimated value in \code{fit_ML})
#' @param theta value used for \eqn{\theta} (default: estimated value in \code{fit_ML})
#' @param rho value used for \eqn{\rho} (default: estimated value in \code{fit_ML});
#'  not needed if both \code{bias_estimator = "theorem2"} and \code{bs_inference = "none"}
#' @param rho0 value of theta under null hypothesis
#'  (default: zero/test significance)
#' @param bs_iterations number of bootstrap iterations
#' @param alpha alpha used to compute bootstrap critical values
#' @param bias_estimator specify estimator of bias term
#' \itemize{
#'   \item \code{"theorem1"} based on theorem in paper, needs estimate of rho
#'   \item \code{"theorem2"} based on theorem in paper, doesn't need estimate of rho
#'   \item \code{"boostrap"} boostrap the bias
#'   }
#' @param bs_inference specify what kind of bootstrap is used to determine whether the test rejects
#' \itemize{
#'   \item \code{"none"} no bootstrap inference
#'   \item \code{"double"} double bootstrap
#'   }
#' @param trimming_bound specify trimming of probabilities used to compute statistics
#' @param bs_trimming_bound specify trimming of probabilities used in
#'  simulation of bootstrap distribution
#' @param force_robust logical, don't use speedglm
#' @param bs_force_robust logical, don't use speedglm in simulation of bootstrap distribution
#' @param compute_other_variances logical, force computation of standard errors for infeasible
#'   triangle tests
#'
#' @return list with test statistic, rejection and auxilary output
#'
#' @export
#' @import data.table
#' @import Matrix
#' @importFrom speedglm speedlm.wfit
triangle_test <- function(fit_ML, fit_ttest_theta, gammaS = fit_ML$gammaS_hat,
                          gammaR = fit_ML$gammaR_hat,
                          theta = fit_ML$theta_hat, rho = fit_ML$rho_hat_MLE, rho0 = 0,
                          bs_iterations = 200,
                          alpha = 0.1,
                          bias_estimator = "theorem2", bs_inference = "none",
                          trimming_bound = 0, bs_trimming_bound = trimming_bound,
                          force_robust = FALSE, bs_force_robust = force_robust,
                          compute_other_variances = FALSE) {

  if (bias_estimator == "bootstrap" & bs_iterations <= 0)
    stop("Cannot bootstrap bias without bootstrap iterations (bs_iterations has to be a positive integer).")
  if (bs_inference == "double" & bs_iterations <= 0)
    stop("Cannot bootstrap bias without bootstrap iterations (bs_iterations has to be a positive integer).")

  N <- fit_ML$model$N

  links <- copy(fit_ML$model$links)
  add_ystar(fit_ML$model, theta, gammaS, gammaR, links, trimming_bound)

  links[, p := compute_p(Ystar)]
  links[, `:=`(p1 = compute_p1(Ystar, p),
               omega = compute_omega(Ystar, p),
               H = compute_H(Ystar, p),
               p_ystar = partial_p_ystar(Ystar),
               p_ystar_ystar = partial_p_ystar_ystar(Ystar) )]

  add_xtilde(fit_ML$model, links, force_robust)

  trans_stats <- fast_excess_transitivity(links, fit_ML$model$col_outcome, "p")

  links <- long_to_wide(links, fit_ML$model$col_sender, fit_ML$model$col_receiver)

  if (bias_estimator == "theorem1") {
    links[, `:=`(r_ij = compute_r(Ystar_ij, Ystar_ji, rho))]
    links <- copy_var_to_new_suffix(links, "r", "_ij", "_ji")
    links <- upper_triangular_to_full(links)
    links[, `:=`(rho_tilde_ij = (r_ij - p_ij * p_ji)/sqrt(p1_ij * p1_ji))]
  } else {
    links <- upper_triangular_to_full(links)
  }

  links_new_vars <- data.table(triangles_biascorr_new_vars(links))
  setkey(links_new_vars, i, j)
  links[links_new_vars, `:=`(boldbeta_ij = boldbeta, bss_ij = bss, bsr_ij = bsr,
                             bssr_ij = bssr)]

  add_predict(fit_ML$model, links, "boldbeta_ij", "boldbeta_tilde_ij", "omega_ij",
              force_robust = force_robust)

  if (length(fit_ML$model$X_names) > 0) {
    links[, (sprintf("U_%s_ij", fit_ML$model$X_names)) :=
            lapply(.SD, function(x, boldbeta_ij, omega_ij) { boldbeta_ij * omega_ij * x },
                   boldbeta_ij = boldbeta_ij, omega_ij = omega_ij),
          .SDcols = sprintf("Xtilde_%s_ij", fit_ML$model$X_names)]

    U_N <- unlist(links[, setNames(lapply(.SD, mean), paste0("U_N_", fit_ML$model$X_names)),
                        .SDcols = paste0("U_", fit_ML$model$X_names, "_ij")])
    U_times_W1_inv <- U_N %*% as.matrix(fit_ttest_theta$W1_inv)

    links[, u_ij := as.vector(U_times_W1_inv %*% t(as.matrix(.SD))),
          .SDcols = sprintf("Xtilde_%s_ij", fit_ML$model$X_names)]
    links[, u_ji := as.vector(U_times_W1_inv %*% t(as.matrix(.SD))),
          .SDcols = sprintf("Xtilde_%s_ji", fit_ML$model$X_names)]

    bias_theta_contribution <- as.vector(U_times_W1_inv %*% fit_ttest_theta$B_theta)
  } else {
    links[, `:=`(u_ij = 0, u_ji = 0)]
    bias_theta_contribution <- 0
  }

  if (bias_estimator %in% c("theorem1", "theorem2")) {

    B_SS <-
      links[, .(S = sum(H_ij * p_ystar_ystar_ij * boldbeta_tilde_ij +
                          bss_ij)/sum(omega_ij)), by = i][, mean(S)]/2
    B_SR <-
      links[, .(R = sum(H_ij * p_ystar_ystar_ij * boldbeta_tilde_ij + bsr_ij)/sum(omega_ij)),
            by = j][, mean(R)]/2

    if (bias_estimator == "theorem1") {
      B_SSR <-
        links[, .(SR = sum(rho_tilde_ij * sqrt(omega_ij * omega_ji)) *
                    sum(bssr_ij)/(sum(omega_ij) * sum(omega_ji))), by = i][, mean(SR)]
    } else if (bias_estimator == "theorem2") {
      B_SSR <-
        links[, .(SR = sum(H_ij * H_ji * (.SD[[1]] - p_ij) * (.SD[[2]] - p_ji)) *
                    sum(bssr_ij)/(sum(omega_ij) * sum(omega_ji))),
              by = i, .SDcols = paste0(fit_ML$model$col_outcome, c("_ij", "_ji"))][, mean(SR)]
    }

    bias <- - (B_SS + B_SR + B_SSR + bias_theta_contribution)/N
  } else if (bias_estimator == "bootstrap") {
    links[, `:=`(Ystar_bs_inner_ij = Ystar_ij, Ystar_bs_inner_ji = Ystar_ji)]
    bias <- mean(double_bootstrap_inner(links, rho, fit_ML$model, bs_iterations,
                                        trimming_bound))
  }

  copy_ij_value_to_ji(links, "boldbeta_tilde_ij")
  # links[, `:=`(v2_spec_test_theorem_ij = (boldbeta_tilde_ij - u_ij)^2 * omega_ij +
  #                (boldbeta_tilde_ij - u_ij) * (boldbeta_tilde_ji - u_ji) * rho_tilde_ij *
  #                sqrt(omega_ij * omega_ji))]
  links[, `:=`(v_spec_test_theorem_ij = (boldbeta_tilde_ij - u_ij)^2 * H_ij^2 * (.SD[[1]] - p_ij)^2 +
                 (boldbeta_tilde_ij - u_ij) * (boldbeta_tilde_ji - u_ji) *
                 H_ij * H_ji * (.SD[[1]] - p_ij) * (.SD[[2]] - p_ji)),
        .SDcols = paste0(fit_ML$model$col_outcome, c("_ij", "_ji"))]

  se_analytic <- sqrt(links[, mean(v_spec_test_theorem_ij)] / num_pairs(N))

  excess_trans_bias_corr <- trans_stats$excess_trans - bias
  teststat <- excess_trans_bias_corr/se_analytic
  pval_analytic <- 2 * pnorm(abs(teststat), lower.tail = FALSE)

  out_list <- list(teststat = teststat, pval_analytic = pval_analytic,
                   se_analytic = se_analytic, bias = bias,
                   obs_trans = trans_stats$num_triangles,
                   expected_trans = trans_stats$expected_triangles,
                   excess_trans = trans_stats$excess_trans,
                   excess_trans_bias_corr = excess_trans_bias_corr)

  if (bs_inference == "double") {
    bs <- double_bootstrap(links, rho, fit_ML$model, bs_iterations, alpha, trimming_bound,
                           bias_estimator, bs_force_robust)
    rej_bs_perc_t <-
      (teststat < bs$bs_quantiles_t[[1]]) | (teststat > bs$bs_quantiles_t[[2]])

    rej_bs_perc_excess <-
      (excess_trans_bias_corr < bs$bs_quantiles_excess[[1]]) | (excess_trans_bias_corr > bs$bs_quantiles_excess[[2]])

    out_list <- c(out_list, list(rej_bs_perc_t = rej_bs_perc_t,
                                 rej_bs_perc_excess = rej_bs_perc_excess))
  }

  if (compute_other_variances) {
    add_predict(fit_ML$model, links, "boldbeta_ij", "Psi_ij", "omega_ij",
                force_robust = force_robust, residuals = FALSE )
    copy_ij_value_to_ji(links, "Psi_ij")
    copy_ij_value_to_ji(links, "boldbeta_ij")

    links[, `:=`(v_triangles_ttest_ij = (Psi_ij - u_ij)^2 * H_ij^2 * (.SD[[1]] - p_ij)^2 +
                   (Psi_ij - u_ij) * (Psi_ji - u_ji) *
                   H_ij * H_ji * (.SD[[1]] - p_ij) * (.SD[[2]] - p_ji),
                 v_oracle_ij = boldbeta_ij^2 * H_ij^2 * (.SD[[1]] - p_ij)^2 +
                   boldbeta_ij * boldbeta_ji * H_ij * H_ji * (.SD[[1]] - p_ij) * (.SD[[2]] - p_ji)),
          .SDcols = paste0(fit_ML$model$col_outcome, c("_ij", "_ji"))]

    se_predicted_triangles <- sqrt(links[, mean(v_triangles_ttest_ij)] / num_pairs(N))
    se_oracle <- sqrt(links[, mean(v_oracle_ij)] / num_pairs(N))

    naive_teststat <- trans_stats$excess_trans / se_oracle
    naive_pval <- 2 * pnorm(abs(naive_teststat), lower.tail = FALSE)

    naive_bc_teststat <- excess_trans_bias_corr / se_oracle
    naive_bc_pval <- 2 * pnorm(abs(naive_bc_teststat), lower.tail = FALSE)

    out_list <- c(out_list, list(se_predicted_triangles = se_predicted_triangles,
                                 se_oracle = se_oracle,
                                 naive_teststat = naive_teststat,
                                 naive_pval = naive_pval,
                                 naive_bc_teststat = naive_bc_teststat,
                                 naive_bc_pval = naive_bc_pval))
  }

  out_list
}

#' Double bootstrap procedure
#'
#' @param wide link data in wide format
#' @param rho_bs0 within-dyad correlation of bootstrap distribution
#' @param model model description in special list format
#' @param bs_iterations number of bootstrap replications
#' @param alpha significance level
#' @param trimming_bound bound for trimming small probabilities
#' @param bias_estimator estimator of bias (see triangle test)
#' @param force_robust logical, force using numerically stable (but possibly slow)
#'  algorithm to fit probit
#'
#' @return list with bootstrapped SE and quantiles
double_bootstrap <- function(wide, rho_bs0, model, bs_iterations = 200, alpha = 0.1,
                             trimming_bound = 0, bias_estimator = "theorem2", force_robust = FALSE) {

  wide <- wide[i < j, c("i", "j", "Ystar_ij", "Ystar_ji", sprintf("%s_ij", model$X_names)),
              with = FALSE]
  wide <- wide[, `:=`(Ystar_bs_ij = Ystar_ij, Ystar_bs_ji = Ystar_ji)]
  model$col_outcome <- "Y_bs"
  model$col_sender <- "i"
  model$col_receiver <- "j"

  bs_tstat <- numeric(bs_iterations)
  bs_excess <- numeric(bs_iterations)

  for (b in seq_len(bs_iterations)) {
    model$links <-
      draw_bootstrap_sample(wide, rho_bs0, model, wide$Ystar_bs_ij, wide$Ystar_bs_ji, "Y_bs")

    fitted_ML_bs <- MLE_stage1(model)

    ttest_theta_bs <-
      ttest_theta(fitted_ML_bs, trimming_bound = trimming_bound, force_robust = force_robust)

    ttest_bs <-
      triangle_test(fitted_ML_bs, ttest_theta_bs, rho = rho_bs0, bs_iterations = bs_iterations,
                    trimming_bound = trimming_bound, force_robust = force_robust,
                    bias_estimator = bias_estimator, compute_other_variances = FALSE)

    bs_tstat[[b]] <- ttest_bs$teststat
    bs_excess[[b]] <- ttest_bs$excess_trans_bias_corr
  }

  # allow a small number of bootstrap iterations to be thrown out due to numerical problems
  if (sum(!is.finite(bs_tstat)) < 0.1 * bs_iterations) {
    bs_tstat <- bs_tstat[is.finite(bs_tstat)]
    bs_excess <- bs_excess[is.finite(bs_excess)]
  }

  list(se_bs = sd(bs_excess),
       bs_quantiles_t = quantile(bs_tstat, c(alpha/2, 1 - alpha/2)),
       bs_quantiles_excess = quantile(bs_excess, c(alpha/2, 1 - alpha/2)))
}

#' Fast computation of bootstrap bias
#'
#' @param wide link data in wide format
#' @param rho_bs0 within-dyad correlation of bootstrap distribution
#' @param model model description in special list format
#' @param bs_iterations number of bootstrap replications
#' @param trimming_bound bound for trimming small probabilities
#'
#' @return bootstrap sample (vector)
double_bootstrap_inner <- function(wide, rho_bs0, model, bs_iterations = 200, trimming_bound = 0) {

  wide <- wide[i < j, c("i", "j", "Ystar_bs_inner_ij", "Ystar_bs_inner_ji",
                        sprintf("%s_ij", model$X_names)),
              with = FALSE]
  model$col_outcome <- "Y_bs_inner"

  excess_trans_bs <- numeric(bs_iterations)

  for (b in seq_len(bs_iterations)) {
    model$links <-
      draw_bootstrap_sample(wide, rho_bs0, model, wide$Ystar_bs_inner_ij,
                            wide$Ystar_bs_inner_ji, "Y_bs_inner")

    fitted_ML_bs <- MLE_stage1(model)

    add_ystar(model, fitted_ML_bs$theta_hat_MLE, fitted_ML_bs$gammaS_hat,
              fitted_ML_bs$gammaR_hat, model$links, trimming_bound = trimming_bound)
    model$links[, p_bs_inner := compute_p(Ystar)]

    excess_trans_bs[[b]] <-
      fast_excess_transitivity(model$links, "Y_bs_inner", "p_bs_inner")$excess_trans
  }

  excess_trans_bs
}

#' Draw from bootstrap distribution
#'
#' @param wide link data wide format (has to contain i and j, canonically ordered)
#' @param rho_bs0 within-dyad correlation of bootstrap distribution
#' @param model model description in special list format
#' @param values_ystar_ij of bootstrap process
#' @param values_ystar_ji of bootstrap process
#' @param col_y_bs column name for bootstrapped Y
#'
#' @return data.table bootstrap sample
draw_bootstrap_sample <-
  function(wide, rho_bs0, model, values_ystar_ij, values_ystar_ji, col_y_bs = "Y_bs") {

  discarded_samples <- 0

  repeat {
    shocks_bs <- mvrnorm(nrow(wide), mu = c(0,0),
                         Sigma = matrix( c(1, rho_bs0, rho_bs0, 1), ncol = 2) )
    cols <- paste0(col_y_bs, c("_ij", "_ji"))
    wide[, (cols) := .(as.integer(values_ystar_ij >= shocks_bs[, 1]),
                     as.integer(values_ystar_ji >= shocks_bs[, 2]))]

    bs_sample <-
      upper_triangular_to_full(wide)[, c(.(i = i, j = j),
                                         setNames(.SD, c(col_y_bs, model$X_names))),
                                     .SDcols = paste0(c(col_y_bs, model$X_names), "_ij")]

    if (!(all(bs_sample[[col_y_bs]] == 0) | all(bs_sample[[col_y_bs]] == 1)))
      break

    warning("Discarding bootstrap sample with no variation.")

    discarded_samples <- discarded_samples + 1

    if (discarded_samples > 1)
      stop("Bootstrap samples have no variation")
  }

  bs_sample
}

#' Compute excess transitivity
#'
#' count transitive triangles and compute excess transitivity at the same time
#' for efficiency gain
#'
#' @param links data frame with link data
#' @param col_link_indicator name of column with binary link indicator
#' @param col_link_probability name of column with link probability
#' @param col_sender name of column that gives sender id
#' @param col_receiver name of column that gives receiver id
#'
#' @return list with number of transitive triangles, expected triangles, excess transitivity
#' @export
fast_excess_transitivity <-
  function(links, col_link_indicator, col_link_probability, col_sender = "i",
           col_receiver = "j") {

    # setkeyv(links, c(col_sender, col_receiver))
    out <- excess_trans(links, col_link_indicator, col_link_probability)
    if (out$error) stop("C++ code threw exception.")

    out[c("num_triangles", "expected_triangles", "excess_trans")]
  }
