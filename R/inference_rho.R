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

#' Compute bias adjusted t-statistic for parameter rho
#'
#' @param fit_ML returned object from ML estimation
#' @param fit_ttest_theta fitted ttest_theta object
#' @param gammaS value used for gammaS (default: estimated value)
#' @param gammaR value used for gammaR (default: estimated value)
#' @param theta value used for theta (default: estimated value)
#' @param rho value used for rho (default: estimated value)
#' @param rho0 value of theta under null hypothesis
#'  (default: zero/test significance)
#' @param var_meat_estimator method to compute "meat" of the variance sandwich
#'  * "theorem" implement exactly as in theorem (put hats on everything)
#'  * "clustering" clustering based on first-order expansion (will not depend on rho_hat! )
#' @param trimming_bound specify trimming of probabilities
#' @param force_robust don't use speedglm (for robustness in sparse networks)
#'
#' @return list with t statistic and auxilary output
#' @export
#' @import data.table
#' @import Matrix
#' @importFrom speedglm speedlm.wfit
ttest_rho <- function(fit_ML, fit_ttest_theta, gammaS = fit_ML$gammaS_hat,
                      gammaR = fit_ML$gammaR_hat,
                      theta = fit_ML$theta_hat, rho = fit_ML$rho_hat, rho0 = 0,
                      var_meat_estimator = "theorem",
                      trimming_bound = 0,
                      force_robust = FALSE) {

  N <- fit_ML$model$N

  long <- copy(fit_ML$model$links)
  add_ystar(fit_ML$model, theta, gammaS, gammaR, long, trimming_bound)

  long[, p := compute_p(Ystar)]
  long[, `:=`(p1 = compute_p1(Ystar, p),
              omega = compute_omega(Ystar, p),
              H = compute_H(Ystar, p),
              p_ystar = partial_p_ystar(Ystar),
              p_ystar_ystar = partial_p_ystar_ystar(Ystar) )]

  add_xtilde(fit_ML$model, long, force_robust)

  wide <- long_to_wide(long, fit_ML$model$col_sender, fit_ML$model$col_receiver)

  wide[, `:=`(r_ij = compute_r(Ystar_ij, Ystar_ji, rho),
              r_rho_ij = partial_r_rho(Ystar_ij, Ystar_ji, rho))]
  wide[, r1_ij := r_ij * (1 - r_ij)]
  wide[, J_ij := r_rho_ij/r1_ij]

  wide <- copy_var_to_new_suffix(wide, c("r", "r_rho", "r1", "J"), "_ij", "_ji")

  wide <- upper_triangular_to_full(wide)

  wide[, `:=`(r_ystar1_ij = partial_r_ystar1(Ystar_ij, Ystar_ji, rho),
              r_ystar1_ji = partial_r_ystar1(Ystar_ji, Ystar_ij, rho),
              r_ystar1_ystar1_ij = partial_r_ystar1_ystar1(Ystar_ij, Ystar_ji, rho),
              r_ystar1_ystar2_ij = partial_r_ystar1_ystar2(Ystar_ij, Ystar_ji, rho))]

  wide[, `:=`(J_ystar1_ij = partial_J_ystar1(Ystar_ij, Ystar_ji, rho, r_ij, r_rho_ij),
              J_ystar1_ji = partial_J_ystar1(Ystar_ji, Ystar_ij, rho, r_ji, r_rho_ji),
              rho_tilde_ij = (r_ij - p_ij * p_ji) / sqrt(p1_ij * p1_ji),
              Omega_prep_ij = J_ij * r_ystar1_ij / omega_ij)]

  wide$Omega_prep_ij[is.infinite(wide$Omega_prep_ij)] <- 0

  if (length(fit_ML$model$X_names) > 0) {
    wide[, (sprintf("T_%s_ij", fit_ML$model$X_names)) :=
           lapply(.SD, function(x, J_ij, r_ystar1_ij) { - J_ij * r_ystar1_ij * x},
                  J_ij = J_ij, r_ystar1_ij = r_ystar1_ij),
         .SDcols = sprintf("Xtilde_%s_ij", fit_ML$model$X_names)]
    T_N <- unlist(wide[, setNames(lapply(.SD, mean), paste0("T_N_", fit_ML$model$X_names)),
                       .SDcols = paste0("T_", fit_ML$model$X_names, "_ij")])
    T_times_W1_inv <- T_N %*% as.matrix(fit_ttest_theta$W1_inv)

    wide[, t_ij := as.vector(T_times_W1_inv %*% t(as.matrix(.SD))),
         .SDcols = sprintf("Xtilde_%s_ij", fit_ML$model$X_names)]
    wide[, t_ji := as.vector(T_times_W1_inv %*% t(as.matrix(.SD))),
         .SDcols = sprintf("Xtilde_%s_ji", fit_ML$model$X_names)]

    bias_theta_contribution <- as.vector(T_times_W1_inv %*% fit_ttest_theta$B_theta)
  } else {
    wide[, `:=`(t_ij = 0, t_ji = 0)]
    bias_theta_contribution <- 0
  }

  # proj <- speedlm.wfit(wide$Omega_prep_ij, fit_ML$model$design_matrix_FE, wide$omega_ij,
  #                      sparse = TRUE)
  # proj_val <- as.vector(fit_ML$model$design_matrix_FE %*% coef(proj))
  # wide[, `:=`(Omega_ij = proj_val)]
  add_predict(fit_ML$model, wide, "Omega_prep_ij", "Omega_ij", "omega_ij", force_robust,
              residuals = FALSE)

  wide_reverse <- wide[, .(j = i, i = j, Omega_ji = Omega_ij)]
  wide <- merge(wide, wide_reverse)
  setkey(wide, i, j)

  wide[, `:=`(num_diag_ij = p_ystar_ij * J_ystar1_ij * r_ij/p_ij +
                0.5 * Omega_ij * H_ij * p_ystar_ystar_ij - J_ystar1_ij * r_ystar1_ij -
                0.5 * J_ij * r_ystar1_ystar1_ij,
              num_off_ij = - (J_ystar1_ij * r_ystar1_ji + J_ystar1_ji * r_ystar1_ij +
                                J_ij * r_ystar1_ystar2_ij),
              num_corr_ij = rho_tilde_ij * sqrt(omega_ij * omega_ji),
              v22_ij = 4 * (t_ij - Omega_ij) * J_ij * p_ystar_ij * r_ij/p_ij +
                2 * (t_ij - Omega_ij)^2 * omega_ij +
                2 * (t_ij - Omega_ij) * (t_ji - Omega_ji) *
                sqrt(omega_ij * omega_ji) * rho_tilde_ij,
              v_meat_clustering_ij = J_ij * (Y_ij * Y_ji - r_ij) +
                (t_ij - Omega_ij) * H_ij * (Y_ij - p_ij) +
                (t_ji - Omega_ji) * H_ji * (Y_ji - p_ji))]

  B_rho_overi <- wide[, .(S = sum(num_diag_ij)/sum(omega_ij),
                          SR = sum(num_corr_ij) * sum(num_off_ij) / (sum(omega_ij) * sum(omega_ji))),
                      by = i]
  B_rho_S <- B_rho_overi[, mean(S)]
  B_rho_SR <- B_rho_overi[, mean(SR)]

  B_rho_overj <- wide[, .(R = sum(num_diag_ij)/sum(omega_ij)), by = j]
  B_rho_R <- B_rho_overj[, mean(R)]

  v1 <- wide[, mean(J_ij * r_rho_ij)]

  if (var_meat_estimator == "theorem") {
    v22 <- wide[, mean(v22_ij)]
    v2 <- v1 + v22
  } else if (var_meat_estimator == "clustering") {
    v2 <- wide[, mean(v_meat_clustering_ij^2)]
  } else
    stop(sprintf("Unknown method to compute variance. You provide %s", var_meat_estimator))

  bias <- 2 * (B_rho_S + B_rho_R + B_rho_SR + bias_theta_contribution) / (N  * v1)
  rho_hat_bias_corr <- fit_ML$rho_hat_MLE - bias

  se_rho <- sqrt(2 * v2)/(v1 * N)
  tstat <- (rho_hat_bias_corr - rho0) / se_rho
  pval <- 2 * pnorm(abs(tstat), lower.tail = FALSE)

  # WriteXLS::WriteXLS(ds2, ExcelFileName = 'newcode.xls')
  c(fit_ML, list(bias_rho = bias, se_rho_hat = se_rho,
                 rho_hat_bias_corr = rho_hat_bias_corr, tstat_rho = tstat,
                 pval_rho = pval))
}
