#' Compute bias adjusted t-statistic for parameter theta
#'
#' @param fit_ML returned object from ML estimation
#' @param gammaS value used for gammaS (default: estimated value)
#' @param gammaR value used for gammaR (default: estimated value)
#' @param theta value used for theta (default: estimated value)
#' @param rho value used for rho (default: estimated value)
#' @param theta0 value of theta under null hypothesis
#'  (default: zero vector/test significance)
#' @param var_meat_estimator  method to compute "meat" of the variance sandwich
#'  * "theorem" implement exactly as in theorem (put hats on everything)
#'  * "clustering" clustering based on first-order expansion (will not depend on rho_hat! )
#' @param trimming_bound specify trimming of probabilities
#' @param force_robust don't attempt to use speedglm
#'
#' @return list with t statistic and auxilary output
#' @export
#' @import data.table
#' @import Matrix
#' @importFrom speedglm speedlm.wfit
ttest_theta <- function(fit_ML, gammaS = fit_ML$gammaS_hat, gammaR = fit_ML$gammaR_hat,
                        theta = fit_ML$theta_hat, rho = fit_ML$rho_hat,
                        theta0 = rep.int(0, length(theta)),
                        var_meat_estimator = "clustering",
                        trimming_bound = 0,
                        force_robust = FALSE) {

  Xtilde_names <- sprintf('Xtilde_%s', fit_ML$model$X_names)
  N <- fit_ML$model$N
  N2 <- num_pairs(N)

  links <- copy(fit_ML$model$links)
  links[, ("tmp_Y") := .SD, .SDcols = fit_ML$model$col_outcome]
  add_ystar(fit_ML$model, theta, gammaS, gammaR, links, trimming_bound)
  links[, `:=`(omega = compute_omega(Ystar))]
  add_xtilde(fit_ML$model, links, force_robust)

  compute_B_or_D <- function(bycolumn) {

    process_node <- function(node_data) {
      Xtilde_node <- as.matrix(node_data[, Xtilde_names, with = FALSE])

      # numerator inner loop in matrix notation
      numerator <- weighted_inner_prod(Xtilde_node, Xtilde_node, node_data$omega)
      denom <- sum(node_data$omega)

      # 'abuse' data.table's functionality of having list typed columns
      list(list(numerator/denom))
    }

    lst_B_or_D <- links[, .(B_or_D_i = process_node(.SD)), by = bycolumn][['B_or_D_i']]

    Reduce(function(a, b) a[[1]] + b[[1]], lst_B_or_D) %*% theta / (2 * fit_ML$model$N)
  }

  links <- data.table(links, key = c(fit_ML$model$col_sender, fit_ML$model$col_receiver))

  B_theta_S <- compute_B_or_D(bycolumn = fit_ML$model$col_sender)
  B_theta_R <- compute_B_or_D(bycolumn = fit_ML$model$col_receiver)
  B_theta <- B_theta_S + B_theta_R

  Xtilde <- as.matrix(links[, Xtilde_names, with = FALSE])
  W1 <- weighted_inner_prod(Xtilde, Xtilde, links[, omega])/N2
  W1_inv <- solve(W1)

  bias <- setNames(as.vector(W1_inv %*% (B_theta)/N), names(fit_ML$theta_hat_MLE))

  theta_hat_bias_corr <- fit_ML$theta_hat_MLE - bias

  if (tolower(var_meat_estimator) == "clustering") {
    links[, `:=`(H = omega / partial_p_ystar(Ystar), p = compute_p(Ystar))]

    wide <- long_to_wide(links, fit_ML$model$col_sender, fit_ML$model$col_receiver,
                         ij_variables = c('H', Xtilde_names, 'tmp_Y', 'p'))

    var_ij_names <- sprintf("var_ij_%s", Xtilde_names)

    m_ij <- as.matrix(wide[, lapply(.SD, function(x) x * H_ij * (tmp_Y_ij - p_ij)),
                            .SDcols = paste0(Xtilde_names, '_ij')])
    m_ji <- as.matrix(wide[, lapply(.SD, function(x) x * H_ij * (tmp_Y_ij - p_ij)),
                            .SDcols = paste0(Xtilde_names, '_ji')])
    m <- m_ij + m_ji

    W2 <- t(m) %*% m / (2 * N2)
  } else if (tolower(var_meat_estimator) == "theorem") {
    links[, `:=`(p = compute_p(Ystar))]
    links[, `:=`(p1 = compute_p1(Ystar, p))]

    wide <- long_to_wide(links, fit_ML$model$col_sender, fit_ML$model$col_receiver,
                         ij_variables = c(Xtilde_names, "tmp_Y", "Ystar", "p", "p1", "omega"))
    wide[, r_ij := compute_r(Ystar_ij, Ystar_ji, rho)]
    wide[, `:=`(rho_tilde_ij = (r_ij - p_ij*p_ji) / sqrt(p1_ij * p1_ji))]
    wide[, `:=`(w_ij = sqrt(omega_ij * omega_ji) * rho_tilde_ij)]

    var_ij_names <- sprintf("var_ij_%s", Xtilde_names)

    m_ij <- as.matrix(wide[, paste0(Xtilde_names, '_ij'), with = FALSE])
    m_ji <- as.matrix(wide[, paste0(Xtilde_names, '_ji'), with = FALSE])

    W2 <- W1 + weighted_inner_prod(m_ij, m_ji, wide$w_ij) / N2
  } else
    stop("Incorrect `estimator' option for variance of theta.")

  variance <- W1_inv %*% W2 %*% W1_inv / N2

  tstat <- (theta_hat_bias_corr - theta0) / sqrt(diag(variance))
  pval <- 2 * pnorm(abs(tstat), lower.tail = FALSE)

  c(fit_ML, list(bias_theta = bias, var_theta = variance,
                 theta_hat_bias_corr = theta_hat_bias_corr, tstat_theta = tstat,
                 pval_theta = pval, W1_inv = W1_inv, B_theta = B_theta))
}
