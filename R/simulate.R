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

#' Run simulations
#'
#' @param sim_design Simulation design in list format
#' @param draw_network function that draws a network (takes argument sim_design)
#' @param bs_iterations number of bootstrap iterations for triangle test
#' @param stage2_method either "mle" or "gradient"
#' @param num_simulations number of simulations
#' @param file_name file name for saving results (relative to path)
#'  if not provided then results are not saved to disk
#' @param verbose print more diagnostic output
#' @param bs_trimming_bound specify trimming of probabilities in bootstrap Monte Carlo
#' @param force_robust don't use speedglm (for robustness in sparse networks)
#' @param bs_force_robust don't use speedglm in bootstrap Monte Carlo
#'
#' @return data frame with simulation results
#' @export
simulate_methods <- function(sim_design, draw_network, num_simulations = 100,
                             bs_iterations = 100,
                             stage2_method = "gradient",
                             triangles_bias_est = "theorem2",
                             triangles_bs_inf = "none",
                             file_name = NULL,
                             verbose = TRUE, bs_trimming_bound = 0,
                             force_robust = FALSE, bs_force_robust = force_robust) {

  options(nwarnings = 1e4)

  if (is.null(names(sim_design$theta)))
    names(sim_design$theta) <- sprintf("X%d", seq_along(sim_design$theta))

  simulate_once <- function(s) {

    if (verbose) cat(sprintf("Starting simulation %d ... ", s))

    simulated_links <- draw_network(sim_design)$links
    simulated_links[, `:=`(Ystar0 = Ystar, Ystar = NULL)]

    df_degree <- data.frame(density = simulated_links[, mean(Y)])

    fit_stage1 <- MLE_stage1(define_model(simulated_links, col_sender = "i", col_receiver = "j",
                                          X_names = names(sim_design$theta)))

    fit_stage2 <- MLE_stage2(fit_stage1, method = stage2_method, rho_start = sim_design$rho)

    fit_ttest_theta <- ttest_theta(fit_stage2, theta0 = sim_design$theta,
                                   force_robust = force_robust)

    fit_ttest_rho <- ttest_rho(fit_stage2, fit_ttest_theta, rho0 = sim_design$rho,
                               force_robust = force_robust)

    df_theta <- as.data.frame(c(setNames(as.list(fit_stage1$theta_hat_MLE), sprintf("theta_%s_ML", names(sim_design$theta))),
    setNames(as.list(fit_ttest_theta$theta_hat_bias_corr), sprintf("theta_%s_bc", names(sim_design$theta))),
    setNames(as.list(fit_ttest_theta$pval_theta), sprintf("theta_%s_pval", names(sim_design$theta))) ))

    df_rho <- data.frame(rho_hat_MLE = fit_stage2$rho_hat_MLE,
                         rho_hat_bias_corr = fit_ttest_rho$rho_hat_bias_corr,
                         rho_pval = fit_ttest_rho$pval_rho)
    fit_triangles <-
      triangle_test(fit_stage2, fit_ttest_theta,
                    bs_iterations = bs_iterations,
                    bias_estimator = triangles_bias_est, bs_inference = triangles_bs_inf,
                    bs_trimming_bound = bs_trimming_bound,
                    force_robust = force_robust, bs_force_robust = bs_force_robust,
                    compute_other_variances = TRUE)

    simulated_links[, p := compute_p(Ystar0)]
    expected_trans_oracle <- fast_excess_transitivity(simulated_links, "Y", "p")$expected_triangles

    # fit_triangles$naive_teststat <- with(fit_triangles, excess_trans / se_oracle)
    # fit_triangles$naive_pval <- 2 * pnorm(abs(fit_triangles$naive_teststat), lower.tail = FALSE)
    #
    # fit_triangles$naive_bc_teststat <- with(fit_triangles, excess_trans_bias_corr / se_oracle)
    # fit_triangles$naive_bc_pval <- 2 * pnorm(abs(fit_triangles$naive_bc_teststat), lower.tail = FALSE)

    fit_triangles$excess_trans_oracle <-
      (fit_triangles$obs_trans - expected_trans_oracle)/(sim_design$N^3)

    fit_triangles$oracle_teststat <- fit_triangles$excess_trans_oracle/fit_triangles$se_oracle
    fit_triangles$oracle_pval <- 2 * pnorm(abs(fit_triangles$oracle_teststat), lower.tail = FALSE)

    fit_triangles$ttest_teststat <- with(fit_triangles,
      (expected_trans - bias - expected_trans_oracle)/
      (sim_design$N^3)/se_predicted_triangles)

    fit_triangles$ttest_pval <- 2 * pnorm(abs(fit_triangles$ttest_teststat), lower.tail = FALSE)

    triangles_show_variables <- c("teststat", "pval_analytic", "se_analytic", "bias",
                                  "excess_trans", "excess_trans_bias_corr",
                                  "rej_bs_perc_t", "rej_bs_perc_excess",
                                  "excess_trans_oracle", "oracle_teststat", "oracle_pval",
                                  "ttest_teststat", "ttest_pval",
                                  "naive_teststat", "naive_pval", "naive_bc_teststat", "naive_bc_pval")

    df_triangles <- fit_triangles[triangles_show_variables]
    names(df_triangles) <- paste0("triangles_", triangles_show_variables)
    df_triangles <- lapply(df_triangles, function(x) if (is.null(x)) NA_real_ else x)
    df_triangles <- data.frame(df_triangles)

    if (verbose) cat("done.\n")

    cbind(df_degree, df_theta, df_rho, df_triangles)
  }

  save_simulate_once <- function(s) {
    tryCatch(simulate_once(s),
             error = function(e) {
               cat("\n")
               warning(sprintf("Error in simulation %s (%s)", s, e$message))
               e$origin <- sprintf("simulation %s", s)
               e
             })
  }

  all_sim <- lapply(seq_len(num_simulations), save_simulate_once)

  results <-
    do.call("rbind", lapply(all_sim, function(x) if (is.data.frame(x)) x else NULL))

  try( # call can go wrong if there are no results (only errors)
    rownames(results) <-
      sprintf("sim %s", seq_len(num_simulations))[sapply(all_sim, function(x) is.data.frame(x))]
  )

  errors <-
    lapply(all_sim, function(x) if (!is.data.frame(x)) x else NULL)

  if (is.character(file_name)) {
    save(results, errors, file = file_name)
    if (verbose) cat(sprintf("Saved simulation results to %s.\n", file_name))
  }

  list(results = results, errors = errors)
}

#' Sample a random network following modified Jochmans design
#'
#' @param design network parameters
#' @import data.table
#' @importFrom MASS mvrnorm
#' @importFrom utils combn
#' @return list containing the generated data and fixed effects
#' @export
draw_network_jochmans_2018 <- function(design) {

  units <- data.frame(unit_index = seq_len(design$N))

  # X = unit index is even
  units$X1 <- as.integer(1 - 2 * units$unit_index %% 2)

  units$gammaS <- - (design$N - units$unit_index)/(design$N - 1) * design$C_N
  units$gammaR <- units$gammaS

  dyads <- combn(design$N, 2)
  shocks <-
    mvrnorm(design$N * (design$N - 1) / 2, mu = c(0,0),
            Sigma = matrix( c(1, design$rho, design$rho, 1), ncol = 2) )

  wide <- data.table(i = dyads[1, ], j = dyads[2, ], U_ij = shocks[, 1], U_ji = shocks[, 2])
  wide[, X1_ij := units$X1[i] * units$X1[j]]
  wide[, X1_ji := X1_ij]

  # translate to long format
  long <- wide[, .(i = c(i, j), j = c(j, i), X1 = c(X1_ij, X1_ji), U = c(U_ij, U_ji))]
  setkey(long, i, j)

  model <- list(col_sender = "i", col_receiver = "j", X_names = "X1")
  add_ystar(model, design$theta, units$gammaS, units$gammaR, long)
  # long[, U := stats::rlogis(nrow(long))]
  long[, Y := as.integer(U <= Ystar)]

  list(links = long)
}

#' Sample from alternative with dynacmic behavior
#'
#' @param design network parameters
#' @import data.table
#' @return list containing the generated data and fixed effects
#' @export
draw_j18_dynamic_meeting <- function(design) {

  if (!("n_iterations" %in% names(design))) design$n_iterations = 1

  links <- draw_network_jochmans_2018(design)$links

  for (iter in seq_len(design$n_iterations)) {

    links[, support := num_trans_closures(links, "Y")]
    links[, support := support > 0]
    links[, missing_link := (support > 0) & !Y]
    links[, in_triangle := num_in_triangles(links, "Y")]
    links[, weak_link := (in_triangle == 0) & Y]
    links[, U_new := rnorm(nrow(links))]
    links[, U_new := missing_link * pmin(U_new, U) + weak_link * pmax(U_new, U)]

    links[, Y := as.integer(U_new <= Ystar)]
  }

  list(links = links)
}

#' Define a simulation design (Jochmans)
#'
#' @param N number of nodes
#' @param theta coefficient for parametric part
#' @param rho within dyad correlation
#' @param type_C_N sparsity parameter
#'  * "0"
#'  * "loglog"
#'  * "sqrt_log"
#'  * "log"
#' @param S number of simulations
#'
#' @return simulation design as list
#' @export
#'
define_sim_design_jochmans <- function(N, theta, rho, type_C_N = "0") {

  get_C_N <- function(N, type) {
    if (type == "0")
      return(0)
    if (type == "loglog")
      return(log(log(N)))
    if (type == "sqrt_log")
      return(sqrt(log(N)))
    if (type == "sqrt_2log")
      return(sqrt(2 * log(N)))
    if (type == "log")
      return(log(N))

    stop("Invalid type_C_N")
  }

  if (is.null(names(theta))) names(theta) <- sprintf("X%d", seq_along(theta))

  list(N = N, theta = theta, rho = rho, C_N = get_C_N(N, type_C_N))
}

