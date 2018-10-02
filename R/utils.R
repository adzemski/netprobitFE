#' Compute a weighted inner product
#'
#' @param A matrix (n x k)
#' @param B matrix (n x k)
#' @param weights vector of weights (1 x n)
#'
#' @return inner product weighted by w
weighted_inner_prod <- function(A, B, weights) {
  B <- matrix(rep(weights, times = ncol(B)), ncol = ncol(B))  * B
  t(A) %*% B
}

#' compute number of (ordered) pairs from 1, ..., N
#'
#' @param N number of elements
#'
#' @return number of ordered pairs
num_pairs <- function(N) {
  N*(N - 1)
}

#' Center a vector
#'
#' @param x vector
#'
#' @return centered vector
center <- function(x) {
  x - mean(x)
}

#' Return vector unchanged
#'
#' @param x vector
#'
#' @return x
identity <- function(x) {
  x
}

#' Add latent index to dataset for given structural parameters.
#'
#' @param model list that describes the model specification
#' @param gammaS vector of sender effects
#' @param gammaR vector of receiver effects
#' @param theta homophily parameter
#' @param long data.table
#' @param trimming_bound specify trimming sequence for probabilities
#'
#' @return ystar (=latent index) is added by reference, a view on the data.table is returned
#' @export
#'
#' @import data.table
add_ystar <- function(model, theta, gammaS, gammaR, long = model$links, trimming_bound = 0) {
  long[, Ystar := compute_ystar(long, theta, gammaS, gammaR, model$X_names,
                                model$col_sender, model$col_receiver, model$N, trimming_bound)]
}

#' Compute ystar
#'
#' @param links data.table with link data
#' @param theta theta
#' @param gammaS gammaS
#' @param gammaR gammaR
#' @param X_names column(s) of regressors
#' @param col_sender column that gives sender
#' @param col_receiver column that gives receiver
#' @param num_nodes number of nodes
#' @param trimming_bound specify trimming sequence for probabilities
#'
#' @return ystar vector
compute_ystar <- function(links, theta, gammaS, gammaR, X_names,
                          col_sender, col_receiver, num_nodes, trimming_bound = 0) {

  ystar <- as.vector(as.matrix(links[ , X_names, with = FALSE])%*% theta) +
    gammaS[links[[col_sender]]] + gammaR[links[[col_receiver]]]

  if (trimming_bound > 0) {
    qn <- qnorm(trimming_bound/num_nodes)
    lower_bound <- min(qn, -qn)
    upper_bound <- max(qn, -qn)
    ystar <- pmin(pmax(ystar, lower_bound), upper_bound)
  }

  ystar
}

#' add xtilde's by reference
#'
#' @param model list that describes the model specification
#' @param long data.table
#'
#' @import data.table
add_xtilde <- function(model, long, force_robust = FALSE) {

  for (x_name in model$X_names) {
    add_predict(model, long, x_name, sprintf("Xtilde_%s", x_name), force_robust = force_robust)
  }
}

#' add residuals (projected out of fixed effects) of a variable
#'
#' @param model list that describes the model specification
#' @param long data.table
#' @param var_name name of variable
#' @param new_var_name name of new variable with residuals
#'
#' @import data.table
add_predict <- function(model, long, var_name, new_var_name,
                          var_name_weights = "omega", force_robust = FALSE, residuals = TRUE) {

  add_predict_robust <-
    function(model, long, var_name, new_var_name, var_name_weights, residuals) {

      # if (any(is.nan(long[[var_name_weights]])) | any(is.infinite(long[[var_name_weights]])) |
      #     any(is.nan(long[[var_name]])) | any(is.infinite(long[[var_name]])) |
      #     any(long[[var_name_weights]] < 0))
      #   browser()

      if (residuals)
        long[, (new_var_name) :=
               residuals(lm.wfit(as.matrix(model$design_matrix_FE), long[[var_name]],
                                 long[[var_name_weights]]))]
      else
        long[, (new_var_name) :=
               fitted.values(lm.wfit(as.matrix(model$design_matrix_FE), long[[var_name]],
                                     long[[var_name_weights]]))]
    }

  if (!(var_name_weights %in% colnames(long))) stop("data.table lacking column with weights.")

  if (force_robust) {
    add_predict_robust(model, long, var_name, new_var_name, var_name_weights, residuals)
  } else {
    tryCatch({
      proj <-
        speedlm.wfit(long[[var_name]], model$design_matrix_FE, long[[var_name_weights]],
                     sparse = TRUE)
      fitted_values <- as.vector(model$design_matrix_FE %*% coef(proj))
      if (residuals)
        long[, (new_var_name) := (long[[var_name]] - fitted_values)]
      else
        long[, (new_var_name) := fitted_values]
    }, error = function(e) {
      warning(sprintf("Speedlm aborted (error: %s). Using fallback method.", e$message))
      add_predict_robust(model, long, var_name, new_var_name, var_name_weights, residuals)
    })
  }
}
