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

#' Compute probability of reciprocated link
#'
#' @param ystar1 ystar for ij link
#' @param ystar2 ystar for ji link
#' @param rho correlation (reciprocity parameter)
#'
#' @return vector of probabilities
#'
#' @importFrom mnormt pmnorm
compute_r <- function(ystar1, ystar2, rho) {
  evaluation_points <- matrix(c(ystar1, ystar2), nrow = length(ystar1), ncol = 2)
  pmnorm(evaluation_points, varcov = matrix(c(1, rho, rho, 1), 2, 2))
}

#' Compute J (see paper for definition)
#'
#' @param ystar1 ystar for ij link
#' @param ystar2 ystar for ji link
#' @param rho correlation (reciprocity parameter)
#' @param r12 precomputed dyad probability can be supplied for efficiency gains
#'
#' @return vector
compute_J <- Vectorize(function(ystar1, ystar2, rho, r12 = NULL) {

    if (is.null(r12)) r12 <- compute_r(ystar1, ystar2, rho)

    partial_r_rho(ystar1, ystar2, rho)/(r12*(1 - r12))
  }
  , vectorize.args = c('ystar1', 'ystar2', 'r12')
)

#' Compute derivative of r wrt rho (see paper for definition)
#'
#' @param ystar1 ystar for ij link
#' @param ystar2 ystar for ji link
#' @param rho correlation (reciprocity parameter)
#'
#' @return vector
partial_r_rho <- Vectorize(
  function(ystar1, ystar2, rho) {

    if (rho == 0)
      return(dnorm(ystar1) * dnorm(ystar2))

    rho1 <- sqrt(1 - rho^2)
    rho3 <- rho1^3

    integrand <- function(t, rho, ystar2) {
      (rho * ystar2 - t) / rho3 * dnorm((ystar2 - rho * t)/rho1) * dnorm(t)
    }

    integrate(function(t) integrand(t, rho, ystar2), lower = -Inf, upper = ystar1, rel.tol = .Machine$double.eps^0.6)$value
  }
  , vectorize.args = c('ystar1', 'ystar2')
)

#' Compute derivative of r wrt rho
#'
#' @param ystar1 ystar for ij link
#' @param ystar2 ystar for ji link
#' @param rho correlation (reciprocity parameter)
#'
#' @return vector
partial_r_ystar1 <- Vectorize(
  function(ystar1, ystar2, rho) {

    pnorm((ystar2 - rho*ystar1)/sqrt(1 - rho^2)) * dnorm(ystar1)
  }
  , vectorize.args = c('ystar1', 'ystar2')
)

#' Compute 2nd order derivative of r wrt ystar1 and rho
#'
#' @param ystar1 ystar for ij link
#' @param ystar2 ystar for ji link
#' @param rho correlation (reciprocity parameter)
#'
#' @return vector
partial_r_ystar1_rho <- Vectorize(
  function(ystar1, ystar2, rho) {

    - (ystar1 - rho * ystar2)/((1 - rho^2)^(1.5)) *
      dnorm((ystar2 - rho*ystar1)/sqrt(1 - rho^2)) * dnorm(ystar1)

  }
  , vectorize.args = c('ystar1', 'ystar2')
)

#' Compute 2nd order derivative of r wrt ystar1 and ystar2
#'
#' @param ystar1 ystar for ij link
#' @param ystar2 ystar for ji link
#' @param rho correlation (reciprocity parameter)
#'
#' @return vector
partial_r_ystar1_ystar2 <- Vectorize(
  function(ystar1, ystar2, rho) {

    (1 - rho^2)^(-0.5) * dnorm((ystar2 - rho*ystar1)/sqrt(1 - rho^2)) * dnorm(ystar1)

  }
  , vectorize.args = c('ystar1', 'ystar2')
)

#' Compute 2nd order derivative of r wrt ystar1 and ystar1
#'
#' @param ystar1 ystar for ij link
#' @param ystar2 ystar for ji link
#' @param rho correlation (reciprocity parameter)
#'
#' @return vector
partial_r_ystar1_ystar1 <- Vectorize(
  function(ystar1, ystar2, rho) {

    - rho/sqrt(1 - rho^2) * dnorm((ystar2 - rho*ystar1)/sqrt(1 - rho^2)) * dnorm(ystar1) -
      ystar1 * pnorm((ystar2 - rho*ystar1)/sqrt(1 - rho^2)) * dnorm(ystar1)
  }
  , vectorize.args = c('ystar1', 'ystar2')
)

#' Compute derivative of J (see paper) wrt ystar1
#'
#' @param ystar1 ystar for ij link
#' @param ystar2 ystar for ji link
#' @param rho correlation (reciprocity parameter)
#'
#' @return vector
partial_J_ystar1 <-  Vectorize(
  function(ystar1, ystar2, rho, r = NULL, r_rho = NULL, r_ystar1 = NULL, r_ystar1_rho = NULL) {

    # first compute all terms that have not been pre-computed
    if (is.null(r)) r <- fcn_r(ystar1, ystar2, rho)
    if (is.null(r_ystar1_rho)) r_ystar1_rho <- partial_r_ystar1_rho(ystar1, ystar2, rho)
    if (is.null(r_rho)) r_rho <- partial_r_rho(ystar1, ystar2, rho)
    if (is.null(r_ystar1)) r_ystar1 <- partial_r_ystar1(ystar1, ystar2, rho)

    1/(r * (1 - r)) * ( r_ystar1_rho - r_rho * 1/(r * (1 - r)) * (1 - 2 * r) * r_ystar1)
  }
  , vectorize.args = c('ystar1', 'ystar2', 'r', 'r_rho', 'r_ystar1', 'r_ystar1_rho')
)

#' Compute omega for ij link (definition see paper)
#'
#' @param ystar ystar for ij link
#' @param p pre-computed linking probability
#'
#' @return vector
compute_omega <- function(ystar, p = compute_p(ystar)) {
  omega <- (dnorm(ystar))^2 / compute_p1(ystar, p)
  omega[!is.finite(omega)] <- 0
  omega
}

#' Compute H for ij link (definition see paper)
#'
#' @param ystar ystar for ij link
#' @param p pre-computed linking probability
#'
#' @return vector
compute_H <- function(ystar, p = compute_p(ystar)) {
  dnorm(ystar) / compute_p1(ystar, p)
}

#' Compute linking probability p for ij link (definition see paper)
#'
#' @param ystar ystar for ij link
#'
#' @return vector
compute_p <- function(ystar) {
  pnorm(ystar)
}

#' Compute p(1 - p) for ij link (definition see paper)
#'
#' @param ystar ystar for ij link
#' @param p pre-computed linking probability
#'
#' @return vector
compute_p1 <- function(ystar, p = compute_p(ystar)) {
  p*(1 - p)
}

#' Compute derivative of p (linking probability) wrt to ystar
#'
#' @param ystar index
#'
#' @return vector of same length as ystar
partial_p_ystar <- function(ystar) {
  dnorm(ystar)
}

#' Compute 2nd order derivative of p (linking probability) wrt to ystar
#'
#' @param ystar index
#'
#' @return vector of same length as ystar
partial_p_ystar_ystar <- function(ystar) {
  - ystar * dnorm(ystar)
}
