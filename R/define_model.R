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

#' Define a network model
#'
#' @param d data.frame with link data
#' @param X_names covariate names (NULL = no covariates)
#' @param col_sender name of column that identifies sender
#' @param col_receiver name of column that identifies receiver
#'
#' @return list describing the network model and data
#' @export
#'
#' @import data.table
#' @importFrom Matrix Matrix
define_model <- function(d, X_names = NULL, col_sender = "i", col_receiver = "j",
                         col_outcome = "Y") {

  d <- as.data.frame(d)

  # does data.frame contain necessary information
  var_not_in_data <-
    c(X_names, col_sender, col_receiver, col_outcome)[!(c(X_names, col_sender, col_receiver, col_outcome) %in% names(d))]

  if (length(var_not_in_data) > 0) {
    stop(sprintf("data is missing some variables (%s)", paste(var_not_in_data, collapse = ", ")))
  }

  # which variables store sender / receiver of link
  model <- list(col_sender = col_sender, col_receiver = col_receiver, X_names = X_names,
                col_outcome = col_outcome)

  # get node names
  model$node_names <- as.character(unique(d[[model$col_sender]]))

  # get number of nodes
  model$N <- length(model$node_names)

  # save original names
  d[[sprintf("%s_original", col_sender)]] <- as.character(d[[col_sender]])
  d[[sprintf("%s_original", col_receiver)]] <- as.character(d[[col_receiver]])

  # enumerate sender/receiver
  translation_table <- setNames(seq_len(model$N), model$node_names)
  d[[col_sender]] <- translation_table[d[[col_sender]]]
  d[[col_receiver]] <- translation_table[d[[col_receiver]]]

  # order by new sender/receiver index
  d <- d[order(d[[col_sender]], d[[col_receiver]]), ]

  # set-up the design matrix for the projections
  FE <- cbind(do.call(cbind, lapply(seq_len(model$N),
                                    FUN = function(x) as.numeric(d[[model$col_sender]]  == x))),
              do.call(cbind, lapply(seq_len(model$N),
                                    FUN = function(x) as.numeric(d[[model$col_receiver]]  == x))))

  FE[d[[model$col_sender]] == 1, ] <-
    FE[d[[model$col_sender]] == 1, ] +
    matrix(rep.int(c(rep.int(-1, times = model$N), rep.int(1, times = model$N)),
                   times = model$N - 1), ncol = 2 * model$N, byrow = TRUE)

  colnames(FE) <- c(sprintf("gammaS_%d", seq_len(model$N)),
                    sprintf("gammaR_%d", seq_len(model$N)))

  FE <- FE[ , -1]

  model$design_matrix_FE <- Matrix(FE, sparse = TRUE)

  X_matrix <- as.matrix(d[, model$X_names])
  colnames(X_matrix) <- model$X_names

  model$design_matrix_MLE <- Matrix(cbind2(X_matrix, model$design_matrix_FE), sparse = TRUE)

  # add link data and return model
  model$links <- data.table(d)
  model
}
