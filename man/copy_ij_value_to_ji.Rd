% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/reshape_dyadic_data.R
\name{copy_ij_value_to_ji}
\alias{copy_ij_value_to_ji}
\title{Make new variables that for entry i,j contains the j,i values of given variables}
\usage{
copy_ij_value_to_ji(wide, var_names)
}
\arguments{
\item{wide}{data table with link data}

\item{var_names}{names of variables to copy}
}
\value{
view on update data table
}
\description{
data table is updated by reference,
new variable names are generated by adding "_ji" subscript
}
