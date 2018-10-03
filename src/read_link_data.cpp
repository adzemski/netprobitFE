// This R package implements the methods proposed in
// Dzemski, Andreas : An empirical model of dyadic link formation in
// a network with unobserved heterogeneity, Review of Economics and Statistics, forthcoming
//
// Copyright(C) 2018 Andreas Dzemski

// This program is free software : you can redistribute it and / or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <https://www.gnu.org/licenses/>.

#include <Rcpp.h>
#include "LinkVariable.h"

using namespace Rcpp;

//' use C++ class LinkVariable to establish link to R data frame
//' 
//' this function is mainly for testing to see if the link is 
//' established correctly 
//' 
//' @param df data frame with link data (ordered canonically)
//' @param i Sender index
//' @param j Receiver index
//' @param var_name name of variabe in data frame 
//' 
//' @return var_name at i, j
//' @export
// [[Rcpp::export]]
double read_link_data(DataFrame df, int i, int j, String var_name) {
  LinkVariable link_var(df, var_name);
  return link_var.get(i, j);
}
