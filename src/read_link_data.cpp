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
