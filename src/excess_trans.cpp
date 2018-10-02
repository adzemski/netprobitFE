#include <Rcpp.h>
#include "LinkVariable.h"
#include <math.h>

using namespace Rcpp;

//' compute observed excess transitivity
//'
//' assumes that link observations are ordered canonically a
//'
//' @param df data frame with link data (ordered canonically)
//' @param y_variable string that contains name of column that stores link indicator
//' @param p_variable string that contains name of column that stores link probability
//'
//' @return number of expected transitive triangles
// [[Rcpp::export]]
List excess_trans (DataFrame df, std::string y_variable,
  std::string p_variable)
{

  double total_p = 0;
  unsigned long total_y = 0;
  double excess_trans = -1;

  try
  {
    LinkVariable y(df, y_variable);
    LinkVariable p(df, p_variable);
    unsigned int n = y.nnodes();
    int y_ij, y_ji;
    double p_ij, p_ji, p_ik, p_ki, p_jk, p_kj;

    // remember that R indexing starts at 1
    for (int i = 1; i <= n; ++i)
    {
      for (int j = i + 1; j <= n; ++j)
      {
        y_ij = y.get(i, j);
        p_ij = p.get(i, j);
        y_ji = y.get(j, i);
        p_ji = p.get(j, i);
        for (int k = j + 1; k <= n; ++k)
        {
          total_y += (y_ij && y.get(i, k) && y.get(j, k)) + \
                     (y_ij && y.get(i, k) && y.get(k, j)) + \
                     (y_ji && y.get(j, k) && y.get(i, k)) + \
                     (y_ji && y.get(j, k) && y.get(k, i)) + \
                     (y_ij && y.get(k, i) && y.get(k, j)) + \
                     (y_ji && y.get(k, i) && y.get(k, j));
          p_ik = p.get (i, k);
          p_ki = p.get (k, i);
          p_jk = p.get (j, k);
          p_kj = p.get (k, j);
          total_p += p_ij * p_ik * p_jk + \
                   p_ij * p_ik * p_kj + \
                   p_ji * p_jk * p_ik + \
                   p_ji * p_jk * p_ki + \
                   p_ij * p_ki * p_kj + \
                   p_ji * p_ki * p_kj;
        }
      }
    }
    excess_trans = (total_y - total_p) / pow(n, 3);
  }
  catch (const char *exception)
  {
    Rcerr << "Error: " << exception << std::endl;
    return List::create(_["error"] = LogicalVector::create(1));
  }

  return List::create(_["num_triangles"] = total_y, _["expected_triangles"] = total_p, \
  _["excess_trans"] = excess_trans, _["error"] = LogicalVector::create(0));
  ;
}
