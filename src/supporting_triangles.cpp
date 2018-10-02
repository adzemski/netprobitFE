#include <Rcpp.h>
#include "LinkVariable.h"
#include <math.h>

using namespace Rcpp;

//' Detect number of transitive relationships closed by this link
//'
//' assumes that link observations are ordered canonically a
//'
//' @param df data frame with link data (ordered canonically)
//' @param y_variable string that contains name of column that stores link indicator
//'
//' @export
//' @return number of supporting transitive triangles for this link
// [[Rcpp::export]]
NumericVector num_trans_closures (DataFrame df, std::string y_variable)
{

    LinkVariable y(df, y_variable);
    int n = y.nnodes();
    NumericVector link(n * (n - 1));

    for (int i = 1; i <= n; ++i)
    {
        for (int j = i + 1; j <= n; ++j)
        {
            for (int k = j + 1; k <= n; ++k)
            {
                link[y.get_index (i, k)] += static_cast<int>(y.get (i, j)) &&
                  static_cast<int>(y.get (j, k));
                link[y.get_index (j, k)] += static_cast<int>(y.get (j, i)) &&
                  static_cast<int>(y.get (i, k));
                link[y.get_index (j, i)] += static_cast<int>(y.get (j, k)) &&
                  static_cast<int>(y.get (k, i));
                link[y.get_index (k, i)] += static_cast<int>(y.get (k, j)) &&
                  static_cast<int>(y.get (j, i));
                link[y.get_index (i, j)] += static_cast<int>(y.get (i, k)) &&
                  static_cast<int>(y.get (k, j));
                link[y.get_index (k, j)] += static_cast<int>(y.get (k, i)) &&
                  static_cast<int>(y.get (i, j));
            }
        }
    }

    return link;
}

//' Detect number of triangles that this link is a part of
//'
//' assumes that link observations are ordered canonically a
//'
//' @param df data frame with link data (ordered canonically)
//' @param y_variable string that contains name of column that stores link indicator
//'
//' @export
//' @return number of supporting transitive triangles for this link
// [[Rcpp::export]]
NumericVector num_in_triangles (DataFrame df, std::string y_variable)
{

  LinkVariable y(df, y_variable);
  int n = y.nnodes();
  NumericVector link(n * (n - 1));

  for (int i = 1; i <= n; ++i)
  {
    for (int j = 1; j <= n; ++j)
    {
      if (i == j) continue;
      for (int k = 1; k <= n; ++k)
      {
        if (i == k) continue;
        if (j == k) continue;
        
        link[y.get_index(i, j)] += 
          (static_cast<int>(y.get(i, k)) && static_cast<int>(y.get(k, j))) + 
          (static_cast<int>(y.get(i, k)) && static_cast<int>(y.get(j, k))) + 
          (static_cast<int>(y.get(k, i)) && static_cast<int>(y.get(k, j)));
      }
    }
  }

  return link;
}
