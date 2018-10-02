#include <Rcpp.h>
#include "LinkVariable.h"

using namespace Rcpp;

//' add new link variables that are required for bias correction
//' 
//' assumes that link observations are ordered canonically and
//' that data frame contains the variables "p", "H", "partial_p"
//' 
//' @param df data frame with link data (ordered canonically)
//' 
//' @return data frame with new variables
// [[Rcpp::export]]
DataFrame triangles_biascorr_new_vars(DataFrame df)
{
    LinkVariable p (df, "p_ij");
    LinkVariable partial_p (df, "p_ystar_ij");
    LinkVariable H (df, "H_ij");
    unsigned long N = df.nrows (); 
    NumericVector boldbeta (N);
    NumericVector bss (N), bsr (N), bssr (N);
    unsigned long index = 0;
    unsigned int n = p.nnodes ();

    double accumulated_boldbeta;
    double accumulated_bss, accumulated_bsr, accumulated_bssr;
    double p_ik, p_jk, p_ki, p_kj, partial_p_ij; 

    for (int i = 1; i <= n; ++i) 
    {
        for (int j = 1; j <= n; ++j)
        {
            if (i == j) continue;
            accumulated_boldbeta = 0; accumulated_bss = 0; 
            accumulated_bsr = 0; accumulated_bssr = 0; 
            for (int k = 1; k <= n; ++k)
            {
                if (k == i) continue;
                if (k == j) continue;
                p_ik = p.get(i, k);
                p_jk = p.get(j, k);
                p_ki = p.get(k, i);
                p_kj = p.get(k, j);
                accumulated_boldbeta += p_ik * p_jk +
                                        p_ik * p_kj + p_ki * p_kj;
                accumulated_bss += partial_p.get (i, k) * (p_jk + p_kj); 
                accumulated_bsr += partial_p.get (k, j) * (p_ik + p_ki);  
                accumulated_bssr += partial_p.get (k, i) * p_kj;  
            }
            partial_p_ij = partial_p.get (i, j); 
            boldbeta[index] = accumulated_boldbeta / H.get(i, j) / n;
            bss[index] = accumulated_bss * partial_p_ij / n;
            bsr[index] = accumulated_bsr * partial_p_ij / n;
            bssr[index] = accumulated_bssr * partial_p_ij / n;
            ++index;
        }
    }
    return DataFrame::create(_["i"] = df["i"], _["j"] = df["j"], _["boldbeta"] = boldbeta,
                             _["bss"] = bss, _["bsr"] = bsr, _["bssr"] = bssr);
}