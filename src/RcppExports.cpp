// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// excess_trans
List excess_trans(DataFrame df, std::string y_variable, std::string p_variable);
RcppExport SEXP _netprobitFE_excess_trans(SEXP dfSEXP, SEXP y_variableSEXP, SEXP p_variableSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< DataFrame >::type df(dfSEXP);
    Rcpp::traits::input_parameter< std::string >::type y_variable(y_variableSEXP);
    Rcpp::traits::input_parameter< std::string >::type p_variable(p_variableSEXP);
    rcpp_result_gen = Rcpp::wrap(excess_trans(df, y_variable, p_variable));
    return rcpp_result_gen;
END_RCPP
}
// read_link_data
double read_link_data(DataFrame df, int i, int j, String var_name);
RcppExport SEXP _netprobitFE_read_link_data(SEXP dfSEXP, SEXP iSEXP, SEXP jSEXP, SEXP var_nameSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< DataFrame >::type df(dfSEXP);
    Rcpp::traits::input_parameter< int >::type i(iSEXP);
    Rcpp::traits::input_parameter< int >::type j(jSEXP);
    Rcpp::traits::input_parameter< String >::type var_name(var_nameSEXP);
    rcpp_result_gen = Rcpp::wrap(read_link_data(df, i, j, var_name));
    return rcpp_result_gen;
END_RCPP
}
// num_trans_closures
NumericVector num_trans_closures(DataFrame df, std::string y_variable);
RcppExport SEXP _netprobitFE_num_trans_closures(SEXP dfSEXP, SEXP y_variableSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< DataFrame >::type df(dfSEXP);
    Rcpp::traits::input_parameter< std::string >::type y_variable(y_variableSEXP);
    rcpp_result_gen = Rcpp::wrap(num_trans_closures(df, y_variable));
    return rcpp_result_gen;
END_RCPP
}
// num_in_triangles
NumericVector num_in_triangles(DataFrame df, std::string y_variable);
RcppExport SEXP _netprobitFE_num_in_triangles(SEXP dfSEXP, SEXP y_variableSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< DataFrame >::type df(dfSEXP);
    Rcpp::traits::input_parameter< std::string >::type y_variable(y_variableSEXP);
    rcpp_result_gen = Rcpp::wrap(num_in_triangles(df, y_variable));
    return rcpp_result_gen;
END_RCPP
}
// triangles_biascorr_new_vars
DataFrame triangles_biascorr_new_vars(DataFrame df);
RcppExport SEXP _netprobitFE_triangles_biascorr_new_vars(SEXP dfSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< DataFrame >::type df(dfSEXP);
    rcpp_result_gen = Rcpp::wrap(triangles_biascorr_new_vars(df));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_netprobitFE_excess_trans", (DL_FUNC) &_netprobitFE_excess_trans, 3},
    {"_netprobitFE_read_link_data", (DL_FUNC) &_netprobitFE_read_link_data, 4},
    {"_netprobitFE_num_trans_closures", (DL_FUNC) &_netprobitFE_num_trans_closures, 2},
    {"_netprobitFE_num_in_triangles", (DL_FUNC) &_netprobitFE_num_in_triangles, 2},
    {"_netprobitFE_triangles_biascorr_new_vars", (DL_FUNC) &_netprobitFE_triangles_biascorr_new_vars, 1},
    {NULL, NULL, 0}
};

RcppExport void R_init_netprobitFE(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}