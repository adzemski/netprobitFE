#ifndef LINK_VARIABLE_H
#define LINK_VARIABLE_H

#include <Rcpp.h>
#include <string>

class LinkVariable
{
private:
    Rcpp::NumericVector column; 
    int n;
public:
    LinkVariable (Rcpp::DataFrame df, std::string column_name);
    double  get (int i, int j);
    long  get_index (int i, int j);
    int nnodes (); 
};

#endif