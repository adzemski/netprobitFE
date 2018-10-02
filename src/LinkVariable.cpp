#include <Rcpp.h>
#include <math.h>
#include "LinkVariable.h"

using namespace Rcpp;

LinkVariable::LinkVariable (DataFrame df, std::string column_name) {
    if (!df.containsElementNamed (column_name.c_str ()))
        throw "requested column not in dataframe";
    int N = df.nrows ();
    this->column = df [column_name.c_str ()];
    this->n = 0.5 + sqrt (N + 0.25);
}

int LinkVariable::nnodes ()
{
    return n;
}

long LinkVariable::get_index(int i, int j)
{
    try
    {
        if (i == j)
            throw "invalid index (need i != j)";
        if (i < 1 || i > n)
            throw "i index out of bounds";

        // int index = (i - 1) * (n - (i / 2)) + j - 1;
        int index = (i - 1) * (n - 1) + j - 1;
        if (i < j)
            --index;
        return index;
    }
    catch (const char *exception)
    {
        Rcerr << "Error: " << exception << std::endl;
        return -1;
    }
}

double LinkVariable::get (int i, int j)
{
    long index; 
    index = get_index (i, j); 

    if (index < 0) 
        throw "computed invalid index";
    
    return column[index];
}
