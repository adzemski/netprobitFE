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
