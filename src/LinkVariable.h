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