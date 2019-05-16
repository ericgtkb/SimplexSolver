//
// Created by Eric Chen on 2019-05-16.
//

#ifndef SIMPLEXSOLVER_SIMPLEX_UTILS_H
#define SIMPLEXSOLVER_SIMPLEX_UTILS_H

#include "simplex_algorithm.h"

void read_and_solve();
void print_problem(int n, int m, const matrix& A, const std::vector<double>& b, const std::vector<double>&c);


#endif //SIMPLEXSOLVER_SIMPLEX_UTILS_H
