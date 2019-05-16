//
// Created by Eric Chen on 2019-05-16.
//

#include <iostream>
#include "simplex_utils.h"

void read_and_solve() {
    int n{};
    int m{};
    std::cin >> n >> m;
    matrix A(n, std::vector<double>(m));
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < m; ++j) {
            std::cin >>A[i][j];
        }
    }
    std::vector<double> b(n);
    for (int i = 0; i < n; ++i) {
        std::cin >> b[i];
    }
    std::vector<double> c(m);
    for (int i = 0; i < m; ++i) {
        std::cin >> c[i];
    }

    print_problem(n, m, A, b, c);

    SimplexAlgorithm simplex_algorithm{n, m, A, b, c};
    simplex_algorithm.print_solution();


}

void print_problem(int n, int m, const matrix &A, const std::vector<double> &b, const std::vector<double> &c) {
    std::cout << "The linear program in standard form is:\n";
    std::cout << "Maximize\n";
    for (int i = 0; i < m; ++i) {
        std::cout << ((i == 0) ? std::noshowpos : std::showpos) << c[i] << " * x" << std::noshowpos << i
                  << ((i == m - 1) ? "\n" : " ");
    }
    std::cout << "Subject to\n";
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < m; ++j) {
            std::cout << ((j == 0) ? std::noshowpos : std::showpos) << A[i][j] << " * x"
                      << std::noshowpos << j << " ";
        }
        std::cout << "<= " << b[i] << "\n";
    }
    for (int i = 0; i < m; ++i) {
        std::cout << "x" << i << ((i != m - 1) ? ", " : " >= 0\n");
    }

    std::cout << "------------------------------------------------------------------\n";
}
