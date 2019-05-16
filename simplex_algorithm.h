//
// Created by Eric Chen on 2019-05-04.
//

#ifndef SIMPLEXSOLVER_SIMPLEX_ALGORITHM_H
#define SIMPLEXSOLVER_SIMPLEX_ALGORITHM_H

#include <vector>
#include <unordered_map>

typedef std::vector<std::vector<double>> matrix;
const double epsilon = 1e-6;


class SimplexAlgorithm {
public:
    SimplexAlgorithm(int n, int m, const matrix& A, const std::vector<double>& b, const std::vector<double>& c);
    void print_solution(int precision = 4);

private:
    int num_basics{}, num_nonbasics{};
    matrix solution_table{};
    std::vector<double> objective_function{};
    std::unordered_map<int, int> basic_to_row{};  // For keeping track of basic variables and their corresponding rows
    std::unordered_map<int, int> row_to_basic{};
    bool has_solution{false};
    bool infinite_solution{false};
    std::vector<double> solution{};

    void initialize(const matrix& A, const std::vector<double>& b, const std::vector<double>& c);
    void pivot(int pivot_row, int pivot_col);
    void solve();
    void print_solution_table() const;  // For debugging only

};



#endif //SIMPLEXSOLVER_SIMPLEX_ALGORITHM_H
