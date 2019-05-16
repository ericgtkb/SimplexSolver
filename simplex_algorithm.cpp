//
// Created by Eric Chen on 2019-05-04.
//

#include <iostream>
#include <iomanip>
#include "simplex_algorithm.h"

SimplexAlgorithm::SimplexAlgorithm(int n, int m, const matrix &A, const std::vector<double> &b, const std::vector<double> &c)
        : num_basics(n), num_nonbasics(m) {
    initialize(A, b, c);
}


void SimplexAlgorithm::initialize(const matrix &A, const std::vector<double> &b, const std::vector<double> &c) {
    int min_b_index = min_element(b.begin(), b.end()) - b.begin();
    double min_b = b[min_b_index];
    basic_to_row = std::unordered_map<int, int>();
    row_to_basic = std::unordered_map<int, int>();
    for (int i = 0; i < num_basics; ++i) {
        basic_to_row.emplace(i + num_nonbasics, i);
        row_to_basic.emplace(i, i + num_nonbasics);
    }
    if (min_b > 0) {
        solution_table = A;
        objective_function = c;
        for (size_t i = 0; i < num_basics + 1; ++i) {
            // num_basics + 1, last entry corresponds to the negative of the current basic solution
            objective_function.emplace_back(0);
        }
        for (size_t i = 0; i < num_basics; ++i) {
            for (size_t j = 0; j < num_basics; ++j) {
                if (i == j) {
                    solution_table[i].emplace_back(1);
                } else {
                    solution_table[i].emplace_back(0);
                }
            }
            solution_table[i].emplace_back(b[i]);
        }
        has_solution = true;

    } else {
        // For this case we solve the auxiliary linear program, with one additional nonbasic variable
        solution_table = A;
        // second last element is the new variable, last element is the current basic solution
        objective_function = std::vector<double>(num_nonbasics + num_basics + 2, 0);
        // Column index of the auxiliary variable
        int aux_index = num_nonbasics + num_basics;
        objective_function[aux_index] = -1;

        for (size_t i = 0; i < num_basics; ++i) {
            for (size_t j = 0; j < num_basics; ++j) {
                if (i == j) {
                    solution_table[i].emplace_back(1);
                } else {
                    solution_table[i].emplace_back(0);
                }
            }
            // The additional nonbasic variable
            solution_table[i].emplace_back(-1);
            solution_table[i].emplace_back(b[i]);
        }

        pivot(min_b_index, aux_index);

        solve();

        if (fabs(*objective_function.rbegin() - 0) < epsilon) {
            if (basic_to_row.find(aux_index) != basic_to_row.end()) {
                // auxiliary variable is basic, need an additional pivot
                int aux_row = basic_to_row[aux_index];
                int pivot_col = -1;
                for (int col = 0; col < num_nonbasics + num_basics; ++col) {
                    if (basic_to_row.find(col) == basic_to_row.end() && fabs(solution_table[aux_row][col] - 0) > epsilon) {
                        pivot_col = col;
                        break;
                    }
                }
                pivot(aux_row, pivot_col);
            }

            // Construct the solution_table for the actual problem
            for (std::vector<double> &row : solution_table) {
                row.erase(row.end() - 2);  // Remove the column for the auxiliary variable
            }
            // Initialize the new objective_function
            objective_function = c;
            for (size_t i = 0; i < num_basics + 1; ++i) {
                // num_basics + 1, last entry corresponds to the negative of the current basic solution
                objective_function.emplace_back(0);
            }
            // Express objective_function with only nonbasic variables
            for (int col = 0; col < num_nonbasics; ++col) {
                if (basic_to_row.find(col) != basic_to_row.end()) {
                    int row = basic_to_row[col];
                    double multiple = objective_function[col];
                    for (int col_in_OF = 0; col_in_OF < objective_function.size(); ++col_in_OF) {
                        objective_function[col_in_OF] -= multiple * solution_table[row][col_in_OF];
                    }
                }
            }
            infinite_solution = false;
            has_solution = true;
        } else {
            has_solution = false;
        }


    }

}



void SimplexAlgorithm::solve() {
    std::vector<double> constraints(num_basics, 0);
    // Loop while there is at least one non-negative nonbasic variable
    while (any_of(objective_function.begin(), objective_function.end() - 1, [](double c) { return c > epsilon; })) {
        // Choose the non basic variable with smallest index
        int pivot_col =
                find_if(objective_function.begin(), objective_function.end() - 1, [](double c) { return c > epsilon; })
                - objective_function.begin();

        for (int row = 0; row < num_basics; ++row) {
            if (solution_table[row][pivot_col] > epsilon) {
                constraints[row] = solution_table[row][solution_table[row].size() - 1] / solution_table[row][pivot_col];
            } else {
                constraints[row] = std::numeric_limits<double>::infinity();
            }
        }
        // Get the pivot row
        int pivot_row = min_element(constraints.begin(), constraints.end()) - constraints.begin();

        if (constraints[pivot_row] == std::numeric_limits<double>::infinity()) {
            infinite_solution = true;
            return;
        }

        // Gaussian elimination, or pivot
        pivot(pivot_row, pivot_col);
    }
}

void SimplexAlgorithm::pivot(int pivot_row, int pivot_col) {
    // Make the coefficient at pivot_row pivot_col 1
    double factor = solution_table[pivot_row][pivot_col];
    for (int col = 0; col < solution_table[pivot_row].size(); ++col) {
        solution_table[pivot_row][col] /= factor;
    }
    // Eliminate pivot_col from all other rows
    for (int row = 0; row < solution_table.size(); ++row) {
        if (row == pivot_row) {
            continue;
        }
        double multiple = solution_table[row][pivot_col];
        for (int col = 0; col < solution_table[row].size(); ++col) {
            solution_table[row][col] -= multiple * solution_table[pivot_row][col];
        }
    }
    // Dealing with objective_function
    double multiple = objective_function[pivot_col];
    for (int col = 0; col < objective_function.size(); ++col) {
        objective_function[col] -= multiple * solution_table[pivot_row][col];
    }
    // Update basic variables and their corresponding rows
    int entering_col = pivot_col;
    int leaving_col = row_to_basic[pivot_row];
    basic_to_row.erase(leaving_col);
    basic_to_row.emplace(entering_col, pivot_row);
    row_to_basic[pivot_row] = entering_col;
}


void SimplexAlgorithm::print_solution_table() const {
    std::cout << "============================================================\n";
    std::cout << "Maximizing: \n";
    for (size_t i = 0; i < objective_function.size() - 1; ++i) {
        std::cout << std::setw(7) << std::right << "x" << i;
    }
    std::cout << std::endl;
    std::cout << std::fixed << std::setprecision(2);
    for (size_t i = 0; i < objective_function.size(); ++i) {
        std::cout << std::setw(8) << std::right << objective_function[i];
    }
    std::cout << std::endl;
    std::cout << "Subject to: \n";
    for (size_t i = 0; i < solution_table.size(); ++i) {
        for (size_t j = 0; j < solution_table[0].size(); ++j) {
            std::cout << std::setw(8) << std::right << solution_table[i][j];
        }
        std::cout << std::endl;
    }
    std::cout << "Basics:" << std::endl;
    for (std::pair<int, int> p : basic_to_row) {
        std::cout << p.first << ": " << p.second << ", ";

    }
    std::cout << std::endl;
    std::cout << "Rows:" << std::endl;
    for (std::pair<int, int> p : row_to_basic) {
        std::cout << p.first << ": " << p.second << ", ";

    }
    std::cout << std::endl;
    std::cout << "============================================================\n";

}

void SimplexAlgorithm::print_solution(int precision) {
    if (!has_solution) {
        std::cout << "Infeasible\n";
    } else {
        solve();
        if (infinite_solution) {
            std::cout << "Unbounded\n";
        } else {
            solution = std::vector<double>(num_nonbasics, 0);
            for (int col = 0; col < num_nonbasics; ++col) {
                if (basic_to_row.find(col) != basic_to_row.end()) {
                    solution[col] = *solution_table[basic_to_row[col]].rbegin();
                }
            }
            std::cout << "The optimal solution is: " << std::endl;
            std::cout << std::fixed << std::setprecision(precision);
            for (int i = 0; i < num_nonbasics; ++i) {
                std::cout << "x" << i << " = " << solution[i] << ((i != num_nonbasics - 1) ? ", " : "\n");
            }
            std::cout << "Optimal = " << -*objective_function.rbegin() << std::endl;
        }
    }
}