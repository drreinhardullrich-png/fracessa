#include <fracessa/fracessa.hpp>
#include <fracessa/bitset64.hpp>
#include <rational_linalg/linear_solver.hpp>
#include <cstdint>
#include <cstdlib>
#include <type_traits>

// Templated function for all types (double, small_rational, rational)
template<typename T>
bool fracessa::find_candidate(const rational_linalg::Matrix<T>& game_matrix, const bitset64& support, size_t support_size, rational_linalg::Matrix<T>& Ab)
{
    // Build KKT system as augmented matrix [A|b] from support
    // Resize if needed (constructor zero-initializes, so no separate zeroing needed)
    if (Ab.rows() != support_size + 1 || Ab.cols() != support_size + 2) 
        Ab = rational_linalg::Matrix<T>(support_size + 1, support_size + 2);

    // Fill all rows in one pass
    size_t ab_row = 0;
    for (size_t i = 0; i < dimension_; ++i) {
        if (support.test(i)) {
            // Fill submatrix columns (0 to support_size-1) from game_matrix
            size_t ab_col = 0;
            for (size_t j = 0; j < dimension_; ++j) {
                if (support.test(j)) {
                    Ab(ab_row, ab_col) = game_matrix(i, j);
                    ab_col++;
                }
            }
            // Column support_size: -1 (last column of A)
            Ab(ab_row, support_size) = T(-1);
            // Column support_size + 1: 0 (b vector, initially zero)
            Ab(ab_row, support_size + 1) = T(0);
            ab_row++;
        }
    }
    
    // Fill last row: all 1s in columns 0..support_size-1, 0 in column support_size, 1 in column n
    for (size_t i = 0; i < support_size; ++i) {
        Ab(support_size, i) = T(1);
    }
    Ab(support_size, support_size) = T(0);  // Column support_size (last column of A)
    Ab(support_size, support_size + 1) = T(1);  // Column support_size + 1 (b vector, last element is 1)

    // Use LinearSolver (automatically selects GaussDouble for double, GaussRational<T> for rational types)
    rational_linalg::LinearSolver<T> solver(Ab);   
    rational_linalg::Matrix<T> solution;

    if (!solver.solve(solution)) 
        return false;    
    //solutions with zeros in it (meaning not full support for this matrix are already eliminated by the solver!!!!!!  

    // Build full solution vector from support by padding with zeros###################################################################
    rational_linalg::Matrix<T> solution_full_n = rational_linalg::Matrix<T>(dimension_, 1);    
    size_t tracker = 0;
    for (size_t i = 0; i < dimension_; i++) {
        if (support.test(i)) {
            solution_full_n(i, 0) = solution(tracker, 0);
            tracker += 1;
        } else 
            solution_full_n(i, 0) = T(0);
    }

    // check (Ap)_i <= v for all rows i not in the support ###################################################################################
    const T payoff = solution(support_size, 0);
    candidate_.extended_support = support; //gets copied!
    
    for (size_t i = 0; i < dimension_; i++) {
        if (!support.test(i)) { //not in the support - rows
            T rowsum = T(0);
            for (size_t j = 0; j < dimension_; j++) {
                if (support.test(j)) { // is in the support - columns
                    rowsum += game_matrix(i,j) * solution_full_n(j, 0);
                }
            }
            if constexpr (std::is_same_v<T, double>) {
                if (rowsum > payoff + 1e-4 * dimension_) // huge margin, false positives eliminated by rational check
                    return false;
            } else {
                if (rowsum > payoff)
                    return false;
                if (rowsum == payoff)
                    candidate_.extended_support.set(i);
            }
        }
    }   
    // Convert solution_full_n to candidate's vector (always rational) - only for rational types
    if constexpr (!std::is_same_v<T, double>) {

        if constexpr (std::is_same_v<T, rational>) {
            candidate_.vector = solution_full_n;
            candidate_.payoff = payoff;
        } else {
            candidate_.vector = rational_linalg::convert_small_to_rational(solution_full_n);
            candidate_.payoff = small_to_rational(payoff);
        }
        candidate_.payoff_double = rational_to_double(payoff);        
        candidate_.extended_support_size = candidate_.extended_support.count();
    }
    return true;
}

// Explicit template instantiations
template bool fracessa::find_candidate<double>(const rational_linalg::Matrix<double>&, const bitset64&, size_t, rational_linalg::Matrix<double>&);
template bool fracessa::find_candidate<small_rational>(const rational_linalg::Matrix<small_rational>&, const bitset64&, size_t, rational_linalg::Matrix<small_rational>&);
template bool fracessa::find_candidate<rational>(const rational_linalg::Matrix<rational>&, const bitset64&, size_t, rational_linalg::Matrix<rational>&);

