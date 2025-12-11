#include <fracessa/fracessa.hpp>
#include <fracessa/bitset64.hpp>
#include <rational_linalg/linear_solver.hpp>
#include <cstdlib>
#include <type_traits>

// Templated function for all types (double, small_rational, rational)
template<typename T>
bool fracessa::find_candidate(const bitset64& support, size_t support_size)
{
    const auto& game_matrix = matrix_server_.get_game_matrix<T>();

    // Use LinearSolver for all types
    auto& Ab = matrix_server_.get_linear_system<T>(support, support_size);
    rational_linalg::LinearSolver<T> solver(Ab);
    rational_linalg::Matrix<T> solution;
    bool solved = solver.solve(solution);
    
    if (!solved)
        return false;    
    //solutions with zeros in it (meaning not full support for this matrix are already eliminated by the solver!!!!!!  

    // Build full solution vector from support by padding with zeros###################################################################
    rational_linalg::Matrix<T> solution_full_n = rational_linalg::Matrix<T>(dimension_, 1);    
    size_t tracker = 0;
    for (size_t i = 0; i < dimension_; i++) {
        if (bs64::test(support, i)) {
            solution_full_n(i, 0) = solution(tracker, 0);
            tracker += 1;
        } else 
            solution_full_n(i, 0) = T(0);
    }

    // check (Ap)_i <= v for all rows i not in the support ###################################################################################
    const T payoff = solution(support_size, 0);
    candidate_.extended_support = support; //gets copied!
    
    for (size_t i = 0; i < dimension_; i++) {
        if (!bs64::test(support, i)) { //not in the support - rows
            T rowsum = T(0);
            for (size_t j = 0; j < dimension_; j++) {
                if (bs64::test(support, j)) { // is in the support - columns
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
                    bs64::set(candidate_.extended_support, i);
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
        candidate_.extended_support_size = bs64::count(candidate_.extended_support);
    }
    return true;
}

// Explicit template instantiations
template bool fracessa::find_candidate<double>(const bitset64&, size_t);
template bool fracessa::find_candidate<small_rational>(const bitset64&, size_t);
template bool fracessa::find_candidate<rational>(const bitset64&, size_t);
