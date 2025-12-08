#include <fracessa/fracessa.hpp>
#include <fracessa/bitset64.hpp>
#include <rational_linalg/bareiss_gauss.hpp>
#include <Eigen/LU>
#include <limits>

// Helper function to build full solution vector from support (double version)
bool fracessa::build_solution_vector(const DoubleVector& solution, DoubleVector& solution_full_n)
{
    solution_full_n = DoubleVector::Zero(dimension_);
    
    size_t tracker = 0;
    for (size_t i = 0; i < dimension_; i++) {
        if (c_.support.test(i)) {
            double x = solution(tracker);
            if (x > -1e-5) { //large margin to be on the safe side, if it is a false positive, it will be eliminated by the rational check!
                solution_full_n(i) = x;
            } else {
                return false;
            }
            tracker += 1;
        } else {
            solution_full_n(i) = 0.;
        }
    }
    
    return true;
}

// Helper function to build full solution vector from support (rational version)
bool fracessa::build_solution_vector(const RationalVector& solution, RationalVector& solution_full_n)
{
    solution_full_n = RationalVector::Zero(dimension_);
    
    size_t tracker = 0;
    for (size_t i = 0; i < dimension_; i++) {
        if (c_.support.test(i)) {
            rational x = solution(tracker, 0);
            if (x > rational(0)) {
                solution_full_n(i, 0) = x;
            } else {
                return false;
            }
            tracker += 1;
        } else {
            solution_full_n(i, 0) = rational(0);
        }
    }
    
    return true;
}

// Helper function to check constraints p'Ap<=v for rows not in support (double version)
bool fracessa::check_constraints(const DoubleVector& solution, const DoubleVector& solution_full_n)
{
    double errorbound_rowsum = 5e-5 * dimension_; // here use a wide margin. if it is a false positive, it will be eliminated by the rational check!

    for (size_t i = 0; i < dimension_; i++) {
        if (!c_.support.test(i)) { //not in the support - rows
            double rowsum = 0.;
            for (size_t j = 0; j < dimension_; j++)
                if (c_.support.test(j)) // is in the support - columns
                    rowsum += game_matrix_double_(i,j) * solution_full_n(j);

            if (!(rowsum <= solution(c_.support_size) + errorbound_rowsum)) //result
                return false;
        }
    }
    return true;
}

// Helper function to check constraints p'Ap<=v for rows not in support (rational version)
bool fracessa::check_constraints(const RationalVector& solution, const RationalVector& solution_full_n, bitset64& extended_support)
{
    for (size_t i = 0; i < dimension_; i++) {
            if (!c_.support.test(i)) //not in the support - rows
        {
            rational rowsum = rational(0);
            for (size_t j = 0; j < dimension_; j++)
                if (c_.support.test(j)) // is in the support - columns
                    rowsum += game_matrix_(i,j) * solution_full_n(j, 0);

            if (rowsum > solution(c_.support_size, 0))
                return false;
            if (rowsum == solution(c_.support_size, 0))
                extended_support.set(i);
        }
    }
    return true;
}

bool fracessa::find_candidate_double_optimized(DoubleMatrix& A_SS)
{
    // Note: Principal submatrices of circulant matrices are NOT necessarily circulant,
    // so we cannot use FFT solver here. We always use LDLT for principal submatrices.
    
    // Use LDLT decomposition to solve A_{S,S} * v = 1_k
    Eigen::LDLT<DoubleMatrix> ldlt(A_SS);
    
    // Check if LDLT succeeded (matrix must be positive or negative definite)
    if (ldlt.info() != Eigen::Success) {
        return false; // Matrix is indefinite, fall back to standard method
    }
    
    // Create vector of ones (1_k)
    DoubleVector ones_k = DoubleVector::Ones(c_.support_size);
    
    // Solve A_{S,S} * v = 1_k
    DoubleVector v = ldlt.solve(ones_k);
    
    // Calculate s = 1_k^T * v
    double s = DoubleVector::Ones(c_.support_size).dot(v);
    
    // Check if s is too close to zero (would cause division by zero)
    const double tol = std::numeric_limits<double>::epsilon();
    if (std::abs(s) < tol) {
        return false;
    }
    
    // Calculate x_S = v / s
    DoubleVector x_S = v / s;
    
    // Calculate lambda = 1 / s
    double lambda = 1.0 / s;
    
    // Construct full solution vector [x_S, lambda] of size (support_size + 1)
    DoubleVector solution(c_.support_size + 1);
    solution.head(c_.support_size) = x_S;
    solution(c_.support_size) = lambda;
    
    DoubleVector solution_full_n;
    
    // Build large vector, if none of the elements is too negative
    if (!build_solution_vector(solution, solution_full_n)) {
        return false;
    }

    // check p'Ap<=v for all rows not in the support
    if (!check_constraints(solution, solution_full_n)) {
        return false;
    }
    
    return true;
}

bool fracessa::find_candidate_double(DoubleMatrix& A, DoubleVector& b)
{
    Eigen::PartialPivLU<DoubleMatrix> lu(A);
    
    auto diag = lu.matrixLU().diagonal().cwiseAbs();
    double minPivot = diag.minCoeff();
    const double tol = std::numeric_limits<double>::epsilon();
   
    if (minPivot < tol) { // use the smallest possible tolerance, false positives will get detected by the rational check!
        return false; // Matrix is singular
    }    
    DoubleVector solution = lu.solve(b);
    DoubleVector solution_full_n;
    
    // Build large vector, if none of the elements is too negative
    if (!build_solution_vector(solution, solution_full_n)) {
        return false;
    }

    // check p'Ap<=v for all rows not in the support
    if (!check_constraints(solution, solution_full_n)) {
        return false;
    }
    
    return true;
}


bool fracessa::find_candidate_rational_optimized(RationalMatrix& A_SS)
{
    // Use BareissLU to solve A_{S,S} * v = 1_k (A_SS is symmetric)
    rational_linalg::BareissGauss<rational> bareiss(A_SS);
    
    // Create vector of ones (1_k)
    RationalVector ones_k = RationalVector::Ones(c_.support_size);
    
    // Solve A_{S,S} * v = 1_k
    RationalVector v;
    if (!bareiss.solve(ones_k, v)) {
        return false; // Matrix is singular, fall back to standard method
    }
    
    // Calculate s = 1_k^T * v
    rational s = v.sum();
    
    // Check if s is zero (would cause division by zero)
    if (s == rational(0)) {
        return false;
    }
    
    // Calculate x_S = v / s
    RationalVector x_S = v / s;
    
    // Calculate lambda = 1 / s
    rational lambda = rational(1) / s;
    
    // Construct full solution vector [x_S, lambda] of size (support_size + 1)
    RationalVector solution(c_.support_size + 1, 1);
    RationalVector x_S_head = x_S.head(c_.support_size);
    for (size_t i = 0; i < c_.support_size; ++i) {
        solution(i, 0) = x_S_head(i, 0);
    }
    solution(c_.support_size, 0) = lambda;
    
    RationalVector solution_full_n;
    
    // Build large vector, if all elements are not too negative
    if (!build_solution_vector(solution, solution_full_n)) {
        return false;
    }

    // check p'Ap<=v for all rows not in the support
    bitset64 extended_support = c_.support;
    if (!check_constraints(solution, solution_full_n, extended_support)) {
        return false;
    }
    
    c_.vector = solution_full_n;
    c_.extended_support = extended_support;
    c_.extended_support_size = extended_support.count();
    c_.payoff = solution(c_.support_size, 0);
    c_.payoff_double = static_cast<double>(c_.payoff);

    return true;
}

bool fracessa::find_candidate_rational(RationalMatrix& A, RationalVector& b)
{
    // Use BareissLU solver (fraction-free for exact arithmetic)
    rational_linalg::BareissGauss<rational> bareiss(A);
    
    RationalVector solution;
    if (!bareiss.solve(b, solution)) {
        return false;
    }
    
    RationalVector solution_full_n;
    
    // Build large vector, if all elements are not too negative
    if (!build_solution_vector(solution, solution_full_n)) {
        return false;
    }

    // check p'Ap<=v for all rows not in the support
    bitset64 extended_support = c_.support;
    if (!check_constraints(solution, solution_full_n, extended_support)) {
        return false;
    }
    
    c_.vector = solution_full_n;
    c_.extended_support = extended_support;
    c_.extended_support_size = extended_support.count();
    c_.payoff = solution(c_.support_size, 0);
    c_.payoff_double = static_cast<double>(c_.payoff);

    return true;
}

