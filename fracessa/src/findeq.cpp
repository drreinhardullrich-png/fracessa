#include <fracessa/fracessa.hpp>
#include <fracessa/bitset64.hpp>
#include <fracessa/eigen_extensions.hpp>
#include <fracessa/helper.hpp>
#include <Eigen/LU>
#include <limits>

/*
 * ============================================================================
 * Why PartialPivLU is the right choice for your constraints
 * ============================================================================
 * 
 * Your requirements were:
 * 
 * 1. Fastest possible algorithm
 *    PartialPivLU is significantly faster than ColPivHouseholderQR
 *    (~2× faster on average, sometimes more).
 * 
 * 2. Never return "no solution" unless 100% certain
 *    To meet this requirement, you must use a very conservative singularity
 *    check (machine epsilon).
 * 
 *    PartialPivLU works perfectly for this:
 *    Pivot = 0 → matrix is mathematically singular → "no solution"
 *    Pivot > epsilon → treat as invertible
 * 
 * 3. Prefer a wrong solution over a false "no solution"
 *    PartialPivLU is more "optimistic" and will attempt a solution unless
 *    pivot is truly zero.
 * 
 * 4. A is square
 *    So QR (which normally handles rectangular) provides no special benefit.
 * 
 * 5. You do NOT want least-norm or least-squares
 *    ColPivHouseholderQR is a least-squares solver and will always return a
 *    solution even when A is singular — unless you enable thresholding, which
 *    becomes harder to tune for your requirement ("100% sure").
 * 
 *    LU is simpler and clearer:
 *    A pivot of 0 = singular.
 *    QR has no such clear criterion.
 * 
 * ============================================================================
 * Why NOT ColPivHouseholderQR in your case
 * ============================================================================
 * 
 * QR with column pivoting is:
 *     more stable, but
 *     slower, and
 *     its rank detection always depends on a tolerance.
 * 
 * There is no "100% certain" singularity detection with QR unless the matrix
 * has an exact zero column after pivoting — a rare event due to floating-point
 * operations.
 * 
 * To avoid false "no solution", you would need a very small tolerance, but
 * then QR loses its stability advantage.
 * 
 * In short:
 *     ✔ LU lets you reliably detect a pivot=0
 *     ✖ QR's rank detection is always tolerance-based
 * 
 * ============================================================================
 * FINAL RECOMMENDATION
 * ============================================================================
 * 
 * 👉 Choose PartialPivLU
 * 
 * with a pivot threshold:
 *     const double tol = std::numeric_limits<double>::epsilon();
 * 
 * This gives you:
 *     • Fastest speed
 *     • Clear singularity detection (pivot=0)
 *     • Zero risk of false "no solution"
 *     • A solution in all borderline cases
 *     • Minimal complexity
 * 
 * ============================================================================
 */

// Helper function to build full solution vector from support (double version)
bool fracessa::build_solution_vector(const DoubleVector& solution, DoubleVector& solution_full_n)
{
    solution_full_n = DoubleVector::Zero(dimension);
    
    size_t tracker = 0;
    for (size_t i = 0; i < dimension; i++) {
        if (_c.support.test(i)) {
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
    solution_full_n = RationalVector::Zero(dimension);
    
    size_t tracker = 0;
    for (size_t i = 0; i < dimension; i++) {
        if (_c.support.test(i)) {
            rational x = solution(tracker);
            if (x > rational(0)) {
                solution_full_n(i) = x;
            } else {
                return false;
            }
            tracker += 1;
        } else {
            solution_full_n(i) = rational(0);
        }
    }
    
    return true;
}

// Helper function to check constraints p'Ap<=v for rows not in support (double version)
bool fracessa::check_constraints(const DoubleVector& solution, const DoubleVector& solution_full_n)
{
    double errorbound_rowsum = 5e-5 * dimension; // here use a wide margin. if it is a false positive, it will be eliminated by the rational check!

    for (size_t i = 0; i < dimension; i++) {
        if (!_c.support.test(i)) { //not in the support - rows
            double rowsum = 0.;
            for (size_t j = 0; j < dimension; j++)
                if (_c.support.test(j)) // is in the support - columns
                    rowsum += game_matrix_double(i,j) * solution_full_n(j);

            if (!(rowsum <= solution(_c.support_size) + errorbound_rowsum)) //result
                return false;
        }
    }
    return true;
}

// Helper function to check constraints p'Ap<=v for rows not in support (rational version)
bool fracessa::check_constraints(const RationalVector& solution, const RationalVector& solution_full_n, bitset64& extended_support)
{
    for (size_t i = 0; i < dimension; i++) {
        if (!_c.support.test(i)) //not in the support - rows
        {
            rational rowsum = rational(0);
            for (size_t j = 0; j < dimension; j++)
                if (_c.support.test(j)) // is in the support - columns
                    rowsum += game_matrix(i,j) * solution_full_n(j);

            if (rowsum > solution(_c.support_size))
                return false;
            if (rowsum == solution(_c.support_size))
                extended_support.set(i);
        }
    }
    return true;
}

bool fracessa::find_candidate_double(DoubleMatrix& A, DoubleVector& b)
{
    Eigen::PartialPivLU<DoubleMatrix> lu(A);
    
    auto diag = lu.matrixLU().diagonal().cwiseAbs();
    double minPivot = diag.minCoeff();
    const double tol = std::numeric_limits<double>::epsilon();
   
    if (minPivot < tol) { // use the smallest possible tolerance, in case of false positives they get eliminated by the rational check!
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


bool fracessa::find_candidate_rational(RationalMatrix& A, RationalVector& b)
{
    // Use Eigen's BareissLU solver (fraction-free for exact arithmetic)
    Eigen::Ext::BareissLU<RationalMatrix> bareiss(A);
    
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
    bitset64 extended_support = _c.support;
    if (!check_constraints(solution, solution_full_n, extended_support)) {
        return false;
    }
    
    _c.vector = solution_full_n;
    _c.extended_support = extended_support;
    _c.extended_support_size = extended_support.count();
    _c.payoff = solution(_c.support_size);
    _c.payoff_double = static_cast<double>(_c.payoff);

    return true;
}

