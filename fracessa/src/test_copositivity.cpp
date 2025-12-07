#include <iostream>
#include <iomanip>
#include <Eigen/Dense>
#include <fracessa/matrix.hpp>
#include <eigen_extensions/copositivity.hpp>

using namespace matrix_ops;
using RationalMatrix = Eigen::Matrix<rational, Eigen::Dynamic, Eigen::Dynamic>;

void test_matrix(const RationalMatrix& A, const std::string& name, bool expected) {
    std::cout << "\n" << std::string(60, '=') << "\n";
    std::cout << "Test: " << name << "\n";
    std::cout << std::string(60, '=') << "\n";
    std::cout << "Matrix:\n" << A << "\n";
    
    // Test with copositivity.hpp (expects Eigen::MatrixXd, so convert)
    Eigen::MatrixXd A_double = matrix_ops::to_double(A);
    bool result_memoized = isStrictlyCopositiveMemoized(A_double);
    
    // Test with existing implementation
    bool result_existing = is_strictly_copositive(A);
    
    std::cout << "\nResults:\n";
    std::cout << "  Memoized algorithm: " << (result_memoized ? "TRUE" : "FALSE") << "\n";
    std::cout << "  Existing algorithm: " << (result_existing ? "TRUE" : "FALSE") << "\n";
    std::cout << "  Expected:          " << (expected ? "TRUE" : "FALSE") << "\n";
    
    if (result_memoized == result_existing) {
        std::cout << "  ✓ Algorithms AGREE\n";
    } else {
        std::cout << "  ✗ Algorithms DISAGREE!\n";
    }
    
    if (result_memoized == expected) {
        std::cout << "  ✓ Memoized matches expected\n";
    } else {
        std::cout << "  ✗ Memoized does NOT match expected\n";
    }
    
    // Print determinant
    rational det = matrix_ops::determinant(A);
    std::cout << "Determinant: " << det << "\n";
}

int main() {
    std::cout << "Copositivity Algorithm Comparison Test\n";
    std::cout << "======================================\n";
    
    // Test 1: Clearly copositive (positive definite)
    {
        RationalMatrix A(2, 2);
        A << rational(2), rational(1),
             rational(1), rational(2);
        test_matrix(A, "Test 1: Positive Definite 2x2", true);
    }
    
    // Test 2: Non-copositive (negative determinant, positive adjugate)
    {
        RationalMatrix A(2, 2);
        A << rational(1), rational(-2),
             rational(-2), rational(1);
        test_matrix(A, "Test 2: Negative det, positive adjugate", false);
    }
    
    // Test 3: Identity matrix (copositive)
    {
        RationalMatrix A = RationalMatrix::Identity(3, 3);
        test_matrix(A, "Test 3: Identity 3x3", true);
    }
    
    // Test 4: Matrix with det = 0 (borderline case)
    {
        RationalMatrix A(2, 2);
        A << rational(1), rational(1),
             rational(1), rational(1);
        test_matrix(A, "Test 4: Singular matrix (det=0)", false);
    }
    
    // Test 5: Negative diagonal (not copositive)
    {
        RationalMatrix A(2, 2);
        A << rational(-1), rational(0),
             rational(0), rational(1);
        test_matrix(A, "Test 5: Negative diagonal element", false);
    }
    
    // Test 6: Another non-copositive case (det=0, positive adjugate)
    {
        RationalMatrix A(3, 3);
        A << rational(1), rational(-1), rational(0),
             rational(-1), rational(1), rational(0),
             rational(0), rational(0), rational(1);
        test_matrix(A, "Test 6: 3x3 with det=0, positive adjugate in 2x2 block", false);
    }
    
    std::cout << "\n" << std::string(60, '=') << "\n";
    std::cout << "Test Summary Complete\n";
    std::cout << std::string(60, '=') << "\n";
    
    return 0;
}

