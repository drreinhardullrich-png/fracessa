#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <Eigen/Dense>
#include <Eigen/Cholesky>
#include <fracessa/matrix.hpp>
#include <fracessa/types.hpp>

using namespace matrix_ops;
using DoubleMatrix = Eigen::MatrixXd;
using RationalMatrix = Eigen::Matrix<rational, Eigen::Dynamic, Eigen::Dynamic>;

// Helper function to check with Eigen's LLT
bool is_positive_definite_eigen_llt(const DoubleMatrix& A) {
    Eigen::LLT<DoubleMatrix> llt(A);
    return (llt.info() == Eigen::Success);
}

// Helper function to get eigenvalues
void print_eigenvalues(const DoubleMatrix& A, const std::string& label) {
    Eigen::SelfAdjointEigenSolver<DoubleMatrix> solver(A);
    if (solver.info() != Eigen::Success) {
        std::cout << "  " << label << " eigenvalues: Failed to compute\n";
        return;
    }
    auto eigenvals = solver.eigenvalues();
    std::cout << "  " << label << " eigenvalues: [";
    for (int i = 0; i < eigenvals.size(); ++i) {
        std::cout << std::scientific << std::setprecision(6) << eigenvals(i);
        if (i < eigenvals.size() - 1) std::cout << ", ";
    }
    std::cout << "]\n";
    std::cout << "  Min eigenvalue: " << std::scientific << std::setprecision(6) << eigenvals.minCoeff() << "\n";
}

// Helper function to print matrix
void print_matrix(const DoubleMatrix& A, const std::string& label) {
    std::cout << "  " << label << ":\n";
    std::cout << std::fixed << std::setprecision(6);
    for (int i = 0; i < A.rows(); ++i) {
        std::cout << "    [";
        for (int j = 0; j < A.cols(); ++j) {
            std::cout << std::setw(10) << A(i, j);
            if (j < A.cols() - 1) std::cout << ", ";
        }
        std::cout << "]\n";
    }
}

// Helper function to print matrix (rational)
void print_matrix_rational(const RationalMatrix& A, const std::string& label) {
    std::cout << "  " << label << ":\n";
    for (int i = 0; i < A.rows(); ++i) {
        std::cout << "    [";
        for (int j = 0; j < A.cols(); ++j) {
            std::cout << A(i, j);
            if (j < A.cols() - 1) std::cout << ", ";
        }
        std::cout << "]\n";
    }
}

// Test a matrix with all three methods
void test_matrix(const DoubleMatrix& A_double, const std::string& test_name) {
    std::cout << "\n" << std::string(80, '=') << "\n";
    std::cout << "Test: " << test_name << "\n";
    std::cout << std::string(80, '=') << "\n";
    
    // Convert to rational for ground truth
    RationalMatrix A_rational = RationalMatrix::Zero(A_double.rows(), A_double.cols());
    for (int i = 0; i < A_double.rows(); ++i) {
        for (int j = 0; j < A_double.cols(); ++j) {
            A_rational(i, j) = rational(A_double(i, j));
        }
    }
    
    // Print matrix
    print_matrix(A_double, "Matrix (double)");
    print_matrix_rational(A_rational, "Matrix (rational)");
    
    // Print eigenvalues
    print_eigenvalues(A_double, "Double");
    
    // Test all three methods
    bool result_custom = is_positive_definite_double(A_double);
    bool result_eigen = is_positive_definite_eigen_llt(A_double);
    bool result_rational = is_positive_definite_rational(A_rational);
    
    std::cout << "\n  Results:\n";
    std::cout << "    Custom is_positive_definite_double: " << (result_custom ? "TRUE" : "FALSE") << "\n";
    std::cout << "    Eigen LLT:                        " << (result_eigen ? "TRUE" : "FALSE") << "\n";
    std::cout << "    Rational (ground truth):          " << (result_rational ? "TRUE" : "FALSE") << "\n";
    
    // Analyze differences
    std::cout << "\n  Analysis:\n";
    if (result_custom != result_eigen) {
        std::cout << "    ⚠️  Custom and Eigen LLT DISAGREE!\n";
    }
    if (result_custom != result_rational) {
        std::cout << "    ❌ Custom method DISAGREES with ground truth (rational)!\n";
        if (result_custom && !result_rational) {
            std::cout << "       → FALSE POSITIVE: Custom says PD but rational says not PD\n";
        } else {
            std::cout << "       → FALSE NEGATIVE: Custom says not PD but rational says PD\n";
        }
    }
    if (result_eigen != result_rational) {
        std::cout << "    ❌ Eigen LLT DISAGREES with ground truth (rational)!\n";
        if (result_eigen && !result_rational) {
            std::cout << "       → FALSE POSITIVE: Eigen says PD but rational says not PD\n";
        } else {
            std::cout << "       → FALSE NEGATIVE: Eigen says not PD but rational says PD\n";
        }
    }
    if (result_custom == result_eigen && result_custom == result_rational) {
        std::cout << "    ✓ All methods agree\n";
    }
    
    // Detailed Cholesky analysis for custom method
    std::cout << "\n  Custom method Cholesky details:\n";
    int n = A_double.rows();
    DoubleMatrix L = DoubleMatrix::Zero(n, n);
    bool failed = false;
    
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < i; j++) {
            double aSum = 0.;
            for (int k = 0; k < j; k++)
                aSum += L(i,k) * L(j,k);
            L(i,j) = (1 / L(j,j)) * (A_double(i,j) - aSum);
        }
        double bSum = 0.;
        for (int k = 0; k < i; k++)
            bSum += L(i,k) * L(i,k);
        double x = A_double(i,i) - bSum;
        std::cout << "    Step " << i << ": diagonal value = " << std::scientific << std::setprecision(10) << x;
        if (x < 0.01) {
            std::cout << " (FAILED: < 0.01 tolerance)";
            if (!failed) {
                failed = true;
            }
        }
        std::cout << "\n";
        if (!failed) {
            L(i,i) = sqrt(x);
        }
    }
    
    // Eigen LLT details
    std::cout << "\n  Eigen LLT details:\n";
    Eigen::LLT<DoubleMatrix> llt(A_double);
    if (llt.info() == Eigen::Success) {
        std::cout << "    Status: Success\n";
        DoubleMatrix L_eigen = llt.matrixL();
        std::cout << "    Diagonal values of L:\n";
        for (int i = 0; i < n; ++i) {
            double diag_val = L_eigen(i, i);
            std::cout << "      L(" << i << "," << i << ") = " << std::scientific << std::setprecision(10) << diag_val << "\n";
        }
    } else {
        std::cout << "    Status: Failed\n";
        if (llt.info() == Eigen::NumericalIssue) {
            std::cout << "    Reason: NumericalIssue\n";
        } else if (llt.info() == Eigen::InvalidInput) {
            std::cout << "    Reason: InvalidInput\n";
        }
    }
}

int main() {
    std::cout << "Positive Definite Check Comparison Test\n";
    std::cout << "========================================\n";
    std::cout << "Comparing:\n";
    std::cout << "  1. Custom is_positive_definite_double (tolerance: 0.01)\n";
    std::cout << "  2. Eigen LLT (adaptive tolerance)\n";
    std::cout << "  3. Rational is_positive_definite_rational (ground truth)\n";
    
    // Test 1: Clearly positive definite matrix
    {
        DoubleMatrix A(3, 3);
        A << 4.0, 1.0, 0.0,
             1.0, 5.0, 2.0,
             0.0, 2.0, 6.0;
        test_matrix(A, "Test 1: Clearly Positive Definite (3x3)");
    }
    
    // Test 2: Identity matrix
    {
        DoubleMatrix A = DoubleMatrix::Identity(4, 4);
        test_matrix(A, "Test 2: Identity Matrix (4x4)");
    }
    
    // Test 3: Matrix with very small eigenvalues (borderline case)
    {
        DoubleMatrix A(2, 2);
        A << 1.0, 0.0,
             0.0, 0.001;  // Very small eigenvalue
        test_matrix(A, "Test 3: Borderline Case - Very Small Eigenvalue");
    }
    
    // Test 4: Matrix with eigenvalue just above 0.01
    {
        DoubleMatrix A(2, 2);
        A << 1.0, 0.0,
             0.0, 0.011;  // Just above custom tolerance
        test_matrix(A, "Test 4: Eigenvalue Just Above Custom Tolerance (0.011)");
    }
    
    // Test 5: Matrix with eigenvalue just below 0.01
    {
        DoubleMatrix A(2, 2);
        A << 1.0, 0.0,
             0.0, 0.009;  // Just below custom tolerance
        test_matrix(A, "Test 5: Eigenvalue Just Below Custom Tolerance (0.009)");
    }
    
    // Test 6: Matrix with eigenvalue very close to 0.01
    {
        DoubleMatrix A(2, 2);
        A << 1.0, 0.0,
             0.0, 0.0100001;  // Very close to tolerance
        test_matrix(A, "Test 6: Eigenvalue Very Close to Custom Tolerance (0.0100001)");
    }
    
    // Test 7: Matrix with eigenvalue exactly 0.01
    {
        DoubleMatrix A(2, 2);
        A << 1.0, 0.0,
             0.0, 0.01;  // Exactly at tolerance
        test_matrix(A, "Test 7: Eigenvalue Exactly at Custom Tolerance (0.01)");
    }
    
    // Test 8: Matrix that should fail (negative eigenvalue)
    {
        DoubleMatrix A(2, 2);
        A << 1.0, 0.0,
             0.0, -1.0;
        test_matrix(A, "Test 8: Negative Eigenvalue (Should Fail)");
    }
    
    // Test 9: Matrix that should fail (indefinite)
    {
        DoubleMatrix A(2, 2);
        A << 1.0, 2.0,
             2.0, 1.0;
        test_matrix(A, "Test 9: Indefinite Matrix (Should Fail)");
    }
    
    // Test 10: Near-singular matrix (very small but positive eigenvalue)
    {
        DoubleMatrix A(2, 2);
        A << 1.0, 0.0,
             0.0, 1e-10;  // Extremely small but positive
        test_matrix(A, "Test 10: Near-Singular (Very Small Positive Eigenvalue: 1e-10)");
    }
    
    // Test 11: Matrix with small off-diagonal elements
    {
        DoubleMatrix A(3, 3);
        A << 0.1, 0.01, 0.0,
             0.01, 0.1, 0.01,
             0.0, 0.01, 0.1;
        test_matrix(A, "Test 11: Small Values Matrix (diagonal ~0.1)");
    }
    
    // Test 12: Matrix that becomes borderline during Cholesky
    {
        // This matrix might have issues during Cholesky decomposition
        DoubleMatrix A(3, 3);
        A << 0.02, 0.01, 0.0,
             0.01, 0.02, 0.01,
             0.0, 0.01, 0.02;
        test_matrix(A, "Test 12: Borderline During Cholesky (diagonal ~0.02)");
    }
    
    // Test 13: Large matrix
    {
        DoubleMatrix A = DoubleMatrix::Random(5, 5);
        A = A * A.transpose();  // Make it positive definite
        A += DoubleMatrix::Identity(5, 5) * 0.1;  // Ensure it's well-conditioned
        test_matrix(A, "Test 13: Random Positive Definite (5x5)");
    }
    
    std::cout << "\n" << std::string(80, '=') << "\n";
    std::cout << "Test Summary Complete\n";
    std::cout << std::string(80, '=') << "\n";
    
    return 0;
}

