#include <iostream>
#include <iomanip>
#include <Eigen/Dense>
#include <Eigen/Cholesky>

using DoubleMatrix = Eigen::MatrixXd;

int main() {
    // Matrix from log file: bee matrix
    // 2, -2,
    // -2, 2
    DoubleMatrix A(2, 2);
    A << 2.0, -2.0,
         -2.0, 2.0;
    
    std::cout << "Matrix A:\n" << A << "\n\n";
    
    // Compute eigenvalues
    Eigen::SelfAdjointEigenSolver<DoubleMatrix> eigen_solver(A);
    if (eigen_solver.info() == Eigen::Success) {
        std::cout << "Eigenvalues: " << eigen_solver.eigenvalues().transpose() << "\n";
        std::cout << "Min eigenvalue: " << eigen_solver.eigenvalues().minCoeff() << "\n\n";
    }
    
    // Try LLT (Cholesky) decomposition
    std::cout << "=== Eigen LLT (Cholesky) Decomposition ===\n";
    Eigen::LLT<DoubleMatrix> llt(A);
    
    std::cout << "LLT info: ";
    if (llt.info() == Eigen::Success) {
        std::cout << "Success\n";
    } else if (llt.info() == Eigen::NumericalIssue) {
        std::cout << "NumericalIssue\n";
    } else if (llt.info() == Eigen::InvalidInput) {
        std::cout << "InvalidInput\n";
    } else {
        std::cout << "Unknown (" << static_cast<int>(llt.info()) << ")\n";
    }
    
    if (llt.info() == Eigen::Success) {
        DoubleMatrix L = llt.matrixL();
        std::cout << "\nLower triangular matrix L:\n" << L << "\n";
        std::cout << "\nL * L^T (should equal A):\n" << L * L.transpose() << "\n";
        std::cout << "\nDiagonal of L:\n" << L.diagonal().transpose() << "\n";
    }
    
    // Try LDLT decomposition to get D matrix
    std::cout << "\n=== Eigen LDLT Decomposition ===\n";
    Eigen::LDLT<DoubleMatrix> ldlt(A);
    
    std::cout << "LDLT info: ";
    if (ldlt.info() == Eigen::Success) {
        std::cout << "Success\n";
    } else if (ldlt.info() == Eigen::NumericalIssue) {
        std::cout << "NumericalIssue\n";
    } else if (ldlt.info() == Eigen::InvalidInput) {
        std::cout << "InvalidInput\n";
    } else {
        std::cout << "Unknown (" << static_cast<int>(ldlt.info()) << ")\n";
    }
    
    if (ldlt.info() == Eigen::Success) {
        DoubleMatrix L = ldlt.matrixL();
        DoubleMatrix D = ldlt.vectorD().asDiagonal();
        auto transpositions = ldlt.transpositionsP();
        
        // Create permutation matrix from transpositions
        DoubleMatrix P = DoubleMatrix::Identity(A.rows(), A.cols());
        for (int i = 0; i < transpositions.size(); ++i) {
            P.row(i).swap(P.row(transpositions.indices()(i)));
        }
        
        std::cout << "\nLower triangular matrix L:\n" << L << "\n";
        std::cout << "\nDiagonal matrix D:\n" << D << "\n";
        std::cout << "\nDiagonal values (vectorD):\n" << ldlt.vectorD().transpose() << "\n";
        std::cout << "\nPermutation indices:\n" << transpositions.indices().transpose() << "\n";
        std::cout << "\nPermutation matrix P:\n" << P << "\n";
        std::cout << "\nL * D * L^T (should equal P * A * P^T):\n" << L * D * L.transpose() << "\n";
        std::cout << "\nP * A * P^T:\n" << P * A * P.transpose() << "\n";
        
        // Check if positive definite (all diagonal elements of D > 0)
        std::cout << "\nDiagonal elements of D:\n";
        for (int i = 0; i < D.rows(); ++i) {
            std::cout << "  D(" << i << "," << i << ") = " << std::scientific << std::setprecision(10) << D(i, i);
            if (D(i, i) > 0) {
                std::cout << " > 0 ✓\n";
            } else if (D(i, i) < 0) {
                std::cout << " < 0 ✗\n";
            } else {
                std::cout << " = 0 ✗\n";
            }
        }
        
        bool is_positive_definite = (ldlt.vectorD().array() > 0).all();
        std::cout << "\nIs positive definite (all D diagonal > 0): " << (is_positive_definite ? "YES" : "NO") << "\n";
    }
    
    // Also try the custom method for comparison
    std::cout << "\n=== Custom is_positive_definite_double method ===\n";
    int n = A.rows();
    DoubleMatrix L_custom = DoubleMatrix::Zero(n, n);
    bool failed = false;
    
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < i; j++) {
            double aSum = 0.;
            for (int k = 0; k < j; k++)
                aSum += L_custom(i,k) * L_custom(j,k);
            L_custom(i,j) = (1 / L_custom(j,j)) * (A(i,j) - aSum);
        }
        double bSum = 0.;
        for (int k = 0; k < i; k++)
            bSum += L_custom(i,k) * L_custom(i,k);
        double x = A(i,i) - bSum;
        std::cout << "Step " << i << ": diagonal value = " << std::scientific << std::setprecision(10) << x;
        if (x < 0.01) {
            std::cout << " (FAILED: < 0.01 tolerance)";
            failed = true;
        }
        std::cout << "\n";
        if (!failed) {
            L_custom(i,i) = sqrt(x);
        }
    }
    
    std::cout << "\nCustom method result: " << (failed ? "NOT POSITIVE DEFINITE" : "POSITIVE DEFINITE") << "\n";
    if (!failed) {
        std::cout << "Custom L matrix:\n" << L_custom << "\n";
    }
    
    return 0;
}

