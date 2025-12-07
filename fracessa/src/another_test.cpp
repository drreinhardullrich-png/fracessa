#include <Eigen/Dense>
#include <iostream>
using namespace Eigen;
using namespace std;

int main() {
    MatrixXd A(2,2);
    A << 2.0, -2.0,
         -2.0,  2.0;

    // 0) Ensure exact symmetry (no tiny asymmetry)
    MatrixXd symA = 0.5 * (A + A.transpose());

    // 1) Try LLT
    LLT<MatrixXd> llt(symA);
    cout << "llt.info() = " << llt.info() << "\n";

    // 2) Print L diagonal (pivots)
    if (llt.info() == Success) {
        MatrixXd L = llt.matrixL();
        cout << "diag(L): " << L.diagonal().transpose() << "\n";
        cout << "reconstruction error ||A - L*L^T||_inf = "
             << (symA - L * L.transpose()).cwiseAbs().maxCoeff() << "\n";
    }

    // 3) Eigenvalue check (definitive)
    SelfAdjointEigenSolver<MatrixXd> es(symA);
    cout << "eigenvalues: " << es.eigenvalues().transpose() << "\n";
    cout << "min eigenvalue = " << es.eigenvalues().minCoeff() << "\n";

    // 4) LDLT check: inspect D
    LDLT<MatrixXd> ldlt(symA);
    if (ldlt.info() == Success) {
        cout << "LDLT succeeded. D diag: " << ldlt.vectorD().transpose() << "\n";
    } else {
        cout << "LDLT failed\n";
    }
}
