#ifndef MATRIX_H
#define MATRIX_H

#include <vector>
#include <sstream>
#include <cmath>

#include <fracessa/helper.hpp>
#include <fracessa/bitset64.hpp>
#include <fracessa/rational_eigen.hpp>
#include <Eigen/Dense>
#include <eigen_extensions/bareiss_ldlt.hpp>

// Type aliases for Eigen matrices and vectors
using RationalMatrix = Eigen::Matrix<rational, Eigen::Dynamic, Eigen::Dynamic>;
using RationalVector = Eigen::Matrix<rational, Eigen::Dynamic, 1>;
using DoubleMatrix = Eigen::MatrixXd;
using DoubleVector = Eigen::VectorXd;

// For backward compatibility, create a namespace or use free functions
namespace matrix_ops {


//*************************************static construction methods********************************************

//circular nxn-matrices given by the first half-row starting at the second element
// Only accepts std::vector<rational>
inline RationalMatrix create_circular_symmetric(size_t n, const std::vector<rational> &half_row)
{
  std::vector<rational> first_row(n);

  first_row[0] = 0;
  if (n%2 == 0) {//even n
    first_row[n/2] = half_row[half_row.size()-1];
    for (size_t i=0; i<n/2-1; i++) {
      first_row[i+1] = half_row[i];
      first_row[n-i-1] = half_row[i];
    }
  } else { //odd n
    for (size_t i=0; i<n/2; i++) {
      first_row[i+1] = half_row[i];
      first_row[n-i-1] = half_row[i];
    }
  } 
  // Create circular matrix from first_row
  RationalMatrix A = RationalMatrix::Zero(n, n);
  for (size_t i=0; i<n; i++)
    for (size_t j=0; j<n; j++)
      A(i, j) = first_row[(j-i+n)%n];
  
  return A;
}

// Symmetric matrix from upper triangular values (n*(n+1)/2 values in row-major order)
// Rational version - only accepts std::vector<rational>
inline RationalMatrix create_symmetric(size_t n, const std::vector<rational> &upper_triangular)
{
  RationalMatrix A = RationalMatrix::Zero(n, n);
  
  size_t idx = 0;
  for (size_t i = 0; i < n; i++) {
    for (size_t j = i; j < n; j++) {
      rational val = upper_triangular[idx];
      A(i, j) = val;
      A(j, i) = val;  // Make symmetric
      idx++;
    }
  }
  
  return A;
}



// Convert rational matrix to double matrix
inline DoubleMatrix to_double(const RationalMatrix& A)
{
  DoubleMatrix result = DoubleMatrix::Zero(A.rows(), A.cols());
  for (Eigen::Index i=0; i<static_cast<Eigen::Index>(A.rows()); i++)
    for (Eigen::Index j=0; j<static_cast<Eigen::Index>(A.cols()); j++) {
      result(i,j) = static_cast<double>(A(i,j));
    }
  return result;
}

// String representation
template<typename MatrixType>
inline std::string to_string(const MatrixType& A)
{
  std::ostringstream oss;
  for (Eigen::Index i=0; i<static_cast<Eigen::Index>(A.rows()); i++) {
    for (Eigen::Index j=0; j<static_cast<Eigen::Index>(A.cols()); j++)
      oss << A(i,j) << ",";
    oss << std::endl;
  }
  return oss.str();
}


// Get bordered matrix A from support (split version for optimization)
inline void get_kkt_bordering(const DoubleMatrix& game_matrix, const bitset64& support, size_t support_size, DoubleMatrix& A)
{
  size_t n = support_size + 1;
  // Only resize if needed, and matrix is square, so check for rows is sufficient
  if (A.rows() != static_cast<Eigen::Index>(n)) {
    A = DoubleMatrix::Zero(n, n);
  } else {
    A.setZero();
  }
  
  size_t row = 0;
  Eigen::Index rows = static_cast<Eigen::Index>(game_matrix.rows());
  Eigen::Index cols = static_cast<Eigen::Index>(game_matrix.cols());

  // Fill rows 0 to support_size-1: submatrix from game_matrix, then -1 in last column
  for (Eigen::Index i=0; i<rows; i++)
    if (support.test(static_cast<size_t>(i)))
    {
      size_t column = 0;
      for (Eigen::Index j=0; j<cols; j++)
        if (support.test(static_cast<size_t>(j)))
        {
          A(static_cast<Eigen::Index>(row), static_cast<Eigen::Index>(column)) = game_matrix(i, j);
          column++;
        }
      A(static_cast<Eigen::Index>(row), static_cast<Eigen::Index>(support_size)) = -1.0;
      row++;
    }
  
  // Fill last row: all 1s, then 0 in last column
  for (size_t i=0; i<support_size; i++)
    A(static_cast<Eigen::Index>(support_size), static_cast<Eigen::Index>(i)) = 1.0;
  A(static_cast<Eigen::Index>(support_size), static_cast<Eigen::Index>(support_size)) = 0.0;
}

// Get RHS vector b (split version for optimization)
inline void get_kkt_rhs(size_t support_size, DoubleVector& b)
{
  size_t n = support_size + 1;
  // Only resize if needed
  if (b.size() != static_cast<Eigen::Index>(n)) {
    b = DoubleVector::Zero(n);
  }
  // Set all elements to 0 first
  b.setZero();
  // Last element is 1
  b(static_cast<Eigen::Index>(support_size)) = 1.0;
}

// Get bordered matrix A from support for rational matrices (split version for optimization)
inline void get_kkt_bordering(const RationalMatrix& game_matrix, const bitset64& support, size_t support_size, RationalMatrix& A)
{
  size_t n = support_size + 1;
  // Only resize if needed, and matrix is square, so check for rows is sufficient
  if (A.rows() != static_cast<Eigen::Index>(n)) {
    A = RationalMatrix::Zero(n, n);
  } else {
    A.setZero();
  }
  
  size_t row = 0;
  Eigen::Index rows = static_cast<Eigen::Index>(game_matrix.rows());
  Eigen::Index cols = static_cast<Eigen::Index>(game_matrix.cols());

  // Fill rows 0 to support_size-1: submatrix from game_matrix, then -1 in last column
  for (Eigen::Index i=0; i<rows; i++)
    if (support.test(static_cast<size_t>(i)))
    {
      size_t column = 0;
      for (Eigen::Index j=0; j<cols; j++)
        if (support.test(static_cast<size_t>(j)))
        {
          A(static_cast<Eigen::Index>(row), static_cast<Eigen::Index>(column)) = game_matrix(i, j);
          column++;
        }
      A(static_cast<Eigen::Index>(row), static_cast<Eigen::Index>(support_size)) = -1;
      row++;
    }
  
  // Fill last row: all 1s, then 0 in last column
  for (size_t i=0; i<support_size; i++)
    A(static_cast<Eigen::Index>(support_size), static_cast<Eigen::Index>(i)) = 1;
  A(static_cast<Eigen::Index>(support_size), static_cast<Eigen::Index>(support_size)) = 0;
}

// Get RHS vector b for rational matrices (split version for optimization)
inline void get_kkt_rhs(size_t support_size, RationalVector& b)
{
  size_t n = support_size + 1;
  // Only resize if needed
  if (b.size() != static_cast<Eigen::Index>(n)) {
    b = RationalVector::Zero(n);
  }
  // Set all elements to 0 first
  b.setZero();
  // Last element is 1
  b(static_cast<Eigen::Index>(support_size)) = 1;
}

// Extract principal submatrix from matrix using bitset64 mask (reuse version for optimization)
template<typename MatrixType>
inline void principal_submatrix(const MatrixType& A, size_t dimension, const bitset64& support, size_t support_size, MatrixType& submatrix)
{
  // Only resize if needed, and matrix is square, so check for rows is sufficient
  if (submatrix.rows() != static_cast<Eigen::Index>(support_size)) {
    submatrix = MatrixType::Zero(support_size, support_size);
  } else {
    submatrix.setZero();
  }
  
  size_t row = 0;
  for (size_t i = 0; i < dimension; ++i) {
    if (support.test(i)) {
      size_t col = 0;
      for (size_t j = 0; j < dimension; ++j) {
        if (support.test(j)) {
          submatrix(row, col) = A(static_cast<Eigen::Index>(i), static_cast<Eigen::Index>(j));
          ++col;
        }
      }
      ++row;
    }
  }
}




// Check if matrix is positive definite (LDLT decomposition for rational)
inline bool is_positive_definite_rational(const RationalMatrix& A)
{
  //ldlt-decomposition for rational numbers
  int n = A.rows();
  RationalMatrix D = RationalMatrix::Zero(n, n);
  RationalMatrix L = RationalMatrix::Identity(n, n);

  for (int i = 0; i < n; i++) {
    for (int j = 0; j < i; j++) {
      rational aSum = rational(0);
      for (int k = 0; k < j; k++)
        aSum += L(i,k) * L(j,k) * D(k,k);
      L(i,j) = (1 / D(j,j)) * (A(i,j) - aSum);
    }
    rational bSum = rational(0);
    for (int k = 0; k < i; k++)
      bSum += L(i,k) * L(i,k) * D(k,k);
    D(i,i) = A(i,i) - bSum;
    if (D(i,i)<=rational(0))
      return false;
  }
  return true;
}

// // Bareiss fraction-free elimination algorithm to check for Positive Definiteness
// inline bool is_positive_definite_rational(const RationalMatrix& A) {
//   int n = A.rows();
//   RationalMatrix B = A;                    // Working copy (will become Δ)
//   RationalMatrix L = RationalMatrix::Identity(n, n);  // Unit lower triangular
//   rational prevPivot = rational(1);
  
//   for (int k = 0; k < n; ++k) {
//       // Extract pivot from transformed matrix
//       rational pivot = B(k, k);
      
//       // Sylvester's criterion: all pivots must be positive
//       if (pivot <= rational(0)) {
//           return false;  // NOT positive definite
//       }
      
//       // Store multipliers for L (BEFORE updating B!)
//       for (int i = k + 1; i < n; ++i) {
//           L(i, k) = B(i, k);  // Store numerator, not quotient
//       }
      
//       // Bareiss update for remaining submatrix
//       // Only update if not last iteration
//       if (k < n - 1) {
//           for (int i = k + 1; i < n; ++i) {
//               // Update diagonal element
//               B(i, i) = (B(i, i) * pivot - B(i, k) * B(i, k)) / prevPivot;
              
//               // Update off-diagonal elements (upper triangle only)
//               for (int j = i + 1; j < n; ++j) {
//                   B(i, j) = (B(i, j) * pivot - B(i, k) * B(k, j)) / prevPivot;
//               }
//           }
          
//           // Zero out column k (except diagonal)
//           for (int i = k + 1; i < n; ++i) {
//               B(i, k) = rational(0);
//           }
          
//           prevPivot = pivot;
//       }
//   }
  
//   return true;  // All pivots positive
// }

// Check if matrix is positive definite (Cholesky for double)
inline bool is_positive_definite_double(const DoubleMatrix& A)
{
  // cholesky-decomposition for double!
  int n = A.rows();
  DoubleMatrix L = DoubleMatrix::Zero(n, n);

  for (int i = 0; i < n; i++) {
    for (int j = 0; j < i; j++) {
      double aSum = 0.;
      for (int k = 0; k < j; k++)
        aSum += L(i,k) * L(j,k);
      L(i,j) = (1 / L(j,j)) * (A(i,j) - aSum);
    }
    double bSum = 0.;
    for (int k = 0; k < i; k++)
      bSum += L(i,k) * L(i,k);
    double x = A(i,i) - bSum;
    if (x< 0.01)
      return false;
    L(i,i) = sqrt(x);
  }
  return true;
}

// Compute determinant (LU decomposition, Crout algorithm)
inline rational determinant(const RationalMatrix& A)
{
  //lu-decomposition, crout-algorithm
  RationalMatrix a = A;
  size_t n = A.rows();
  std::vector<size_t> indx = std::vector<size_t>(n);
  std::vector<rational> vv = std::vector<rational>(n);
  int d = 1;

  size_t i, imax, j, k;
  rational big, sum, temp;

  /* search for the largest element in each row; save the scaling in the
  temporary array vv and return zero if the matrix is singular */
  for(i=0; i<n; i++) {
    big=0.;
    for(j=0;j<n;j++) if((temp=abs(a(i,j)))>big) big=temp;
    if(big==0) return(0);
    vv[i]=big;
   }
   /* the main loop for the Crout's algorithm */
   for(j=0;j<n;j++) {
    /* this is the part a) of the algorithm except for i==j */
    for(i=0;i<j;i++) {
      sum=a(i,j);
      for(k=0;k<i;k++) sum-=a(i,k)*a(k,j);
      a(i,j)=sum;
    }
    /* initialize for the search for the largest pivot element */
    big=0;imax=j;
    /* this is the part a) for i==j and part b) for i>j + pivot search */
    for(i=j;i<n;i++) {
      sum=a(i,j);
      for(k=0;k<j;k++) sum-=a(i,k)*a(k,j);
      a(i,j)=sum;
      /* is the figure of merit for the pivot better than the best so far? */
      if((temp=vv[i]*abs(sum))>=big) {big=temp;imax=i;}
    }
    /* interchange rows, if needed, change parity and the scale factor */
    if(imax!=j) {
      for(k=0;k<n;k++) {temp=a(imax,k);a(imax,k)=a(j,k);a(j,k)=temp;}
      d=-d;vv[imax]=vv[j];
    }
    /* store the index */
    indx[j]=imax;
    if(a(j,j)==0) return(0);
    /* finally, divide by the pivot element */
    if(j<n-1) {
      temp=1/a(j,j);
      for(i=j+1;i<n;i++) a(i,j)*=temp;
    }
  }
  rational res = d;
  for(j=0; j<n; j++) res *= a(j,j);
  return(res);
}

// Transpose
inline RationalMatrix transpose(const RationalMatrix& A)
{
  return A.transpose();
}

// Develop by (remove row and column)
inline RationalMatrix develop_by(const RationalMatrix& A, size_t row, size_t col)
{
  Eigen::Index rows = static_cast<Eigen::Index>(A.rows());
  Eigen::Index cols = static_cast<Eigen::Index>(A.cols());
  RationalMatrix result = RationalMatrix::Zero(rows-1, cols-1);
  for (Eigen::Index i=0; i<rows-1; i++)
    for (Eigen::Index j=0; j<cols-1; j++) {
      Eigen::Index row_index_real = (i>=static_cast<Eigen::Index>(row)) ? i+1 : i;
      Eigen::Index col_index_real = (j>=static_cast<Eigen::Index>(col)) ? j+1 : j;
      result(i,j) = A(row_index_real, col_index_real);
    }
  return result;
}

// Adjugate
inline RationalMatrix adjugate(const RationalMatrix& A)
{
  Eigen::Index rows = static_cast<Eigen::Index>(A.rows());
  Eigen::Index cols = static_cast<Eigen::Index>(A.cols());
  RationalMatrix result = RationalMatrix::Zero(rows, cols);

  if (rows == 1) {
    result(0,0) = 1;
    return result;
  }

  for (Eigen::Index i=0; i<rows; i++) {
    for (Eigen::Index j=0; j<cols; j++) {
      result(i,j) = std::pow(-1, static_cast<int>(i+j)) * determinant(develop_by(A, static_cast<size_t>(i), static_cast<size_t>(j)));
    }
  }
  return transpose(result);
}

// Check if all elements are greater than zero
inline bool greater_zero(const RationalMatrix& A)
{
  for (Eigen::Index i=0; i<static_cast<Eigen::Index>(A.rows()); i++)
    for (Eigen::Index j=0; j<static_cast<Eigen::Index>(A.cols()); j++)
      if (A(i,j) <= 0)
        return false;
  return true;
}

// Check if matrix is copositive
inline bool is_strictly_copositive(const RationalMatrix& A)
{
  bool result = true;
  size_t n_rows = static_cast<size_t>(A.rows());
  bitset64::iterate_all_supports(n_rows, [&](const bitset64& support, unsigned) {
    if (!result) return;  // Early exit optimization
    size_t n = support.count();
    RationalMatrix subA = RationalMatrix::Zero(n, n);

    size_t row = 0;
    size_t column = 0;
    for (size_t i=0; i<n_rows; i++) {
      if (support.test(i)) {
        column = 0;
        for (size_t j=0; j<static_cast<size_t>(A.cols()); j++)
          if (support.test(j)) {
            subA(row,column) = A(static_cast<Eigen::Index>(i), static_cast<Eigen::Index>(j));
            column++;
          }
        row++;
      }
    }
    if (determinant(subA) <= 0 && greater_zero(adjugate(subA))) {
      result = false;  // Found counterexample
    }
  });
  return result;
}

} // namespace matrix_ops

#endif // MATRIX_H
