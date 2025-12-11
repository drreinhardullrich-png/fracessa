#ifndef RATIONAL_LINALG_MATRIX_HPP
#define RATIONAL_LINALG_MATRIX_HPP

#include <rational_linalg/types_rational.hpp>
#include <fracessa/bitset64.hpp>
#include <vector>
#include <cstring>
#include <sstream>
#include <string>
#include <type_traits>

namespace rational_linalg {

template<typename T>
class Matrix {
public:
    // Constructors
    inline Matrix() : rows_(0), cols_(0), data_() {}
    
    inline Matrix(size_t rows, size_t cols) : rows_(rows), cols_(cols), data_(rows * cols) {
        // data_ now contains exactly rows*cols elements, all default-initialized to 0
    }
    
    inline Matrix(const Matrix& other) : rows_(other.rows_), cols_(other.cols_), data_(other.data_) {}
    
    inline Matrix(Matrix&& other) noexcept : rows_(other.rows_), cols_(other.cols_), data_(std::move(other.data_)) {
        other.rows_ = 0;
        other.cols_ = 0;
    }
    
    inline Matrix& operator=(const Matrix& other) {
        rows_ = other.rows_;
        cols_ = other.cols_;
        data_ = other.data_;
        return *this;
    }
    
    inline Matrix& operator=(Matrix&& other) noexcept {
        rows_ = other.rows_;
        cols_ = other.cols_;
        data_ = std::move(other.data_);
        other.rows_ = 0;
        other.cols_ = 0;
        return *this;
    }

    // Accessors - NO bounds checking
    inline T& operator()(size_t i, size_t j) {
        return data_[i * cols_ + j];
    }
    
    inline const T& operator()(size_t i, size_t j) const {
        return data_[i * cols_ + j];
    }
    
    inline size_t rows() const { return rows_; }
    inline size_t cols() const { return cols_; }
    
    inline T* data() { return data_.data(); }
    inline const T* data() const { return data_.data(); }

    // Arithmetic operators
    inline Matrix operator+(const Matrix& other) const {
        Matrix result(rows_, cols_);
        const size_t n = rows_ * cols_;
        for (size_t i = 0; i < n; ++i) {
            result.data_[i] = data_[i] + other.data_[i];
        }
        return result;
    }
    
    inline Matrix operator-(const Matrix& other) const {
        Matrix result(rows_, cols_);
        const size_t n = rows_ * cols_;
        for (size_t i = 0; i < n; ++i) {
            result.data_[i] = data_[i] - other.data_[i];
        }
        return result;
    }
    
    inline Matrix operator*(const Matrix& other) const {
        Matrix result(rows_, other.cols_);
        for (size_t i = 0; i < rows_; ++i) {
            for (size_t j = 0; j < other.cols_; ++j) {
                T sum = T(0);
                for (size_t k = 0; k < cols_; ++k) {
                    sum += data_[i * cols_ + k] * other.data_[k * other.cols_ + j];
                }
                result.data_[i * result.cols_ + j] = sum;
            }
        }
        return result;
    }
    
    inline Matrix operator*(const T& scalar) const {
        Matrix result(rows_, cols_);
        const size_t n = rows_ * cols_;
        for (size_t i = 0; i < n; ++i) {
            result.data_[i] = data_[i] * scalar;
        }
        return result;
    }
    
    inline Matrix operator/(const T& scalar) const {
        Matrix result(rows_, cols_);
        const size_t n = rows_ * cols_;
        for (size_t i = 0; i < n; ++i) {
            result.data_[i] = data_[i] / scalar;
        }
        return result;
    }
    
    inline Matrix& operator+=(const Matrix& other) {
        const size_t n = rows_ * cols_;
        for (size_t i = 0; i < n; ++i) {
            data_[i] += other.data_[i];
        }
        return *this;
    }
    
    inline Matrix& operator-=(const Matrix& other) {
        const size_t n = rows_ * cols_;
        for (size_t i = 0; i < n; ++i) {
            data_[i] -= other.data_[i];
        }
        return *this;
    }
    
    inline Matrix& operator*=(const T& scalar) {
        const size_t n = rows_ * cols_;
        for (size_t i = 0; i < n; ++i) {
            data_[i] *= scalar;
        }
        return *this;
    }
    
    inline Matrix& operator/=(const T& scalar) {
        const size_t n = rows_ * cols_;
        for (size_t i = 0; i < n; ++i) {
            data_[i] /= scalar;
        }
        return *this;
    }

    // Comparison operators
    inline bool operator==(const Matrix& other) const {
        if (rows_ != other.rows_ || cols_ != other.cols_) return false;
        const size_t n = rows_ * cols_;
        for (size_t i = 0; i < n; ++i) {
            if (data_[i] != other.data_[i]) return false;
        }
        return true;
    }
    
    inline bool operator!=(const Matrix& other) const {
        return !(*this == other);
    }

    // Matrix operations
    inline Matrix transpose() const {
        Matrix result(cols_, rows_);
        for (size_t i = 0; i < rows_; ++i) {
            for (size_t j = 0; j < cols_; ++j) {
                result.data_[j * result.cols_ + i] = data_[i * cols_ + j];
            }
        }
        return result;
    }
    
    inline bool is_square() const {
        return rows_ == cols_;
    }
    
    inline T determinant() const {
        // Bareiss algorithm for determinant - optimized for small matrices
        if (rows_ != cols_) return T(0);
        if (rows_ == 0) return T(1);
        if (rows_ == 1) return data_[0];
        
        Matrix work = *this;
        T div_prev = T(1);
        int swap_count = 0;
        
        for (size_t k = 0; k < rows_ - 1; ++k) {
            // Partial pivoting
            size_t max_row = k;
            T max_val = work(k, k);
            if (max_val < T(0)) max_val = -max_val;
            
            for (size_t i = k + 1; i < rows_; ++i) {
                T val = work(i, k);
                if (val < T(0)) val = -val;
                if (val > max_val) {
                    max_val = val;
                    max_row = i;
                }
            }
            
            if (max_row != k) {
                // Swap rows
                for (size_t j = 0; j < cols_; ++j) {
                    T tmp = work(k, j);
                    work(k, j) = work(max_row, j);
                    work(max_row, j) = tmp;
                }
                ++swap_count;
            }
            
            const T pivot = work(k, k);
            if (pivot == T(0)) return T(0);
            
            // Bareiss update
            for (size_t i = k + 1; i < rows_; ++i) {
                for (size_t j = k + 1; j < cols_; ++j) {
                    work(i, j) = (work(i, j) * pivot - work(i, k) * work(k, j)) / div_prev;
                }
                work(i, k) = T(0);
            }
            
            div_prev = pivot;
        }
        
        T det = T(1);
        if (swap_count % 2 == 1) det = T(-1);
        
        for (size_t i = 0; i < rows_; ++i) {
            det *= work(i, i);
        }
        
        return det;
    }
    
    inline Matrix inverse() const {
        // Gauss-Jordan elimination for inverse - assumes non-singular
        Matrix inv(rows_, cols_);
        
        // Initialize augmented matrix [A | I]
        Matrix aug(rows_, cols_ * 2);
        for (size_t i = 0; i < rows_; ++i) {
            for (size_t j = 0; j < cols_; ++j) {
                aug(i, j) = data_[i * cols_ + j];
            }
            for (size_t j = 0; j < cols_; ++j) {
                aug(i, cols_ + j) = (i == j) ? T(1) : T(0);
            }
        }
        
        // Forward elimination with partial pivoting
        for (size_t k = 0; k < rows_; ++k) {
            // Partial pivoting
            size_t max_row = k;
            T max_val = aug(k, k);
            if (max_val < T(0)) max_val = -max_val;
            
            for (size_t i = k + 1; i < rows_; ++i) {
                T val = aug(i, k);
                if (val < T(0)) val = -val;
                if (val > max_val) {
                    max_val = val;
                    max_row = i;
                }
            }
            
            if (max_row != k) {
                // Swap rows
                for (size_t j = 0; j < aug.cols_; ++j) {
                    T tmp = aug(k, j);
                    aug(k, j) = aug(max_row, j);
                    aug(max_row, j) = tmp;
                }
            }
            
            const T pivot = aug(k, k);
            
            // Normalize pivot row
            for (size_t j = k + 1; j < aug.cols_; ++j) {
                aug(k, j) = aug(k, j) / pivot;
            }
            aug(k, k) = T(1);
            
            // Eliminate column k
            for (size_t i = 0; i < rows_; ++i) {
                if (i == k) continue;
                const T factor = aug(i, k);
                for (size_t j = k + 1; j < aug.cols_; ++j) {
                    aug(i, j) = aug(i, j) - aug(k, j) * factor;
                }
                aug(i, k) = T(0);
            }
        }
        
        // Extract inverse from right half
        for (size_t i = 0; i < rows_; ++i) {
            for (size_t j = 0; j < cols_; ++j) {
                inv(i, j) = aug(i, cols_ + j);
            }
        }
        
        return inv;
    }

    // Utility functions
    inline void swap_rows(size_t i, size_t j) {
        for (size_t k = 0; k < cols_; ++k) {
            T tmp = data_[i * cols_ + k];
            data_[i * cols_ + k] = data_[j * cols_ + k];
            data_[j * cols_ + k] = tmp;
        }
    }
    
    inline void swap_cols(size_t i, size_t j) {
        for (size_t k = 0; k < rows_; ++k) {
            T tmp = data_[k * cols_ + i];
            data_[k * cols_ + i] = data_[k * cols_ + j];
            data_[k * cols_ + j] = tmp;
        }
    }
    
    static inline Matrix zero(size_t rows, size_t cols) {
        Matrix result(rows, cols);
        // Constructor already zero-initializes all elements
        return result;
    }
    
    static inline Matrix identity(size_t n) {
        Matrix result(n, n);
        // Constructor already zero-initializes all elements, only need to set diagonal
        for (size_t i = 0; i < n; ++i) {
            result(i, i) = T(1);
        }
        return result;
    }

    // Vector helper functions (for Matrix used as vector)
    // head() - return first n elements as column vector
    inline Matrix head(size_t n) const {
        // Assume this is a column vector (cols_ == 1)
        Matrix result(n, 1);
        for (size_t i = 0; i < n && i < rows_; ++i) {
            result(i, 0) = data_[i * cols_];
        }
        return result;
    }
    
    // sum() - sum all elements
    inline T sum() const {
        T result = T(0);
        const size_t n = rows_ * cols_;
        for (size_t i = 0; i < n; ++i) {
            result += data_[i];
        }
        return result;
    }
    
    // Static Ones() - create column vector of ones
    static inline Matrix Ones(size_t n) {
        Matrix result(n, 1);
        for (size_t i = 0; i < n; ++i) {
            result(i, 0) = T(1);
        }
        return result;
    }
    
    // Static Zero() - create zero column vector (already zero from constructor, but explicit)
    static inline Matrix Zero(size_t n) {
        Matrix result(n, 1);
        // Constructor already zero-initializes
        return result;
    }
    
    // Dot product (for column vectors)
    inline T dot(const Matrix& other) const {
        // Assume both are column vectors (cols_ == 1)
        T result = T(0);
        const size_t n = (rows_ < other.rows_) ? rows_ : other.rows_;
        for (size_t i = 0; i < n; ++i) {
            result += data_[i * cols_] * other.data_[i * other.cols_];
        }
        return result;
    }

private:
    size_t rows_;
    size_t cols_;
    std::vector<T> data_;
};

// Factory functions (used in main.cpp)
template<typename T>
inline Matrix<T> create_circular_symmetric(size_t n, const std::vector<T>& half_row) {
    Matrix<T> result(n, n);
    // Constructor already zero-initializes all elements
    
    // Build first row
    std::vector<T> first_row(n);
    first_row[0] = T(0);
    
    if (n % 2 == 0) {
        // Even n
        first_row[n / 2] = half_row[half_row.size() - 1];
        for (size_t i = 0; i < n / 2 - 1; ++i) {
            first_row[i + 1] = half_row[i];
            first_row[n - i - 1] = half_row[i];
        }
    } else {
        // Odd n
        for (size_t i = 0; i < n / 2; ++i) {
            first_row[i + 1] = half_row[i];
            first_row[n - i - 1] = half_row[i];
        }
    }
    
    // Create circular matrix from first_row
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            result(i, j) = first_row[(j - i + n) % n];
        }
    }
    
    return result;
}

template<typename T>
inline Matrix<T> create_symmetric(size_t n, const std::vector<T>& upper_triangular) {
    Matrix<T> result(n, n);
    // Constructor already zero-initializes all elements
    
    size_t idx = 0;
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = i; j < n; ++j) {
            T val = upper_triangular[idx];
            result(i, j) = val;
            result(j, i) = val;  // Make symmetric
            ++idx;
        }
    }
    
    return result;
}

// Scalar * Matrix operator (free function)
template<typename T>
inline Matrix<T> operator*(const T& scalar, const Matrix<T>& m) {
    return m * scalar;  // Scalar multiplication is commutative
}


// Convert Matrix<T> to Matrix<double>
template<typename T>
inline Matrix<double> convert_t_to_double(const Matrix<T>& A)
{
    Matrix<double> result(A.rows(), A.cols());
    const size_t n = A.rows() * A.cols();
    const T* src = A.data();
    double* dst = result.data();
    for (size_t i = 0; i < n; ++i) {
        dst[i] = rational_to_double(src[i]);
    }
    return result;
}

// Convert Matrix<small_rational> to Matrix<rational>
inline Matrix<rational> convert_small_to_rational(const Matrix<small_rational>& A)
{
    Matrix<rational> result(A.rows(), A.cols());
    const size_t n = A.rows() * A.cols();
    const small_rational* src = A.data();
    rational* dst = result.data();
    for (size_t i = 0; i < n; ++i) {
        dst[i] = small_to_rational(src[i]);
    }
    return result;
}

// String representation
template<typename T>
inline std::string matrix_to_log(const Matrix<T>& A)
{
    std::ostringstream oss;
    for (size_t i = 0; i < A.rows(); ++i) {
        oss << "\t\t\t";  // Add tabs before each line
        for (size_t j = 0; j < A.cols(); ++j) {
            oss << A(i, j) << ",";
        }
        if (i < A.rows() - 1) {
            oss << std::endl;
        }
    }
    return oss.str();
}

// Specialized string representation for Matrix<small_rational>
// Uses to_string overload for boost::rational which handles denominator==1 logic
template<>
inline std::string matrix_to_log(const Matrix<small_rational>& A)
{
    std::ostringstream oss;
    for (size_t i = 0; i < A.rows(); ++i) {
        oss << "\t\t\t";  // Add tabs before each line
        for (size_t j = 0; j < A.cols(); ++j) {
            oss << to_string(A(i, j)) << ",";
        }
        if (i < A.rows() - 1) {
            oss << std::endl;
        }
    }
    return oss.str();
}

// Extract principal submatrix from matrix using bitset64 mask
// Optimized: only iterates over SET bits using for_each_set_bit_no_exit for better performance
// NOTE: submatrix must already be correctly sized (support_size x support_size) before calling!
template<typename T>
inline void principal_submatrix(const Matrix<T>& A, size_t /*dimension*/, const bitset64& support, size_t /*support_size*/, Matrix<T>& submatrix)
{
    // Only iterate over SET bits for efficiency
    size_t row = 0;
    bs64::for_each_set_bit_no_exit(support, [&](unsigned i) {
        size_t col = 0;
        bs64::for_each_set_bit_no_exit(support, [&](unsigned j) {
            submatrix(row, col) = A(static_cast<size_t>(i), static_cast<size_t>(j));
            ++col;
        });
        ++row;
    });
}

// Check if matrix is positive definite (LDLT decomposition for rational)
template<typename T>
inline bool is_positive_definite_rational(const Matrix<T>& A)
{
    // LDLT decomposition for rational numbers
    size_t n = A.rows();
    Matrix<T> D(n, n);
    Matrix<T> L(n, n);
    
    // Initialize L as identity
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            L(i, j) = (i == j) ? T(1) : T(0);
        }
    }

    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < i; ++j) {
            T aSum = T(0);
            for (size_t k = 0; k < j; ++k) {
                aSum += L(i, k) * L(j, k) * D(k, k);
            }
            L(i, j) = (T(1) / D(j, j)) * (A(i, j) - aSum);
        }
        T bSum = T(0);
        for (size_t k = 0; k < i; ++k) {
            bSum += L(i, k) * L(i, k) * D(k, k);
        }
        D(i, i) = A(i, i) - bSum;
        if (D(i, i) <= T(0)) {
            return false;
        }
    }
    return true;
}

// Check if all entries of a matrix are greater than zero
template<typename T>
inline bool all_entries_greater_zero(const Matrix<T>& A) {
    const size_t n = A.rows() * A.cols();
    const T* data = A.data();
    for (size_t i = 0; i < n; ++i) {
        if (data[i] <= T(0)) {
            return false;
        }
    }
    return true;
}

} // namespace rational_linalg

#endif // RATIONAL_LINALG_MATRIX_HPP

