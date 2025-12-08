#ifndef RATIONAL_LINALG_MATRIX_HPP
#define RATIONAL_LINALG_MATRIX_HPP

#include <fracessa/rational.hpp>
#include <vector>
#include <cstring>

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

// Type aliases for convenience
using RationalMatrix = Matrix<rational>;
using SmallRationalMatrix = Matrix<small_rational>;

} // namespace rational_linalg

#endif // RATIONAL_LINALG_MATRIX_HPP

