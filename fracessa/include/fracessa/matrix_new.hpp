#ifndef MATRIX_REFACTORED_HPP
#define MATRIX_REFACTORED_HPP

#include <vector>
#include <stdexcept>
#include <iostream>
#include <algorithm>
#include <cmath>
#include <string>
#include <functional>
#include <numeric>
#include <type_traits>
#include <concepts>

// Forward declaration for helper types
namespace fracessa {
    // Assuming rational type exists
    class rational {
    public:
        rational() = default;
        rational(int num) : numerator(num), denominator(1) {}
        rational(int num, int den) : numerator(num), denominator(den) {}
        
        double to_double() const { return static_cast<double>(numerator) / denominator; }
        std::string to_string() const;
        
        rational operator+(const rational&) const;
        rational operator*(const rational&) const;
        // ... other operators
        
        int numerator = 0;
        int denominator = 1;
    };
    
    // Helper functions
    size_t support_size_from_int(uint64_t support);
    uint64_t shift_right(uint64_t x, size_t dimension);
}

// ========================================
// MODERN MATRIX CLASS IMPLEMENTATION
// ========================================

// Concepts for template constraints
template <typename T>
concept Numeric = std::is_arithmetic_v<T> || 
    std::is_same_v<T, fracessa::rational> ||
    requires(T t) {
        typename T::value_type;
        t = T{0};
        t + t; t * t; t / t;
        static_cast<double>(t);
    };

// Custom exception for matrix operations
class matrix_error : public std::runtime_error {
public:
    explicit matrix_error(const std::string& msg) 
        : std::runtime_error(msg) {}
    
    explicit matrix_error(const char* msg) 
        : std::runtime_error(msg) {}
};

// Main matrix class template
template <Numeric T>
class matrix {
public:
    // ========================================
    // TYPE DEFINITIONS
    // ========================================
    using value_type = T;
    using size_type = size_t;
    using reference = T&;
    using const_reference = const T&;
    using pointer = T*;
    using const_pointer = const T*;
    using iterator = typename std::vector<T>::iterator;
    using const_iterator = typename std::vector<T>::const_iterator;
    using reverse_iterator = typename std::vector<T>::reverse_iterator;
    using const_reverse_iterator = typename std::vector<T>::const_reverse_iterator;

private:
    // ========================================
    // MEMBER VARIABLES
    // ========================================
    std::vector<T> _data;
    size_type _rows = 0;
    size_type _cols = 0;
    bool _is_circular_symmetric = false;

public:
    // ========================================
    // CONSTRUCTORS & DESTRUCTORS
    // ========================================
    
    // Default constructor
    matrix() noexcept = default;
    
    // Size constructor
    matrix(size_type rows, size_type cols) {
        resize(rows, cols);
    }
    
    // Size + value constructor
    matrix(size_type rows, size_type cols, const T& value) {
        resize(rows, cols);
        fill(value);
    }
    
    // Initializer list constructor (2D)
    matrix(std::initializer_list<std::initializer_list<T>> list) {
        if (list.size() == 0) return;
        
        _rows = list.size();
        _cols = list.begin()->size();
        
        _data.reserve(_rows * _cols);
        for (const auto& row : list) {
            if (row.size() != _cols) {
                throw matrix_error("Inconsistent row sizes in initializer list");
            }
            _data.insert(_data.end(), row.begin(), row.end());
        }
    }
    
    // Copy constructor
    matrix(const matrix&) = default;
    
    // Move constructor
    matrix(matrix&& other) noexcept 
        : _data(std::move(other._data))
        , _rows(other._rows)
        , _cols(other._cols)
        , _is_circular_symmetric(other._is_circular_symmetric) {
        other._rows = 0;
        other._cols = 0;
        other._is_circular_symmetric = false;
    }
    
    // Destructor
    ~matrix() = default;

    // ========================================
    // ASSIGNMENT OPERATORS
    // ========================================
    
    // Copy assignment
    matrix& operator=(const matrix&) = default;
    
    // Move assignment
    matrix& operator=(matrix&& other) noexcept {
        if (this != &other) {
            _data = std::move(other._data);
            _rows = other._rows;
            _cols = other._cols;
            _is_circular_symmetric = other._is_circular_symmetric;
            
            other._rows = 0;
            other._cols = 0;
            other._is_circular_symmetric = false;
        }
        return *this;
    }
    
    // Initializer list assignment
    matrix& operator=(std::initializer_list<std::initializer_list<T>> list) {
        matrix temp(list);
        *this = std::move(temp);
        return *this;
    }

    // ========================================
    // BASIC ACCESSORS
    // ========================================
    
    size_type rows() const noexcept { return _rows; }
    size_type cols() const noexcept { return _cols; }
    size_type size() const noexcept { return _data.size(); }
    bool empty() const noexcept { return _data.empty(); }
    bool is_square() const noexcept { return _rows == _cols; }
    bool is_circular_symmetric() const noexcept { return _is_circular_symmetric; }
    
    // Mark as circular symmetric
    void set_circular_symmetric(bool value = true) noexcept {
        _is_circular_symmetric = value;
    }

    // ========================================
    // ELEMENT ACCESS
    // ========================================
    
    // Bounds-checked access (safe)
    reference at(size_type row, size_type col) {
        check_bounds(row, col);
        return _data[row * _cols + col];
    }
    
    const_reference at(size_type row, size_type col) const {
        check_bounds(row, col);
        return _data[row * _cols + col];
    }
    
    // Unchecked access (fast, for performance-critical code)
    reference operator()(size_type row, size_type col) noexcept {
        return _data[row * _cols + col];
    }
    
    const_reference operator()(size_type row, size_type col) const noexcept {
        return _data[row * _cols + col];
    }
    
    // Row access (returns pointer to row start)
    T* row_data(size_type row) noexcept {
        return _data.data() + row * _cols;
    }
    
    const T* row_data(size_type row) const noexcept {
        return _data.data() + row * _cols;
    }

    // ========================================
    // ITERATORS
    // ========================================
    
    iterator begin() noexcept { return _data.begin(); }
    const_iterator begin() const noexcept { return _data.begin(); }
    const_iterator cbegin() const noexcept { return _data.cbegin(); }
    
    iterator end() noexcept { return _data.end(); }
    const_iterator end() const noexcept { return _data.end(); }
    const_iterator cend() const noexcept { return _data.cend(); }
    
    reverse_iterator rbegin() noexcept { return _data.rbegin(); }
    const_reverse_iterator rbegin() const noexcept { return _data.rbegin(); }
    
    reverse_iterator rend() noexcept { return _data.rend(); }
    const_reverse_iterator rend() const noexcept { return _data.rend(); }

    // ========================================
    // RESIZING & CAPACITY
    // ========================================
    
    void resize(size_type rows, size_type cols) {
        if (rows == 0 || cols == 0) {
            throw matrix_error("Matrix dimensions must be positive");
        }
        
        // Prevent unreasonably large matrices
        if (rows > 100000 || cols > 100000 || 
            rows * cols > 100000000ULL) {  // 100M elements max
            throw matrix_error("Matrix dimensions too large");
        }
        
        _rows = rows;
        _cols = cols;
        
        try {
            _data.resize(rows * cols);
        } catch (const std::bad_alloc&) {
            throw matrix_error("Failed to allocate memory for matrix");
        }
    }
    
    void reserve(size_type capacity) {
        _data.reserve(capacity);
    }
    
    void shrink_to_fit() {
        _data.shrink_to_fit();
    }
    
    size_type capacity() const noexcept {
        return _data.capacity();
    }

    // ========================================
    // ELEMENT MODIFICATION
    // ========================================
    
    void fill(const T& value) noexcept {
        std::fill(_data.begin(), _data.end(), value);
    }
    
    void fill_diagonal(const T& value) {
        if (!is_square()) {
            throw matrix_error("Matrix must be square for diagonal operations");
        }
        for (size_type i = 0; i < _rows; ++i) {
            (*this)(i, i) = value;
        }
    }
    
    void swap(matrix& other) noexcept {
        _data.swap(other._data);
        std::swap(_rows, other._rows);
        std::swap(_cols, other._cols);
        std::swap(_is_circular_symmetric, other._is_circular_symmetric);
    }

    // ========================================
    // ARITHMETIC OPERATIONS
    // ========================================
    
    // Compound assignment operators
    matrix& operator+=(const matrix& other) {
        check_dimensions_match(other, "addition");
        std::transform(_data.begin(), _data.end(), other._data.begin(),
                      _data.begin(), std::plus<T>{});
        return *this;
    }
    
    matrix& operator-=(const matrix& other) {
        check_dimensions_match(other, "subtraction");
        std::transform(_data.begin(), _data.end(), other._data.begin(),
                      _data.begin(), std::minus<T>{});
        return *this;
    }
    
    matrix& operator*=(const T& scalar) noexcept {
        std::transform(_data.begin(), _data.end(), _data.begin(),
                      [scalar](const T& val) { return val * scalar; });
        return *this;
    }
    
    matrix& operator/=(const T& scalar) {
        if constexpr (std::is_integral_v<T>) {
            if (scalar == T{0}) throw matrix_error("Division by zero");
        } else {
            if (std::abs(static_cast<double>(scalar)) < 1e-15) {
                throw matrix_error("Division by near-zero scalar");
            }
        }
        
        std::transform(_data.begin(), _data.end(), _data.begin(),
                      [scalar](const T& val) { return val / scalar; });
        return *this;
    }
    
    // Binary operators
    matrix operator+(const matrix& other) const {
        matrix result = *this;
        result += other;
        return result;
    }
    
    matrix operator-(const matrix& other) const {
        matrix result = *this;
        result -= other;
        return result;
    }
    
    matrix operator*(const T& scalar) const {
        matrix result = *this;
        result *= scalar;
        return result;
    }
    
    matrix operator/(const T& scalar) const {
        matrix result = *this;
        result /= scalar;
        return result;
    }
    
    // Matrix multiplication
    matrix operator*(const matrix& other) const {
        if (_cols != other._rows) {
            throw matrix_error("Matrix dimensions incompatible for multiplication: " +
                             std::to_string(_rows) + "x" + std::to_string(_cols) + " * " +
                             std::to_string(other._rows) + "x" + std::to_string(other._cols));
        }
        
        matrix result(_rows, other._cols, T{0});
        
        for (size_type i = 0; i < _rows; ++i) {
            for (size_type j = 0; j < other._cols; ++j) {
                for (size_type k = 0; k < _cols; ++k) {
                    result(i, j) += (*this)(i, k) * other(k, j);
                }
            }
        }
        
        return result;
    }

    // ========================================
    // COMPARISON OPERATORS
    // ========================================
    
    bool operator==(const matrix& other) const noexcept {
        return _rows == other._rows && _cols == other._cols && 
               _data == other._data;
    }
    
    bool operator!=(const matrix& other) const noexcept {
        return !(*this == other);
    }

    // ========================================
    // MATRIX TRANSFORMATIONS
    // ========================================
    
    matrix transpose() const {
        matrix result(_cols, _rows);
        for (size_type i = 0; i < _rows; ++i) {
            for (size_type j = 0; j < _cols; ++j) {
                result(j, i) = (*this)(i, j);
            }
        }
        return result;
    }
    
    matrix submatrix(size_type start_row, size_type start_col,
                     size_type end_row, size_type end_col) const {
        if (start_row > end_row || start_col > end_col ||
            end_row > _rows || end_col > _cols) {
            throw matrix_error("Invalid submatrix bounds");
        }
        
        size_type new_rows = end_row - start_row;
        size_type new_cols = end_col - start_col;
        
        matrix result(new_rows, new_cols);
        for (size_type i = 0; i < new_rows; ++i) {
            for (size_type j = 0; j < new_cols; ++j) {
                result(i, j) = (*this)(start_row + i, start_col + j);
            }
        }
        return result;
    }
    
    // Extract submatrix based on support bitmask (for ESS finding)
    matrix extract_support(uint64_t support_mask) const {
        size_type support_size = fracessa::support_size_from_int(support_mask);
        matrix result(support_size, support_size);
        
        size_type row_idx = 0;
        for (size_type i = 0; i < _rows; ++i) {
            if ((support_mask & (1ULL << i)) == 0) continue;
            
            size_type col_idx = 0;
            for (size_type j = 0; j < _cols; ++j) {
                if ((support_mask & (1ULL << j)) == 0) continue;
                
                result(row_idx, col_idx) = (*this)(i, j);
                ++col_idx;
            }
            ++row_idx;
        }
        
        return result;
    }

    // ========================================
    // ELEMENT-WISE OPERATIONS
    // ========================================
    
    // Apply function to each element
    matrix apply(std::function<T(T)> func) const {
        matrix result(_rows, _cols);
        std::transform(_data.begin(), _data.end(), result._data.begin(), func);
        return result;
    }
    
    // Apply function with element indices
    matrix apply(std::function<T(T, size_type, size_type)> func) const {
        matrix result(_rows, _cols);
        for (size_type i = 0; i < _rows; ++i) {
            for (size_type j = 0; j < _cols; ++j) {
                result(i, j) = func((*this)(i, j), i, j);
            }
        }
        return result;
    }

    // ========================================
    // MATHEMATICAL PROPERTIES
    // ========================================

    // Unified positive definite check (works for both rational and floating-point types)
    bool is_positive_definite() const {
        if (!is_square()) {
            throw matrix_error("Positive definite check requires square matrix");
        }

        // Use template specialization to choose appropriate algorithm
        return is_positive_definite_impl(*this);
    }
    
    // Determinant (only for rational matrices to avoid precision issues)
    template <typename U = T>
    std::enable_if_t<std::is_same_v<U, fracessa::rational>, fracessa::rational> 
    determinant() const {
        if (!is_square()) {
            throw matrix_error("Determinant requires square matrix");
        }
        
        // Use LU decomposition for numerical stability
        return lu_determinant();
    }
    
    // Trace (sum of diagonal elements)
    T trace() const {
        if (!is_square()) {
            throw matrix_error("Trace requires square matrix");
        }
        
        T result = T{0};
        for (size_type i = 0; i < _rows; ++i) {
            result += (*this)(i, i);
        }
        return result;
    }
    
    // Check if matrix contains only positive elements
    bool all_positive() const noexcept {
        return std::all_of(_data.begin(), _data.end(), 
                          [](const T& val) { return val > T{0}; });
    }
    
    // Frobenius norm (for numerical stability checks)
    double frobenius_norm() const {
        double sum_squares = 0.0;
        for (const auto& val : _data) {
            double real_val = static_cast<double>(val);
            sum_squares += real_val * real_val;
        }
        return std::sqrt(sum_squares);
    }

    // ========================================
    // CONVERSION METHODS
    // ========================================
    
    // Convert to double matrix (for numerical computations)
    matrix<double> to_double() const {
        matrix<double> result(_rows, _cols);
        std::transform(_data.begin(), _data.end(), result._data.begin(),
                      [](const T& val) { return static_cast<double>(val); });
        return result;
    }
    
    // Convert to rational matrix (if applicable)
    template <typename U = T>
    std::enable_if_t<!std::is_same_v<U, fracessa::rational>, matrix<fracessa::rational>>
    to_rational() const {
        matrix<fracessa::rational> result(_rows, _cols);
        std::transform(_data.begin(), _data.end(), result._data.begin(),
                      [](const T& val) { return fracessa::rational(static_cast<int>(val)); });
        return result;
    }

    // ========================================
    // OUTPUT METHODS
    // ========================================
    
    void print(std::ostream& os = std::cout, 
               const std::string& separator = ", ") const {
        for (size_type i = 0; i < _rows; ++i) {
            for (size_type j = 0; j < _cols; ++j) {
                if (j > 0) os << separator;
                print_element(os, (*this)(i, j));
            }
            os << '\n';
        }
    }
    
    void print_csv(std::ostream& os = std::cout) const {
        print(os, ",");
    }
    
    std::string to_string() const {
        std::ostringstream oss;
        print(oss);
        return oss.str();
    }
    
    // Friend function for stream output
    friend std::ostream& operator<<(std::ostream& os, const matrix& mat) {
        mat.print(os);
        return os;
    }

    // ========================================
    // FACTORY METHODS
    // ========================================
    
    static matrix zeros(size_type rows, size_type cols) {
        return matrix(rows, cols, T{0});
    }
    
    static matrix ones(size_type rows, size_type cols) {
        return matrix(rows, cols, T{1});
    }
    
    static matrix identity(size_type n) {
        matrix result = zeros(n, n);
        result.fill_diagonal(T{1});
        return result;
    }
    
    static matrix diagonal(size_type n, const T& value) {
        matrix result = zeros(n, n);
        result.fill_diagonal(value);
        return result;
    }
    
    // Random matrix (for testing)
    static matrix random(size_type rows, size_type cols, T min_val = T{0}, T max_val = T{1}) {
        matrix result(rows, cols);
        std::random_device rd;
        std::mt19937 gen(rd());
        
        if constexpr (std::is_floating_point_v<T>) {
            std::uniform_real_distribution<T> dist(min_val, max_val);
            for (auto& val : result._data) {
                val = dist(gen);
            }
        } else if constexpr (std::is_integral_v<T>) {
            std::uniform_int_distribution<T> dist(min_val, max_val);
            for (auto& val : result._data) {
                val = dist(gen);
            }
        }
        
        return result;
    }

private:
    // ========================================
    // PRIVATE HELPER METHODS
    // ========================================

    // Single unified positive definite check using optimized algorithm
    template <typename MatrixType>
    static bool is_positive_definite_impl(const MatrixType& mat) {
        using ValueType = typename MatrixType::value_type;

        // For rational types: Use LDLT decomposition (exact arithmetic)
        if constexpr (std::is_same_v<ValueType, fracessa::rational>) {
            return is_positive_definite_unified<ValueType, true>(mat);
        }
        // For floating-point types: Use Cholesky decomposition (numerical stability)
        else if constexpr (std::is_floating_point_v<ValueType>) {
            return is_positive_definite_unified<ValueType, false>(mat);
        }
        // For other types: Convert to double and use Cholesky
        else {
            size_t n = mat.rows();
            matrix<double> double_mat(n, n);
            for (size_t i = 0; i < n; ++i) {
                for (size_t j = 0; j < n; ++j) {
                    double_mat(i, j) = static_cast<double>(mat(i, j));
                }
            }
            return is_positive_definite_unified<double, false>(double_mat);
        }
    }

    // SINGLE UNIFIED ALGORITHM for both LDLT and Cholesky decomposition
    // UseLDLT template parameter controls the variant:
    // - true: LDLT decomposition (rational numbers, exact)
    // - false: Cholesky decomposition (floating-point, numerical)
    template <typename ValueType, bool UseLDLT>
    static bool is_positive_definite_unified(const matrix<ValueType>& mat) {
        const size_t n = mat.rows();

        // Working matrix L (lower triangular)
        matrix<ValueType> L = matrix<ValueType>::zeros(n, n);

        // For LDLT: D matrix stores diagonal scaling factors
        // For Cholesky: D is not used (we store sqrt directly in L)
        matrix<ValueType> D = UseLDLT ? matrix<ValueType>::zeros(n, n) :
                                       matrix<ValueType>::zeros(0, 0);

        // Initialize D to identity for LDLT
        if constexpr (UseLDLT) {
            for (size_t i = 0; i < n; ++i) {
                D(i, i) = ValueType{1};
            }
        }

        // Main decomposition loop - SAME for both algorithms!
        for (size_t i = 0; i < n; ++i) {
            // Step 1: Compute L[i][j] for j < i (off-diagonal elements)
            for (size_t j = 0; j < i; ++j) {
                ValueType sum = ValueType{0};

                // Inner sum: Σ_{k=0}^{j-1} L[i][k] * L[j][k] * D[k][k] (LDLT) or just L[i][k] * L[j][k] (Cholesky)
                for (size_t k = 0; k < j; ++k) {
                    if constexpr (UseLDLT) {
                        sum += L(i, k) * L(j, k) * D(k, k);
                    } else {
                        sum += L(i, k) * L(j, k);
                    }
                }

                // L[i][j] = (A[i][j] - sum) / D[j][j] (LDLT) or L[j][j] (Cholesky)
                ValueType denominator = UseLDLT ? D(j, j) : L(j, j);
                L(i, j) = (ValueType{1} / denominator) * (mat(i, j) - sum);
            }

            // Step 2: Compute diagonal element
            ValueType sum = ValueType{0};

            // Sum: Σ_{k=0}^{i-1} L[i][k]² * D[k][k] (LDLT) or just L[i][k]² (Cholesky)
            for (size_t k = 0; k < i; ++k) {
                if constexpr (UseLDLT) {
                    sum += L(i, k) * L(i, k) * D(k, k);
                } else {
                    sum += L(i, k) * L(i, k);
                }
            }

            ValueType diagonal_element = mat(i, i) - sum;

            // Step 3: Check positive definiteness and store result
            if constexpr (UseLDLT) {
                // LDLT: Check if D[i][i] > 0 (exact check for rational numbers)
                if (diagonal_element <= ValueType{0}) {
                    return false;
                }
                D(i, i) = diagonal_element;
            } else {
                // Cholesky: Check if diagonal_element > 0, then store sqrt
                if constexpr (std::is_floating_point_v<ValueType>) {
                    if (diagonal_element < ValueType{1e-12}) {  // Robust numerical tolerance
                        return false;
                    }
                    L(i, i) = std::sqrt(diagonal_element);
                } else {
                    // For non-floating types in Cholesky mode
                    if (diagonal_element <= ValueType{0}) {
                        return false;
                    }
                    L(i, i) = diagonal_element;  // Approximate for non-floating types
                }
            }
        }

        return true;  // Matrix passed all positive definiteness checks
    }
    
    void check_bounds(size_type row, size_type col) const {
        if (row >= _rows || col >= _cols) {
            throw std::out_of_range(
                "Matrix access out of bounds: [" + std::to_string(row) + 
                "," + std::to_string(col) + "] in " + std::to_string(_rows) + 
                "x" + std::to_string(_cols) + " matrix"
            );
        }
    }
    
    void check_dimensions_match(const matrix& other, const std::string& operation) const {
        if (_rows != other._rows || _cols != other._cols) {
            throw matrix_error(
                "Matrix dimensions must match for " + operation + ": " +
                std::to_string(_rows) + "x" + std::to_string(_cols) + " vs " +
                std::to_string(other._rows) + "x" + std::to_string(other._cols)
            );
        }
    }
    
    void print_element(std::ostream& os, const T& val) const {
        if constexpr (std::is_same_v<T, fracessa::rational>) {
            os << val.to_string();
        } else {
            os << val;
        }
    }
    
    // LU decomposition for determinant calculation
    fracessa::rational lu_determinant() const {
        matrix<fracessa::rational> lu = *this;  // Copy for decomposition
        fracessa::rational det = fracessa::rational(1);
        int sign = 1;
        
        // Forward elimination
        for (size_type i = 0; i < _rows - 1; ++i) {
            // Find pivot
            size_type pivot_row = i;
            for (size_type k = i + 1; k < _rows; ++k) {
                if (std::abs(lu(k, i).to_double()) > std::abs(lu(pivot_row, i).to_double())) {
                    pivot_row = k;
                }
            }
            
            // Swap rows if needed
            if (pivot_row != i) {
                for (size_type j = 0; j < _cols; ++j) {
                    std::swap(lu(i, j), lu(pivot_row, j));
                }
                sign = -sign;
            }
            
            // Eliminate
            for (size_type k = i + 1; k < _rows; ++k) {
                fracessa::rational factor = lu(k, i) / lu(i, i);
                for (size_type j = i; j < _cols; ++j) {
                    if (i == j) {
                        lu(k, j) = fracessa::rational(0);
                    } else {
                        lu(k, j) = lu(k, j) - factor * lu(i, j);
                    }
                }
            }
        }
        
        // Calculate determinant from diagonal
        for (size_type i = 0; i < _rows; ++i) {
            det = det * lu(i, i);
        }
        
        return (sign == 1) ? det : -det;
    }
};

// ========================================
// NON-MEMBER FUNCTIONS
// ========================================

// Swap specialization
template <Numeric T>
void swap(matrix<T>& lhs, matrix<T>& rhs) noexcept {
    lhs.swap(rhs);
}

// Scalar multiplication (commutative)
template <Numeric T>
matrix<T> operator*(const T& scalar, const matrix<T>& mat) {
    return mat * scalar;
}

// Utility functions
template <Numeric T>
bool is_symmetric(const matrix<T>& mat) {
    if (!mat.is_square()) return false;
    
    for (size_t i = 0; i < mat.rows(); ++i) {
        for (size_t j = i + 1; j < mat.cols(); ++j) {
            if (mat(i, j) != mat(j, i)) return false;
        }
    }
    return true;
}

template <Numeric T>
matrix<T> absolute(const matrix<T>& mat) {
    return mat.apply([](const T& val) {
        using std::abs;
        return abs(val);
    });
}

// ========================================
// TYPE ALIASES FOR COMMON USE CASES
// ========================================

using rational_matrix = matrix<fracessa::rational>;
using double_matrix = matrix<double>;
using int_matrix = matrix<int>;

// ========================================
// BACKWARD COMPATIBILITY HELPERS
// ========================================

// Legacy exception name (for compatibility)
using matrix_exception = matrix_error;

#endif // MATRIX_REFACTORED_HPP