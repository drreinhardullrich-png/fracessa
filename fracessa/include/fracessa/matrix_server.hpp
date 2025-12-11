#ifndef MATRIX_SERVER_HPP
#define MATRIX_SERVER_HPP

#include <rational_linalg/matrix.hpp>
#include <rational_linalg/types_rational.hpp>
#include <fracessa/bitset64.hpp>
#include <utility>

/**
 * MatrixServer - Centralized matrix storage and operations for fracessa
 * 
 * This class manages all matrix representations (small_rational, rational, double)
 * and provides the matrix-building functionality needed by find_candidate and check_stability.
 */
class MatrixServer {
public:
    /**
     * Constructor - initializes all matrix representations from the input game matrix
     * @param game_matrix The input game matrix (small_rational)
     */
    MatrixServer(const rational_linalg::Matrix<small_rational>& game_matrix)
        : game_small_(game_matrix)
    {
        // Always create both game matrices (faster than checking every time)
        game_double_ = rational_linalg::convert_t_to_double(game_small_);
        game_rational_ = rational_linalg::convert_small_to_rational(game_small_);
        dimension_ = game_matrix.rows();
    }

    // =========================================================================
    // Public API - Get matrices as needed
    // =========================================================================


    template<typename T>
    rational_linalg::Matrix<T>& get_linear_system(const bitset64& support, size_t support_size) {
        auto& Ab = get_augmented<T>();
        build_kkt_augmented<T>(support, support_size, Ab);
        return Ab;
    }

    template<typename T>
    rational_linalg::Matrix<T>& get_bee_matrix(const bitset64& extended_support_reduced, size_t m) {
        auto& Bee = get_bee<T>();
        const auto& game = get_game<T>();
        
        size_t extended_support_size_reduced = bs64::count(extended_support_reduced);
        
        // Resize if needed (constructor zero-initializes)
        if (Bee.rows() != extended_support_size_reduced) {
            Bee = rational_linalg::Matrix<T>(extended_support_size_reduced, extended_support_size_reduced);
        }

        size_t row = 0;
        size_t column = 0;
        for (size_t i = 0; i < dimension_; i++) {
            if (bs64::test(extended_support_reduced, i)) {
                column = 0;
                for (size_t j = 0; j < i + 1; j++) {
                    if (bs64::test(extended_support_reduced, j)) {
                        // Bee formula from Bomze 1992
                        Bee(row, column) = Bee(column, row) = 
                            game(m, j) + game(j, m) + game(i, m) + game(m, i) -
                            game(i, j) - game(j, i) - T(2) * game(m, m);
                        column += 1;
                    }
                }
                row += 1;
            }
        }
        
        return Bee;
    }

    template<typename T>
    const rational_linalg::Matrix<T>& get_game_matrix() const {
        return get_game<T>();
    }

private:
    // =========================================================================
    // Internal Matrix Builders (private)
    // =========================================================================
    
    template<typename T>
    void build_kkt_augmented(const bitset64& support, size_t support_size, rational_linalg::Matrix<T>& Ab) const {
        const auto& game = get_game<T>();
        
        // Resize if needed (constructor zero-initializes)
        if (Ab.rows() != support_size + 1) {
            Ab = rational_linalg::Matrix<T>(support_size + 1, support_size + 2);
        }

        // Fill all rows in one pass
        size_t ab_row = 0;
        for (size_t i = 0; i < dimension_; ++i) {
            if (bs64::test(support, i)) {
                // Fill submatrix columns (0 to support_size-1) from game_matrix
                size_t ab_col = 0;
                for (size_t j = 0; j < dimension_; ++j) {
                    if (bs64::test(support, j)) {
                        Ab(ab_row, ab_col) = game(i, j);
                        ab_col++;
                    }
                }
                // Column support_size: -1 (last column of A)
                Ab(ab_row, support_size) = T(-1);
                // Column support_size + 1: 0 (b vector, initially zero)
                Ab(ab_row, support_size + 1) = T(0);
                ab_row++;
            }
        }
        
        // Fill last row: all 1s in columns 0..support_size-1, 0 in column support_size, 1 in column n
        for (size_t i = 0; i < support_size; ++i) {
            Ab(support_size, i) = T(1);
        }
        Ab(support_size, support_size) = T(0);      // Column support_size (last column of A)
        Ab(support_size, support_size + 1) = T(1);  // Column support_size + 1 (b vector, last element is 1)
    }

    // =========================================================================
    // Internal Helpers
    // =========================================================================
    
    /**
     * Get the appropriate game matrix for the template type
     */
    template<typename T>
    const rational_linalg::Matrix<T>& get_game() const {
        if constexpr (std::is_same_v<T, double>) {
            return game_double_;
        } else if constexpr (std::is_same_v<T, small_rational>) {
            return game_small_;
        } else {
            return game_rational_;
        }
    }
    
    /**
     * Get the appropriate augmented matrix reference
     */
    template<typename T>
    rational_linalg::Matrix<T>& get_augmented() {
        if constexpr (std::is_same_v<T, double>) {
            return subgame_augmented_double_;
        } else if constexpr (std::is_same_v<T, small_rational>) {
            return subgame_augmented_small_;
        } else {
            return subgame_augmented_rational_;
        }
    }
    
    /**
     * Get the appropriate Bee matrix reference
     */
    template<typename T>
    rational_linalg::Matrix<T>& get_bee() {
        if constexpr (std::is_same_v<T, small_rational>) {
            return bee_small_;
        } else {
            return bee_rational_;
        }
    }
    

    // =========================================================================
    // Member Variables
    // =========================================================================
    
    // Game matrices
    rational_linalg::Matrix<small_rational> game_small_;
    rational_linalg::Matrix<rational> game_rational_;
    rational_linalg::Matrix<double> game_double_;
    
    // Augmented matrices for linear solver (KKT system)
    rational_linalg::Matrix<double> subgame_augmented_double_;
    rational_linalg::Matrix<small_rational> subgame_augmented_small_;
    rational_linalg::Matrix<rational> subgame_augmented_rational_;
    
    // Bee matrices for copositivity check
    rational_linalg::Matrix<small_rational> bee_small_;
    rational_linalg::Matrix<rational> bee_rational_;
    
    // State

    size_t dimension_;
};

#endif // MATRIX_SERVER_HPP

