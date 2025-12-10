// Migrated header for modern include path
#include <vector>
#include <string>
#include <memory>

#include <spdlog/spdlog.h>
#include <spdlog/sinks/rotating_file_sink.h>

#include <rational_linalg/matrix.hpp>
#include <rational_linalg/types_rational.hpp>
#include <fracessa/candidate.hpp>
#include <fracessa/bitset64.hpp>
#include <fracessa/supports.hpp>
#include <boost/math/special_functions/binomial.hpp>
#include <boost/integer/common_factor.hpp>


class fracessa
{
    public:

        fracessa(const rational_linalg::Matrix<small_rational>& matrix, bool is_cs, bool with_candidates = false, bool exact = false, bool full_support = false, bool with_log = false, int matrix_id = -1);

        size_t ess_count_ = 0;
        std::vector<candidate> candidates_;

    private:

        rational_linalg::Matrix<small_rational> game_small_;
        rational_linalg::Matrix<rational> game_rational_;  // Converted to arbitrary precision when overflow occurs
        rational_linalg::Matrix<double> game_double_;
        rational_linalg::Matrix<double> subgame_augmented_double_;      // Augmented matrix [A|b] of size (n x n+1)
        rational_linalg::Matrix<small_rational> subgame_augmented_small_;  // Augmented matrix [A|b] of size (n x n+1)
        rational_linalg::Matrix<rational> subgame_augmented_rational_;     // Augmented matrix [A|b] of size (n x n+1)
        rational_linalg::Matrix<small_rational> bee_small_;  // Matrix for copositivity check (small_rational)
        rational_linalg::Matrix<rational> bee_rational_;     // Matrix for copositivity check (rational)
        bool use_small_ = true;  // Flag to track if we're using small_rational (true) or rational (false)

        size_t dimension_;
        bool is_cs_;
        int matrix_id_;

        bool conf_with_candidates_;
        bool conf_exact_;
        bool conf_full_support_;
        bool conf_with_log_;

        
        candidate candidate_;

        Supports supports_;

        std::shared_ptr<spdlog::logger> logger_;

        void search_one_support(const bitset64& support, size_t support_size, bool is_cs_and_coprime = false);
        // COMMENTED OUT: Optimized and find_candidate functions
        // bool find_candidate_double_optimized(rational_linalg::Matrix<double>& A_SS);
        
        // Templated function for all types (double, small_rational, rational)
        template<typename T>
        bool find_candidate(const rational_linalg::Matrix<T>& game_matrix, const bitset64& support, size_t support_size, rational_linalg::Matrix<T>& Ab);
        
        // Templated check_stability function for both rational and small_rational
        template<typename T>
        void check_stability(const rational_linalg::Matrix<T>& game_matrix, rational_linalg::Matrix<T>& Bee);

};

