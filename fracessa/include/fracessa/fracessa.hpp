// Migrated header for modern include path
#include <vector>
#include <string>
#include <bitset>
#include <memory>

#include <spdlog/spdlog.h>
#include <spdlog/sinks/rotating_file_sink.h>

#include <rational_linalg/matrix.hpp>
#include <fracessa/rational.hpp>
#include <fracessa/candidate.hpp>
#include <fracessa/bitset64.hpp>
#include <fracessa/supports.hpp>
#include <Eigen/Dense>
#include <boost/math/special_functions/binomial.hpp>
#include <boost/integer/common_factor.hpp>

// Type aliases for Eigen double matrices (still used for double operations)
using DoubleMatrix = Eigen::MatrixXd;
using DoubleVector = Eigen::VectorXd;

// Type aliases for rational matrices
using RationalMatrix = rational_linalg::Matrix<rational>;
using RationalVector = rational_linalg::Matrix<rational>;  // Column vector: Matrix<rational>(n, 1)

// Named constants
static constexpr size_t CANDIDATE_RESERVE_MULTIPLIER = 100;

class fracessa
{
    public:

        fracessa(const RationalMatrix& matrix, bool is_cs, bool with_candidates = false, bool exact = false, bool full_support = false, bool with_log = false, int matrix_id = -1);

        size_t ess_count_ = 0;
        std::vector<candidate> candidates_;

    private:

        RationalMatrix game_matrix_;
        DoubleMatrix game_matrix_double_;

        size_t dimension_;
        bool is_cs_;
        int matrix_id_;

        bool conf_with_candidates_;
        bool conf_exact_;
        bool conf_full_support_;
        bool conf_with_log_;

        
        candidate c_;

        std::vector<bitset64> candidates_supports_;
        std::vector<bool> coprime_sizes_;
        Supports supports_;

        std::shared_ptr<spdlog::logger> logger_;

        void search_one_support(const bitset64& support, size_t support_size, bool is_cs_and_coprime = false);
        bool find_candidate_double_optimized(DoubleMatrix& A_SS);
        bool find_candidate_double(DoubleMatrix& A, DoubleVector& b);
        bool find_candidate_rational_optimized(RationalMatrix& A_SS);
        bool find_candidate_rational(RationalMatrix& A, RationalVector& b);
        void check_stability();
        
        // Helper functions to build full solution vector from support
        bool build_solution_vector(const DoubleVector& solution, DoubleVector& solution_full_n);
        bool build_solution_vector(const RationalVector& solution, RationalVector& solution_full_n);
        
        // Helper functions to check constraints p'Ap<=v for rows not in support
        bool check_constraints(const DoubleVector& solution, const DoubleVector& solution_full_n);
        bool check_constraints(const RationalVector& solution, const RationalVector& solution_full_n, bitset64& extended_support);
};

