// Migrated header for modern include path
#include <vector>
#include <string>
#include <bitset>
#include <memory>

#include <spdlog/spdlog.h>
#include <spdlog/sinks/rotating_file_sink.h>

#include <fracessa/matrix.hpp>
#include <fracessa/candidate.hpp>
#include <fracessa/bitset64.hpp>
#include <fracessa/rational_eigen.hpp>


class fracessa
{
    public:

        fracessa(const RationalMatrix& matrix, bool is_cs, bool with_candidates = false, bool exact = false, bool full_support = false, bool with_log = false);

        size_t ess_count = 0;
        std::vector<candidate> candidates;

        RationalMatrix game_matrix;
        DoubleMatrix game_matrix_double;

        size_t dimension;
        bool is_cs;

        bool conf_with_candidates;
        bool conf_exact;
        bool conf_full_support;
        bool conf_with_log;

    private:

        candidate _c;

        std::vector< std::vector<bitset64>> _supports;
        std::vector<bool> _coprime_sizes;

        std::shared_ptr<spdlog::logger> _logger;

        void search_support_size(size_t support_size);
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

