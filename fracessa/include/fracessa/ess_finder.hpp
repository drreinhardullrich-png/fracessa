// Migrated header for modern include path
#include <vector>
#include <string>
#include <fstream>
#include <bitset>
#include <optional>

#include <fracessa/matrix.hpp>
#include <fracessa/candidate.hpp>
#include <fracessa/logger.hpp>
#include <memory>


class ess_finder
{
    public:

        ess_finder(const matrix<rational>& matrix, bool with_candidates = false, bool exact = false, bool full_support = false, bool with_log = false);

        size_t get_ess_count() const { return ess_count; }
        const std::vector<candidate>& get_candidates() const { return candidates; }

    private:
        size_t ess_count = 0;
        std::vector<candidate> candidates;

        const matrix<rational>& game_matrix;
        matrix<double> game_matrix_double;

        size_t dimension;

        bool conf_with_candidates;
        bool conf_exact;
        bool conf_full_support;
        bool conf_with_log;

        candidate _c;

        std::vector< std::vector<uint64_t>> _supports;
        std::vector<bool> _coprime_sizes;

        std_unique_ptr<Logger> _log;

        void search_support_size(size_t support_size);
        bool find_candidate_double(matrix<double> &le_matrix);
        bool find_candidate_rational(matrix<rational> &le_matrix);
        void check_stability();
};
