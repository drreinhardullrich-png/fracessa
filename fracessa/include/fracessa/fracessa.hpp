// Migrated header for modern include path
#include <vector>
#include <string>
#include <bitset>
#include <memory>

#include <spdlog/spdlog.h>
#include <spdlog/sinks/rotating_file_sink.h>

#include <fracessa/matrix.hpp>
#include <fracessa/candidate.hpp>


class fracessa
{
    public:

        fracessa(const matrix<rational>& matrix, bool with_candidates = false, bool exact = false, bool full_support = false, bool with_log = false);

        size_t ess_count = 0;
        std::vector<candidate> candidates;

        matrix<rational> game_matrix;
        matrix<double> game_matrix_double;

        size_t dimension;

        bool conf_with_candidates;
        bool conf_exact;
        bool conf_full_support;
        bool conf_with_log;

    private:

        candidate _c;

        std::vector< std::vector<uint64_t>> _supports;
        std::vector<bool> _coprime_sizes;

        std::shared_ptr<spdlog::logger> _logger;

        void search_support_size(size_t support_size);
        bool find_candidate_double(matrix<double> &le_matrix);
        bool find_candidate_rational(matrix<rational> &le_matrix);
        void check_stability();
};

