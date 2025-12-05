#ifndef ESSFINDER_H
#define ESSFINDER_H

#include <vector>
#include <string>
#include <fstream>
#include <bitset>

#include "../include/Matrix.hpp"
#include "../include/Candidate.hpp"


class EssFinder
{
    public:

        EssFinder(Matrix<rational> &matrix, bool with_candidates = false, bool exact = false, bool full_support = false, bool with_log = false);

        size_t ess_count = 0;
        std::vector<Candidate> candidates;

        Matrix<rational> game_matrix;
        Matrix<double> game_matrix_double;

        size_t dimension;

        bool conf_with_candidates;
        bool conf_exact;
        bool conf_full_support;
        bool conf_with_log;

    private:

        Candidate _c;

        std::vector< std::vector<uint64_t>> _supports;
		std::vector<bool> _coprime_sizes;

		std::ofstream* _log;

		void search_support_size(size_t support_size);
        bool find_candidate_double(Matrix<double> &le_matrix);
        bool find_candidate_rational(Matrix<rational> &le_matrix);
        void check_stability();
};

#endif // ESSFINDER_H
