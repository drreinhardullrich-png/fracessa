#include <fracessa/fracessa.hpp>
#include <fracessa/bitset64.hpp>

fracessa::fracessa(const RationalMatrix& matrix, bool is_cs, bool with_candidates, bool exact, bool full_support, bool with_log)
{

    game_matrix = matrix;
    this->is_cs = is_cs;
    conf_with_candidates = with_candidates;
    conf_exact = exact;
    conf_full_support = full_support;
    conf_with_log = with_log;

    if (!conf_exact)
        game_matrix_double = matrix_ops::to_double(game_matrix);

    dimension = matrix.rows();
    _supports = std::vector< std::vector<bitset64>>(dimension);

    if (conf_with_candidates)
        candidates.reserve(CANDIDATE_RESERVE_MULTIPLIER * dimension);

    if (conf_with_log) {
        auto rotating_sink = std::make_shared<spdlog::sinks::rotating_file_sink_mt>(
            "log/fracessa.log", 20*1024*1024, 5);  // 20MB, 5 files
        _logger = std::make_shared<spdlog::logger>("fracessa", rotating_sink);
        _logger->set_pattern("[%Y-%m-%d %H:%M:%S.%e] [%l] %v");
        _logger->set_level(spdlog::level::info);
        
        _logger->info("n={}", dimension);
        _logger->info("game matrix:\n{}", matrix_ops::to_string(game_matrix));
    }

	for (size_t i = 0; i < dimension; i++) {
		_supports[i].reserve(static_cast<uint64_t>(boost::math::binomial_coefficient<double>(dimension, i + 1)));
	}

	if (is_cs) {
        _coprime_sizes.resize(dimension);
        for (size_t i=0; i<dimension; i++)
            _coprime_sizes[i] = (boost::integer::gcd(i+1, dimension) == 1); //support size and dimension are coprime

        bitset64::iterate_all_supports(dimension, [&](const bitset64& support, unsigned nbits) {
            size_t support_size_minus_one = support.count() - 1;
            if (_coprime_sizes[support_size_minus_one]) {
                if (support.smallest_representation(nbits).to_uint64() == support.to_uint64())
                    (_supports[support_size_minus_one]).push_back(support);
            } else {
                (_supports[support_size_minus_one]).push_back(support);
            }
        });
	} else {
        bitset64::iterate_all_supports(dimension, [&](const bitset64& support, unsigned) {
            (_supports[support.count() - 1]).push_back(support);
        });
	}

    if (conf_full_support) {
        search_support_size(1);
        size_t temp_count = ess_count;
        search_support_size(dimension);
        if (ess_count == temp_count) //no ess in full support found, i.e. we have to search all supports, otherwise we are done
            for (size_t i=2; i<dimension;i++)
                search_support_size(i);
    } else {
        for (size_t i=1; i<dimension+1;i++)
            search_support_size(i);
    }

}


void fracessa::search_support_size(size_t support_size) //uses real supportsize not c-array-style!
{
    _c.support_size = support_size;

    // Initialize matrices once per support_size for reuse
    DoubleMatrix A_double;
    DoubleVector b_double;
    DoubleMatrix A_SS;  // Principal submatrix for optimized path (double)
    RationalMatrix A_rational;
    RationalMatrix A_SS_rational;  // Principal submatrix for optimized path (rational)
    RationalVector b_rational;
    if (!conf_exact) {
        matrix_ops::get_kkt_rhs(support_size, b_double);
    }
    matrix_ops::get_kkt_rhs(support_size, b_rational);

    if (conf_with_log && _logger)
        _logger->info("Searching support size {}", support_size);

    for (auto support : _supports[support_size-1]) {

        _c.support = support;

        if (!conf_exact) {
            // Extract principal submatrix A_{S,S} for optimized path
            matrix_ops::principal_submatrix(game_matrix_double, dimension, _c.support, _c.support_size, A_SS);
            
            // Try optimized block inversion approach first
            bool optimized_success = find_candidate_double_optimized(A_SS);
            
            if (!optimized_success) {
                // Fall back to standard method if matrix is indefinite
                matrix_ops::get_kkt_bordering(game_matrix_double, _c.support, _c.support_size, A_double);
                if (!find_candidate_double(A_double, b_double))
                    continue;
            }
        }

        if (conf_with_log && _logger) {
            _logger->info("[rational: {}]", _c.support.to_string());
        }

        // Extract principal submatrix A_{S,S} for optimized path (rational)
        matrix_ops::principal_submatrix(game_matrix, dimension, _c.support, _c.support_size, A_SS_rational);
        
        // Try optimized block inversion approach first
        bool optimized_success = find_candidate_rational_optimized(A_SS_rational);
        
        if (!optimized_success) {
            // Fall back to standard method if matrix is singular
            matrix_ops::get_kkt_bordering(game_matrix, _c.support, _c.support_size, A_rational);
            if (!find_candidate_rational(A_rational, b_rational))
                continue;
        }

        _c.candidate_id++;

        if (conf_with_log && _logger)
            _logger->info("Found candidate! Check stability:");

        check_stability();

        if (_c.is_ess)
            ess_count++;

        if (is_cs && _coprime_sizes[support_size-1])
            _c.shift_reference = _c.candidate_id;
        else
            _c.shift_reference = 0;

        if (conf_with_candidates)
            candidates.push_back(_c);

        if (conf_with_log && _logger) {
            _logger->info("{}", candidate::header());
            _logger->info("{}", _c.to_string());
        }

        //remove all supersets for this support
        for (size_t i=support_size; i<dimension; i++) {
            _supports[i].erase(std::remove_if(
                _supports[i].begin(),
                _supports[i].end(),
                [=](const bitset64& x) {return _c.support.is_subset_of(x);}),
                _supports[i].end());
        }

        if (is_cs && _coprime_sizes[support_size-1]) {

            for (size_t i=0; i<dimension-1;i++) {

                _c.support = _c.support.rot_r(1, dimension);

                _c.candidate_id++;

                if (_c.is_ess)
                    ess_count++;

                if (conf_with_candidates) {
                    // Rotate Eigen vector: move first element to end
                    rational first = _c.vector(0);
                    for (Eigen::Index j = 0; j < _c.vector.size() - 1; j++) {
                        _c.vector(j) = _c.vector(j + 1);
                    }
                    _c.vector(_c.vector.size() - 1) = first;
                    _c.extended_support = _c.extended_support.rot_r(1, dimension);
                    candidates.push_back(_c);

                    if (conf_with_log && _logger) {
                        _logger->info("{}", candidate::header());
                        _logger->info("{}", _c.to_string());
                    }
                }

                //remove all supersets for shifted support
                for (size_t i=support_size; i<dimension; i++) {
                    _supports[i].erase(std::remove_if(
                        _supports[i].begin(),
                        _supports[i].end(),
                        [=](const bitset64& x) {return _c.support.is_subset_of(x);}),
                        _supports[i].end());
                }
            }
        }
    }
}
