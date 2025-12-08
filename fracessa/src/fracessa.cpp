#include <fracessa/fracessa.hpp>
#include <fracessa/bitset64.hpp>
#include <rational_linalg/matrix_ops.hpp>

fracessa::fracessa(const RationalMatrix& matrix, bool is_cs, bool with_candidates, bool exact, bool full_support, bool with_log, int matrix_id)
{

    game_matrix_ = matrix;
    is_cs_ = is_cs;
    matrix_id_ = matrix_id;
    conf_with_candidates_ = with_candidates;
    conf_exact_ = exact;
    conf_full_support_ = full_support;
    conf_with_log_ = with_log;

    if (!conf_exact_)
        game_matrix_double_ = rational_linalg::to_double(game_matrix_);

    dimension_ = matrix.rows();
    supports_ = std::vector< std::vector<bitset64>>(dimension_);

    if (conf_with_candidates_)
        candidates_.reserve(CANDIDATE_RESERVE_MULTIPLIER * dimension_);

    if (conf_with_log_) {
        auto rotating_sink = std::make_shared<spdlog::sinks::rotating_file_sink_mt>(
            "log/fracessa.log", 20*1024*1024, 5);  // 20MB, 5 files
        logger_ = std::make_shared<spdlog::logger>("fracessa", rotating_sink);
        logger_->set_pattern("[%Y-%m-%d %H:%M:%S.%e] [%l] %v");
        logger_->set_level(spdlog::level::info);
        
        // Write empty line and 3 lines of asterisks and hash symbols as first lines in log
        logger_->info("");
        logger_->info("*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*");
        logger_->info("#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#");
        logger_->info("*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*");
        
        // Write matrix_id as first line in log
        if (matrix_id >= 0) {
            logger_->info("matrix_id={}", matrix_id);
        }
        
        logger_->info("n={}", dimension_);
        logger_->info("game matrix:\n{}", rational_linalg::to_string(game_matrix_));
    }

	for (size_t i = 0; i < dimension_; i++) {
		supports_[i].reserve(static_cast<uint64_t>(boost::math::binomial_coefficient<double>(dimension_, i + 1)));
	}

	if (is_cs_) {
        coprime_sizes_.resize(dimension_);
        for (size_t i=0; i<dimension_; i++)
            coprime_sizes_[i] = (boost::integer::gcd(i+1, dimension_) == 1); //support size and dimension are coprime

        bitset64::iterate_all_supports(dimension_, [&](const bitset64& support) {
            size_t this_size = support.count() - 1;
            if (coprime_sizes_[this_size]) {
                if (support.smallest_representation(dimension_) == support)
                    (supports_[this_size]).push_back(support);
            } else {
                (supports_[this_size]).push_back(support);
            }
        });
	} else {
        bitset64::iterate_all_supports(dimension_, [&](const bitset64& support) {
            (supports_[support.count() - 1]).push_back(support);
        });
	}

    if (conf_full_support_) {
        search_support_size(1);
        size_t temp_count = ess_count_;
        search_support_size(dimension_);
        if (ess_count_ == temp_count) //no ess in full support found, i.e. we have to search all supports, otherwise we are done
            for (size_t i=2; i<dimension_;i++)
                search_support_size(i);
    } else {
        for (size_t i=1; i<dimension_+1;i++)
            search_support_size(i);
    }

}


void fracessa::search_support_size(size_t support_size) //uses real supportsize not c-array-style!
{
    c_.support_size = support_size;

    // Initialize matrices once per support_size for reuse
    DoubleMatrix A_double;
    DoubleVector b_double;
    DoubleMatrix A_SS;  // Principal submatrix for optimized path (double)
    RationalMatrix A_rational;
    RationalMatrix A_SS_rational;  // Principal submatrix for optimized path (rational)
    RationalVector b_rational;
    if (!conf_exact_) {
        matrix_ops::get_kkt_rhs(support_size, b_double);
    }
    rational_linalg::get_kkt_rhs(support_size, b_rational);

    if (conf_with_log_ && logger_)
        logger_->info("Searching support size {}:", support_size);

    for (auto support : supports_[support_size-1]) {

        c_.support = support;

        if (!conf_exact_) {
            // Extract principal submatrix A_{S,S} for optimized path
            // matrix_ops::principal_submatrix(game_matrix_double_, dimension_, c_.support, c_.support_size, A_SS);
            
            // Try optimized block inversion approach first
            // bool optimized_success = find_candidate_double_optimized(A_SS);
            
            // if (!optimized_success) {
                // Fall back to standard method if matrix is indefinite
                matrix_ops::get_kkt_bordering(game_matrix_double_, c_.support, c_.support_size, A_double);
                if (!find_candidate_double(A_double, b_double))
                    continue;
            // }
        }

        // Extract principal submatrix A_{S,S} for optimized path (rational)
        // rational_linalg::principal_submatrix(game_matrix_, dimension_, c_.support, c_.support_size, A_SS_rational);
        
        // Try optimized block inversion approach first
        // bool optimized_success = find_candidate_rational_optimized(A_SS_rational);
        
        // if (!optimized_success) {
            // Fall back to standard method if matrix is singular
            rational_linalg::get_kkt_bordering(game_matrix_, c_.support, c_.support_size, A_rational);
            if (!find_candidate_rational(A_rational, b_rational))
                continue;
        // }

        c_.candidate_id++;

        if (conf_with_log_ && logger_)
            logger_->info("Found candidate! Check stability:");

        check_stability();

        if (c_.is_ess)
            ess_count_++;

        if (is_cs_ && coprime_sizes_[support_size-1])
            c_.shift_reference = c_.candidate_id;
        else
            c_.shift_reference = 0;

        if (conf_with_candidates_)
            candidates_.push_back(c_);

        if (conf_with_log_ && logger_) {
            logger_->info("{}", candidate::header());
            logger_->info("{}", c_.to_string());
        }

        //remove all supersets for this support
        for (size_t i=support_size; i<dimension_; i++) {
            supports_[i].erase(std::remove_if(
                supports_[i].begin(),
                supports_[i].end(),
                [=](const bitset64& x) {return c_.support.is_subset_of(x);}),
                supports_[i].end());
        }

        if (is_cs_ && coprime_sizes_[support_size-1]) {

            for (size_t i=0; i<dimension_-1;i++) {

                c_.support.rot_one_right(dimension_);

                c_.candidate_id++;

                if (c_.is_ess)
                    ess_count_++;

                if (conf_with_candidates_) {
                    // Rotate vector: move first element to end
                    rational first = c_.vector(0, 0);
                    size_t vec_size = c_.vector.rows();
                    for (size_t j = 0; j < vec_size - 1; j++) {
                        c_.vector(j, 0) = c_.vector(j + 1, 0);
                    }
                    c_.vector(vec_size - 1, 0) = first;
                    c_.extended_support.rot_one_right(dimension_);
                    candidates_.push_back(c_);

                    if (conf_with_log_ && logger_) {
                        logger_->info("{}", candidate::header());
                        logger_->info("{}", c_.to_string());
                    }
                }

                //remove all supersets for shifted support
                for (size_t i=support_size; i<dimension_; i++) {
                    supports_[i].erase(std::remove_if(
                        supports_[i].begin(),
                        supports_[i].end(),
                        [=](const bitset64& x) {return c_.support.is_subset_of(x);}),
                        supports_[i].end());
                }
            }
        }
    }
}