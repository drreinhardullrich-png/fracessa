#include <fracessa/fracessa.hpp>
#include <fracessa/bitset64.hpp>
#include <rational_linalg/matrix_ops.hpp>

fracessa::fracessa(const RationalMatrix& matrix, bool is_cs, bool with_candidates, bool exact, bool full_support, bool with_log, int matrix_id)
    : game_matrix_(matrix)
    , dimension_(matrix.rows())
    , is_cs_(is_cs)
    , matrix_id_(matrix_id)
    , conf_with_candidates_(with_candidates)
    , conf_exact_(exact)
    , conf_full_support_(full_support)
    , conf_with_log_(with_log)
    , c_()
    , candidates_supports_()
    , supports_(dimension_, is_cs_)
    , logger_()
{

    if (!conf_exact_)
        game_matrix_double_ = rational_linalg::to_double(game_matrix_);

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
        
    // Initialize supports
    supports_.initialize();

    if (conf_full_support_) {      
        // Search full support (dimension_)
        for (const auto& support : supports_.get_supports(dimension_)) {
            search_one_support(support, dimension_);
        }
        if (ess_count_ > 0) 
            return;                
    }
    // Search all support sizes
    for (size_t i = 1; i <= (conf_full_support_ ? dimension_-1: dimension_) ; i++) {
        if (conf_with_log_ && logger_)
            logger_->info("Searching support size {}:", i);

        bool is_cs_and_coprime = false;
        if (is_cs_) {
            is_cs_and_coprime = (boost::integer::gcd(i, dimension_) == 1) && is_cs_;
        }
        for (const auto& support : supports_.get_supports(i)) {
            search_one_support(support, i, is_cs_and_coprime);
        }
    }
}


void fracessa::search_one_support(const bitset64& support, size_t support_size, bool is_cs_and_coprime)
{
    c_.support_size = support_size;
    c_.support = support;

    // Initialize matrices once per support for reuse
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

    if (!conf_exact_) {
        // Extract principal submatrix A_{S,S} for optimized path
        // matrix_ops::principal_submatrix(game_matrix_double_, dimension_, c_.support, c_.support_size, A_SS);
        
        // Try optimized block inversion approach first
        // bool optimized_success = find_candidate_double_optimized(A_SS);
        
        // if (!optimized_success) {
            // Fall back to standard method if matrix is indefinite
            matrix_ops::get_kkt_bordering(game_matrix_double_, c_.support, c_.support_size, A_double);
            if (!find_candidate_double(A_double, b_double))
                return;
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
            return;
    // }

    c_.candidate_id++;

    if (conf_with_log_ && logger_)
        logger_->info("Found candidate! Check stability:");

    check_stability();

    if (c_.is_ess)
        ess_count_++;

    if (is_cs_and_coprime)
        c_.shift_reference = c_.candidate_id;
    else
        c_.shift_reference = 0;

    if (conf_with_candidates_)
        candidates_.push_back(c_);

    if (conf_with_log_ && logger_) {
        logger_->info("{}", candidate::header());
        logger_->info("{}", c_.to_string());
    }
    // Remove all supersets for this support using Supports class
    supports_.remove_supersets(c_.support, support_size);

    if (is_cs_and_coprime) { //do the same for the n-1 more candidates we get for free because of the coprime property!

        for (size_t i = 0; i < dimension_ - 1; i++) {

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
            // Remove all supersets for shifted support using Supports class
            supports_.remove_supersets(c_.support, support_size);
        }
    }
}

