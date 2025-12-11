#include <fracessa/fracessa.hpp>
#include <fracessa/bitset64.hpp>
#include <rational_linalg/matrix.hpp>
#include <exception>

fracessa::fracessa(const rational_linalg::Matrix<small_rational>& matrix, bool is_cs, bool with_candidates, bool exact, bool full_support, bool with_log, int matrix_id)
    : matrix_server_(matrix)
    , dimension_(matrix.rows())
    , is_cs_(is_cs)
    , matrix_id_(matrix_id)
    , conf_with_candidates_(with_candidates)
    , conf_exact_(exact)
    , conf_full_support_(full_support)
    , conf_with_log_(with_log)
    , candidate_()
    , supports_(dimension_, is_cs_)
    , logger_()
{
    if (conf_with_candidates_)
        candidates_.reserve(250 * dimension_);

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
        logger_->info("game matrix:\n{}", rational_linalg::matrix_to_log(matrix_server_.get_game_matrix<small_rational>()));
    }
        
    // Initialize supports
    supports_.initialize();

    if (conf_full_support_) {      
        // Search full support (dimension_)
        // BATCH: supports_to_remove_.clear();
        for (const auto& support : supports_.get_supports(dimension_)) {
            search_one_support(support, dimension_);
        }
        // BATCH: Batch remove all collected supersets
        // supports_.remove_supersets_batch(supports_to_remove_, dimension_);
        // supports_to_remove_.clear();
        if (ess_count_ > 0) 
            return;                
    }
    // Search all support sizes
    for (size_t i = 1; i <= (conf_full_support_ ? dimension_-1: dimension_) ; i++) {
        if (conf_with_log_)
            logger_->info("Searching support size {}:", i);

        // BATCH: supports_to_remove_.clear();
        const bool is_cs_and_coprime = is_cs_ && (boost::integer::gcd(i, dimension_) == 1);
        for (const auto& support : supports_.get_supports(i)) {
            search_one_support(support, i, is_cs_and_coprime);
        }
        // BATCH: Batch remove all collected supersets for this support size
        // supports_.remove_supersets_batch(supports_to_remove_, i);
        // supports_to_remove_.clear();
    }
}


void fracessa::search_one_support(const bitset64& support, size_t support_size, bool is_cs_and_coprime)
{
    if (!conf_exact_) 
            if (!find_candidate<double>(support, support_size))
                return;

    // Try small_rational first (fast path)
        try {
        if (!find_candidate<small_rational>(support, support_size))
                return;
        } catch (const std::exception& e) {
            // Check if it's an overflow error
            std::string error_msg = e.what();
            if (error_msg.find("overflow") != std::string::npos) {
            // Fall back to rational precision
            if (!find_candidate<rational>(support, support_size))
                    return;
            } else {              
                throw;  // Re-throw if it's not an overflow error
            }
    }

    candidate_.support_size = support_size;
    candidate_.support = support;
    candidate_.candidate_id++;
    
    if (conf_with_log_)
        logger_->info("Found candidate! Check stability:");

    // Try small_rational first (fast path)
        try {
        check_stability<small_rational>();
        } catch (const std::exception& e) {
            // Check if it's an overflow error
            std::string error_msg = e.what();
            if (error_msg.find("overflow") != std::string::npos) {
            // Fall back to rational precision
            check_stability<rational>();
            } else {
                // Re-throw if it's not an overflow error
                throw;
            }
    }

    if (candidate_.is_ess)
        ess_count_++;

    if (is_cs_and_coprime)
        candidate_.shift_reference = candidate_.candidate_id;
    else
        candidate_.shift_reference = 0;

    if (conf_with_candidates_)
        candidates_.push_back(candidate_);

    if (conf_with_log_) {
        logger_->info("{}", candidate::header());
        logger_->info("{}", candidate_.to_string());
    }
    // BATCH: Collect support for batch removal
    // supports_to_remove_.push_back(candidate_.support);
    // Remove all supersets for this support using Supports class
    supports_.remove_supersets(candidate_.support, support_size);

    if (is_cs_and_coprime) { //do the same for the n-1 more candidates we get for free because of the coprime property!

        for (size_t i = 0; i < dimension_ - 1; i++) {

            bs64::rot_one_right(candidate_.support, dimension_);

            candidate_.candidate_id++;

            if (candidate_.is_ess)
                ess_count_++;

            if (conf_with_candidates_) {
                // Rotate vector: move first element to end
                rational first = candidate_.vector(0, 0);
                size_t vec_size = candidate_.vector.rows();
                for (size_t j = 0; j < vec_size - 1; j++) {
                    candidate_.vector(j, 0) = candidate_.vector(j + 1, 0);
                }
                candidate_.vector(vec_size - 1, 0) = first;
                bs64::rot_one_right(candidate_.extended_support, dimension_);
                candidates_.push_back(candidate_);

                if (conf_with_log_) {
                    logger_->info("{}", candidate::header());
                    logger_->info("{}", candidate_.to_string());
                }
            }
            // BATCH: Collect shifted support for batch removal
            // supports_to_remove_.push_back(candidate_.support);
            // Remove all supersets for shifted support using Supports class
            supports_.remove_supersets(candidate_.support, support_size);
        }
    }
}

