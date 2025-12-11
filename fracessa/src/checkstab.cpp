#include <fracessa/fracessa.hpp>
#include <fracessa/bitset64.hpp>
#include <rational_linalg/copositivity.hpp>
#include <rational_linalg/matrix.hpp>

// Templated check_stability function
template<typename T>
void fracessa::check_stability()
{
    bitset64 bitsetm = bs64::lowest_set_bit(candidate_.support); //get lowest set bit as bitfield
    bitset64 extended_support_reduced = bs64::subtract(candidate_.extended_support, bitsetm); //ext support without m
    size_t m = bs64::find_first(candidate_.support);
    size_t extended_support_size_reduced = candidate_.extended_support_size - 1;

    if (conf_with_log_ && logger_) {
        logger_->info("Support: {}", bs64::to_bitstring(candidate_.support, dimension_));
        logger_->info("Support size: {}", bs64::count(candidate_.support));
        logger_->info("Extended support: {}", bs64::to_bitstring(candidate_.extended_support, dimension_));
        logger_->info("Extended support size: {}", candidate_.extended_support_size);
        logger_->info("Extended support reduced: {}", bs64::to_bitstring(extended_support_reduced, dimension_));
        logger_->info("index m: {}", m);
    }

    if (extended_support_size_reduced == 0)
    {
        if (conf_with_log_ && logger_)
            logger_->info("Reason: true_pure_ess");
        candidate_.stability = "T_pure_ess";
        candidate_.is_ess = true;
        return;
    }
    
    // Get Bee matrix from MatrixServer (builds it internally)
    auto& Bee = matrix_server_.get_bee_matrix<T>(extended_support_reduced, m);

    if (conf_with_log_ && logger_) {
        logger_->info("matrix bee:\n{}", rational_linalg::matrix_to_log(Bee));
    }

    if (rational_linalg::is_positive_definite_rational(Bee)) {
        if (conf_with_log_ && logger_)
            logger_->info("Reason: true_posdef_rational");
        candidate_.stability = "T_pd_rat";
        candidate_.is_ess = true;
        return;
    }

    bitset64 kay = bs64::subtract(candidate_.extended_support, candidate_.support); //extended_support without support
    size_t kay_size = bs64::count(kay);

    if (conf_with_log_ && logger_)
        logger_->info("kay: {}", bs64::to_bitstring(kay, dimension_));

    if (kay_size==0 || kay_size==1) {
        if (conf_with_log_ && logger_)
            logger_->info("Reason: false_not_posdef_and_kay_0_1");
        candidate_.stability = "F_not_pd_kay_0_1";
        candidate_.is_ess = false;
        return;
    }

    //do partial copositivity-check as in bomze_1992, p. 321/322
    bitset64 jay = extended_support_reduced;
    bitset64 jay_minus_kay = bs64::subtract(jay, kay);
    size_t r = bs64::count(jay_minus_kay);

    std::vector<bitset64> kay_vee(r+1);
    std::vector<size_t> kay_vee_size(r+1);
    std::vector<bitset64> jay_without_kay_vee(r+1);
    std::vector<rational_linalg::Matrix<T>> bee_vee(r+1);

    kay_vee[0] = jay;
    kay_vee_size[0] = extended_support_size_reduced;
    jay_without_kay_vee[0] = jay_minus_kay;
    bee_vee[0] = Bee;

    if (conf_with_log_ && logger_) {
        logger_->info("Partial Copositivity Check:");
        logger_->info("v=0:");
        logger_->info("kay_vee[0]: {}", bs64::to_bitstring(kay_vee[0], dimension_));
        logger_->info("kay_vee_size[0]: {}", kay_vee_size[0]);
        logger_->info("jay_without_kay_vee[0]: {}", bs64::to_bitstring(jay_without_kay_vee[0], dimension_));
        logger_->info("r: {}", r);
        logger_->info("bee_vee[0]:\n{}", rational_linalg::matrix_to_log(bee_vee[0]));
    }

    for (size_t v=1; v<=r; v++) {
        bitset64 iv = bs64::lowest_set_bit(jay_without_kay_vee[v-1]); //iv is lowest set bit!
        unsigned iv_pos = bs64::find_first(iv);
        jay_without_kay_vee[v] = bs64::subtract(jay_without_kay_vee[v-1], iv); //remove iv from jay\kay
        kay_vee[v] = bs64::subtract(kay_vee[v-1], iv); //build kay_vee
        kay_vee_size[v] = kay_vee_size[v-1]-1; //kay_vee_size
        // Resize bee_vee[v] only if size changed (reuse existing matrix if size matches)
        bee_vee[v] = rational_linalg::Matrix<T>(kay_vee_size[v], kay_vee_size[v]);

        // Find the position of iv within kay_vee[v-1] by counting set bits before it
        size_t pivot_pos = 0;
        for (unsigned i = 0; i < iv_pos; i++) {
            if (bs64::test(kay_vee[v-1], i)) pivot_pos++;
        }

        if (conf_with_log_ && logger_) {
            logger_->info("v={}:", v);
            logger_->info("kay_vee: {}", bs64::to_bitstring(kay_vee[v], dimension_));
            logger_->info("kay_vee_size: {}", kay_vee_size[v]);
            logger_->info("jay_without_kay_vee: {}", bs64::to_bitstring(jay_without_kay_vee[v], dimension_));
            logger_->info("iv: {}", bs64::to_bitstring(iv, dimension_));
            logger_->info("Real index (distance) to remove: {}", pivot_pos);
        }

        T pivot = bee_vee[v-1](pivot_pos, pivot_pos);
        if (pivot <= T(0)) {
            if (conf_with_log_ && logger_)
                logger_->info("Reason: false_not_partial_copositive");
            candidate_.stability = "F_not_part_copos";
            candidate_.is_ess = false;
            return;
        }

        // Equation (20) in bomze 1992, p. 321: Apply rank-1 update, develops matrix by pivot
        // Compute rank-1 product from original matrix (before scaling)
        rational_linalg::Matrix<T> rank1(kay_vee_size[v-1], kay_vee_size[v-1]);
        for (size_t i = 0; i < kay_vee_size[v-1]; ++i) {
            for (size_t j = 0; j < kay_vee_size[v-1]; ++j) {
                rank1(i, j) = -bee_vee[v-1](i, pivot_pos) * bee_vee[v-1](pivot_pos, j);
            }
        }
        // Scale matrix: pivot * B (element-wise)
        for (size_t i = 0; i < kay_vee_size[v-1]; ++i) {
            for (size_t j = 0; j < kay_vee_size[v-1]; ++j) {
                bee_vee[v-1](i, j) = bee_vee[v-1](i, j) * pivot;
            }
        }
        // Add rank-1 product: pivot*B - col*row
        for (size_t i = 0; i < kay_vee_size[v-1]; ++i) {
            for (size_t j = 0; j < kay_vee_size[v-1]; ++j) {
                bee_vee[v-1](i, j) += rank1(i, j);
            }
        }
        // Remove pivot row/column using principal_submatrix, get bitset for this dimension, and remove pivot bitset from it
        bitset64 keep_mask = 0ULL;
        bs64::set_all(keep_mask, static_cast<unsigned>(kay_vee_size[v-1]));
        bs64::reset(keep_mask, static_cast<unsigned>(pivot_pos));
        rational_linalg::principal_submatrix(bee_vee[v-1], kay_vee_size[v-1], keep_mask, kay_vee_size[v], bee_vee[v]);

        if (conf_with_log_ && logger_) {
            logger_->info("bee_vee:\n{}", rational_linalg::matrix_to_log(bee_vee[v]));
        }
    }

    //copositivity check as in hadeler_1983
    if (conf_with_log_ && logger_)
        logger_->info("Copositivity Check:");

    if (rational_linalg::isStrictlyCopositiveMemoized(bee_vee[r])) {
        if (conf_with_log_ && logger_)
            logger_->info("Reason: true_copositive");
        candidate_.stability = "T_copos";
        candidate_.is_ess = true;
        return;
    } else {
        if (conf_with_log_ && logger_)
            logger_->info("Reason: false_not_copositive");
        candidate_.stability = "F_not_copos";
        candidate_.is_ess = false;
        return;
    }
}

// Explicit template instantiations
template void fracessa::check_stability<rational>();
template void fracessa::check_stability<small_rational>();
