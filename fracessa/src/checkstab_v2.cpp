#include <fracessa/fracessa.hpp>
#include <fracessa/bitset64.hpp>
#include <rational_linalg/copositivity.hpp>
#include <rational_linalg/matrix_ops.hpp>
#include <Eigen/Cholesky>

void fracessa::check_stability()
{
    bitset64 bitsetm = c_.support.lowest_set_bit(); //get lowest set bit as bitfield
    bitset64 extended_support_reduced = c_.extended_support.subtract(bitsetm); //ext support without m
    size_t m = c_.support.find_first();
    size_t extended_support_size_reduced = c_.extended_support_size - 1;

    if (conf_with_log_ && logger_) {
        logger_->info("Support: {}", c_.support.to_bitstring(dimension_));
        logger_->info("Support size: {}", c_.support.count());
        logger_->info("Extended support: {}", c_.extended_support.to_bitstring(dimension_));
        logger_->info("Extended support size: {}", c_.extended_support_size);
        logger_->info("Extended support reduced: {}", extended_support_reduced.to_bitstring(dimension_));
        logger_->info("index m: {}", m);
    }

    if (extended_support_size_reduced == 0)
    {

        if (conf_with_log_ && logger_)
            logger_->info("Reason: true_pure_ess");
        c_.stability = "T_pure_ess";
        c_.is_ess = true;
        return;
    }

    RationalMatrix bee(extended_support_size_reduced, extended_support_size_reduced);

    size_t row = 0;
    size_t column = 0;
    for (size_t i = 0;i<dimension_;i++)
        if (extended_support_reduced.test(i)) {
            column = 0;
            for (size_t j = 0; j < i+1; j++)
                if (extended_support_reduced.test(j)) {
                    bee(row,column) = bee(column,row) = game_matrix_(m, j) + game_matrix_(j, m) + game_matrix_(i, m) + game_matrix_(m, i) -
                        game_matrix_(i, j) - game_matrix_(j, i) - 2 * game_matrix_(m, m);
                    column += 1;
                }
            row += 1;
        }

    if (conf_with_log_ && logger_) {
        logger_->info("matrix bee:\n{}", rational_linalg::to_string(bee));
    }

    if (rational_linalg::is_positive_definite_rational(bee)) {

        if (conf_with_log_ && logger_)
            logger_->info("Reason: true_posdef_rational");
        c_.stability = "T_pd_rat";
        c_.is_ess = true;
        return;
    }

    bitset64 kay = c_.extended_support.subtract(c_.support); //extended_support without support
    size_t kay_size = kay.count();

    if (conf_with_log_ && logger_)
        logger_->info("kay: {}", kay.to_bitstring(dimension_));

    if (kay_size==0 || kay_size==1) {
        if (conf_with_log_ && logger_)
            logger_->info("Reason: false_not_posdef_and_kay_0_1");
        c_.stability = "F_not_pd_kay_0_1";
        c_.is_ess = false;
        return;
    }

    //do partial copositivity-check as in bomze_1992, p. 321/322
    bitset64 jay = extended_support_reduced;
    bitset64 jay_minus_kay = jay.subtract(kay);
    size_t r = jay_minus_kay.count();

    std::vector<bitset64> kay_vee(r+1);
    std::vector<size_t> kay_vee_size(r+1);
    std::vector<bitset64> jay_without_kay_vee(r+1);
    std::vector< RationalMatrix> bee_vee(r+1);


    kay_vee[0] = jay;
    kay_vee_size[0] = extended_support_size_reduced;
    jay_without_kay_vee[0] = jay_minus_kay;
    bee_vee[0] = bee;

    if (conf_with_log_ && logger_) {
        logger_->info("Partial Copositivity Check:");
        logger_->info("v=0:");
        logger_->info("kay_vee[0]: {}", kay_vee[0].to_bitstring(dimension_));
        logger_->info("kay_vee_size[0]: {}", kay_vee_size[0]);
        logger_->info("jay_without_kay_vee[0]: {}", jay_without_kay_vee[0].to_bitstring(dimension_));
        logger_->info("r: {}", r);
        logger_->info("bee_vee[0]:\n{}", rational_linalg::to_string(bee_vee[0]));
    }

    for (size_t v=1; v<=r; v++) {

        bitset64 iv = jay_without_kay_vee[v-1].lowest_set_bit(); //iv is lowest set bit!
        unsigned iv_pos = iv.find_first();
        jay_without_kay_vee[v] = jay_without_kay_vee[v-1].subtract(iv); //remove iv from jay\kay
        kay_vee[v] = kay_vee[v-1].subtract(iv); //build kay_vee
        kay_vee_size[v] = kay_vee_size[v-1]-1; //kay_vee_size
        // Resize bee_vee[v] only if size changed (reuse existing matrix if size matches)
        bee_vee[v] = RationalMatrix(kay_vee_size[v], kay_vee_size[v]);

        // Find the position of iv within kay_vee[v-1] by counting set bits before it
        size_t pivot_pos = 0;
        for (unsigned i = 0; i < iv_pos; i++) {
            if (kay_vee[v-1].test(i)) pivot_pos++;
        }

        if (conf_with_log_ && logger_) {
            logger_->info("v={}:", v);
            logger_->info("kay_vee: {}", kay_vee[v].to_bitstring(dimension_));
            logger_->info("kay_vee_size: {}", kay_vee_size[v]);
            logger_->info("jay_without_kay_vee: {}", jay_without_kay_vee[v].to_bitstring(dimension_));
            logger_->info("iv: {}", iv.to_bitstring(dimension_));
            logger_->info("Real index (distance) to remove: {}", pivot_pos);
        }

        rational pivot = bee_vee[v-1](pivot_pos, pivot_pos);
        if (pivot <= rational(0)) {
            if (conf_with_log_ && logger_)
                logger_->info("Reason: false_not_partial_copositive");
            c_.stability = "F_not_part_copos";
            c_.is_ess = false;
            return;
        }

        // Equation (20) in bomze 1992, p. 321: Apply rank-1 update, develops matrix by pivot
        // Compute rank-1 product from original matrix (before scaling)
        RationalMatrix rank1(kay_vee_size[v-1], kay_vee_size[v-1]);
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
        bitset64 keep_mask;
        keep_mask.set_all(static_cast<unsigned>(kay_vee_size[v-1]));
        keep_mask.reset(static_cast<unsigned>(pivot_pos));
        rational_linalg::principal_submatrix(bee_vee[v-1], kay_vee_size[v-1], keep_mask, kay_vee_size[v], bee_vee[v]);

        if (conf_with_log_ && logger_) {
            logger_->info("bee_vee:\n{}", rational_linalg::to_string(bee_vee[v]));
        }
    }

    //copositivity check as in hadeler_1983
    if (conf_with_log_ && logger_)
        logger_->info("Copositivity Check:");

    if (rational_linalg::isStrictlyCopositiveMemoized(bee_vee[r])) {
        if (conf_with_log_ && logger_)
            logger_->info("Reason: true_copositive");
        c_.stability = "T_copos";
        c_.is_ess = true;
        return;
    } else {
        if (conf_with_log_ && logger_)
            logger_->info("Reason: false_not_copositive");
        c_.stability = "F_not_copos";
        c_.is_ess = false;
        return;
    }
}

