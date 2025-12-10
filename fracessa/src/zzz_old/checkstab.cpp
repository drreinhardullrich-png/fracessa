#include <fracessa/fracessa.hpp>
#include <fracessa/bitset64.hpp>
#include <rational_linalg/copositivity.hpp>
#include <rational_linalg/matrix.hpp>
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

    // //dont need to care too much about precision/tolerance here, since if true we will do rational check, otherwise we find it later anyways....
    // if (!conf_exact_) {       
    //     DoubleMatrix bee_double = matrix_ops::to_double(bee);
        
    //     // // Use is_positive_definite_double to check positive definiteness
    //     // if (matrix_ops::is_positive_definite_double(bee_double)) {
    //     //     if (conf_with_log_ && logger_)
    //     //         logger_->info("Reason: true_posdef_double");
    //     //     c_.stability = "T_pd_double";
    //     //     c_.is_ess = true;
    //     //     return;
    //     // }
    //     // Use Eigen's LDLT decomposition with tolerance check on D diagonal
    //     Eigen::LDLT<DoubleMatrix> ldlt(bee_double);
    //     if (ldlt.info() == Eigen::Success) {
    //         // Check that all diagonal elements of D are strictly positive
    //         auto D_diag = ldlt.vectorD();
    //         double tolerance = std::max(1e-10, 
    //             std::numeric_limits<double>::epsilon() * bee_double.diagonal().cwiseAbs().maxCoeff());
            
    //         if ((D_diag.array() > tolerance).all()) {
    //             // Strictly positive definite
    //             if (conf_with_log_ && logger_)
    //                 logger_->info("Reason: true_posdef_double");
    //             c_.stability = "T_pd_double";
    //             c_.is_ess = true;
    //             return;
    //         }
    //     }
    // }

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

/*
    for (size_t v = 1; v <= r; v++) {
        bitset64 iv = jay_without_kay_vee[v-1].lowest_set_bit();
        unsigned iv_pos = iv.find_first();
        kay_vee[v] = kay_vee[v-1].subtract(iv);
        kay_vee_size[v] = kay_vee_size[v-1] - 1;
        
        size_t pivot_pos = 0;
        for (unsigned i = 0; i < iv_pos; i++) {
            if (kay_vee[v-1].test(i)) pivot_pos++;
        }
        
        rational pivot = bee_vee[v-1](pivot_pos, pivot_pos);
        if (pivot <= rational(0)) {
            if (conf_with_log && _logger)
                _logger->info("Reason: false_not_partial_copositive");
            _c.stability = "F_not_part_copos";
            _c.is_ess = false;
            return;
        }
        
        // Equation (20): Apply rank-1 update
        RationalMatrix temp = bee_vee[v-1];
        temp.noalias() -= (temp.col(pivot_pos) * temp.row(pivot_pos)) / pivot;
        
        // Remove pivot row/column using principal_submatrix
        bitset64 keep_mask;
        keep_mask.set_all(temp.rows());
        keep_mask.reset(pivot_pos);
        matrix_ops::principal_submatrix(temp, temp.rows(), keep_mask, kay_vee_size[v], bee_vee[v]);
    }*/




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



        // for (size_t i = 0;i< kay_vee_size[v];i++) { //build matrix bee_vee[v]
        //     for (size_t j = 0; j< kay_vee_size[v]; j++) {

        //         size_t row_index_kv_minus_one = (i>=index_position_to_remove) ? i+1 : i;
        //         size_t col_index_kv_minus_one = (j>=index_position_to_remove) ? j+1 : j;

        //         bee_vee[v](i,j) = bee_vee[v-1](index_position_to_remove,index_position_to_remove) * bee_vee[v-1](row_index_kv_minus_one,col_index_kv_minus_one) -
        //             bee_vee[v-1](index_position_to_remove,row_index_kv_minus_one) * bee_vee[v-1](index_position_to_remove,col_index_kv_minus_one);
        //     }
        // }

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

