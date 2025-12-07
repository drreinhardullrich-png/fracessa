#include <fracessa/fracessa.hpp>
#include <fracessa/bitset64.hpp>
#include <eigen_extensions/copositivity.hpp>
#include <Eigen/Cholesky>

void fracessa::check_stability()
{
    bitset64 bitsetm = _c.support.lowest_set_bit(); //get lowest set bit as bitfield
    bitset64 extended_support_reduced = _c.extended_support.subtract(bitsetm); //ext support without m
    size_t m = _c.support.find_first();
    size_t extended_support_size_reduced = _c.extended_support_size - 1;

    if (conf_with_log && _logger) {
        _logger->info("Support: {}", _c.support.to_bitstring(dimension));
        _logger->info("Support size: {}", _c.support.count());
        _logger->info("Extended support: {}", _c.extended_support.to_bitstring(dimension));
        _logger->info("Extended support size: {}", _c.extended_support_size);
        _logger->info("Extended support reduced: {}", extended_support_reduced.to_bitstring(dimension));
        _logger->info("index m: {}", m);
    }

    if (extended_support_size_reduced == 0)
    {

        if (conf_with_log && _logger)
            _logger->info("Reason: true_pure_ess");
        _c.stability = "T_pure_ess";
        _c.is_ess = true;
        return;
    }

    RationalMatrix bee = RationalMatrix::Zero(extended_support_size_reduced,extended_support_size_reduced);

    size_t row = 0;
    size_t column = 0;
    for (size_t i = 0;i<dimension;i++)
        if (extended_support_reduced.test(i)) {
            column = 0;
            for (size_t j = 0; j < i+1; j++)
                if (extended_support_reduced.test(j)) {
                    bee(row,column) = bee(column,row) = game_matrix(m, j) + game_matrix(j, m) + game_matrix(i, m) + game_matrix(m, i) -
                        game_matrix(i, j) - game_matrix(j, i) - 2 * game_matrix(m, m);
                    column += 1;
                }
            row += 1;
        }

    if (conf_with_log && _logger) {
        _logger->info("matrix bee:\n{}", matrix_ops::to_string(bee));
    }

    // //dont need to care too much about precision/tolerance here, since if true we will do rational check, otherwise we find it later anyways....
    // if (!conf_exact) {       
    //     DoubleMatrix bee_double = matrix_ops::to_double(bee);
        
    //     // // Use is_positive_definite_double to check positive definiteness
    //     // if (matrix_ops::is_positive_definite_double(bee_double)) {
    //     //     if (conf_with_log && _logger)
    //     //         _logger->info("Reason: true_posdef_double");
    //     //     _c.stability = "T_pd_double";
    //     //     _c.is_ess = true;
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
    //             if (conf_with_log && _logger)
    //                 _logger->info("Reason: true_posdef_double");
    //             _c.stability = "T_pd_double";
    //             _c.is_ess = true;
    //             return;
    //         }
    //     }
    // }

    if (matrix_ops::is_positive_definite_rational(bee)) {

        if (conf_with_log && _logger)
            _logger->info("Reason: true_posdef_rational");
        _c.stability = "T_pd_rat";
        _c.is_ess = true;
        return;
    }

    bitset64 kay = _c.extended_support.subtract(_c.support); //extended_support without support
    size_t kay_size = kay.count();

    if (conf_with_log && _logger)
        _logger->info("kay: {}", kay.to_bitstring(dimension));

    if (kay_size==0 || kay_size==1) {
        if (conf_with_log && _logger)
            _logger->info("Reason: false_not_posdef_and_kay_0_1");
        _c.stability = "F_not_pd_kay_0_1";
        _c.is_ess = false;
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

    if (conf_with_log && _logger) {
        _logger->info("Partial Copositivity Check:");
        _logger->info("v=0:");
        _logger->info("kay_vee[0]: {}", kay_vee[0].to_bitstring(dimension));
        _logger->info("kay_vee_size[0]: {}", kay_vee_size[0]);
        _logger->info("jay_without_kay_vee[0]: {}", jay_without_kay_vee[0].to_bitstring(dimension));
        _logger->info("r: {}", r);
        _logger->info("bee_vee[0]:\n{}", matrix_ops::to_string(bee_vee[0]));
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
        bee_vee[v] = RationalMatrix::Zero(kay_vee_size[v],kay_vee_size[v]);

        // Find the position of iv within kay_vee[v-1] by counting set bits before it
        size_t pivot_pos = 0;
        for (unsigned i = 0; i < iv_pos; i++) {
            if (kay_vee[v-1].test(i)) pivot_pos++;
        }

        if (conf_with_log && _logger) {
            _logger->info("v={}:", v);
            _logger->info("kay_vee: {}", kay_vee[v].to_bitstring(dimension));
            _logger->info("kay_vee_size: {}", kay_vee_size[v]);
            _logger->info("jay_without_kay_vee: {}", jay_without_kay_vee[v].to_bitstring(dimension));
            _logger->info("iv: {}", iv.to_bitstring(dimension));
            _logger->info("Real index (distance) to remove: {}", pivot_pos);
        }

        rational pivot = bee_vee[v-1](pivot_pos, pivot_pos);
        if (pivot <= rational(0)) {
            if (conf_with_log && _logger)
                _logger->info("Reason: false_not_partial_copositive");
            _c.stability = "F_not_part_copos";
            _c.is_ess = false;
            return;
        }

        // Equation (20) in bomze 1992, p. 321: Apply rank-1 update, develops matrix by pivot
        // Compute rank-1 product from original matrix (before scaling) - uses Eigen's optimized outer product
        RationalMatrix rank1 = -bee_vee[v-1].col(pivot_pos) * bee_vee[v-1].row(pivot_pos);
        // Scale matrix: pivot * B (element-wise to avoid Boost operator interception)
        for (Eigen::Index i = 0; i < bee_vee[v-1].rows(); ++i) {
            for (Eigen::Index j = 0; j < bee_vee[v-1].cols(); ++j) {
                bee_vee[v-1](i, j) = bee_vee[v-1](i, j) * pivot;
            }
        }
        // Add rank-1 product: pivot*B - col*row (uses Eigen's optimized addition)
        bee_vee[v-1].noalias() += rank1;
        // Remove pivot row/column using principal_submatrix, get bitset for this dimension, and remove pivot bitset from it
        bitset64 keep_mask;
        keep_mask.set_all(kay_vee_size[v-1]);
        keep_mask.reset(pivot_pos);
        matrix_ops::principal_submatrix(bee_vee[v-1], kay_vee_size[v-1], keep_mask, kay_vee_size[v], bee_vee[v]);



        // for (size_t i = 0;i< kay_vee_size[v];i++) { //build matrix bee_vee[v]
        //     for (size_t j = 0; j< kay_vee_size[v]; j++) {

        //         size_t row_index_kv_minus_one = (i>=index_position_to_remove) ? i+1 : i;
        //         size_t col_index_kv_minus_one = (j>=index_position_to_remove) ? j+1 : j;

        //         bee_vee[v](i,j) = bee_vee[v-1](index_position_to_remove,index_position_to_remove) * bee_vee[v-1](row_index_kv_minus_one,col_index_kv_minus_one) -
        //             bee_vee[v-1](index_position_to_remove,row_index_kv_minus_one) * bee_vee[v-1](index_position_to_remove,col_index_kv_minus_one);
        //     }
        // }

        if (conf_with_log && _logger) {
            _logger->info("bee_vee:\n{}", matrix_ops::to_string(bee_vee[v]));
        }
    }

    //copositivity check as in hadeler_1983
    if (conf_with_log && _logger)
        _logger->info("Copositivity Check:");

    if (isStrictlyCopositiveMemoized(bee_vee[r])) {
        if (conf_with_log && _logger)
            _logger->info("Reason: true_copositive");
        _c.stability = "T_copos";
        _c.is_ess = true;
        return;
    } else {
        if (conf_with_log && _logger)
            _logger->info("Reason: false_not_copositive");
        _c.stability = "F_not_copos";
        _c.is_ess = false;
        return;
    }
}

