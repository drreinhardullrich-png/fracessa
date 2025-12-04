#include <fracessa/fracessa.hpp>

fracessa::fracessa(const matrix<rational>& matrix, bool with_candidates, bool exact, bool full_support, bool with_log)
{

    game_matrix = matrix;
    conf_with_candidates = with_candidates;
    conf_exact = exact;
    conf_full_support = full_support;
    conf_with_log = with_log;

    if (!conf_exact)
        game_matrix_double = game_matrix.to_double();

    dimension = matrix.rows();
    _supports = std::vector< std::vector<uint64_t>>(dimension);

    if (conf_with_candidates)
        candidates.reserve(CANDIDATE_RESERVE_MULTIPLIER * dimension);

    if (conf_with_log) {
        auto rotating_sink = std::make_shared<spdlog::sinks::rotating_file_sink_mt>(
            "log/fracessa.log", 20*1024*1024, 5);  // 20MB, 5 files
        _logger = std::make_shared<spdlog::logger>("fracessa", rotating_sink);
        _logger->set_pattern("[%Y-%m-%d %H:%M:%S.%e] [%l] %v");
        _logger->set_level(spdlog::level::info);
        
        _logger->info("n={}", dimension);
        _logger->info("game matrix:\n{}", game_matrix.to_string());
    }

	for (size_t i = 0; i < dimension; i++) {
		_supports[i].reserve(static_cast<uint64_t>(boost::math::binomial_coefficient<double>(dimension, i + 1)));
	}

	if (game_matrix.is_cs()) {
        _coprime_sizes.resize(dimension);
        for (size_t i=0; i<dimension; i++)
            _coprime_sizes[i] = (boost::integer::gcd(i+1, dimension) == 1); //support size and dimension are coprime

        for (uint64_t i=1ull; i<(1ull<<dimension); i++) {
            size_t support_size_minus_one = support_size_from_int(i)-1;
            if (_coprime_sizes[support_size_minus_one]) {
                if (smallest_representation(i, dimension) == i)
                    (_supports[support_size_minus_one]).push_back(i);
            } else {
                (_supports[support_size_minus_one]).push_back(i);
            }
        }
	} else {
        for (uint64_t i=1ull; i<(1ull<<dimension); i++) {
            (_supports[support_size_from_int(i)-1]).push_back(i);
        }
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
    matrix<double> le_matrix_double = matrix<double>(support_size + 1, support_size + 2);
    matrix<rational> le_matrix_rational = matrix<rational>(support_size + 1, support_size + 2);

    if (conf_with_log && _logger)
        _logger->info("Searching support size {}", support_size);

    for (auto support : _supports[support_size-1]) {

        _c.support = support;

        std::string support_log_msg;
        if (!conf_exact) {
            if (conf_with_log && _logger) {
                support_log_msg = "[" + std::to_string(_c.support) + "]";
            }
            if (!find_candidate_double(le_matrix_double))
                continue;
        }

        if (conf_with_log && _logger) {
            // Combine support and rational messages like the old format
            if (!conf_exact) {
                _logger->info("{}[rational: {}]", support_log_msg, _c.support);
            } else {
                _logger->info("[rational: {}]", _c.support);
            }
        }

        if (!find_candidate_rational(le_matrix_rational))
            continue;

        _c.candidate_id++;

        if (conf_with_log && _logger)
            _logger->info("Found candidate! Check stability:");

        check_stability();

        if (_c.is_ess)
            ess_count++;

        if (game_matrix.is_cs() && _coprime_sizes[support_size-1])
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
                [=](const uint64_t& x) {return ((_c.support & x) == _c.support);}),
                _supports[i].end());
        }

        if (game_matrix.is_cs() && _coprime_sizes[support_size-1]) {

            for (size_t i=0; i<dimension-1;i++) {

                _c.support = shift_right(_c.support,dimension);

                _c.candidate_id++;

                if (_c.is_ess)
                    ess_count++;

                if (conf_with_candidates) {

                    std::rotate(_c.vector.begin(), _c.vector.begin() + 1, _c.vector.end());
                    _c.extended_support = shift_right(_c.extended_support,dimension);
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
                        [=](const uint64_t& x) {return ((_c.support & x) == _c.support);}),
                        _supports[i].end());
                }
            }
        }
    }
}



bool fracessa::find_candidate_double(matrix<double> &le_matrix)
{
    int n = _c.support_size + 1;

    std::vector<double> result_le = std::vector<double>(n, 0.);
    std::vector<double> result_vector = std::vector<double>(dimension);

    game_matrix_double.get_le_matrix(_c.support, _c.support_size, le_matrix);

    ////////////////////////////////////////////////////////////////////////// gauss with partial pivoting
    for (int i=0; i<n; i++) {
        // Search for maximum in this column
        double max_element = std::abs(le_matrix(i,i));
        int max_row = i;
        for (int k=i+1; k<n; k++) {
            if (std::abs(le_matrix(k,i)) > max_element) {
                max_element = std::abs(le_matrix(k,i));
                max_row = k;
            }
        }
        // Swap maximum row with current row (column by column)
        for (int k=i; k<n+1;k++) {
            double tmp = le_matrix(max_row,k);
            le_matrix(max_row,k) = le_matrix(i,k);
            le_matrix(i,k) = tmp;
        }

        if (std::abs(le_matrix(i,i)) < 1e-15) {
            return false;
        }
        // Make all rows below this one 0 in current column
        for (int k=i+1; k<n; k++) {
            double c = -le_matrix(k,i)/le_matrix(i,i);
            for (int j=i; j<n+1; j++) {
                if (i==j) {
                    le_matrix(k,j)= 0;
                } else {
                    le_matrix(k,j) += c * le_matrix(i,j);
                }
            }
        }
    }
    // Solve equation Ax=b for an upper triangular matrix A
    for (int i=n-1; i>=0; i--) {
        result_le[i] = le_matrix(i,n)/le_matrix(i,i);
        for (int k=i-1; k>=0; k--) {
            le_matrix(k,n) -= le_matrix(k,i) * result_le[i];
        }
    }

    //build large vector, if none of the elements is too negative
    size_t tracker = 0;
    for (size_t i = 0; i < dimension; i++) {
        if ((_c.support & (1ull << i)) != 0) {
            double x = result_le[tracker];
            if (x > -1e-5)
                result_vector[i] = x;
            else
                return false;
            tracker += 1;
        }
        else
            result_vector[i] = 0.;
    }

    // check p'Ap<=v for all rows not in the support
    double errorbound_rowsum=5e-5*dimension;

    for (size_t i = 0; i < dimension; i++) {
        if ((_c.support & (1ull << i)) == 0) { //not in the support - rows
            double rowsum = 0.;
            for (size_t j = 0; j < dimension; j++)
                if ((_c.support & (1ull << j)) != 0) // is in the support - columns
                    rowsum += game_matrix_double(i,j) * result_vector[j];

            if (!(rowsum <= result_le[_c.support_size] + errorbound_rowsum))
                return false;
        }
    }
    return true;
}


bool fracessa::find_candidate_rational(matrix<rational> &le_matrix)
{
    int n = _c.support_size + 1;

    std::vector<rational> result_le= std::vector<rational>(n, 0);
    std::vector<rational> result_vector= std::vector<rational>(dimension);

    game_matrix.get_le_matrix(_c.support,_c.support_size, le_matrix);

    ////////////////////////////////////////////////////////////////////////// gauss with partial pivoting
    for (int i=0; i<n; i++) {
        // Search for maximum in this column
        rational max_element = abs(le_matrix(i,i));
        int max_row = i;
        for (int k=i+1; k<n; k++) {
            if (abs(le_matrix(k,i)) > max_element) {
                max_element = abs(le_matrix(k,i));
                max_row = k;
            }
        }
        // Swap maximum row with current row (column by column)
        for (int k=i; k<n+1;k++) {
            rational tmp = le_matrix(max_row,k);
            le_matrix(max_row,k) = le_matrix(i,k);
            le_matrix(i,k) = tmp;
        }

        if (le_matrix(i,i)==rational(0)) {
            return false;
        }
        // Make all rows below this one 0 in current column
        for (int k=i+1; k<n; k++) {
            rational c = -le_matrix(k,i)/le_matrix(i,i);
            for (int j=i; j<n+1; j++) {
                if (i==j) {
                    le_matrix(k,j)= rational(0);
                } else {
                    le_matrix(k,j) += c * le_matrix(i,j);
                }
            }
        }

    }
    // Solve equation Ax=b for an upper triangular matrix A
    for (int i=n-1; i>=0; i--) {
        result_le[i] = le_matrix(i,n)/le_matrix(i,i);
        for (int k=i-1;k>=0; k--) {
            le_matrix(k,n) -= le_matrix(k,i) * result_le[i];
        }
    }

    //build large vector, if all elements are not too negative
    size_t tracker = 0;
    for (size_t i = 0; i < dimension; i++)
    {
        if ((_c.support & (1ull << i)) != 0)
        {
            rational x = result_le[tracker];
            if (x > rational(0))
                result_vector[i] = x;
            else
                return false;
            tracker += 1;
        }
        else
            result_vector[i] = rational(0);
    }

    // check p'Ap<=v for all rows not in the support
    uint64_t extended_support = _c.support;
    for (size_t i = 0; i < dimension; i++)
    {
        if ((_c.support & (1ull << i)) == 0) //not in the support - rows
        {
            rational rowsum = rational(0);
            for (size_t j = 0; j < dimension; j++)
                if ((_c.support & (1ull << j)) != 0) // is in the support - columns
                    rowsum += game_matrix(i,j) * result_vector[j];

            if (rowsum > result_le[_c.support_size])
                return false;
            if (rowsum==result_le[_c.support_size])
                extended_support = (extended_support | (1ull << i));
        }

    }
    _c.vector = result_vector;
    _c.extended_support = extended_support;
    _c.extended_support_size = support_size_from_int(extended_support);
    _c.payoff = result_le[_c.support_size];
    _c.payoff_double = static_cast<double>(_c.payoff);

    return true;
}


void fracessa::check_stability()
{
    uint64_t bitsetm = _c.support & (~(_c.support - 1)); //get lowest set bit as bitfield
    int extended_support_reduced = _c.extended_support & (~bitsetm); //ext support without m
    size_t m = postion_of_lowest_setbit(bitsetm);
    size_t extended_support_size_reduced = _c.extended_support_size - 1;

    if (conf_with_log && _logger) {
        _logger->info("Support: {}", std::bitset<64>(_c.support).to_string());
        _logger->info("Support size: {}", support_size_from_int(_c.support));
        _logger->info("Extended support: {}", std::bitset<64>(_c.extended_support).to_string());
        _logger->info("Extended support size: {}", _c.extended_support_size);
        _logger->info("Extended support reduced: {}", std::bitset<64>(extended_support_reduced).to_string());
        _logger->info("index m: {}", m);
    }

    if (extended_support_size_reduced == 0)
    {

        if (conf_with_log && _logger)
            _logger->info("Reason: true_pure_ess");
        _c.reason_ess = ReasonEss::true_pure_ess;
        _c.is_ess = true;
        return;
    }

    matrix<rational> bee = matrix<rational>(extended_support_size_reduced,extended_support_size_reduced);

    size_t row = 0;
    size_t column = 0;
    for (size_t i = 0;i<dimension;i++)
        if ((extended_support_reduced & (1ull << i)) != 0) {
            column = 0;
            for (size_t j = 0; j < i+1; j++)
                if ((extended_support_reduced & (1ull << j)) != 0) {
                    bee(row,column) = bee(column,row) = game_matrix(m, j) + game_matrix(j, m) + game_matrix(i, m) + game_matrix(m, i) -
                        game_matrix(i, j) - game_matrix(j, i) - 2 * game_matrix(m, m);
                    column += 1;
                }
            row += 1;
        }

    if (conf_with_log && _logger) {
        _logger->info("matrix bee:\n{}", bee.to_string());
    }

    if (!conf_exact) {
        if (bee.to_double().is_positive_definite_double()) {

            if (conf_with_log && _logger)
                _logger->info("Reason: true_posdef_double");
            _c.reason_ess = ReasonEss::true_posdef_double;
            _c.is_ess = true;
            return;
        }
    }

    if (bee.is_positive_definite()) {

        if (conf_with_log && _logger)
            _logger->info("Reason: true_posdef_rational");
        _c.reason_ess = ReasonEss::true_posdef_rational;
        _c.is_ess = true;
        return;
    }

    uint64_t kay = (_c.extended_support & (~_c.support)); //extended_support without support
    size_t kay_size = support_size_from_int(kay);

    if (conf_with_log && _logger)
        _logger->info("kay: {}", std::bitset<64>(kay).to_string());

    if (kay_size==0 || kay_size==1) {
        if (conf_with_log && _logger)
            _logger->info("Reason: false_not_posdef_and_kay_0_1");
        _c.reason_ess = ReasonEss::false_not_posdef_and_kay_0_1;
        _c.is_ess = false;
        return;
    }

    //do partial copositivity-check as in bomze_1992, p. 321/322
    uint64_t jay = extended_support_reduced;
    size_t r = support_size_from_int(jay & (~kay));

    std::vector<uint64_t> kay_vee(r+1);
    std::vector<size_t> kay_vee_size(r+1);
    std::vector<uint64_t> jay_without_kay_vee(r+1);
    std::vector< matrix<rational>> bee_vee(r+1);


    kay_vee[0] = jay;
    kay_vee_size[0] = extended_support_size_reduced;
    jay_without_kay_vee[0] = jay & (~kay);
    bee_vee[0] = bee;

    if (conf_with_log && _logger) {
        _logger->info("Partial Copositivity Check:");
        _logger->info("v=0:");
        _logger->info("kay_vee[0]: {}", std::bitset<64>(kay_vee[0]).to_string());
        _logger->info("kay_vee_size[0]: {}", kay_vee_size[0]);
        _logger->info("jay_without_kay_vee[0]: {}", std::bitset<64>(jay_without_kay_vee[0]).to_string());
        _logger->info("r: {}", r);
        _logger->info("bee_vee[0]:\n{}", bee_vee[0].to_string());
    }

    for (size_t v=1; v<=r; v++) {

        uint64_t iv = jay_without_kay_vee[v-1] & (~(jay_without_kay_vee[v-1]-1)); //iv is lowest set bit!
        jay_without_kay_vee[v] = jay_without_kay_vee[v-1] & (~iv); //remove iv from jay\kay
        kay_vee[v] = kay_vee[v-1] & (~iv); //build kay_vee
        kay_vee_size[v] = kay_vee_size[v-1]-1; //kay_vee_size
        bee_vee[v] = matrix<rational>(kay_vee_size[v],kay_vee_size[v]);

        uint64_t lowest_kv_minus_one = kay_vee[v-1] & (~(kay_vee[v-1] - 1)); //get lowest set bit of kay_vee[v-1]
        size_t index_position_to_remove = 0; //get the real distance of iv and lowest set bit of kay_vee[v-1]
        size_t iterater = 0;
        while (true) {
            uint64_t actual_bit_kv_minus_one = (lowest_kv_minus_one << iterater);
            if ((kay_vee[v-1] & actual_bit_kv_minus_one) != 0) {
                if ((iv & actual_bit_kv_minus_one) != 0)
                    break;
                else
                    index_position_to_remove++;
            }
            iterater++;
        }

        if (conf_with_log && _logger) {
            _logger->info("v={}:", v);
            _logger->info("kay_vee: {}", std::bitset<64>(kay_vee[v]).to_string());
            _logger->info("kay_vee_size: {}", kay_vee_size[v]);
            _logger->info("jay_without_kay_vee: {}", std::bitset<64>(jay_without_kay_vee[v]).to_string());
            _logger->info("iv: {}", std::bitset<64>(iv).to_string());
            _logger->info("Real index (distance) to remove: {}", index_position_to_remove);
        }

        if (bee_vee[v-1](index_position_to_remove,index_position_to_remove) <=0) {

            if (conf_with_log && _logger)
                _logger->info("Reason: false_not_partial_copositive");
            _c.reason_ess = ReasonEss::false_not_partial_copositive;
            _c.is_ess = false;
            return;
        }

        for (size_t i = 0;i< kay_vee_size[v];i++) { //build matrix bee_vee[v]
            for (size_t j = 0; j< kay_vee_size[v]; j++) {

                size_t row_index_kv_minus_one = (i>=index_position_to_remove) ? i+1 : i;
                size_t col_index_kv_minus_one = (j>=index_position_to_remove) ? j+1 : j;

                bee_vee[v](i,j) = bee_vee[v-1](index_position_to_remove,index_position_to_remove) * bee_vee[v-1](row_index_kv_minus_one,col_index_kv_minus_one) -
                    bee_vee[v-1](index_position_to_remove,row_index_kv_minus_one) * bee_vee[v-1](index_position_to_remove,col_index_kv_minus_one);
            }
        }

        if (conf_with_log && _logger) {
            _logger->info("bee_vee:\n{}", bee_vee[v].to_string());
        }
    }

    //copositivity check as in hadeler_1983
    if (conf_with_log && _logger)
        _logger->info("Copositivity Check:");

    if (bee_vee[r].is_copositive()) {
        if (conf_with_log && _logger)
            _logger->info("Reason: true_copositive");
        _c.reason_ess = ReasonEss::true_copositive;
        _c.is_ess = true;
        return;
    } else {
        if (conf_with_log && _logger)
            _logger->info("Reason: false_not_copositive");
        _c.reason_ess = ReasonEss::false_not_copositive;
        _c.is_ess = false;
        return;
    }
}
