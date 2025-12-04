#include <fracessa/fracessa.hpp>
#include <fracessa/bitset64.hpp>

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
        if (_c.support.test(i)) {
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
        if (!_c.support.test(i)) { //not in the support - rows
            double rowsum = 0.;
            for (size_t j = 0; j < dimension; j++)
                if (_c.support.test(j)) // is in the support - columns
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

    game_matrix.get_le_matrix(_c.support, _c.support_size, le_matrix);

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
        if (_c.support.test(i))
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
    bitset64 extended_support = _c.support;
    for (size_t i = 0; i < dimension; i++)
    {
        if (!_c.support.test(i)) //not in the support - rows
        {
            rational rowsum = rational(0);
            for (size_t j = 0; j < dimension; j++)
                if (_c.support.test(j)) // is in the support - columns
                    rowsum += game_matrix(i,j) * result_vector[j];

            if (rowsum > result_le[_c.support_size])
                return false;
            if (rowsum==result_le[_c.support_size])
                extended_support.set(i);
        }

    }
    _c.vector = result_vector;
    _c.extended_support = extended_support;
    _c.extended_support_size = extended_support.count();
    _c.payoff = result_le[_c.support_size];
    _c.payoff_double = static_cast<double>(_c.payoff);

    return true;
}

