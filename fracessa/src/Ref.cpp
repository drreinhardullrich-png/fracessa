#include "../include/Ref.hpp"

EssFinder ref_from_cli(const std::string &matrix, bool with_candidates, bool exact, bool full_support, bool with_log)
{
    Matrix<rational> A;
    try {
        A = Matrix<rational>::create_from_cli_string(matrix);
    } catch (std::exception& e) {
        throw e;
    }
    EssFinder x = EssFinder(A, with_candidates, exact, full_support, with_log);
    return x;
}
