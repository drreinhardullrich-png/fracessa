#include <fracessa/ref.hpp>
#include <fracessa/helper.hpp>

ess_finder ref_from_cli(const std::string &matrix, bool with_candidates, bool exact, bool full_support, bool with_log)
{
    ::matrix<rational> A;
    try {
        A = ::matrix<rational>::create_from_cli_string(matrix);
    } catch (std::exception& e) {
        throw e;
    }
    ess_finder x = ess_finder(A, with_candidates, exact, full_support, with_log);
    return x;
}
