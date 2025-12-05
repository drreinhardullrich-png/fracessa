#include <iostream>
#include <vector>
#include <cassert>
#include <string>

#include <fracessa/fracessa.hpp>
#include <fracessa/matrix.hpp>
#include <argparse/argparse.hpp>
#include <boost/algorithm/string.hpp>

int main(int argc, char *argv[])
{
    argparse::ArgumentParser program("fracessa", "3.0.0");

    program.add_description("FRACESSA - Fractional ESS Analyzer - A solver for Standard Quadratic Problems");

    program.add_argument("-c", "--candidates")
        .help("include the found candidates for ESS/solutions in the output")
        .flag();

    program.add_argument("-l", "--log")
        .help("output a detailed log file named 'fracessa.log' in the directory of the program, for diagnostic of learning purposes only")
        .flag();

    program.add_argument("-e", "--exact")
        .help("only uses rational numbers, for matrices with extreme differences in the input, is much much slower!")
        .flag();

    program.add_argument("-f", "--fullsupport")
        .help("searches the full support directly after searching support size one. Enable if you expect the matrix to have exactly one ess in the interior of the simplex!")
        .flag();

    program.add_argument("matrix")
        .help("the matrix to compute");

    try {
        program.parse_args(argc, argv);
    }
    catch (const std::runtime_error& err) {
        std::cerr << err.what() << std::endl;
        std::cerr << program;
        return EXIT_FAILURE;
    }

    auto matrix_str = program.get<std::string>("matrix");
    auto candidates = program.get<bool>("--candidates");
    auto logger = program.get<bool>("--log");
    auto exact = program.get<bool>("--exact");
    auto fullsupport = program.get<bool>("--fullsupport");

    // Parse CLI string format: "n#values"
    std::vector<std::string> first_split;
    boost::split(first_split, matrix_str, boost::is_any_of("#"));
    if (first_split.size() != 2) {
        std::cerr << "Error: String for the matrix does not include '#' as a separator between dimension and matrix, or several '#' were found!" << std::endl;
        return EXIT_FAILURE;
    }
    
    size_t n;
    try {
        n = std::stoull(first_split[0]);
    } catch (std::exception& e) {
        std::cerr << "Error: The given dimension could not be converted into an integer number!" << std::endl;
        return EXIT_FAILURE;
    }
    
    std::vector<std::string> second_split;
    boost::split(second_split, first_split[1], boost::is_any_of(","));
    
    // Convert string values to rational
    std::vector<rational> rational_values;
    try {
        for (const auto& str_val : second_split) {
            rational_values.push_back(rational(str_val));
        }
    } catch (std::exception& e) {
        std::cerr << "Error: Could not convert matrix values to rational numbers!" << std::endl;
        std::cerr << "  " << e.what() << std::endl;
        return EXIT_FAILURE;
    }
    
    RationalMatrix A;
    bool is_cs;
    
    if (rational_values.size() == n/2) {
        // Circular symmetric matrix
        A = matrix_ops::create_circular_symmetric(n, rational_values);
        is_cs = true;
    } else if (rational_values.size() == n*(n+1)/2) {
        // Symmetric matrix (upper triangular)
        A = matrix_ops::create_symmetric(n, rational_values);
        is_cs = false;
    } else {
        std::cerr << "Error: The number of matrix-elements must either be floor(dimension/2) (for a circular symmetric matrix) or dimension*(dimension+1)/2 (for a symmetric matrix)!" << std::endl;
        std::cerr << "  Got " << rational_values.size() << " values, but expected " << n/2 << " (circular symmetric) or " << n*(n+1)/2 << " (symmetric)." << std::endl;
        return EXIT_FAILURE;
    }
    
    ::fracessa x = ::fracessa(A, is_cs, candidates, exact, fullsupport, logger);

    std::cout << x.ess_count << std::endl;

    if (candidates) {
        std::cout << candidate::header() << std::endl;
        for (auto c: x.candidates) {
            std::cout << c.to_string() << std::endl;
        }
    }

    return EXIT_SUCCESS;
}
