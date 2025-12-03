#include <iostream>
#include <vector>
#include <cstdint>
#include <cassert>
#include <string>

#include <fracessa/ess_finder.hpp>
#include <fracessa/matrix.hpp>
#include <argparse/argparse.hpp>

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

    try {
        ::matrix<rational> A = ::matrix<rational>::create_from_cli_string(matrix_str);
        ess_finder x = ess_finder(A, candidates, exact, fullsupport, logger);

        std::cout << x.get_ess_count() << std::endl;

        if (candidates) {
            std::cout << candidate::header() << std::endl;
            for (const auto& c : x.get_candidates()) {
                std::cout << c.to_string() << std::endl;
            }
        }
    }
    catch (const std::runtime_error& err) {
        std::cerr << "Error: " << err.what() << std::endl;
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
