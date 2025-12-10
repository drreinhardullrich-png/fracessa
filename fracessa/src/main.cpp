#include <iostream>
#include <vector>
#include <string>
#include <chrono>
#include <iomanip>
#include <cstdint>

#include <fracessa/fracessa.hpp>
#include <rational_linalg/types_rational.hpp>
#include <rational_linalg/matrix.hpp>
#include <argparse/argparse.hpp>

// Helper function to parse matrix string format: "n#values"
// Returns true on success, false on error
bool parse_matrix_string(const std::string& matrix_str, rational_linalg::Matrix<small_rational>& A, bool& is_cs)
{
    // Parse CLI string format: "n#values" - optimized manual parsing
    const size_t hash_pos = matrix_str.find('#');
    if (hash_pos == std::string::npos || hash_pos == 0 || hash_pos == matrix_str.length() - 1) {
        std::cerr << "Error: String for the matrix does not include '#' as a separator between dimension and matrix!" << std::endl;
        return false;
    }
    
    // Check for multiple '#' characters
    if (matrix_str.find('#', hash_pos + 1) != std::string::npos) {
        std::cerr << "Error: Multiple '#' characters found in matrix string!" << std::endl;
        return false;
    }
    
    size_t n;
    try {
        n = std::stoull(matrix_str.substr(0, hash_pos));
    } catch (const std::exception& e) {
        std::cerr << "Error: The given dimension could not be converted into an integer number!" << std::endl;
        return false;
    }
    
    // Parse comma-separated values - optimized manual parsing
    const std::string& values_str = matrix_str.substr(hash_pos + 1);
    std::vector<small_rational> rational_values;
    
    // Pre-allocate vector with estimated size (at least n/2 for circular symmetric)
    rational_values.reserve(n / 2);
    
    try {
        size_t start = 0;
        size_t comma_pos;
        
        while (start < values_str.length()) {
            comma_pos = values_str.find(',', start);
            if (comma_pos == std::string::npos) {
                comma_pos = values_str.length();
            }
            
            // Parse value: "numerator/denominator" or just "numerator"
            const size_t slash_pos = values_str.find('/', start);
            if (slash_pos != std::string::npos && slash_pos < comma_pos) {
                // Has denominator
                int64_t num = std::stoll(values_str.substr(start, slash_pos - start));
                int64_t den = std::stoll(values_str.substr(slash_pos + 1, comma_pos - slash_pos - 1));
                rational_values.emplace_back(num, den);
            } else {
                // No denominator, treat as integer
                int64_t num = std::stoll(values_str.substr(start, comma_pos - start));
                rational_values.emplace_back(num, 1);
            }
            
            start = comma_pos + 1;
        }
    } catch (const std::exception& e) {
        std::cerr << "Error: Could not convert matrix values to rational numbers!" << std::endl;
        std::cerr << "  " << e.what() << std::endl;
        return false;
    }
    
    // Determine matrix type and create matrix
    const size_t expected_cs = n / 2;
    const size_t expected_sym = n * (n + 1) / 2;
    const size_t actual_size = rational_values.size();
    
    if (actual_size == expected_cs) {
        // Circular symmetric matrix
        A = rational_linalg::create_circular_symmetric<small_rational>(n, rational_values);
        is_cs = true;
    } else if (actual_size == expected_sym) {
        // Symmetric matrix (upper triangular)
        A = rational_linalg::create_symmetric<small_rational>(n, rational_values);
        is_cs = false;
    } else {
        std::cerr << "Error: The number of matrix-elements must either be floor(dimension/2) (for a circular symmetric matrix) or dimension*(dimension+1)/2 (for a symmetric matrix)!" << std::endl;
        std::cerr << "  Got " << actual_size << " values, but expected " << expected_cs << " (circular symmetric) or " << expected_sym << " (symmetric)." << std::endl;
        return false;
    }
    
    return true;
}

// Unsafe version: parses matrix string without any safety checks for maximum performance
// Assumes valid input format: "n#values" where values are comma-separated rationals
void parse_matrix_string_unsafe(const std::string& matrix_str, rational_linalg::Matrix<small_rational>& A, bool& is_cs)
{
    // Find '#' by direct iteration
    size_t hash_pos = 0;
    while (matrix_str[hash_pos] != '#') {
        ++hash_pos;
    }
    
    // Parse dimension n manually (no stoull, no substr)
    size_t n = 0;
    for (size_t i = 0; i < hash_pos; ++i) {
        n = n * 10 + (matrix_str[i] - '0');
    }
    
    // Pre-allocate vector with estimated size
    std::vector<small_rational> rational_values;
    rational_values.reserve(n / 2);
    
    // Parse values by direct character iteration
    size_t pos = hash_pos + 1;
    const size_t len = matrix_str.length();
    
    while (pos < len) {
        // Parse numerator
        int64_t num = 0;
        bool num_negative = false;
        
        if (matrix_str[pos] == '-') {
            num_negative = true;
            ++pos;
        }
        
        while (pos < len && matrix_str[pos] != '/' && matrix_str[pos] != ',') {
            num = num * 10 + (matrix_str[pos] - '0');
            ++pos;
        }
        
        if (num_negative) {
            num = -num;
        }
        
        // Check if there's a denominator
        if (pos < len && matrix_str[pos] == '/') {
            ++pos; // skip '/'
            
            // Parse denominator
            int64_t den = 0;
            bool den_negative = false;
            
            if (matrix_str[pos] == '-') {
                den_negative = true;
                ++pos;
            }
            
            while (pos < len && matrix_str[pos] != ',') {
                den = den * 10 + (matrix_str[pos] - '0');
                ++pos;
            }
            
            if (den_negative) {
                den = -den;
            }
            
            rational_values.emplace_back(num, den);
        } else {
            // No denominator, treat as integer
            rational_values.emplace_back(num, 1);
        }
        
        // Skip comma if present
        if (pos < len && matrix_str[pos] == ',') {
            ++pos;
        }
    }
    
    // Determine matrix type and create matrix
    const size_t expected_cs = n / 2;
    const size_t actual_size = rational_values.size();
    
    if (actual_size == expected_cs) {
        // Circular symmetric matrix
        A = rational_linalg::create_circular_symmetric<small_rational>(n, rational_values);
        is_cs = true;
    } else {
        // Symmetric matrix (upper triangular)
        A = rational_linalg::create_symmetric<small_rational>(n, rational_values);
        is_cs = false;
    }
}

int main(int argc, char *argv[])
{
    argparse::ArgumentParser program("fracessa", "3.0.0");

    program.add_description("FRACESSA - Fractional ESS Analyzer - A solver for Standard Quadratic Problems");

    program.add_argument("-c", "--candidates")
        .help("include the found candidates for ESS/solutions in the output")
        .implicit_value(true)
        .default_value(false);

    program.add_argument("-l", "--log")
        .help("output a detailed log file named 'fracessa.log' in the directory of the program, for diagnostic of learning purposes only")
        .implicit_value(true)
        .default_value(false);

    program.add_argument("-e", "--exact")
        .help("only uses rational numbers, for matrices with extreme differences in the input, is much much slower!")
        .implicit_value(true)
        .default_value(false);

    program.add_argument("-f", "--fullsupport")
        .help("searches the full support directly after searching support size one. Enable if you expect the matrix to have exactly one ess in the interior of the simplex!")
        .implicit_value(true)
        .default_value(false);

    program.add_argument("-t", "--timing")
        .help("output the computation time in seconds on a new line after the ESS count")
        .implicit_value(true)
        .default_value(false);

    program.add_argument("-m", "--matrixid")
        .help("optional matrix ID to write in the log file")
        .scan<'i', int>()
        .default_value(-1);

    program.add_argument("-u", "--unsafe")
        .help("use unsafe matrix string parsing for maximum performance (no input validation)")
        .implicit_value(true)
        .default_value(false);

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

    const auto& matrix_str = program.get<std::string>("matrix");
    const auto candidates = program.get<bool>("--candidates");
    const auto logger = program.get<bool>("--log");
    const auto exact = program.get<bool>("--exact");
    const auto fullsupport = program.get<bool>("--fullsupport");
    const auto timing = program.get<bool>("--timing");
    const auto matrix_id = program.get<int>("--matrixid");
    const auto unsafe = program.get<bool>("--unsafe");

    // Parse matrix from CLI string format: "n#values"
    bool is_cs;
    rational_linalg::Matrix<small_rational> A;
    if (unsafe) {
        parse_matrix_string_unsafe(matrix_str, A, is_cs);
    } else {
        if (!parse_matrix_string(matrix_str, A, is_cs)) {
            return EXIT_FAILURE;
        }
    }
    
    // Measure computation time only if timing flag is set
    std::chrono::high_resolution_clock::time_point start_time, end_time;
    double elapsed_seconds = 0.0;
    
    if (timing) {
        start_time = std::chrono::high_resolution_clock::now();
    }
    
    ::fracessa x(A, is_cs, candidates, exact, fullsupport, logger, matrix_id);
    
    if (timing) {
        end_time = std::chrono::high_resolution_clock::now();
        // Calculate elapsed time in seconds
        const auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
        elapsed_seconds = duration.count() / 1000000.0;
    }

    std::cout << x.ess_count_ << std::endl;
    
    // Output timing on second line if -t flag is present
    if (timing) {
        std::cout << std::fixed << std::setprecision(6) << elapsed_seconds << std::endl;
    }

    if (candidates) {
        std::cout << candidate::header() << std::endl;
        for (auto& c : x.candidates_) {
            std::cout << c.to_string() << std::endl;
        }
    }

    return EXIT_SUCCESS;
}
