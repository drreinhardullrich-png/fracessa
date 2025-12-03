# FRACESSA - Fractional ESS Analyzer

A solver for Standard Quadratic Problems

Introduction
------------

FRACESSA is a command-line C++ program for calculating all strict local maximizers of a Standard Quadratic Problem (StQP). From the point of view of Evolutionary Game Theory this is equivalent to finding all Evolutionary Stable Strategies (ESSs) of a partnership game. The algorithm for doing so was originally described in I.M. Bomze: "Detecting all evolutionarily stable strategies", Journal of Optimization Theory and Applications 75.2 (1992), and this implementation is a further development and improvement of it. It was most recently used [here](http://www.optimization-online.org/DB_HTML/2016/05/5452.html), although the implementational details are not explained.

FRACESSA takes a rational, not necessarily symmetric n × n-matrix as an argument and outputs the number of ESSs (as well as detailed information about them) in exact arithmetics. It does so by traversing (potentially) all the 2^n support sets of {1,...,n}, and FRACESSA can (theoretically) handle matrices up to n=64. In practice, depending on your CPU, your memory and your time, n up to 30 or 40 is feasible.

Downloads
---------

-   [Windows 64bit](fracessa_win64.zip)
-   [Ubuntu/Debian 64bit](fracessa_ubuntu64.zip)
-   MacOsX - cross-compilation is quite complicated, please compile it yourself!
-   [C++ Source](fracessa_source.zip) (CMake project)

Initial Test
------------

-   Download the appropriate .zip for your platform.
-   Unpack the .zip and navigate to the directory.
-   Initial test: run `fracessa 21\#15,15,7,15,15,7,7,15,7,15`
-   FRACESSA should run for 5-30 seconds (depending on your hardware) and write "4410" as output.

Documentation
-------------

Now you are able to use FRACESSA, this is possible in two ways.

-   You manually type into the command line - as done for the initial test. This is only recommended for experimenting with the command line options or for quick checks.
-   You call the command line from a script/program of your choice. This is the way the program is designed to run.

The program can be seperated in three parts: fracessa \[options\] \[matrix\]

-   fracessa is just the name of the program you call, it is always the same.
-   \[options\] are the so-called command-line-options, see below.
-   \[matrix\] is the matrix you want to search for ESS, see below.

### Input \[options\]

-   **-v or --vectors**
    Usually the program only outputs the number of ESS a matrix has, this is useful if you want to search many matrices fast. The option -v enables you to output the whole information which is stored about the solution of your problem. This is done in a .csv-style manner, see section "Resulting Vectors".
-   **-f or --fullsupport**
    Usually the program starts with the smalles support sets and searches up to the bigger ones. This behaviour is favoured if you expect the matrix to have many ESS. But if you expect the matrix to have only one ESS with full support, you should set this option. Now first the support sets of size one are searched, next the whole support is searched, and afterwards the support sets from size 2 to n-1 are searched in increasing order.
-   **-e or --exact**
    Usually most of the calculations are done as floating points (double), and only if the outcome seems to be a solution then the calculation is repeated with rational numbers. This vastly increases the speed of the algorithm and works in almost all cases. But if the differences of the magnitude of the matrix inputs are huge (say more than 10^10), this could lead to unnoticed errors, since the double-datatyp only has a precision of 15-17 digits. This option disables the use of floating points, all calculations are done entirely as rationals. You can also enable this option to double-check a result.
-   **-l or --log**
    This option generates a detailed log-file named "Log.txt" in the directory where the program runs. For diagnostic- or learning-purposes only!

### Input \[matrix\]

At the moment two different ways of passing a matrix are supported

-   **Whole matrix**
    Example: 2\#5,6/7,-1,-2/3
    First the dimension of the matrix is given, then a hash-symbol (\#), afterwards the numbers as integers or rationals, where rationals are written with a slash (/). All numbers can be positive (nothing) or negative (-) and are seperated by a comma (,). Whitespaces or other characters are not allowed! The numbers must be in row-major-order, that means that the first *row* of the matrix is given as (5,6/7) and the second row as (-1,-2/3).
-   **Cyclically symmetric matrix**
    Example: 7\#2/3,5,9
    Since cyclically symmetric matrices play a major role in the search for hard cases for StQP's, there is an shortcut for inputting them. The input in the example generates a 7x7 matrix where the first line is given as (2/3,5,9,9,5,2/3,0), ie. take the vecor, glue it together with its mirrored version and glue a zero at the end. Now this row is the last row of the matrix, all other rows are generated by rotational left shifts of order n-i, where i is current row index. Note that if n is even then the last number of the input is not doubled. Use this type of input also if you want to enable the optimziations for cyclically symmetric matrices, see section "Cyclically Symmetric Opimization".

### Full Examples

fracessa -v 2\#5,6/7,-1,-2/3
fracessa -l -f -e 19\#1,2,2,2,2,2,1,1,2
fracessa 23\#27478,22664,10976,25676,18552,18552,25676,10976,22664,27478,17939

A script in R which tests some matrices in multi-threading can be found [here](ref_tester.R).

### Resulting Vectors

If the option -v is enabled, then the lines of output of FRACESSA consist of the following:
\[Number of ESSs\]
\[CSV Header\]
\[CSV Data\]
\[CSV Data\]
\[CSV Data\]
...

One easy way to process the data is to seperate the first from the following lines, and parse these lines by a .csv-parser. The field sperator is a semicolon, field delimiters are not used.

See a full example here, the output of `fracessa -v 5\#1,0,2,2,2,0,1,2,2,2,2,2,1,0,0,2,2,0,1,0,2,2,0,0,0`:
4
VectorID;Vector;Support;SupportSize;ExtendedSupport;ExtendedSupportSize;ShiftReference;IsEss;Reason;Payoff;PayoffDecimal
1;1/2,0,1/2,0,0;5;2;5;2;0;1;2;3/2;1.500000
2;0,1/2,1/2,0,0;6;2;6;2;0;1;2;3/2;1.500000
3;1/2,0,0,1/2,0;9;2;9;2;0;1;2;3/2;1.500000
4;0,1/2,0,1/2,0;10;2;10;2;0;1;2;3/2;1.500000
5;2/3,0,0,0,1/3;17;2;29;4;0;0;7;4/3;1.333333
6;0,2/3,0,0,1/3;18;2;30;4;0;0;7;4/3;1.333333

FRACESSA with option -v does not only save all the ESSs, but all the "Candidates". A serious candidate for an ESS is an isolated Nash Equilibrium strategy (equivalently an isolated KKT-point), which is in the relative interior of the considered support, but only if no other candidate's support is contained in its support. Sloppily speaking, only the smaller isolated Nash Equilibrium strategies are candidates. For details see the references in the "Introduction".

Description of the different fields:

-   **VectorID**
    Auto-incrementing integer.
-   **Vector**
    Returns the vector of the ESS of dimension n as a rational number. The encoding is the same as for the input matrices.
-   **Support**
    The support set of the ESS coded in binary, e.g. support 6 means that the support set is {2,3}. Always smaller than 2^n.
-   **SupportSize**
    The amount of numbers in the support set
-   **ExtendedSupport**
    The extended support coded in binary. See the references in the "Introduction" for a definition.
-   **ExtendedSupportSize**
    The amount of numbers in the extended support set.
-   **ShiftReference**
    Only used for cyclically symmetric matrices, see section "Cyclically Symmetric Optimization", otherwise 0.
-   **IsEss**
    Is 1 if the candidate is an ESS, otherwise 0.
-   **Reason**
    Since the procdure to determine if a candidate is an ESS is quite cumbersome, here the reason for the decision is stated as an integer. The coding is the following:
    true\_pure\_ess=1
    true\_posdef\_double = 2
    true\_posdef\_rational=3
    true\_copositive=4
    false\_not\_posdef\_and\_K\_0\_1 = 5
    false\_not\_partial\_copositive = 6
    false\_not\_copositive=7
    Numbers 1 to 4 mean that the candidate is an ESS, 5 to 7 mean that it is not. For details see the references in the "Introduction" and the source code.
-   **Payoff**
    The value of the maximum for this particular maximizer (or equivalently the average population payoff for the ESS) as a rational number.
-   **PayoffDecimal**
    The same as a decimal number.

### Cyclically Symmetric Opimization

Whenever the shortcut for inputting cyclically symmetric matrices is used, a certain optimization is performed, utilizing the symmetrie properties of the matrix. Careful, if a cyclically symmetric matrix is inputted without the shortcut, then this optimization is not performed!
If a candidate/ESS for a cyclically symmetric matrix is found and the support size of the candidate/ESS and n are coprime, then it is known that another n-1 candidates/ESSs exist, they are just shifted (rotated) versions of the original candidate/ESS. The optimization exploits that fact, and the n-1 other candidates/ESSs are just recorded and not calculated. This increases the speed of FRACESSA, depending on n and the size of the numbers of the input matrix. For n being prime the increase is clearly the most.
The field ShiftReference for the vector-output references the VectorID of the originally calculated ESS (and this ESS also references itself). For analyzing this cyclical structure this can be used by querying the .csv for VectorID==ShiftReference.

License
-------

The software FRACESSA is open-source under the GNU General Public License, Version 3. No warranty or liability of any kind is provided.
Other open-source software used by FRACESSA, with special thanks:

-   [Boost C++ libraries](http://www.boost.org/)
-   [argparse](https://github.com/p-ranav/argparse)

Contact
-------

For questions, bug reports and suggestions please write to: reinhard.ullrich (at) univie.ac.at.
## Project Structure

This project follows modern C++ project organization:

```
fracessa/
├── include/                 # Public headers
│   └── fracessa/           # Main public API headers
│       ├── candidate.hpp
│       ├── ess_finder.hpp
│       ├── helper.hpp
│       ├── matrix.hpp
│       ├── ref.hpp
│       └── argparse.hpp
├── src/                    # Source files
│   ├── ess_finder.cpp
│   ├── main.cpp
│   └── ref.cpp
├── third_party/            # External dependencies
│   ├── argparse/          # Command-line argument parsing
│   └── boost_1_89_0_linux/ # Boost C++ libraries
├── CMakeLists.txt          # Build configuration
├── README.md              # This file
└── .gitignore            # Git ignore rules
```

## Building

### Prerequisites

- CMake 3.10.2 or later
- C++17 compatible compiler
- Boost libraries (included in third_party/)

### Build Instructions

```bash
# Create build directory
mkdir build
cd build

# Configure with CMake
cmake ..

# Build
make

# Install (optional)
make install
```

### Build Options

- `BUILD_SHARED_LIBS`: Build shared libraries instead of static (default: OFF)
- `ENABLE_TESTING`: Enable testing infrastructure (default: OFF)

## Library Usage

FRACESSA can also be used as a library in your C++ projects:

```cpp
#include <fracessa/EssFinder.hpp>
#include <fracessa/Matrix.hpp>

// Create a matrix from string input
auto matrix = Matrix<rational>::create_from_cli_string("2#5,6/7,-1,-2/3");

// Find ESS with candidates output
auto finder = EssFinder(matrix, true, false, false, false);

std::cout << "Found " << finder.ess_count << " ESS solutions" << std::endl;
```

## Usage

See the original documentation below for command-line usage.

---

Last updated December 2025 - Project renamed to FRACESSA and restructured for modern C++ development.
