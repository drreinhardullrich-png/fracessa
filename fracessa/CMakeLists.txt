cmake_minimum_required(VERSION 3.10.2)

# FORCE Release build type - this project only builds in Release mode
# Set before project() to ensure it's applied from the start
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
    message(STATUS "Build type not specified, forcing Release")
endif()

# Maximum optimization flags for Release builds
# -O3: Highest optimization level
# -DNDEBUG: Disable assertions for performance
# -flto: Link-Time Optimization (cross-module optimizations)
#   Note: For faster compilation with minimal runtime loss, use -flto=thin (requires GCC 6+ or Clang)
#   For older compilers, regular -flto is used (slower compile but same runtime performance)
# Note: -march=native and -mtune=native are commented out for portability
#       Uncomment them if you want CPU-specific optimizations (10-20% faster but not portable)
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG -flto" CACHE STRING "Flags used by the C++ compiler during release builds" FORCE)
set(CMAKE_C_FLAGS_RELEASE "-O3 -DNDEBUG -flto" CACHE STRING "Flags used by the C compiler during release builds" FORCE)
set(CMAKE_EXE_LINKER_FLAGS_RELEASE "-flto" CACHE STRING "Linker flags for Release builds" FORCE)
set(CMAKE_SHARED_LINKER_FLAGS_RELEASE "-flto" CACHE STRING "Shared linker flags for Release builds" FORCE)

# Disable Debug builds - they are not supported
# If someone tries to build Debug, they'll get an error above

project(fracessa
    VERSION 1.0.0
    DESCRIPTION "FRACESSA - Fractional ESS Analyzer - A solver for Standard Quadratic Problems"
    LANGUAGES CXX C
)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Enforce Release build type after project() is called
if(CMAKE_BUILD_TYPE STREQUAL "")
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
    message(STATUS "Build type not specified, setting to Release")
elseif(NOT CMAKE_BUILD_TYPE STREQUAL "Release")
    message(FATAL_ERROR "Only Release builds are supported. Current build type: ${CMAKE_BUILD_TYPE}")
endif()

# Apply maximum optimization flags to all targets
# These are applied in addition to CMAKE_CXX_FLAGS_RELEASE
if(CMAKE_BUILD_TYPE STREQUAL "Release" AND NOT MSVC)
    # Check if compiler supports thin LTO (GCC 6+ or Clang)
    # Thin LTO compiles faster than full LTO with ~95-99% of runtime performance
    include(CheckCXXCompilerFlag)
    check_cxx_compiler_flag("-flto=thin" COMPILER_SUPPORTS_THIN_LTO)
    
    if(COMPILER_SUPPORTS_THIN_LTO)
        add_compile_options(-flto=thin)
        add_link_options(-flto=thin)
        message(STATUS "Release build with maximum optimizations: -O3 -flto=thin (faster compilation)")
    else()
        # Fall back to regular LTO for older compilers
        add_compile_options(-flto)
        add_link_options(-flto)
        # Note: -flto-partition=none is not supported by all compiler versions
        # Removing it to avoid compilation errors
        message(STATUS "Release build with maximum optimizations: -O3 -flto")
    endif()
    
    # Additional aggressive optimizations (optional - uncomment if needed)
    # add_compile_options(-march=native -mtune=native)  # CPU-specific optimizations (not portable)
    # add_compile_options(-ffast-math)  # Aggressive math optimizations (may affect precision)
endif()

# Enable compiler cache for faster rebuilds (if available)
# ccache significantly speeds up recompilation by caching object files
find_program(CCACHE_PROGRAM ccache)
if(CCACHE_PROGRAM)
    set(CMAKE_CXX_COMPILER_LAUNCHER ${CCACHE_PROGRAM})
    set(CMAKE_C_COMPILER_LAUNCHER ${CCACHE_PROGRAM})
    message(STATUS "Using ccache for faster compilation: ${CCACHE_PROGRAM}")
else()
    message(STATUS "ccache not found - install it for faster rebuilds: sudo apt-get install ccache")
endif()

# Generate compile_commands.json for IDE support
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Optimize compilation speed with parallel builds
# This sets the default number of parallel jobs for cmake --build
# Users can still override with: cmake --build . --parallel N
include(ProcessorCount)
ProcessorCount(N)
if(NOT N EQUAL 0)
    set(CMAKE_BUILD_PARALLEL_LEVEL ${N} CACHE STRING "Number of parallel build jobs")
    message(STATUS "Parallel build level set to ${N} cores (override with: cmake --build . --parallel N)")
endif()

# Build configuration options
option(BUILD_SHARED_LIBS "Build shared libraries instead of static" OFF)
option(ENABLE_TESTING "Enable testing" OFF)
option(BUILD_TESTING "Build unit tests" OFF)
option(ENABLE_STATIC_LINKING "Enable static linking of all dependencies" OFF)

# Add compiler warnings
if(MSVC)
    add_compile_options(/W4)
else()
    add_compile_options(-Wall -Wextra -Wpedantic)
endif()

# Use FetchContent to download dependencies
# Set base directory to third_party so dependencies persist across cache clears
set(FETCHCONTENT_BASE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/third_party CACHE PATH "Directory for FetchContent downloads")
include(FetchContent)

# Fetch spdlog
FetchContent_Declare(
    spdlog
    GIT_REPOSITORY https://github.com/gabime/spdlog.git
    GIT_TAG        v1.13.0
)

# Fetch argparse
FetchContent_Declare(
    argparse
    GIT_REPOSITORY https://github.com/p-ranav/argparse.git
    GIT_TAG        v2.9
)

# Fetch Boost libraries (headers-only, individual repos for smaller downloads)
# Core dependencies first
FetchContent_Declare(
    boost_core
    GIT_REPOSITORY https://github.com/boostorg/core.git
    GIT_TAG        boost-1.89.0
    GIT_SHALLOW    TRUE
)

FetchContent_Declare(
    boost_config
    GIT_REPOSITORY https://github.com/boostorg/config.git
    GIT_TAG        boost-1.89.0
    GIT_SHALLOW    TRUE
)

FetchContent_Declare(
    boost_predef
    GIT_REPOSITORY https://github.com/boostorg/predef.git
    GIT_TAG        boost-1.89.0
    GIT_SHALLOW    TRUE
)

FetchContent_Declare(
    boost_assert
    GIT_REPOSITORY https://github.com/boostorg/assert.git
    GIT_TAG        boost-1.89.0
    GIT_SHALLOW    TRUE
)

FetchContent_Declare(
    boost_static_assert
    GIT_REPOSITORY https://github.com/boostorg/static_assert.git
    GIT_TAG        boost-1.89.0
    GIT_SHALLOW    TRUE
)

FetchContent_Declare(
    boost_throw_exception
    GIT_REPOSITORY https://github.com/boostorg/throw_exception.git
    GIT_TAG        boost-1.89.0
    GIT_SHALLOW    TRUE
)

FetchContent_Declare(
    boost_type_traits
    GIT_REPOSITORY https://github.com/boostorg/type_traits.git
    GIT_TAG        boost-1.89.0
    GIT_SHALLOW    TRUE
)

# multiprecision
FetchContent_Declare(
    boost_multiprecision
    GIT_REPOSITORY https://github.com/boostorg/multiprecision.git
    GIT_TAG        boost-1.89.0
    GIT_SHALLOW    TRUE
)

# math
FetchContent_Declare(
    boost_math
    GIT_REPOSITORY https://github.com/boostorg/math.git
    GIT_TAG        boost-1.89.0
    GIT_SHALLOW    TRUE
)

# integer
FetchContent_Declare(
    boost_integer
    GIT_REPOSITORY https://github.com/boostorg/integer.git
    GIT_TAG        boost-1.89.0
    GIT_SHALLOW    TRUE
)

# algorithm
FetchContent_Declare(
    boost_algorithm
    GIT_REPOSITORY https://github.com/boostorg/algorithm.git
    GIT_TAG        boost-1.89.0
    GIT_SHALLOW    TRUE
)

# rational (needed for boost/rational.hpp)
FetchContent_Declare(
    boost_rational
    GIT_REPOSITORY https://github.com/boostorg/rational.git
    GIT_TAG        boost-1.89.0
    GIT_SHALLOW    TRUE
)

# safe_numerics (needed for boost/safe_numerics/safe_integer.hpp)
FetchContent_Declare(
    boost_safe_numerics
    GIT_REPOSITORY https://github.com/boostorg/safe_numerics.git
    GIT_TAG        boost-1.89.0
    GIT_SHALLOW    TRUE
)

# Make all dependencies available (core dependencies first)
FetchContent_MakeAvailable(
    spdlog 
    argparse 
    boost_config
    boost_predef
    boost_assert
    boost_static_assert
    boost_throw_exception
    boost_type_traits
    boost_core
    boost_multiprecision 
    boost_math 
    boost_integer 
    boost_algorithm
    boost_rational
    boost_safe_numerics
)

message(STATUS "Using Boost 1.89.0 headers (from FetchContent)")
    set(BOOST_LOCAL_FOUND TRUE)

# Find GMP library for arbitrary precision rational numbers
# Try pkg-config first, then fall back to find_library
find_package(PkgConfig QUIET)
if(PkgConfig_FOUND)
    pkg_check_modules(GMP QUIET gmp)
    # Also find gmpxx (C++ wrapper) which is needed for stream operators
    find_library(GMPXX_LIBRARIES
        NAMES gmpxx
        PATHS
        /usr/lib
        /usr/lib/x86_64-linux-gnu
        /usr/local/lib
        /opt/local/lib
    )
    if(GMP_FOUND AND GMPXX_LIBRARIES)
        set(GMP_LIBRARIES ${GMP_LIBRARIES} ${GMPXX_LIBRARIES})
    endif()
endif()

if(NOT GMP_FOUND)
    # Fall back to manual search
    find_path(GMP_INCLUDE_DIR NAMES gmpxx.h
        PATHS
        /usr/include
        /usr/local/include
        /opt/local/lib
    )
    
    find_library(GMP_LIBRARIES
        NAMES gmp
        PATHS
        /usr/lib
        /usr/lib/x86_64-linux-gnu
        /usr/local/lib
        /opt/local/lib
    )
    
    find_library(GMPXX_LIBRARIES
        NAMES gmpxx
        PATHS
        /usr/lib
        /usr/lib/x86_64-linux-gnu
        /usr/local/lib
        /opt/local/lib
    )
    
    if(GMP_INCLUDE_DIR AND GMP_LIBRARIES AND GMPXX_LIBRARIES)
        set(GMP_FOUND TRUE)
        set(GMP_LIBRARIES ${GMP_LIBRARIES} ${GMPXX_LIBRARIES})
    endif()
endif()

if(GMP_FOUND)
    message(STATUS "Found GMP: ${GMP_LIBRARIES}")
    message(STATUS "GMP include dir: ${GMP_INCLUDE_DIR}")
    message(STATUS "Using GMP mpq_class backend for arbitrary precision rational numbers")
else()
    message(FATAL_ERROR "GMP library not found. Please install GMP: sudo apt-get install libgmp-dev libgmpxx4ldbl (or libgmpxx-dev)")
endif()

# Static linking support
if(ENABLE_STATIC_LINKING)
    message(STATUS "Static linking enabled - searching for static GMP libraries")
    
    # Find static GMP libraries
    # On Linux, look in standard library paths
    # On macOS, check both Homebrew locations (Intel and Apple Silicon)
    if(APPLE)
        # macOS: Check Homebrew paths (Apple Silicon first, then Intel)
        find_library(GMP_STATIC_LIB
            NAMES libgmp.a gmp
            PATHS
                /opt/homebrew/lib
                /usr/local/lib
                /opt/local/lib
            NO_DEFAULT_PATH
        )
        
        find_library(GMPXX_STATIC_LIB
            NAMES libgmpxx.a gmpxx
            PATHS
                /opt/homebrew/lib
                /usr/local/lib
                /opt/local/lib
            NO_DEFAULT_PATH
        )
    else()
        # Linux: Check standard library paths
        find_library(GMP_STATIC_LIB
            NAMES libgmp.a gmp
            PATHS
                /usr/lib
                /usr/lib/x86_64-linux-gnu
                /usr/lib64
                /usr/local/lib
            NO_DEFAULT_PATH
        )
        
        find_library(GMPXX_STATIC_LIB
            NAMES libgmpxx.a gmpxx
            PATHS
                /usr/lib
                /usr/lib/x86_64-linux-gnu
                /usr/lib64
                /usr/local/lib
            NO_DEFAULT_PATH
        )
    endif()
    
    # Also try with default paths if not found
    if(NOT GMP_STATIC_LIB)
        find_library(GMP_STATIC_LIB NAMES libgmp.a gmp)
    endif()
    if(NOT GMPXX_STATIC_LIB)
        find_library(GMPXX_STATIC_LIB NAMES libgmpxx.a gmpxx)
    endif()
    
    if(GMP_STATIC_LIB AND GMPXX_STATIC_LIB)
        set(GMP_LIBRARIES ${GMP_STATIC_LIB} ${GMPXX_STATIC_LIB})
        message(STATUS "Found static GMP libraries: ${GMP_LIBRARIES}")
        
        # Add platform-specific static linker flags
        if(UNIX AND NOT APPLE)
            # Linux: Only statically link GMP and standard C++ libraries
            # Do NOT use -static as it tries to link glibc statically (often fails)
            # Instead, link GMP statically and use static libgcc/libstdc++
            set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static-libgcc -static-libstdc++")
            message(STATUS "Linux static linking enabled: GMP static, libgcc/libstdc++ static")
        elseif(APPLE)
            # macOS: static libgcc/libstdc++ but dynamic system libraries
            set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static-libgcc -static-libstdc++")
            message(STATUS "macOS static linking enabled: -static-libgcc -static-libstdc++")
        endif()
    else()
        message(WARNING "Static linking requested but static GMP libraries not found.")
        message(WARNING "GMP_STATIC_LIB: ${GMP_STATIC_LIB}")
        message(WARNING "GMPXX_STATIC_LIB: ${GMPXX_STATIC_LIB}")
        message(WARNING "Falling back to dynamic linking. Binary will require GMP to be installed.")
        if(ENABLE_STATIC_LINKING)
            message(WARNING "Consider installing static GMP libraries or building GMP from source with --enable-static")
        endif()
    endif()
endif()


# Build the core library
add_library(fracessa_lib
    src/fracessa.cpp
    src/findeq.cpp
    src/checkstab.cpp
)

# Set library properties
set_target_properties(fracessa_lib PROPERTIES
    VERSION ${PROJECT_VERSION}
    SOVERSION ${PROJECT_VERSION_MAJOR}
)

# Include directories for the library
target_include_directories(fracessa_lib
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/fracessa>
        $<BUILD_INTERFACE:${argparse_SOURCE_DIR}/include>
        $<BUILD_INTERFACE:${spdlog_SOURCE_DIR}/include>
        ${GMP_INCLUDE_DIR}
    PRIVATE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>
)

# Boost libraries as system headers (suppresses deprecation warnings)
target_include_directories(fracessa_lib
    SYSTEM PUBLIC
        # Boost libraries (individual repos, core dependencies first)
        $<BUILD_INTERFACE:${boost_config_SOURCE_DIR}/include>
        $<BUILD_INTERFACE:${boost_predef_SOURCE_DIR}/include>
        $<BUILD_INTERFACE:${boost_assert_SOURCE_DIR}/include>
        $<BUILD_INTERFACE:${boost_static_assert_SOURCE_DIR}/include>
        $<BUILD_INTERFACE:${boost_throw_exception_SOURCE_DIR}/include>
        $<BUILD_INTERFACE:${boost_type_traits_SOURCE_DIR}/include>
        $<BUILD_INTERFACE:${boost_core_SOURCE_DIR}/include>
        $<BUILD_INTERFACE:${boost_multiprecision_SOURCE_DIR}/include>
        $<BUILD_INTERFACE:${boost_math_SOURCE_DIR}/include>
        $<BUILD_INTERFACE:${boost_integer_SOURCE_DIR}/include>
        $<BUILD_INTERFACE:${boost_algorithm_SOURCE_DIR}/include>
        $<BUILD_INTERFACE:${boost_rational_SOURCE_DIR}/include>
        $<BUILD_INTERFACE:${boost_safe_numerics_SOURCE_DIR}/include>
)

# Suppress Boost deprecation pragma messages
target_compile_definitions(fracessa_lib PRIVATE BOOST_ALLOW_DEPRECATED_HEADERS)

# Link GMP library
target_link_libraries(fracessa_lib PRIVATE ${GMP_LIBRARIES})


# Platform-specific configuration
if(UNIX AND NOT APPLE)
    # Linux executable
    add_executable(fracessa
        src/main.cpp
    )

    target_include_directories(fracessa PRIVATE 
        ${argparse_SOURCE_DIR}/include
        ${GMP_INCLUDE_DIR}
    )
    target_include_directories(fracessa SYSTEM PRIVATE 
        ${boost_math_SOURCE_DIR}/include
    )
    target_compile_definitions(fracessa PRIVATE BOOST_ALLOW_DEPRECATED_HEADERS)
    target_link_libraries(fracessa PRIVATE fracessa_lib ${GMP_LIBRARIES})

    # Install targets
    install(TARGETS fracessa fracessa_lib
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib
    )

elseif(WIN32 OR MSVC OR MSYS OR MINGW)
    # Windows-specific linker flags
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static-libgcc -static-libstdc++ -static")

    # Windows executable
    add_executable(fracessa
        src/main.cpp
    )

    target_include_directories(fracessa PRIVATE 
        ${argparse_SOURCE_DIR}/include
        ${GMP_INCLUDE_DIR}
    )
    target_include_directories(fracessa SYSTEM PRIVATE 
        ${boost_math_SOURCE_DIR}/include
    )
    target_compile_definitions(fracessa PRIVATE BOOST_ALLOW_DEPRECATED_HEADERS)
    target_link_libraries(fracessa PRIVATE fracessa_lib ${GMP_LIBRARIES})

    # Install targets
    install(TARGETS fracessa
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib
    )
endif()

# Install headers
install(DIRECTORY include/fracessa/
    DESTINATION include/fracessa
    FILES_MATCHING PATTERN "*.hpp"
)

# Enable testing if requested
if(ENABLE_TESTING)
    enable_testing()
    # Add test subdirectory if tests are added later
    # add_subdirectory(tests)
endif()

# Testing (optional, enabled with -DBUILD_TESTING=ON)
if(BUILD_TESTING)
    enable_testing()
    
    FetchContent_Declare(
        googletest
        GIT_REPOSITORY https://github.com/google/googletest.git
        GIT_TAG v1.14.0
    )
    FetchContent_MakeAvailable(googletest)
    
    add_subdirectory(tests)
endif()
