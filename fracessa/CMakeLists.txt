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

# Add compiler warnings
if(MSVC)
    add_compile_options(/W4)
else()
    add_compile_options(-Wall -Wextra -Wpedantic)
endif()

# Use local Boost headers directly
if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/third_party/boost_1_89_0_linux")
    message(STATUS "Using local Boost 1.89.0 headers")
    set(BOOST_LOCAL_FOUND TRUE)
else()
    message(FATAL_ERROR "Local Boost headers not found!")
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
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/third_party/argparse/include>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/third_party/boost_1_89_0_linux>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/third_party/spdlog/include>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/third_party/eigen>
    PRIVATE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>
)

# Test executable for positive definite check comparison
add_executable(test_posdef
    src/test_posdef.cpp
)
target_include_directories(test_posdef PRIVATE third_party/argparse/include)
target_link_libraries(test_posdef PRIVATE fracessa_lib)

# Test executable for Eigen LLT
add_executable(test_eigen_llt
    src/test_eigen_llt.cpp
)
target_include_directories(test_eigen_llt PRIVATE third_party/argparse/include)
target_link_libraries(test_eigen_llt PRIVATE fracessa_lib)

# Another test executable
add_executable(another_test
    src/another_test.cpp
)
target_include_directories(another_test PRIVATE third_party/argparse/include)
target_link_libraries(another_test PRIVATE fracessa_lib)

# Test executable for copositivity comparison
add_executable(test_copositivity
    src/test_copositivity.cpp
)
target_include_directories(test_copositivity PRIVATE third_party/argparse/include)
target_link_libraries(test_copositivity PRIVATE fracessa_lib)

# Platform-specific configuration
if(UNIX AND NOT APPLE)
    # Linux executable
    add_executable(fracessa
        src/main.cpp
    )

    target_include_directories(fracessa PRIVATE third_party/argparse/include)
    target_link_libraries(fracessa PRIVATE fracessa_lib)

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

    target_include_directories(fracessa PRIVATE third_party/argparse/include)
    target_link_libraries(fracessa PRIVATE fracessa_lib)

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
