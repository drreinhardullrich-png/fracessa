cmake_minimum_required(VERSION 3.10.2)

# FORCE Release build type - this project only builds in Release mode
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
    message(STATUS "Build type not specified, forcing Release")
endif()

# Maximum optimization flags for Release builds
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG" CACHE STRING "Flags used by the C++ compiler during release builds" FORCE)
set(CMAKE_C_FLAGS_RELEASE "-O3 -DNDEBUG" CACHE STRING "Flags used by the C compiler during release builds" FORCE)

project(fracessa
    VERSION 1.0.0
    DESCRIPTION "FRACESSA - Fractional ESS Analyzer - A solver for Standard Quadratic Problems"
    LANGUAGES CXX C
)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Enforce Release build type
if(CMAKE_BUILD_TYPE STREQUAL "")
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
elseif(NOT CMAKE_BUILD_TYPE STREQUAL "Release")
    message(FATAL_ERROR "Only Release builds are supported. Current build type: ${CMAKE_BUILD_TYPE}")
endif()

# Apply LTO only if supported (check after project())
if(CMAKE_BUILD_TYPE STREQUAL "Release" AND NOT MSVC)
    include(CheckIPOSupported)
    check_ipo_supported(RESULT IPO_SUPPORTED OUTPUT IPO_ERROR)
    
    if(IPO_SUPPORTED)
        set(CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)
        message(STATUS "LTO/IPO enabled via CMAKE_INTERPROCEDURAL_OPTIMIZATION")
    else()
        message(STATUS "LTO/IPO not supported: ${IPO_ERROR}")
    endif()
endif()

# Enable compiler cache for faster rebuilds
find_program(CCACHE_PROGRAM ccache)
if(CCACHE_PROGRAM)
    set(CMAKE_CXX_COMPILER_LAUNCHER ${CCACHE_PROGRAM})
    set(CMAKE_C_COMPILER_LAUNCHER ${CCACHE_PROGRAM})
    message(STATUS "Using ccache: ${CCACHE_PROGRAM}")
endif()

# Generate compile_commands.json
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Parallel builds
include(ProcessorCount)
ProcessorCount(N)
if(NOT N EQUAL 0)
    set(CMAKE_BUILD_PARALLEL_LEVEL ${N} CACHE STRING "Number of parallel build jobs")
    message(STATUS "Parallel build level: ${N} cores")
endif()

# Build options
option(BUILD_SHARED_LIBS "Build shared libraries" OFF)
option(ENABLE_TESTING "Enable testing" OFF)
option(BUILD_TESTING "Build unit tests" OFF)
option(ENABLE_STATIC_LINKING "Enable static linking" OFF)

# Compiler warnings
if(MSVC)
    add_compile_options(/W4)
else()
    add_compile_options(-Wall -Wextra -Wpedantic)
endif()

# FetchContent setup
set(FETCHCONTENT_BASE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/third_party CACHE PATH "FetchContent base dir")
include(FetchContent)

# Declare all dependencies
FetchContent_Declare(spdlog
    GIT_REPOSITORY https://github.com/gabime/spdlog.git
    GIT_TAG v1.13.0
)

FetchContent_Declare(argparse
    GIT_REPOSITORY https://github.com/p-ranav/argparse.git
    GIT_TAG v2.9
)

# Boost dependencies
set(BOOST_TAG boost-1.89.0)
set(BOOST_SHALLOW TRUE)

FetchContent_Declare(boost_core
    GIT_REPOSITORY https://github.com/boostorg/core.git
    GIT_TAG ${BOOST_TAG}
    GIT_SHALLOW ${BOOST_SHALLOW}
)

FetchContent_Declare(boost_config
    GIT_REPOSITORY https://github.com/boostorg/config.git
    GIT_TAG ${BOOST_TAG}
    GIT_SHALLOW ${BOOST_SHALLOW}
)

FetchContent_Declare(boost_predef
    GIT_REPOSITORY https://github.com/boostorg/predef.git
    GIT_TAG ${BOOST_TAG}
    GIT_SHALLOW ${BOOST_SHALLOW}
)

FetchContent_Declare(boost_assert
    GIT_REPOSITORY https://github.com/boostorg/assert.git
    GIT_TAG ${BOOST_TAG}
    GIT_SHALLOW ${BOOST_SHALLOW}
)

FetchContent_Declare(boost_static_assert
    GIT_REPOSITORY https://github.com/boostorg/static_assert.git
    GIT_TAG ${BOOST_TAG}
    GIT_SHALLOW ${BOOST_SHALLOW}
)

FetchContent_Declare(boost_throw_exception
    GIT_REPOSITORY https://github.com/boostorg/throw_exception.git
    GIT_TAG ${BOOST_TAG}
    GIT_SHALLOW ${BOOST_SHALLOW}
)

FetchContent_Declare(boost_type_traits
    GIT_REPOSITORY https://github.com/boostorg/type_traits.git
    GIT_TAG ${BOOST_TAG}
    GIT_SHALLOW ${BOOST_SHALLOW}
)

FetchContent_Declare(boost_multiprecision
    GIT_REPOSITORY https://github.com/boostorg/multiprecision.git
    GIT_TAG ${BOOST_TAG}
    GIT_SHALLOW ${BOOST_SHALLOW}
)

FetchContent_Declare(boost_math
    GIT_REPOSITORY https://github.com/boostorg/math.git
    GIT_TAG ${BOOST_TAG}
    GIT_SHALLOW ${BOOST_SHALLOW}
)

FetchContent_Declare(boost_integer
    GIT_REPOSITORY https://github.com/boostorg/integer.git
    GIT_TAG ${BOOST_TAG}
    GIT_SHALLOW ${BOOST_SHALLOW}
)

FetchContent_Declare(boost_algorithm
    GIT_REPOSITORY https://github.com/boostorg/algorithm.git
    GIT_TAG ${BOOST_TAG}
    GIT_SHALLOW ${BOOST_SHALLOW}
)

# Make dependencies available
FetchContent_MakeAvailable(
    spdlog argparse
    boost_config boost_predef boost_assert boost_static_assert
    boost_throw_exception boost_type_traits boost_core
    boost_multiprecision boost_math boost_integer boost_algorithm
)

message(STATUS "Using Boost 1.89.0 headers")

# Find GMP library
find_package(PkgConfig QUIET)
if(PkgConfig_FOUND)
    pkg_check_modules(GMP QUIET gmp)
    if(GMP_FOUND)
        # pkg_check_modules sets GMP_INCLUDE_DIRS (plural), extract to GMP_INCLUDE_DIR (singular)
        if(GMP_INCLUDE_DIRS)
            list(GET GMP_INCLUDE_DIRS 0 GMP_INCLUDE_DIR)
            message(STATUS "Found GMP via pkg-config: include=${GMP_INCLUDE_DIR}, libs=${GMP_LIBRARIES}")
        else()
            # pkg-config found GMP but didn't provide include dirs, unset to trigger manual search
            set(GMP_FOUND FALSE)
            message(STATUS "pkg-config found GMP but no include dirs, falling back to manual search")
        endif()
        
        # Find gmpxx library (C++ wrapper)
        if(GMP_FOUND)
            find_library(GMPXX_LIB NAMES gmpxx
                PATHS ${GMP_LIBRARY_DIRS}
                      /usr/lib /usr/lib/x86_64-linux-gnu
                      /usr/local/lib /opt/local/lib
                      /opt/homebrew/lib
            )
            if(GMPXX_LIB)
                list(APPEND GMP_LIBRARIES ${GMPXX_LIB})
            endif()
        endif()
    endif()
endif()

# Define library search paths (used by both pkg-config fallback and static linking)
set(GMP_LIB_SEARCH_PATHS
    /usr/lib
    /usr/lib/x86_64-linux-gnu
    /usr/local/lib
    /opt/local/lib
    /opt/homebrew/lib
)

# Use environment variables if set
if(DEFINED ENV{GMP_ROOT})
    list(INSERT GMP_LIB_SEARCH_PATHS 0 "$ENV{GMP_ROOT}/lib")
endif()
if(DEFINED ENV{GMP_LIBRARIES_DIR})
    list(INSERT GMP_LIB_SEARCH_PATHS 0 "$ENV{GMP_LIBRARIES_DIR}")
endif()

# Fallback: if GMP_FOUND is TRUE but GMP_INCLUDE_DIR is missing, search for it
if(GMP_FOUND AND (NOT GMP_INCLUDE_DIR OR GMP_INCLUDE_DIR STREQUAL ""))
    message(STATUS "GMP_FOUND is TRUE but GMP_INCLUDE_DIR missing, searching for header...")
    find_path(GMP_INCLUDE_DIR NAMES gmpxx.h
        PATHS
            /usr/include
            /usr/local/include
            /opt/local/include
            /opt/homebrew/include
    )
    if(GMP_INCLUDE_DIR)
        message(STATUS "Found GMP_INCLUDE_DIR via fallback search: ${GMP_INCLUDE_DIR}")
    endif()
endif()

# Manual GMP search if pkg-config failed
if(NOT GMP_FOUND)
    # Search paths
    set(GMP_SEARCH_PATHS
        /usr/include
        /usr/local/include
        /opt/local/include
        /opt/homebrew/include
    )
    
    # Use environment variables if set
    if(DEFINED ENV{GMP_ROOT})
        list(INSERT GMP_SEARCH_PATHS 0 "$ENV{GMP_ROOT}/include")
    endif()
    
    if(DEFINED ENV{GMP_INCLUDE_DIR})
        list(INSERT GMP_SEARCH_PATHS 0 "$ENV{GMP_INCLUDE_DIR}")
    endif()
    
    # Find headers
    find_path(GMP_INCLUDE_DIR
        NAMES gmpxx.h
        PATHS ${GMP_SEARCH_PATHS}
    )
    
    # Find libraries
    find_library(GMP_LIB
        NAMES gmp
        PATHS ${GMP_LIB_SEARCH_PATHS}
    )
    
    find_library(GMPXX_LIB
        NAMES gmpxx
        PATHS ${GMP_LIB_SEARCH_PATHS}
    )
    
    if(GMP_INCLUDE_DIR AND GMP_LIB AND GMPXX_LIB)
        set(GMP_FOUND TRUE)
        set(GMP_LIBRARIES ${GMP_LIB} ${GMPXX_LIB})
        message(STATUS "Found GMP manually: include=${GMP_INCLUDE_DIR}, libs=${GMP_LIBRARIES}")
    endif()
endif()

# Verify GMP was found
if(NOT GMP_FOUND)
    message(FATAL_ERROR "GMP not found. Install with: sudo apt-get install libgmp-dev (Linux) or brew install gmp (macOS)")
endif()

if(NOT GMP_INCLUDE_DIR OR GMP_INCLUDE_DIR STREQUAL "")
    message(FATAL_ERROR "GMP found but GMP_INCLUDE_DIR is not set. This is a configuration error.")
endif()

if(NOT GMP_LIBRARIES OR GMP_LIBRARIES STREQUAL "")
    message(FATAL_ERROR "GMP found but GMP_LIBRARIES is not set. This is a configuration error.")
endif()

message(STATUS "GMP configuration:")
message(STATUS "  Include directory: ${GMP_INCLUDE_DIR}")
message(STATUS "  Libraries: ${GMP_LIBRARIES}")

# Verify header exists
if(NOT EXISTS "${GMP_INCLUDE_DIR}/gmpxx.h")
    message(FATAL_ERROR "GMP header not found: ${GMP_INCLUDE_DIR}/gmpxx.h")
endif()

# Static linking support
if(ENABLE_STATIC_LINKING)
    message(STATUS "Static linking enabled")
    
    find_library(GMP_STATIC NAMES libgmp.a
        PATHS ${GMP_LIB_SEARCH_PATHS}
    )
    find_library(GMPXX_STATIC NAMES libgmpxx.a
        PATHS ${GMP_LIB_SEARCH_PATHS}
    )
    
    if(GMP_STATIC AND GMPXX_STATIC)
        set(GMP_LIBRARIES ${GMP_STATIC} ${GMPXX_STATIC})
        message(STATUS "Using static GMP: ${GMP_LIBRARIES}")
        
        if(UNIX AND NOT APPLE)
            # Linux: static link libgcc and libstdc++ but not glibc
            set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static-libgcc -static-libstdc++")
        elseif(APPLE)
            # macOS: clang doesn't support -static-libgcc/-static-libstdc++
            # Static GMP libraries are already linked, no additional flags needed
            message(STATUS "macOS: Using static GMP libraries (no additional linker flags needed)")
        endif()
    else()
        message(WARNING "Static GMP not found, using dynamic")
    endif()
endif()

# Build library
add_library(fracessa_lib
    src/fracessa.cpp
    src/findeq.cpp
    src/checkstab.cpp
)

set_target_properties(fracessa_lib PROPERTIES
    VERSION ${PROJECT_VERSION}
    SOVERSION ${PROJECT_VERSION_MAJOR}
)

# Include directories
target_include_directories(fracessa_lib
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/fracessa>
        $<BUILD_INTERFACE:${argparse_SOURCE_DIR}/include>
        $<BUILD_INTERFACE:${spdlog_SOURCE_DIR}/include>
        ${GMP_INCLUDE_DIR}
    PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/src
)

# Boost as system headers
target_include_directories(fracessa_lib SYSTEM PUBLIC
    ${boost_config_SOURCE_DIR}/include
    ${boost_predef_SOURCE_DIR}/include
    ${boost_assert_SOURCE_DIR}/include
    ${boost_static_assert_SOURCE_DIR}/include
    ${boost_throw_exception_SOURCE_DIR}/include
    ${boost_type_traits_SOURCE_DIR}/include
    ${boost_core_SOURCE_DIR}/include
    ${boost_multiprecision_SOURCE_DIR}/include
    ${boost_math_SOURCE_DIR}/include
    ${boost_integer_SOURCE_DIR}/include
    ${boost_algorithm_SOURCE_DIR}/include
)

target_compile_definitions(fracessa_lib PRIVATE
    BOOST_ALLOW_DEPRECATED_HEADERS
)

target_link_libraries(fracessa_lib PRIVATE
    ${GMP_LIBRARIES}
)

# Build executable (only on supported platforms)
if(UNIX OR WIN32 OR MSVC OR MSYS OR MINGW)
    add_executable(fracessa src/main.cpp)
    
    target_include_directories(fracessa PRIVATE
        ${argparse_SOURCE_DIR}/include
        ${GMP_INCLUDE_DIR}
    )
    
    target_include_directories(fracessa SYSTEM PRIVATE
        ${boost_math_SOURCE_DIR}/include
    )
    
    target_compile_definitions(fracessa PRIVATE
        BOOST_ALLOW_DEPRECATED_HEADERS
    )
    
    target_link_libraries(fracessa PRIVATE
        fracessa_lib
        ${GMP_LIBRARIES}
    )
    
    # Windows-specific linker flags
    if(WIN32 OR MSVC OR MSYS OR MINGW)
        set_target_properties(fracessa PROPERTIES
            LINK_FLAGS "-static-libgcc -static-libstdc++ -static"
        )
    endif()
    
    # Install
    install(TARGETS fracessa fracessa_lib
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

# Testing
if(BUILD_TESTING)
    enable_testing()
    
    FetchContent_Declare(googletest
        GIT_REPOSITORY https://github.com/google/googletest.git
        GIT_TAG v1.14.0
    )
    FetchContent_MakeAvailable(googletest)
    
    add_subdirectory(tests)
endif()