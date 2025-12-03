cmake_minimum_required(VERSION 3.10.2)

project(fracessa
    VERSION 1.0.0
    DESCRIPTION "FRACESSA - Fractional ESS Analyzer - A solver for Standard Quadratic Problems"
    LANGUAGES CXX C
)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

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

# Python bindings (optional)
# find_package(pybind11 CONFIG QUIET)

# Build the core library
add_library(fracessa_lib
    src/ess_finder.cpp
    src/ref.cpp
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
    PRIVATE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>
)

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

# Python bindings (commented out - using Python wrapper instead)
# if(pybind11_FOUND)
#     pybind11_add_module(fracessa_py src/bindings.cpp)
#     target_link_libraries(fracessa_py PRIVATE fracessa_lib)
#     set_target_properties(fracessa_py PROPERTIES
#         OUTPUT_NAME "fracessa"
#         LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}
#     )
# endif()

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
