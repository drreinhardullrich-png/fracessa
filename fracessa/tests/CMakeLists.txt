# Test executables
add_executable(test_rational64 test_rational64.cpp)
target_link_libraries(test_rational64 PRIVATE fracessa_lib gtest_main ${GMP_LIBRARIES})
target_include_directories(test_rational64 PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../include)

add_executable(test_matrix test_matrix.cpp)
target_link_libraries(test_matrix PRIVATE fracessa_lib gtest_main ${GMP_LIBRARIES})
target_include_directories(test_matrix PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../include)

add_executable(test_linear_solver test_linear_solver.cpp)
target_link_libraries(test_linear_solver PRIVATE fracessa_lib gtest_main ${GMP_LIBRARIES})
target_include_directories(test_linear_solver PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../include)

add_executable(test_bareiss_lu test_bareiss_lu.cpp)
target_link_libraries(test_bareiss_lu PRIVATE fracessa_lib gtest_main ${GMP_LIBRARIES})
target_include_directories(test_bareiss_lu PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../include)

add_executable(test_copositivity test_copositivity.cpp)
target_link_libraries(test_copositivity PRIVATE fracessa_lib gtest_main ${GMP_LIBRARIES})
target_include_directories(test_copositivity PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../include)

add_executable(test_bitset64 test_bitset64.cpp)
target_link_libraries(test_bitset64 PRIVATE fracessa_lib gtest_main)
target_include_directories(test_bitset64 PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../include)

add_executable(test_integration test_integration.cpp)
target_link_libraries(test_integration PRIVATE fracessa_lib gtest_main ${GMP_LIBRARIES})
target_include_directories(test_integration PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../include)
# Copy test data to build directory
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/test_data DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

# Register tests with CTest
add_test(NAME Rational64Tests COMMAND test_rational64)
add_test(NAME MatrixTests COMMAND test_matrix)
add_test(NAME LinearSolverTests COMMAND test_linear_solver)
add_test(NAME BareissLUTests COMMAND test_bareiss_lu)
add_test(NAME CopositivityTests COMMAND test_copositivity)
add_test(NAME Bitset64Tests COMMAND test_bitset64)
add_test(NAME IntegrationTests COMMAND test_integration)

