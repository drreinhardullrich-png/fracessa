# GProf Profiling Summary - Matrices 27-33

## Execution Overview

All matrices were successfully profiled with `-pg` flag enabled. Results are stored in individual report files.

## Matrix Execution Times

| Matrix ID | Dimension | Execution Time | ESS Count |
|-----------|-----------|----------------|-----------|
| 27        | 18        | ~0.24s         | 1152      |
| 28        | 19        | ~0.06s         | 1444      |
| 29        | 19        | ~0.03s         | 19        |
| 30        | 20        | ~1.40s         | 2560      |
| 31        | 21        | ~2.97s         | 4410      |
| 32        | 22        | ~16.52s        | 5632      |
| 33        | 23        | ~1.12s         | 2507      |

## Top Performance Bottlenecks

### 1. Vector Iterator Operations (Highest Impact)
- **Total Time**: ~13.87s across all matrices
- **Most affected**: Matrix 32 (71.49% of time, 11.81s)
- **Location**: `__gnu_cxx::__normal_iterator<bitset64*>` operations
- **Impact**: This is the single largest time consumer, particularly in larger matrices
- **Recommendation**: Consider optimizing vector iteration patterns, possibly using raw pointers or more efficient iteration methods

### 2. search_one_support Function
- **Total Time**: ~6.81s across all matrices
- **Most affected**: Matrix 32 (24.82%, 4.10s), Matrix 31 (40.74%, 1.21s)
- **Location**: `fracessa::search_one_support()`
- **Impact**: Core search logic is the second largest bottleneck
- **Recommendation**: Profile this function specifically to identify hot paths within it

### 3. Boost Multiprecision Operations
- **Total Time**: ~0.78s across all matrices
- **Components**:
  - `boost::multiprecision::cpp_int` operations: ~0.25s
  - `boost::safe_numerics::safe_base` operations: ~0.18s
  - Rational arithmetic: ~0.35s
- **Impact**: Moderate overhead from arbitrary precision arithmetic
- **Note**: This is expected when using arbitrary precision for overflow handling

### 4. Matrix Operations
- **Total Time**: ~0.17s
- **Components**:
  - `is_positive_definite_rational`: ~0.06s (Matrix 31)
  - Matrix conversions: ~0.11s
- **Impact**: Relatively minor compared to iteration overhead

## Key Findings

1. **Vector Iteration is the Primary Bottleneck**: The iterator overhead (13.87s) exceeds the actual computation time in `search_one_support` (6.81s). This suggests that the iteration pattern itself is inefficient.

2. **Matrix 32 is Exceptionally Slow**: At 16.52s, Matrix 32 takes significantly longer than others. The vector iterator operations consume 71.49% of its time, indicating a scalability issue.

3. **search_one_support Performance**: While this function is called frequently, the per-call overhead is relatively low. The total time is high due to the large number of calls (hundreds of thousands to millions).

4. **Boost Overhead is Acceptable**: The arbitrary precision arithmetic overhead (~0.78s total) is reasonable given the safety it provides for overflow handling.

## Recommendations

1. **Optimize Vector Iteration**:
   - Consider using raw pointers or range-based for loops where possible
   - Profile the specific iterator operations to identify the exact bottleneck
   - Consider using `std::array` for fixed-size collections if applicable

2. **Profile search_one_support Internally**:
   - Use more granular profiling (e.g., callgrind) to identify hot paths within this function
   - Look for opportunities to reduce the number of calls or optimize the algorithm

3. **Matrix 32 Specific Investigation**:
   - Matrix 32 shows unusual behavior with 71.49% time in iterators
   - Investigate why this matrix triggers so many iterator operations
   - Consider matrix-specific optimizations if patterns are identified

4. **Consider Algorithmic Optimizations**:
   - The high number of `search_one_support` calls suggests potential for pruning or early termination
   - Review the support search strategy for optimization opportunities

## Report Files

- Individual reports: `gprof_reports/report_*.txt` (one per matrix)
- Raw profiling data: `gprof_reports/gmon.*` (one per matrix)
- This summary: `gprof_reports/SUMMARY.md`

## Next Steps

1. Use callgrind/kcachegrind for more detailed call graph analysis
2. Focus profiling efforts on Matrix 32 to understand the iterator overhead
3. Consider micro-benchmarking vector iteration patterns
4. Review `search_one_support` implementation for optimization opportunities

