# Positive Definite Check Comparison - Test Results

## Summary

This document summarizes the differences between the custom `is_positive_definite_double` implementation and Eigen's `LLT` method, using `is_positive_definite_rational` as the ground truth.

## Key Findings

### 1. **False Negatives in Custom Method**

The custom `is_positive_definite_double` method incorrectly rejects **3 out of 13 test cases** that are actually positive definite:

- **Test 3**: Matrix with eigenvalue 0.001 → Custom says FALSE, Rational/Eigen say TRUE
- **Test 5**: Matrix with eigenvalue 0.009 → Custom says FALSE, Rational/Eigen say TRUE  
- **Test 10**: Matrix with eigenvalue 1e-10 → Custom says FALSE, Rational/Eigen say TRUE

**Root Cause**: The hardcoded tolerance of `0.01` in line 321 of `matrix.hpp`:
```cpp
if (x < 0.01)
    return false;
```

This rejects any matrix where the Cholesky diagonal value is less than 0.01, even if the matrix is mathematically positive definite.

### 2. **Eigen LLT is More Accurate**

Eigen's `LLT` method:
- **Agrees with rational ground truth in all 13 test cases**
- Uses adaptive tolerance based on machine epsilon and matrix condition number
- Correctly handles borderline cases with very small eigenvalues

### 3. **No False Positives Found**

Neither method incorrectly accepts non-positive-definite matrices in our test cases. Both correctly reject:
- Matrices with negative eigenvalues (Test 8)
- Indefinite matrices (Test 9)

## Detailed Test Results

| Test | Description | Custom | Eigen LLT | Rational | Agreement |
|------|-------------|--------|-----------|----------|-----------|
| 1 | Clearly PD (3x3) | TRUE | TRUE | TRUE | ✓ All agree |
| 2 | Identity (4x4) | TRUE | TRUE | TRUE | ✓ All agree |
| 3 | Small eigenvalue (0.001) | **FALSE** | TRUE | TRUE | ❌ Custom wrong |
| 4 | Eigenvalue 0.011 | TRUE | TRUE | TRUE | ✓ All agree |
| 5 | Eigenvalue 0.009 | **FALSE** | TRUE | TRUE | ❌ Custom wrong |
| 6 | Eigenvalue 0.0100001 | TRUE | TRUE | TRUE | ✓ All agree |
| 7 | Eigenvalue exactly 0.01 | TRUE | TRUE | TRUE | ✓ All agree |
| 8 | Negative eigenvalue | FALSE | FALSE | FALSE | ✓ All agree |
| 9 | Indefinite matrix | FALSE | FALSE | FALSE | ✓ All agree |
| 10 | Very small (1e-10) | **FALSE** | TRUE | TRUE | ❌ Custom wrong |
| 11 | Small values (~0.1) | TRUE | TRUE | TRUE | ✓ All agree |
| 12 | Borderline (~0.02) | TRUE | TRUE | TRUE | ✓ All agree |
| 13 | Random PD (5x5) | TRUE | TRUE | TRUE | ✓ All agree |

## Tolerance Comparison

### Custom Method
- **Tolerance**: Fixed `0.01` (hardcoded)
- **Problem**: Too large for matrices with small eigenvalues
- **Impact**: Rejects valid positive definite matrices

### Eigen LLT
- **Tolerance**: Adaptive, based on:
  - Machine epsilon (`std::numeric_limits<double>::epsilon()`)
  - Maximum diagonal element of the matrix
  - Matrix condition number
- **Formula**: Approximately `epsilon * max_diagonal * condition_number`
- **Result**: Correctly handles all test cases

## Recommendations

1. **Replace custom implementation with Eigen LLT**: Eigen's method is more accurate and handles edge cases correctly.

2. **If keeping custom implementation**: 
   - Replace hardcoded `0.01` tolerance with adaptive tolerance
   - Use: `epsilon * max_diagonal` or similar adaptive approach
   - Consider: `std::numeric_limits<double>::epsilon() * A.diagonal().cwiseAbs().maxCoeff()`

3. **Current behavior impact**: 
   - The custom method may incorrectly classify some ESS candidates as non-ESS
   - This could lead to false negatives in the stability check
   - However, the rational check will catch these later, so the impact is performance-related (slower fallback to rational check)

## Code Locations

- Custom implementation: `fracessa/include/fracessa/matrix.hpp:304-326`
- Old Eigen LLT code (commented): `fracessa/src/checkstab.cpp:63-71`
- Rational ground truth: `fracessa/include/fracessa/matrix.hpp:233-255`

