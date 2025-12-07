# Copositivity Algorithm Analysis

## Mathematical Definition (from Bomze 1992, Hadeler 1983)

A symmetric matrix B is **strictly copositive** if:
- x^T B x > 0 for all x ≥ 0, x ≠ 0

**Equivalent condition** (Hadeler 1983, Theorem):
- For ALL principal submatrices, if det ≤ 0, then the adjugate must NOT have all positive entries
- In other words: If det ≤ 0 AND adjugate has all positive entries → NOT strictly copositive

## What copositivity.hpp Does

### Algorithm Structure:

1. **Base Case (1x1)**:
   - Checks: `A[i,i] > 0`
   - Returns true if diagonal element is positive
   - **This is CORRECT** - a 1x1 matrix [a] is strictly copositive iff a > 0

2. **Recursive Step**:
   - For current mask (set of indices), recursively checks ALL submatrices of size (current_dim - 1)
   - If ANY submatrix fails → returns false
   - **This checks all proper submatrices recursively**

3. **Determinant Check**:
   - If `det < -EPSILON` AND `adjugate has all positive entries` → returns false
   - **This matches the mathematical condition**

4. **Final Result**:
   - If all checks pass → returns true

## Detailed Trace Through Algorithm

### Example: 3x3 Matrix

**Recursion tree:**
```
checkRecursive(A, {0,1,2}, 3)  // Full 3x3 matrix
├─ checkRecursive(A, {1,2}, 3)  // Remove index 0 → 2x2 submatrix
│  ├─ checkRecursive(A, {2}, 3)  // Remove index 1 → 1x1: A[2,2] > 0?
│  └─ checkRecursive(A, {1}, 3)  // Remove index 2 → 1x1: A[1,1] > 0?
│  └─ Check 2x2: if det < 0 and adj all positive → false
├─ checkRecursive(A, {0,2}, 3)  // Remove index 1 → 2x2 submatrix
│  ├─ checkRecursive(A, {2}, 3)  // Remove index 0 → 1x1: A[2,2] > 0?
│  └─ checkRecursive(A, {0}, 3)  // Remove index 2 → 1x1: A[0,0] > 0?
│  └─ Check 2x2: if det < 0 and adj all positive → false
├─ checkRecursive(A, {0,1}, 3)  // Remove index 2 → 2x2 submatrix
│  ├─ checkRecursive(A, {1}, 3)  // Remove index 0 → 1x1: A[1,1] > 0?
│  └─ checkRecursive(A, {0}, 3)  // Remove index 2 → 1x1: A[0,0] > 0?
│  └─ Check 2x2: if det < 0 and adj all positive → false
└─ Check 3x3: if det < 0 and adj all positive → false
```

**This DOES check all principal submatrices:**
- All 1x1: {0}, {1}, {2} ✓
- All 2x2: {0,1}, {0,2}, {1,2} ✓
- Full 3x3: {0,1,2} ✓

## Comparison with Existing Implementation

### Existing `is_strictly_copositive` in matrix.hpp (lines 438-465):

```cpp
bitset64::iterate_all_supports(n_rows, [&](const bitset64& support, unsigned) {
    // For each principal submatrix:
    if (determinant(subA) <= 0 && greater_zero(adjugate(subA))) {
        result = false;  // NOT strictly copositive
    }
});
```

**Checks**: For ALL non-empty principal submatrices, if det ≤ 0 AND adjugate all positive → NOT strictly copositive.

### copositivity.hpp Algorithm:

**Checks**: 
- Recursively verifies all submatrices of size (dim-1) are strictly copositive
- Then checks the current submatrix: if det < 0 AND adjugate all positive → NOT strictly copositive

## Critical Issues Found

### Issue 1: Missing det = 0 Case

**Current code**: `if (det < -EPSILON)`
**Should be**: `if (det <= EPSILON)` (to include det = 0)

**Mathematical definition**: "if det ≤ 0" (includes det = 0)

**Impact**: The algorithm might miss cases where det = 0 and adjugate is all positive, which should return false.

### Issue 2: Logic Completeness

The recursive structure DOES check all principal submatrices (verified by trace above), but:

**Question**: Does the recursive check of (dim-1) submatrices ensure the current matrix is copositive?

**Answer**: The recursion ensures all proper submatrices are copositive, then checks the current matrix's determinant condition. This is correct IF the determinant condition is the only requirement.

However, according to Hadeler 1983, the condition is:
- For ALL principal submatrices: if det ≤ 0, then adjugate must NOT be all positive

The recursive algorithm checks:
1. All (dim-1) submatrices are copositive (recursively)
2. Current matrix: if det < 0 and adj all positive → false

**This seems correct**, but the det = 0 case is missing.

### Issue 3: Base Case Logic

**Current**: For 1x1, checks if diagonal > 0
**Correct**: A 1x1 matrix [a] is strictly copositive iff a > 0 ✓

## Verification Against Mathematical Definition

According to Hadeler 1983 (Theorem 3.1):
> A symmetric matrix B is strictly copositive if and only if for all principal submatrices C:
> - If det(C) ≤ 0, then adj(C) does NOT have all positive entries

**What the algorithm checks**:
- ✓ All principal submatrices (through recursion)
- ✗ Only checks det < 0, missing det = 0
- ✓ Checks if adjugate has all positive entries
- ✓ Returns false if det < 0 AND adj all positive

## Conclusion

### The algorithm structure is CORRECT:
- It checks all principal submatrices through recursion
- It applies the determinant/adjugate condition

### But there is a BUG:
- **Missing det = 0 case**: Should check `det <= EPSILON` not just `det < -EPSILON`
- This could cause false positives (returning true when should return false)

### Recommendation:
Change line 56 from:
```cpp
if (det < -EPSILON) {
```
to:
```cpp
if (det <= EPSILON) {  // Include det = 0 case
```

This will make it match the mathematical definition: "if det ≤ 0"
