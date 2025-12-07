# Detailed Copositivity Algorithm Trace and Analysis

## Test Case: Simple 2x2 Matrix

Let's trace through a 2x2 matrix:
```
A = [2, 1]
    [1, 2]
```

This matrix is strictly copositive (all eigenvalues positive, symmetric positive definite).

## Trace Through copositivity.hpp Algorithm

### Step 1: Initial Call
```
isStrictlyCopositiveMemoized(A) where A is 2x2
- n = 2
- memo.clear()
- full_mask = {0, 1} (bits: 11)
- Call checkRecursive(A, {0,1}, 2)
```

### Step 2: checkRecursive(A, {0,1}, 2) - Full 2x2 Matrix
```
current_dim = mask.count() = 2

Not base case (dim > 1), so continue...

Recursive Step: Check all submatrices of size (2-1) = 1
  Loop i=0: mask.test(0) = true
    sub_mask = {0,1} with bit 0 reset = {1} (bits: 10)
    Call checkRecursive(A, {1}, 2)
    
  Loop i=1: mask.test(1) = true  
    sub_mask = {0,1} with bit 1 reset = {0} (bits: 01)
    Call checkRecursive(A, {0}, 2)
```

### Step 3: checkRecursive(A, {1}, 2) - 1x1 Submatrix
```
current_dim = 1
BASE CASE: 
  idx = mask.find_first() = 1
  result = A(1,1) > 0 = 2 > 0 = true
  memo[{1}] = 1
  return true
```

### Step 4: checkRecursive(A, {0}, 2) - 1x1 Submatrix
```
current_dim = 1
BASE CASE:
  idx = mask.find_first() = 0
  result = A(0,0) > 0 = 2 > 0 = true
  memo[{0}] = 1
  return true
```

### Step 5: Back to checkRecursive(A, {0,1}, 2) - Determinant Check
```
All (dim-1) submatrices passed, so continue...

Extract 2x2 submatrix:
  subMat = [2, 1]
           [1, 2]
  
  det = 2*2 - 1*1 = 4 - 1 = 3 > 0
  
  Condition: if (det < -EPSILON) → 3 < -1e-9 → FALSE
  So skip the adjugate check
  
  memo[{0,1}] = 1
  return true
```

**Result**: Returns true ✓ (correct for this matrix)

## Test Case: Non-Copositive 2x2 Matrix

Let's test with a matrix that should fail:
```
B = [1, -2]
    [-2, 1]
```

This matrix has negative determinant and should NOT be strictly copositive.

### Trace Through copositivity.hpp Algorithm

### Step 1: checkRecursive(B, {0,1}, 2)
```
current_dim = 2

Recursive Step:
  checkRecursive(B, {1}, 2): B(1,1) = 1 > 0 → true ✓
  checkRecursive(B, {0}, 2): B(0,0) = 1 > 0 → true ✓

Determinant Check:
  subMat = [1, -2]
           [-2, 1]
  det = 1*1 - (-2)*(-2) = 1 - 4 = -3 < 0
  
  Condition: if (det < -EPSILON) → -3 < -1e-9 → TRUE
  
  Check adjugate:
    adj = B.inverse() * det
    For 2x2: adj = [1, 2] * (-3) = [-3, -6]
                    [2, 1]          [-6, -3]
    
    Check if all positive:
      adj(0,0) = -3 <= EPSILON → all_positive = false
      
    Since all_positive = false, we DON'T return false here
    Continue...
  
  memo[{0,1}] = 1
  return true
```

**Result**: Returns true ✗ (WRONG! Should return false)

## Comparison with Existing Implementation

### Existing `is_strictly_copositive` in matrix.hpp

For matrix B = [[1, -2], [-2, 1]]:

```cpp
bitset64::iterate_all_supports(2, [&](const bitset64& support, unsigned) {
  // Check support {0}: 1x1 submatrix [1]
    det = 1 > 0 → skip (condition is det <= 0)
  
  // Check support {1}: 1x1 submatrix [1]
    det = 1 > 0 → skip
  
  // Check support {0,1}: 2x2 submatrix [[1,-2],[-2,1]]
    det = -3 <= 0 ✓
    adjugate = [[1,2],[2,1]] (all positive) ✓
    → result = false (NOT strictly copositive)
});
```

**Result**: Returns false ✓ (correct)

## Critical Issues Identified

### Issue 1: Wrong Adjugate Calculation

**Current code (line 58)**:
```cpp
Eigen::MatrixXd adj = subMat.inverse() * det;
```

**Problem**: This calculates `A^(-1) * det`, but the adjugate is `det * A^(-1)` which is the same mathematically, BUT:
- If det < 0, then `adj = inverse * det` will have negative values
- The code then checks if adj has all positive entries
- But if det is negative, the adjugate entries will be negative!

**What should be checked**: According to Hadeler 1983, we need to check if the **adjugate matrix itself** (not multiplied by det) has all positive entries when det ≤ 0.

**Correct calculation**: 
```cpp
Eigen::MatrixXd adj = subMat.adjoint();  // or compute adjugate directly
// adjugate = transpose of cofactor matrix
```

### Issue 2: Logic Error in Adjugate Check

The mathematical condition (Hadeler 1983) states:
- If det ≤ 0 AND adjugate has all positive entries → NOT strictly copositive

But the current code:
1. Calculates `adj = inverse * det` (which will be negative if det < 0)
2. Checks if `adj` has all positive entries
3. If yes → returns false

**Problem**: If det < 0, then `adj = inverse * det` will have negative entries, so the check `adj(r,c) > 0` will always fail, making the condition never trigger!

**Example with B = [[1,-2],[-2,1]]**:
- det = -3
- inverse = [[-1/3, -2/3], [-2/3, -1/3]]
- adj = inverse * (-3) = [[1, 2], [2, 1]] ✓ (this is correct!)
- But wait... let me recalculate:
  - B^(-1) = (1/det) * adjugate(B)
  - So adjugate(B) = B^(-1) * det
  - For B = [[1,-2],[-2,1]]:
    - det = -3
    - B^(-1) = [[-1/3, -2/3], [-2/3, -1/3]]
    - adjugate = B^(-1) * det = [[1, 2], [2, 1]] ✓

Actually, the calculation `inverse * det` IS the adjugate! So that part is correct.

But then checking if adjugate has all positive entries when det < 0... the adjugate itself should be checked, not multiplied by anything.

Wait, let me reconsider. The adjugate of a matrix A is defined as:
- adj(A) = det(A) * A^(-1)

So `A^(-1) * det(A)` = adjugate. That's correct.

But the issue is: when we check if adjugate has all positive entries, we're checking the adjugate matrix itself, which is correct.

However, in the example above, the adjugate [[1,2],[2,1]] DOES have all positive entries, so the condition should trigger and return false. But in my trace, I got `all_positive = false` because I incorrectly calculated the adjugate.

Let me recalculate properly:
- B = [[1, -2], [-2, 1]]
- det(B) = 1*1 - (-2)*(-2) = 1 - 4 = -3
- B^(-1) = (1/-3) * [[1, 2], [2, 1]] = [[-1/3, -2/3], [-2/3, -1/3]]
- adjugate = B^(-1) * det = [[-1/3, -2/3], [-2/3, -1/3]] * (-3) = [[1, 2], [2, 1]]

So adjugate = [[1, 2], [2, 1]] which has all positive entries!
- adj(0,0) = 1 > 0 ✓
- adj(0,1) = 2 > 0 ✓
- adj(1,0) = 2 > 0 ✓
- adj(1,1) = 1 > 0 ✓

So `all_positive = true`, and the code should return false. But in my trace I incorrectly said it would be false. Let me fix the trace.

Actually, I need to check what Eigen's `inverse()` returns and how the multiplication works.

## Corrected Analysis

### Issue: Adjugate Calculation Method

The code uses:
```cpp
Eigen::MatrixXd adj = subMat.inverse() * det;
```

This is mathematically correct (adjugate = inverse * determinant), BUT:
- If the matrix is singular (det = 0), `inverse()` will fail or give incorrect results
- Eigen's `inverse()` may have numerical issues

**Better approach**: Calculate adjugate directly using the definition:
- adjugate[i,j] = (-1)^(i+j) * det(minor(i,j))

### Issue: Missing det = 0 Case

As identified before, the code checks `if (det < -EPSILON)` but should check `if (det <= EPSILON)` to include det = 0.

## Summary of Issues

1. **Missing det = 0 case**: Should check `det <= EPSILON` not `det < -EPSILON`
2. **Potential numerical issues**: Using `inverse() * det` instead of computing adjugate directly
3. **Logic appears correct**: The recursive structure does check all principal submatrices
4. **Adjugate check logic**: The condition "if det < 0 and adjugate all positive → false" is correct according to Hadeler 1983

## Verification Needed

To fully verify, we should:
1. Test with a matrix where det = 0 and adjugate is all positive
2. Test with a matrix where det < 0 and adjugate is all positive (like the B example above)
3. Compare results directly with the existing `is_strictly_copositive` implementation

