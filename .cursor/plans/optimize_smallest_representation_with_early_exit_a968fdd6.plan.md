# Add is_smallest_representation Method to bitset64

## Problem

Currently, when checking if a support is canonical (`support.smallest_representation(dimension_) == support`), we compute all rotations to find the minimum, then compare. This is expensive when called billions of times. We need a dedicated method that does the rotations with early exit when it finds ANY rotation smaller than the original.

## Solution

Create a new `is_smallest_representation` method in `bitset64` that:
1. Performs circular rotations (same logic as `smallest_representation`)
2. Returns `false` immediately when it finds ANY rotation smaller than the original (early exit)
3. Returns `true` only if all rotations are >= original (canonical form)

This avoids computing all rotations for non-canonical supports, which are the majority.

## Implementation

### Step 1: Add is_smallest_representation method to bitset64.hpp

Add a new method that checks if the bitset is in its smallest representation, with early exit optimization.

**File**: `[fracessa/include/fracessa/bitset64.hpp](fracessa/include/fracessa/bitset64.hpp)`

Add after `smallest_representation` method (around line 158):

```cpp
// Check if this bitset is in its smallest representation (canonical form)
// Uses early exit: returns false immediately if any rotation is smaller than original
// Much faster than smallest_representation() == *this for canonical checks
FORCE_INLINE bool is_smallest_representation(unsigned nbits) const noexcept {
    if (nbits == 0 || bits_ == 0ULL) return true;
    uint64_t mask = (1ULL << nbits) - 1ULL;
    uint64_t original = bits_ & mask;
    
    // Early exit: 0 is always canonical (smallest possible)
    if (original == 0ULL) return true;
    
    uint64_t current = original;
    for (unsigned i = 1; i < nbits; i++) {
        // Circular rotate right by 1 (same logic as smallest_representation)
        uint64_t lo = current << (nbits - 1);
        uint64_t hi = current >> 1;
        current = (hi | lo) & mask;
        
        // Early exit: if we find ANY rotation smaller than original, not canonical
        if (current < original) {
            return false;
        }
    }
    // All rotations are >= original, so it's canonical
    return true;
}
```

### Step 2: Update supports.hpp to use is_smallest_representation

Replace the `smallest_representation() == support` check with the faster `is_smallest_representation()` method.

**File**: `[fracessa/include/fracessa/supports.hpp](fracessa/include/fracessa/supports.hpp)`

Change lines 49-53 from:

```cpp
if (is_coprime[current_index]) {
    // Only add if it's the smallest representation
    if (support.smallest_representation(dimension_) == support) {
        supports_[current_index].push_back(support);
    }
}
```

To:

```cpp
if (is_coprime[current_index]) {
    // Only add if it's the smallest representation (canonical form)
    if (support.is_smallest_representation(dimension_)) {
        supports_[current_index].push_back(support);
    }
}
```

## Expected Performance Gain

- **Early exit**: When a support is not canonical, the loop exits as soon as it finds a smaller rotation (often after 1-2 rotations instead of `dimension_` rotations)
- **No overhead for canonical supports**: Still checks all rotations, but this is necessary to confirm canonicity
- **Overall**: Should reduce canonical check time by ~50-80% for non-canonical supports, which are the majority

## Testing

- Verify that the same supports are generated (correctness check)
- Measure initialization time improvement
- Ensure `smallest_representation()` still works correctly for other use cases
