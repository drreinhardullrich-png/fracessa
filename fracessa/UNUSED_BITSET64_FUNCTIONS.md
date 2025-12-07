# Analysis Report: Unused Functions in bitset64.hpp

## Summary
This report identifies functions in `bitset64.hpp` that are never called in the codebase. **Operators are excluded from this analysis** as they may be used implicitly and are kept for API completeness.

## Unused Functions

The following functions are **NOT USED** anywhere in the codebase (except in their definitions):

### 1. `from_mask()` (Line 75)
```cpp
static FORCE_INLINE bitset64 from_mask(unsigned nbits, uint64_t mask) noexcept
```
- **Status**: Unused
- **Purpose**: Construct from raw mask with nbits masking
- **Recommendation**: Can be removed if not needed for future API

### 2. `flip()` (Line 107)
```cpp
FORCE_INLINE void flip(unsigned pos) noexcept
```
- **Status**: Unused
- **Purpose**: Toggle a single bit
- **Recommendation**: Can be removed if not needed

### 3. `inplace_and()` (Line 117)
```cpp
FORCE_INLINE void inplace_and(const bitset64 &o) noexcept
```
- **Status**: Unused
- **Purpose**: In-place bitwise AND operation
- **Recommendation**: Can be removed (use `operator&` instead)

### 4. `inplace_or()` (Line 118)
```cpp
FORCE_INLINE void inplace_or(const bitset64 &o) noexcept
```
- **Status**: Unused
- **Purpose**: In-place bitwise OR operation
- **Recommendation**: Can be removed (use `operator|` instead)

### 5. `inplace_xor()` (Line 119)
```cpp
FORCE_INLINE void inplace_xor(const bitset64 &o) noexcept
```
- **Status**: Unused
- **Purpose**: In-place bitwise XOR operation
- **Recommendation**: Can be removed (use `operator^` instead)

### 6. `inplace_subtract()` (Line 121)
```cpp
FORCE_INLINE void inplace_subtract(const bitset64 &o) noexcept
```
- **Status**: Unused
- **Purpose**: In-place set difference operation
- **Recommendation**: Can be removed (use `subtract()` which returns a new object)

### 7. `complement()` (Line 165)
```cpp
FORCE_INLINE bitset64 complement(unsigned nbits) const noexcept
```
- **Status**: Unused
- **Purpose**: Return complement with nbits masking
- **Recommendation**: Can be removed if not needed

### 8. `any()` (Line 170)
```cpp
FORCE_INLINE bool any() const noexcept
```
- **Status**: Unused
- **Purpose**: Check if any bit is set
- **Recommendation**: Can be removed (use `!none()` or check `count() > 0`)

### 9. `none()` (Line 171)
```cpp
FORCE_INLINE bool none() const noexcept
```
- **Status**: Unused
- **Purpose**: Check if no bits are set
- **Recommendation**: Can be removed (use `count() == 0` or check `bits_ == 0ULL`)

### 10. `for_each_set()` (Line 202)
```cpp
template<typename F>
FORCE_INLINE void for_each_set(F &&f) const
```
- **Status**: Unused
- **Purpose**: Iterate over all set bits with callback
- **Recommendation**: Can be removed (use `find_first()`/`find_next()` pattern instead)

### 11. `shl()` (Line 215)
```cpp
FORCE_INLINE bitset64 shl(unsigned s, unsigned nbits) const noexcept
```
- **Status**: Unused
- **Purpose**: Logical left shift
- **Recommendation**: Can be removed if not needed

### 12. `shr()` (Line 224)
```cpp
FORCE_INLINE bitset64 shr(unsigned s, unsigned nbits) const noexcept
```
- **Status**: Unused
- **Purpose**: Logical right shift
- **Recommendation**: Can be removed if not needed

### 13. `rot_l()` (Line 245)
```cpp
FORCE_INLINE bitset64 rot_l(unsigned s, unsigned nbits) const noexcept
```
- **Status**: Unused
- **Purpose**: Circular rotate left
- **Note**: Only `rot_r()` is used in the codebase
- **Recommendation**: Can be removed (use `rot_r(nbits - s, nbits)` if needed)

### 14. `to_uint64(unsigned nbits)` (Line 254)
```cpp
FORCE_INLINE uint64_t to_uint64(unsigned nbits) const noexcept
```
- **Status**: Unused (overloaded version)
- **Purpose**: Return uint64 with nbits masking
- **Note**: Only the parameterless `to_uint64()` is used
- **Recommendation**: Can be removed if not needed

## Functions That ARE Used

The following functions are actively used in the codebase:

- **Constructors/Destructors**: Default, copy, move constructors and assignment operators
- **`iterate_all_supports()`**: Used in `fracessa.cpp` and `matrix.hpp`
- **`set_all()`**: Used in `copositivity.hpp` and `checkstab.cpp`
- **`reset()`**: Used in `copositivity.hpp` and `checkstab.cpp`
- **`set()`**: Used in `findeq.cpp` and internally in `lowest_set_bit()`
- **`test()`**: Used extensively throughout the codebase
- **`subtract()`**: Used in `checkstab.cpp`
- **`lowest_set_bit()`**: Used in `checkstab.cpp`
- **`smallest_representation()`**: Used in `fracessa.cpp`
- **`is_subset_of()`**: Used in `fracessa.cpp`
- **`count()`**: Used extensively
- **`find_first()`**: Used in `copositivity.hpp`, `checkstab.cpp`, `matrix.hpp`
- **`find_next()`**: Used in `copositivity.hpp`, `matrix.hpp`
- **`rot_r()`**: Used in `fracessa.cpp`, `findeq.cpp`
- **`to_uint64()`**: Used in `fracessa.cpp` (parameterless version)
- **`to_string()`**: Used in `fracessa.cpp`, `candidate.hpp`, `checkstab.cpp`
- **`to_bitstring()`**: Used in `checkstab.cpp`
- **`hash()`**: Used in `copositivity.hpp` (via `bitset64_hash`)
- **`set_mask()`**: Used internally in `iterate_all_supports()`

## Operators (KEPT - Not Analyzed)

The following operators are kept regardless of usage, as they may be used implicitly:

- `operator&()` - Bitwise AND
- `operator|()` - Bitwise OR
- `operator^()` - Bitwise XOR
- `operator~()` - Bitwise NOT
- `operator==()` - Equality (used by `std::unordered_map`)
- `operator!=()` - Inequality

## Potential Bug Found

**Location**: `checkstab.cpp:136`
```cpp
kay_vee[v] = kay_vee[v-1].subtract(iv, dimension);
```

**Issue**: The `subtract()` function only takes one parameter (`const bitset64 &o`), but it's being called with two parameters (`iv, dimension`).

**Fix**: Should be:
```cpp
kay_vee[v] = kay_vee[v-1].subtract(iv);
```

## Recommendations

1. **Safe to Remove**: All 14 unused functions listed above can be safely removed if they're not needed for future API compatibility.

2. **Keep for API Completeness**: Consider keeping some functions if they provide a more intuitive API:
   - `any()` / `none()` - More readable than `count() > 0` / `count() == 0`
   - `for_each_set()` - Convenient iteration pattern
   - `rot_l()` - Completeness with `rot_r()`

3. **Fix Bug**: Fix the `subtract()` call in `checkstab.cpp:136` to use only one parameter.

## Total Unused Functions: 14

