// bitset_perf_test.cpp

#include <iostream>
#include <chrono>
#include <iomanip>
#include <cstdint>
#include <cstdlib>
#include <ctime>

// Force-inline hint
#if defined(_MSC_VER)
#  define FORCE_INLINE __forceinline
#else
#  define FORCE_INLINE __attribute__((always_inline)) inline
#endif

// ==================== YOUR ORIGINAL CLASS ====================

class bitset64_original {
private:
    uint64_t bits_;
    
public:
    FORCE_INLINE bitset64_original() noexcept : bits_(0ULL) {}
    FORCE_INLINE bitset64_original(uint64_t bits) noexcept : bits_(bits) {}
    
    FORCE_INLINE void set(unsigned pos) noexcept { bits_ |= (1ULL << pos); }
    FORCE_INLINE void reset(unsigned pos) noexcept { bits_ &= ~(1ULL << pos); }
    FORCE_INLINE bool test(unsigned pos) const noexcept { return (bits_ >> pos) & 1ULL; }
    FORCE_INLINE unsigned count() const noexcept { 
        #ifdef _MSC_VER
            return (unsigned)__popcnt64(bits_);
        #else
            return (unsigned)__builtin_popcountll(bits_);
        #endif
    }
    
    // Iteration with for_each pattern
    template<typename F>
    FORCE_INLINE bool for_each_set_bit(F&& callback) const noexcept {
        uint64_t temp = bits_;
        while (temp) {
            unsigned pos;
            #ifdef _MSC_VER
                unsigned long index;
                _BitScanForward64(&index, temp);
                pos = (unsigned)index;
            #else
                pos = (unsigned)__builtin_ctzll(temp);
            #endif
            if (!callback(pos)) return false;
            temp &= temp - 1;
        }
        return true;
    }
    
    FORCE_INLINE uint64_t raw() const noexcept { return bits_; }
};

// ==================== NAMESPACE VERSION ====================

namespace bitset64_ns {
    using bits_t = uint64_t;
    
    inline void set(bits_t& bits, unsigned pos) noexcept { 
        bits |= (1ULL << pos); 
    }
    
    inline void reset(bits_t& bits, unsigned pos) noexcept { 
        bits &= ~(1ULL << pos); 
    }
    
    inline bool test(bits_t bits, unsigned pos) noexcept { 
        return (bits >> pos) & 1ULL; 
    }
    
    inline unsigned count(bits_t bits) noexcept { 
        #ifdef _MSC_VER
            return (unsigned)__popcnt64(bits);
        #else
            return (unsigned)__builtin_popcountll(bits);
        #endif
    }
    
    template<typename F>
    inline bool for_each_set_bit(bits_t bits, F&& callback) noexcept {
        while (bits) {
            unsigned pos;
            #ifdef _MSC_VER
                unsigned long index;
                _BitScanForward64(&index, bits);
                pos = (unsigned)index;
            #else
                pos = (unsigned)__builtin_ctzll(bits);
            #endif
            if (!callback(pos)) return false;
            bits &= bits - 1;
        }
        return true;
    }
}

// ==================== RAW UINT64_T ====================

// Inline functions for raw uint64_t operations (no namespace wrapper)
FORCE_INLINE void bitset_raw_set(uint64_t& bits, unsigned pos) noexcept { 
    bits |= (1ULL << pos); 
}

FORCE_INLINE void bitset_raw_reset(uint64_t& bits, unsigned pos) noexcept { 
    bits &= ~(1ULL << pos); 
}

FORCE_INLINE bool bitset_raw_test(uint64_t bits, unsigned pos) noexcept { 
    return (bits >> pos) & 1ULL; 
}

FORCE_INLINE unsigned bitset_raw_count(uint64_t bits) noexcept { 
    #ifdef _MSC_VER
        return (unsigned)__popcnt64(bits);
    #else
        return (unsigned)__builtin_popcountll(bits);
    #endif
}

template<typename F>
FORCE_INLINE bool bitset_raw_for_each_set_bit(uint64_t bits, F&& callback) noexcept {
    while (bits) {
        unsigned pos;
        #ifdef _MSC_VER
            unsigned long index;
            _BitScanForward64(&index, bits);
            pos = (unsigned)index;
        #else
            pos = (unsigned)__builtin_ctzll(bits);
        #endif
        if (!callback(pos)) return false;
        bits &= bits - 1;
    }
    return true;
}

// ==================== MACROS VERSION ====================

#define BS_SET(bits, pos) ((bits) |= (1ULL << (pos)))
#define BS_RESET(bits, pos) ((bits) &= ~(1ULL << (pos)))
#define BS_TEST(bits, pos) (((bits) >> (pos)) & 1ULL)
#define BS_COUNT(bits) __builtin_popcountll(bits)
#define BS_CTZ(bits) __builtin_ctzll(bits)
#define BS_FOREACH(bits, pos_var, code) \
    for (uint64_t _bs_tmp = (bits); _bs_tmp; _bs_tmp &= _bs_tmp - 1) { \
        unsigned pos_var = BS_CTZ(_bs_tmp); \
        code \
    }

// ==================== BENCHMARK HELPERS ====================

class Timer {
    std::chrono::high_resolution_clock::time_point start_;
    
public:
    Timer() : start_(std::chrono::high_resolution_clock::now()) {}
    
    double elapsed_ms() const {
        auto end = std::chrono::high_resolution_clock::now();
        return std::chrono::duration<double, std::milli>(end - start_).count();
    }
    
    double elapsed_ns() const {
        auto end = std::chrono::high_resolution_clock::now();
        return std::chrono::duration<double, std::nano>(end - start_).count();
    }
};

void print_result(const std::string& name, double time_ns, double ops_per_ns, int width = 30) {
    std::cout << std::left << std::setw(width) << name 
              << ": " << std::fixed << std::setprecision(2) << time_ns << " ns"
              << " (" << ops_per_ns << " ops/ns)" << std::endl;
}

// ==================== TEST CASES ====================

// Test 1: Basic set/test/reset operations
void test_basic_operations() {
    std::cout << "\n=== Basic Set/Test/Reset Operations ===\n";
    const int ITERATIONS = 100000000;
    volatile bool dummy = false; // Prevent optimization
    
    // Original class
    {
        Timer timer;
        bitset64_original bs;
        for (int i = 0; i < ITERATIONS; ++i) {
            bs.set(i & 63);
            dummy = bs.test(i & 63);
            bs.reset(i & 63);
        }
        double time_ns = timer.elapsed_ns() / ITERATIONS;
        print_result("Class version", time_ns, 3.0 / time_ns);
    }
    
    // Namespace version
    {
        Timer timer;
        bitset64_ns::bits_t bs = 0;
        for (int i = 0; i < ITERATIONS; ++i) {
            bitset64_ns::set(bs, i & 63);
            dummy = bitset64_ns::test(bs, i & 63);
            bitset64_ns::reset(bs, i & 63);
        }
        double time_ns = timer.elapsed_ns() / ITERATIONS;
        print_result("Namespace version", time_ns, 3.0 / time_ns);
    }
    
    // Raw inline functions
    {
        Timer timer;
        uint64_t bs = 0;
        for (int i = 0; i < ITERATIONS; ++i) {
            bitset_raw_set(bs, i & 63);
            dummy = bitset_raw_test(bs, i & 63);
            bitset_raw_reset(bs, i & 63);
        }
        double time_ns = timer.elapsed_ns() / ITERATIONS;
        print_result("Raw inline functions", time_ns, 3.0 / time_ns);
    }
    
    // Raw uint64_t
    {
        Timer timer;
        uint64_t bs = 0;
        for (int i = 0; i < ITERATIONS; ++i) {
            bs |= (1ULL << (i & 63));
            dummy = (bs >> (i & 63)) & 1ULL;
            bs &= ~(1ULL << (i & 63));
        }
        double time_ns = timer.elapsed_ns() / ITERATIONS;
        print_result("Raw uint64_t", time_ns, 3.0 / time_ns);
    }
    
    // Macro version
    {
        Timer timer;
        uint64_t bs = 0;
        for (int i = 0; i < ITERATIONS; ++i) {
            BS_SET(bs, i & 63);
            dummy = BS_TEST(bs, i & 63);
            BS_RESET(bs, i & 63);
        }
        double time_ns = timer.elapsed_ns() / ITERATIONS;
        print_result("Macro version", time_ns, 3.0 / time_ns);
    }
}

// Test 2: Popcount (count set bits)
void test_popcount() {
    std::cout << "\n=== Popcount Operations ===\n";
    const int ITERATIONS = 500000000;
    volatile unsigned dummy = 0;
    
    uint64_t test_values[] = {
        0xFFFFFFFFFFFFFFFF, 0x5555555555555555, 0xAAAAAAAAAAAAAAAA,
        0x123456789ABCDEF0, 0xFEDCBA9876543210, 0x0F0F0F0F0F0F0F0F
    };
    const int NUM_VALUES = sizeof(test_values) / sizeof(test_values[0]);
    
    // Original class
    {
        Timer timer;
        for (int i = 0; i < ITERATIONS; ++i) {
            bitset64_original bs(test_values[i % NUM_VALUES]);
            dummy += bs.count();
        }
        double time_ns = timer.elapsed_ns() / ITERATIONS;
        print_result("Class popcount", time_ns, 1.0 / time_ns);
    }
    
    // Namespace version
    {
        Timer timer;
        for (int i = 0; i < ITERATIONS; ++i) {
            dummy += bitset64_ns::count(test_values[i % NUM_VALUES]);
        }
        double time_ns = timer.elapsed_ns() / ITERATIONS;
        print_result("Namespace popcount", time_ns, 1.0 / time_ns);
    }
    
    // Raw inline functions
    {
        Timer timer;
        for (int i = 0; i < ITERATIONS; ++i) {
            dummy += bitset_raw_count(test_values[i % NUM_VALUES]);
        }
        double time_ns = timer.elapsed_ns() / ITERATIONS;
        print_result("Raw inline popcount", time_ns, 1.0 / time_ns);
    }
    
    // Raw/Macro (same for popcount)
    {
        Timer timer;
        for (int i = 0; i < ITERATIONS; ++i) {
            dummy += BS_COUNT(test_values[i % NUM_VALUES]);
        }
        double time_ns = timer.elapsed_ns() / ITERATIONS;
        print_result("Builtin popcount", time_ns, 1.0 / time_ns);
    }
}

// Test 3: Iteration over set bits
void test_iteration() {
    std::cout << "\n=== Iteration Over Set Bits ===\n";
    const int ITERATIONS = 10000000;
    volatile unsigned dummy = 0;
    
    // Test with 16 bits set (so iteration does work)
    uint64_t test_value = 0xAAAAAAAAAAAAAAAA; // 1010... pattern
    
    // Original class
    {
        Timer timer;
        bitset64_original bs(test_value);
        for (int i = 0; i < ITERATIONS; ++i) {
            bs.for_each_set_bit([&dummy](unsigned pos) {
                dummy += pos;
                return true;
            });
        }
        double time_ns = timer.elapsed_ns() / ITERATIONS;
        print_result("Class iteration", time_ns, 32.0 / time_ns); // 32 bits set
    }
    
    // Namespace version
    {
        Timer timer;
        for (int i = 0; i < ITERATIONS; ++i) {
            bitset64_ns::for_each_set_bit(test_value, [&dummy](unsigned pos) {
                dummy += pos;
                return true;
            });
        }
        double time_ns = timer.elapsed_ns() / ITERATIONS;
        print_result("Namespace iteration", time_ns, 32.0 / time_ns);
    }
    
    // Raw inline functions
    {
        Timer timer;
        for (int i = 0; i < ITERATIONS; ++i) {
            bitset_raw_for_each_set_bit(test_value, [&dummy](unsigned pos) {
                dummy += pos;
                return true;
            });
        }
        double time_ns = timer.elapsed_ns() / ITERATIONS;
        print_result("Raw inline iteration", time_ns, 32.0 / time_ns);
    }
    
    // Macro iteration
    {
        Timer timer;
        for (int i = 0; i < ITERATIONS; ++i) {
            uint64_t bits = test_value;
            BS_FOREACH(bits, pos, {
                dummy += pos;
            });
        }
        double time_ns = timer.elapsed_ns() / ITERATIONS;
        print_result("Macro iteration", time_ns, 32.0 / time_ns);
    }
    
    // Raw while loop (fastest possible)
    {
        Timer timer;
        for (int i = 0; i < ITERATIONS; ++i) {
            uint64_t bits = test_value;
            while (bits) {
                unsigned pos;
                #ifdef _MSC_VER
                    unsigned long index;
                    _BitScanForward64(&index, bits);
                    pos = (unsigned)index;
                #else
                    pos = (unsigned)__builtin_ctzll(bits);
                #endif
                dummy += pos;
                bits &= bits - 1;
            }
        }
        double time_ns = timer.elapsed_ns() / ITERATIONS;
        print_result("Raw while iteration", time_ns, 32.0 / time_ns);
    }
}

// Test 4: Mixed operations (realistic workload)
void test_mixed_workload() {
    std::cout << "\n=== Mixed Operations (Realistic) ===\n";
    const int ITERATIONS = 50000000;
    volatile unsigned dummy = 0;
    
    srand(time(nullptr));
    
    // Original class
    {
        Timer timer;
        bitset64_original bs1(0x123456789ABCDEF0);
        bitset64_original bs2(0xFEDCBA9876543210);
        
        for (int i = 0; i < ITERATIONS; ++i) {
            // Simulate typical operations
            bs1.set(rand() & 63);
            bs2.reset(rand() & 63);
            dummy += bs1.count();
            bs1.for_each_set_bit([&dummy](unsigned pos) {
                dummy += pos;
                return true;
            });
        }
        double time_ns = timer.elapsed_ns() / ITERATIONS;
        print_result("Class mixed ops", time_ns, 5.0 / time_ns);
    }
    
    // Namespace version
    {
        Timer timer;
        bitset64_ns::bits_t bs1 = 0x123456789ABCDEF0;
        bitset64_ns::bits_t bs2 = 0xFEDCBA9876543210;
        
        for (int i = 0; i < ITERATIONS; ++i) {
            bitset64_ns::set(bs1, rand() & 63);
            bitset64_ns::reset(bs2, rand() & 63);
            dummy += bitset64_ns::count(bs1);
            bitset64_ns::for_each_set_bit(bs1, [&dummy](unsigned pos) {
                dummy += pos;
                return true;
            });
        }
        double time_ns = timer.elapsed_ns() / ITERATIONS;
        print_result("Namespace mixed ops", time_ns, 5.0 / time_ns);
    }
    
    // Raw inline functions
    {
        Timer timer;
        uint64_t bs1 = 0x123456789ABCDEF0;
        uint64_t bs2 = 0xFEDCBA9876543210;
        
        for (int i = 0; i < ITERATIONS; ++i) {
            bitset_raw_set(bs1, rand() & 63);
            bitset_raw_reset(bs2, rand() & 63);
            dummy += bitset_raw_count(bs1);
            bitset_raw_for_each_set_bit(bs1, [&dummy](unsigned pos) {
                dummy += pos;
                return true;
            });
        }
        double time_ns = timer.elapsed_ns() / ITERATIONS;
        print_result("Raw inline mixed ops", time_ns, 5.0 / time_ns);
    }
    
    // Raw operations
    {
        Timer timer;
        uint64_t bs1 = 0x123456789ABCDEF0;
        uint64_t bs2 = 0xFEDCBA9876543210;
        
        for (int i = 0; i < ITERATIONS; ++i) {
            bs1 |= (1ULL << (rand() & 63));
            bs2 &= ~(1ULL << (rand() & 63));
            dummy += __builtin_popcountll(bs1);
            
            uint64_t temp = bs1;
            while (temp) {
                unsigned pos = __builtin_ctzll(temp);
                dummy += pos;
                temp &= temp - 1;
            }
        }
        double time_ns = timer.elapsed_ns() / ITERATIONS;
        print_result("Raw mixed ops", time_ns, 5.0 / time_ns);
    }
}

// Test 5: Memory usage and cache effects
void test_memory_access() {
    std::cout << "\n=== Array Operations (Cache Effects) ===\n";
    const int ARRAY_SIZE = 1000000;
    const int ITERATIONS = 1000;
    
    // Create arrays
    bitset64_original* class_array = new bitset64_original[ARRAY_SIZE];
    bitset64_ns::bits_t* ns_array = new bitset64_ns::bits_t[ARRAY_SIZE];
    uint64_t* raw_array = new uint64_t[ARRAY_SIZE];
    
    // Initialize with random values
    srand(42);
    for (int i = 0; i < ARRAY_SIZE; ++i) {
        uint64_t val = (uint64_t)rand() << 32 | rand();
        class_array[i] = bitset64_original(val);
        ns_array[i] = val;
        raw_array[i] = val;
    }
    
    volatile uint64_t dummy = 0;
    
    // Class array processing
    {
        Timer timer;
        for (int iter = 0; iter < ITERATIONS; ++iter) {
            for (int i = 0; i < ARRAY_SIZE; ++i) {
                class_array[i].set(iter & 63);
                dummy += class_array[i].count();
            }
        }
        double time_ms = timer.elapsed_ms() / (ARRAY_SIZE * ITERATIONS);
        std::cout << "Class array: " << std::fixed << std::setprecision(3) 
                  << time_ms * 1e6 << " ns/element" << std::endl;
    }
    
    // Namespace array processing
    {
        Timer timer;
        for (int iter = 0; iter < ITERATIONS; ++iter) {
            for (int i = 0; i < ARRAY_SIZE; ++i) {
                bitset64_ns::set(ns_array[i], iter & 63);
                dummy += bitset64_ns::count(ns_array[i]);
            }
        }
        double time_ms = timer.elapsed_ms() / (ARRAY_SIZE * ITERATIONS);
        std::cout << "Namespace array: " << time_ms * 1e6 << " ns/element" << std::endl;
    }
    
    // Raw inline functions array processing
    {
        Timer timer;
        for (int iter = 0; iter < ITERATIONS; ++iter) {
            for (int i = 0; i < ARRAY_SIZE; ++i) {
                bitset_raw_set(raw_array[i], iter & 63);
                dummy += bitset_raw_count(raw_array[i]);
            }
        }
        double time_ms = timer.elapsed_ms() / (ARRAY_SIZE * ITERATIONS);
        std::cout << "Raw inline array: " << time_ms * 1e6 << " ns/element" << std::endl;
    }
    
    // Raw array processing
    {
        Timer timer;
        for (int iter = 0; iter < ITERATIONS; ++iter) {
            for (int i = 0; i < ARRAY_SIZE; ++i) {
                raw_array[i] |= (1ULL << (iter & 63));
                dummy += __builtin_popcountll(raw_array[i]);
            }
        }
        double time_ms = timer.elapsed_ms() / (ARRAY_SIZE * ITERATIONS);
        std::cout << "Raw array: " << time_ms * 1e6 << " ns/element" << std::endl;
    }
    
    // Cleanup
    delete[] class_array;
    delete[] ns_array;
    delete[] raw_array;
}

// ==================== MAIN ====================

int main() {
    std::cout << "=============================================\n";
    std::cout << "Bitset Performance Comparison Test\n";
    std::cout << "Compiled with: ";
    #ifdef __clang__
        std::cout << "Clang " << __clang_major__ << "." << __clang_minor__;
    #elif __GNUC__
        std::cout << "GCC " << __GNUC__ << "." << __GNUC_MINOR__;
    #elif _MSC_VER
        std::cout << "MSVC " << _MSC_VER;
    #else
        std::cout << "Unknown compiler";
    #endif
    std::cout << "\n";
    std::cout << "=============================================\n";
    
    // Run tests
    test_basic_operations();
    test_popcount();
    test_iteration();
    test_mixed_workload();
    test_memory_access();

    
    return 0;
}

