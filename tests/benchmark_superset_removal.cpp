#include <iostream>
#include <vector>
#include <algorithm>
#include <chrono>
#include <random>
#include <cstdint>
#include <unordered_set>
#include <numeric>

struct bitset64 {
    uint64_t bits;
    
    bitset64(uint64_t b = 0) : bits(b) {}
    
    bool is_subset_of(const bitset64& other) const {
        return (other.bits & bits) == bits;
    }
    
    int popcount() const {
        return __builtin_popcountll(bits);
    }
};

// ==================== APPROACH 1: Iterative (Current) ====================
class IterativeApproach {
private:
    std::vector<std::vector<bitset64>> supports_;
    size_t dimension_;
    
public:
    IterativeApproach(size_t dim) : dimension_(dim), supports_(dim + 1) {}
    
    void add(const bitset64& bs) {
        supports_[bs.popcount()].push_back(bs);
    }
    
    void remove_supersets(const std::vector<bitset64>& subset_list) {
        // For each good set, do a separate remove_if + erase
        for (const auto& subset : subset_list) {
            size_t support_size = subset.popcount();
            
            for (size_t i = support_size + 1; i <= dimension_; ++i) {
                supports_[i].erase(
                    std::remove_if(
                        supports_[i].begin(),
                        supports_[i].end(),
                        [&](const bitset64& x) { return subset.is_subset_of(x); }
                    ),
                    supports_[i].end()
                );
            }
        }
    }
    
    size_t count_remaining() const {
        size_t total = 0;
        for (const auto& vec : supports_) total += vec.size();
        return total;
    }
    
    uint64_t checksum() const {
        uint64_t sum = 0;
        for (const auto& vec : supports_) {
            for (const auto& bs : vec) {
                sum += bs.bits;
            }
        }
        return sum;
    }
};

// ==================== APPROACH 2: Batch (Your Proposed) ====================
class BatchApproach {
private:
    std::vector<std::vector<bitset64>> supports_;
    size_t dimension_;
    
public:
    BatchApproach(size_t dim) : dimension_(dim), supports_(dim + 1) {}
    
    void add(const bitset64& bs) {
        supports_[bs.popcount()].push_back(bs);
    }
    
    void remove_supersets(const std::vector<bitset64>& subset_list) {
        // Determine minimum popcount to start from
        size_t min_support_size = dimension_;
        for (const auto& subset : subset_list) {
            min_support_size = std::min(min_support_size, (size_t)subset.popcount());
        }
        
        // Single pass for each size class
        for (size_t i = min_support_size + 1; i <= dimension_; ++i) {
            auto& vec = supports_[i];
            
            vec.erase(
                std::remove_if(
                    vec.begin(),
                    vec.end(),
                    [&](const bitset64& x) {
                        return std::any_of(
                            subset_list.begin(),
                            subset_list.end(),
                            [&](const bitset64& s) { return s.is_subset_of(x); }
                        );
                    }
                ),
                vec.end()
            );
        }
    }
    
    size_t count_remaining() const {
        size_t total = 0;
        for (const auto& vec : supports_) total += vec.size();
        return total;
    }
    
    uint64_t checksum() const {
        uint64_t sum = 0;
        for (const auto& vec : supports_) {
            for (const auto& bs : vec) {
                sum += bs.bits;
            }
        }
        return sum;
    }
};

// ==================== APPROACH 3: Batch with Partition ====================
class BatchPartitionApproach {
private:
    std::vector<std::vector<bitset64>> supports_;
    size_t dimension_;
    
public:
    BatchPartitionApproach(size_t dim) : dimension_(dim), supports_(dim + 1) {}
    
    void add(const bitset64& bs) {
        supports_[bs.popcount()].push_back(bs);
    }
    
    void remove_supersets(const std::vector<bitset64>& subset_list) {
        // Determine minimum popcount
        size_t min_support_size = dimension_;
        for (const auto& subset : subset_list) {
            min_support_size = std::min(min_support_size, (size_t)subset.popcount());
        }
        
        // Single partition pass for each size class
        for (size_t i = min_support_size + 1; i <= dimension_; ++i) {
            auto& vec = supports_[i];
            
            auto new_end = std::partition(
                vec.begin(),
                vec.end(),
                [&](const bitset64& x) {
                    return std::none_of(
                        subset_list.begin(),
                        subset_list.end(),
                        [&](const bitset64& s) { return s.is_subset_of(x); }
                    );
                }
            );
            
            vec.erase(new_end, vec.end());
        }
    }
    
    size_t count_remaining() const {
        size_t total = 0;
        for (const auto& vec : supports_) total += vec.size();
        return total;
    }
    
    uint64_t checksum() const {
        uint64_t sum = 0;
        for (const auto& vec : supports_) {
            for (const auto& bs : vec) {
                sum += bs.bits;
            }
        }
        return sum;
    }
};

// ==================== BENCHMARK ====================
template<typename Approach>
void benchmark(const char* name, size_t n, const std::vector<bitset64>& good_sets) {
    std::cout << "\n=== " << name << " ===\n";
    
    Approach approach(n);
    
    // Add ALL 2^n bitsets
    std::cout << "  Adding items...\n";
    auto start = std::chrono::high_resolution_clock::now();
    size_t total_items = 1ULL << n;
    for (uint64_t i = 0; i < total_items; ++i) {
        approach.add(bitset64(i));
    }
    auto after_add = std::chrono::high_resolution_clock::now();
    
    // Remove supersets (BATCH - all good sets at once)
    std::cout << "  Removing supersets (batch of " << good_sets.size() << " good sets)...\n";
    auto before_remove = std::chrono::high_resolution_clock::now();
    
    approach.remove_supersets(good_sets);
    
    auto after_remove = std::chrono::high_resolution_clock::now();
    
    // Get results
    size_t remaining = approach.count_remaining();
    uint64_t check = approach.checksum();
    
    auto add_time = std::chrono::duration_cast<std::chrono::milliseconds>(
        after_add - start);
    auto remove_time = std::chrono::duration_cast<std::chrono::milliseconds>(
        after_remove - before_remove);
    
    std::cout << "\nResults:\n";
    std::cout << "  Add time:    " << add_time.count() << " ms\n";
    std::cout << "  Remove time: " << remove_time.count() << " ms *** KEY METRIC ***\n";
    std::cout << "  Remaining:   " << remaining << " / " << total_items << "\n";
    std::cout << "  Checksum:    " << check << "\n";
}

int main() {
    // Configuration
    const size_t N = 22;
    const size_t NUM_GOOD_SETS = 5000;
    const uint64_t TOTAL_ITEMS = 1ULL << N;
    
    std::cout << "=== BATCH vs ITERATIVE ERASE BENCHMARK ===\n";
    std::cout << "Configuration:\n";
    std::cout << "  n = " << N << "\n";
    std::cout << "  Total items = 2^" << N << " = " << TOTAL_ITEMS << "\n";
    std::cout << "  Number of good sets = " << NUM_GOOD_SETS << "\n\n";
    
    // Generate random good sets
    std::cout << "Generating " << NUM_GOOD_SETS << " random good sets...\n";
    std::random_device rd;
    std::mt19937_64 gen(42);  // Fixed seed for reproducibility
    std::uniform_int_distribution<uint64_t> dist(0, TOTAL_ITEMS - 1);
    
    std::vector<bitset64> good_sets;
    good_sets.reserve(NUM_GOOD_SETS);
    
    std::unordered_set<uint64_t> used;
    while (good_sets.size() < NUM_GOOD_SETS) {
        uint64_t val = dist(gen);
        if (used.insert(val).second) {
            good_sets.push_back(bitset64(val));
        }
    }
    
    std::cout << "Good sets generated!\n";
    std::cout << "Average popcount: " 
              << (std::accumulate(good_sets.begin(), good_sets.end(), 0.0,
                    [](double sum, const bitset64& bs) { return sum + bs.popcount(); })
                  / good_sets.size()) << "\n";
    
    std::cout << "\n" << std::string(70, '=') << "\n";
    
    // Run benchmarks
    benchmark<IterativeApproach>(
        "1. ITERATIVE (loop over each good set separately)", 
        N, good_sets);
    
    std::cout << "\n" << std::string(70, '=') << "\n";
    
    benchmark<BatchApproach>(
        "2. BATCH with remove_if (YOUR PROPOSED)", 
        N, good_sets);
    
    std::cout << "\n" << std::string(70, '=') << "\n";
    
    benchmark<BatchPartitionApproach>(
        "3. BATCH with partition", 
        N, good_sets);
    
    std::cout << "\n" << std::string(70, '=') << "\n";
    std::cout << "\n=== ANALYSIS ===\n";
    std::cout << "ITERATIVE: Multiple passes (5000x), multiple erases (5000x)\n";
    std::cout << "BATCH:     Single pass, single erase, but checks all 5000 subsets per element\n";
    std::cout << "\nWhich is faster? Let's see!\n";
    
    return 0;
}
