#include <iostream>



#include <vector>

#include <unordered_set>

#include <cstdint>

#include <chrono>

#include <algorithm>

#include <random>

struct bitset64 {

    uint64_t bits;

    

    bitset64(uint64_t b = 0) : bits(b) {}

    

    bool is_subset_of(const bitset64& other) const {

        return (other.bits & bits) == bits;

    }

    

    int popcount() const {

        return __builtin_popcountll(bits);

    }

    

    bool operator==(const bitset64& other) const {

        return bits == other.bits;

    }

};

namespace std {

    template<>

    struct hash<bitset64> {

        size_t operator()(const bitset64& bs) const {

            return std::hash<uint64_t>{}(bs.bits);

        }

    };

}

// ==================== APPROACH 1: Vector with remove_if ====================

class VectorRemoveIfApproach {

private:

    std::vector<std::vector<bitset64>> supports_;

    size_t dimension_;

    

public:

    VectorRemoveIfApproach(size_t dim) : dimension_(dim), supports_(dim + 1) {}

    

    void add(const bitset64& bs) {

        supports_[bs.popcount()].push_back(bs);

    }

    

    void remove_supersets(const bitset64& subset) {

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

    

    size_t count_remaining() const {

        size_t total = 0;

        for (const auto& vec : supports_) total += vec.size();

        return total;

    }

    

    uint64_t process_all() const {

        uint64_t checksum = 0;

        for (const auto& vec : supports_) {

            for (const auto& bs : vec) {

                checksum += bs.bits;

            }

        }

        return checksum;

    }

};

// ==================== APPROACH 2: Marking with vector<bool> ====================

class MarkingApproach {

private:

    std::vector<std::vector<bitset64>> supports_;

    std::vector<std::vector<bool>> is_valid_;

    size_t dimension_;

    

public:

    MarkingApproach(size_t dim) : dimension_(dim), supports_(dim + 1), is_valid_(dim + 1) {}

    

    void add(const bitset64& bs) {

        supports_[bs.popcount()].push_back(bs);

    }

    

    void finalize() {

        for (size_t i = 0; i <= dimension_; ++i) {

            is_valid_[i].resize(supports_[i].size(), true);

        }

    }

    

    void remove_supersets(const bitset64& subset) {

        size_t support_size = subset.popcount();

        

        for (size_t i = support_size + 1; i <= dimension_; ++i) {

            for (size_t j = 0; j < supports_[i].size(); ++j) {

                if (is_valid_[i][j] && subset.is_subset_of(supports_[i][j])) {

                    is_valid_[i][j] = false;

                }

            }

        }

    }

    

    size_t count_remaining() const {

        size_t total = 0;

        for (size_t i = 0; i <= dimension_; ++i) {

            for (size_t j = 0; j < is_valid_[i].size(); ++j) {

                if (is_valid_[i][j]) total++;

            }

        }

        return total;

    }

    

    uint64_t process_all() const {

        uint64_t checksum = 0;

        for (size_t i = 0; i <= dimension_; ++i) {

            for (size_t j = 0; j < supports_[i].size(); ++j) {

                if (is_valid_[i][j]) {

                    checksum += supports_[i][j].bits;

                }

            }

        }

        return checksum;

    }

};

// ==================== APPROACH 3: unordered_set ====================

class UnorderedSetApproach {

private:

    std::vector<std::unordered_set<bitset64>> supports_;

    size_t dimension_;

    

public:

    UnorderedSetApproach(size_t dim) : dimension_(dim), supports_(dim + 1) {}

    

    void add(const bitset64& bs) {

        supports_[bs.popcount()].insert(bs);

    }

    

    void remove_supersets(const bitset64& subset) {

        size_t support_size = subset.popcount();

        

        for (size_t i = support_size + 1; i <= dimension_; ++i) {

            std::vector<bitset64> to_remove;

            to_remove.reserve(supports_[i].size() / 10);

            

            for (const auto& x : supports_[i]) {

                if (subset.is_subset_of(x)) {

                    to_remove.push_back(x);

                }

            }

            

            for (const auto& x : to_remove) {

                supports_[i].erase(x);

            }

        }

    }

    

    size_t count_remaining() const {

        size_t total = 0;

        for (const auto& s : supports_) total += s.size();

        return total;

    }

    

    uint64_t process_all() const {

        uint64_t checksum = 0;

        for (const auto& set : supports_) {

            for (const auto& bs : set) {

                checksum += bs.bits;

            }

        }

        return checksum;

    }

};

// ==================== APPROACH 4: Manual Compact ====================

class ManualCompactApproach {

private:

    std::vector<std::vector<bitset64>> supports_;

    size_t dimension_;

    

public:

    ManualCompactApproach(size_t dim) : dimension_(dim), supports_(dim + 1) {}

    

    void add(const bitset64& bs) {

        supports_[bs.popcount()].push_back(bs);

    }

    

    void remove_supersets(const bitset64& subset) {

        size_t support_size = subset.popcount();

        

        for (size_t i = support_size + 1; i <= dimension_; ++i) {

            auto& vec = supports_[i];

            size_t write = 0;

            

            for (size_t read = 0; read < vec.size(); ++read) {

                if (!subset.is_subset_of(vec[read])) {

                    if (write != read) {

                        vec[write] = vec[read];

                    }

                    ++write;

                }

            }

            

            vec.resize(write);

        }

    }

    

    size_t count_remaining() const {

        size_t total = 0;

        for (const auto& vec : supports_) total += vec.size();

        return total;

    }

    

    uint64_t process_all() const {

        uint64_t checksum = 0;

        for (const auto& vec : supports_) {

            for (const auto& bs : vec) {

                checksum += bs.bits;

            }

        }

        return checksum;

    }

};

// ==================== APPROACH 5: Partition + Shrink ====================

class PartitionShrinkApproach {

private:

    std::vector<std::vector<bitset64>> supports_;

    size_t dimension_;

    

public:

    PartitionShrinkApproach(size_t dim) : dimension_(dim), supports_(dim + 1) {}

    

    void add(const bitset64& bs) {

        supports_[bs.popcount()].push_back(bs);

    }

    

    void remove_supersets(const bitset64& subset) {

        size_t support_size = subset.popcount();

        

        for (size_t i = support_size + 1; i <= dimension_; ++i) {

            auto& vec = supports_[i];

            

            auto new_end = std::partition(

                vec.begin(), vec.end(),

                [&](const bitset64& x) { return !subset.is_subset_of(x); }

            );

            

            vec.erase(new_end, vec.end());

        }

    }

    

    size_t count_remaining() const {

        size_t total = 0;

        for (const auto& vec : supports_) total += vec.size();

        return total;

    }

    

    uint64_t process_all() const {

        uint64_t checksum = 0;

        for (const auto& vec : supports_) {

            for (const auto& bs : vec) {

                checksum += bs.bits;

            }

        }

        return checksum;

    }

};

// ==================== BENCHMARK ====================

template<typename Approach>

void benchmark(const char* name, size_t n, size_t num_good_sets, const std::vector<bitset64>& good_sets) {

    std::cout << "\n=== " << name << " ===\n";

    

    Approach approach(n);

    

    // Add ALL 2^n bitsets

    auto start = std::chrono::high_resolution_clock::now();

    size_t total_items = 1ULL << n;

    for (uint64_t i = 0; i < total_items; ++i) {

        approach.add(bitset64(i));

        

        if (i % (total_items / 10) == 0 && i > 0) {

            std::cout << "  Adding: " << (100 * i / total_items) << "%\r" << std::flush;

        }

    }

    std::cout << "  Adding: 100%        \n";

    

    auto after_add = std::chrono::high_resolution_clock::now();

    

    // Finalize if needed (for marking approach)

    auto after_finalize = after_add;

    if constexpr (std::is_same_v<Approach, MarkingApproach>) {

        approach.finalize();

        after_finalize = std::chrono::high_resolution_clock::now();

    }

    

    // Remove supersets for all good sets

    std::cout << "  Removing supersets...\n";

    auto before_remove = std::chrono::high_resolution_clock::now();

    

    for (size_t idx = 0; idx < good_sets.size(); ++idx) {

        approach.remove_supersets(good_sets[idx]);

        

        if ((idx + 1) % 500 == 0) {

            std::cout << "    Processed " << (idx + 1) << " / " << good_sets.size() 

                     << " good sets\r" << std::flush;

        }

    }

    std::cout << "    Processed " << good_sets.size() << " / " << good_sets.size() 

             << " good sets\n";

    

    auto after_remove = std::chrono::high_resolution_clock::now();

    

    // Process all remaining items

    std::cout << "  Processing remaining items...\n";

    uint64_t checksum = approach.process_all();

    auto after_process = std::chrono::high_resolution_clock::now();

    

    size_t remaining = approach.count_remaining();

    

    auto add_time = std::chrono::duration_cast<std::chrono::milliseconds>(

        after_add - start);

    auto finalize_time = std::chrono::duration_cast<std::chrono::milliseconds>(

        after_finalize - after_add);

    auto remove_time = std::chrono::duration_cast<std::chrono::milliseconds>(

        after_remove - before_remove);

    auto process_time = std::chrono::duration_cast<std::chrono::milliseconds>(

        after_process - after_remove);

    

    std::cout << "\nResults:\n";

    std::cout << "  Add time:      " << add_time.count() << " ms\n";

    if (finalize_time.count() > 0) {

        std::cout << "  Finalize time: " << finalize_time.count() << " ms\n";

    }

    std::cout << "  Remove time:   " << remove_time.count() << " ms *** KEY METRIC ***\n";

    std::cout << "  Process time:  " << process_time.count() << " ms\n";

    std::cout << "  Remaining:     " << remaining << " / " << total_items << "\n";

    std::cout << "  Checksum:      " << checksum << "\n";

}

int main() {

    // Configuration

    const size_t N = 22;

    const size_t NUM_GOOD_SETS = 5000;

    const uint64_t TOTAL_ITEMS = 1ULL << N;

    

    std::cout << "=== SUPERSET REMOVAL BENCHMARK ===\n";

    std::cout << "Configuration:\n";

    std::cout << "  n = " << N << "\n";

    std::cout << "  Total items = 2^" << N << " = " << TOTAL_ITEMS << "\n";

    std::cout << "  Number of good sets = " << NUM_GOOD_SETS << "\n\n";

    

    // Generate random good sets

    std::cout << "Generating " << NUM_GOOD_SETS << " random good sets...\n";

    std::random_device rd;

    std::mt19937_64 gen(rd());

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

    std::cout << "Average popcount of good sets: " 

              << (std::accumulate(good_sets.begin(), good_sets.end(), 0.0,

                    [](double sum, const bitset64& bs) { return sum + bs.popcount(); })

                  / good_sets.size()) << "\n";

    

    std::cout << "\n" << std::string(70, '=') << "\n";

    

    // Run benchmarks

    benchmark<MarkingApproach>("1. MARKING ONLY (vector<bool>)", N, NUM_GOOD_SETS, good_sets);

    

    std::cout << "\n" << std::string(70, '=') << "\n";

    

    benchmark<VectorRemoveIfApproach>("2. Vector + remove_if + erase", N, NUM_GOOD_SETS, good_sets);

    

    std::cout << "\n" << std::string(70, '=') << "\n";

    

    benchmark<ManualCompactApproach>("3. Vector + Manual Compact", N, NUM_GOOD_SETS, good_sets);

    

    std::cout << "\n" << std::string(70, '=') << "\n";

    

    benchmark<PartitionShrinkApproach>("4. Vector + Partition + Shrink", N, NUM_GOOD_SETS, good_sets);

    

    std::cout << "\n" << std::string(70, '=') << "\n";

    

    benchmark<UnorderedSetApproach>("5. unordered_set (O(1) erase)", N, NUM_GOOD_SETS, good_sets);

    

    std::cout << "\n" << std::string(70, '=') << "\n";

    std::cout << "\nBenchmark complete!\n";

    std::cout << "Compare the 'Remove time' values to see which is fastest.\n";

    

    return 0;

}
