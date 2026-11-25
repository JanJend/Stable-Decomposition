#pragma once
#include <grlina/general.hpp>
#include <chrono>
#include <functional>
#include "pruning_options.hpp"

namespace stable_decomposition {
using graded_linalg::timed_with_progress;
void print_progress(int iteration, std::size_t current, std::size_t total);

struct PruningTiming {
    std::size_t calls = 0;
    double seconds = 0;
};

// Optional output statistics, filled when the caller passes a pointer to this
// object. Records durations, operation counts, iterations, and map-family sizes.
struct PruningProfile {
    int i_iterations = 0, k_iterations = 0;
    std::size_t original_maps = 0, quotient_dimension = 0, end_dimension = 0;
    std::size_t i_maps = 0, k_maps = 0;
    std::size_t i_restriction_candidates = 0, i_restriction_kept = 0;
    std::size_t k_restriction_candidates = 0, k_restriction_kept = 0;
    std::size_t i_restriction_max = 0, k_restriction_max = 0;
    PruningTiming total, basis, i_phase, k_phase;
    PruningTiming end0, hom_quotient, composition, map_lifts;
    PruningTiming i_preimage, k_preimage, i_image, k_image;
    PruningTiming i_intersection, k_sum, i_reduce, k_reduce, i_equal, k_equal;
    PruningTiming quotient, presentation, minimize;
    PruningTiming i_adaptive_restriction, k_adaptive_restriction, i_target_prepare;
    PruningTiming i_shift, i_target_reduce;
    PruningTiming k_merge;
    // Optional start/completion events, possibly concurrent for the same field.
    // Parallel I/K use worker-local counters; callbacks must pair events per thread.
    std::function<void(PruningTiming PruningProfile::*, bool, double)> trace_step;
    // Stored generator columns, without additional algebra. I/K can report concurrently.
    // Arguments: object (I/K/J/L), stage, pass (0 = initialization), map (0 = no map), count.
    std::function<void(const char*, const char*, int, std::size_t, std::size_t)> trace_generators;
};

inline void record_pruning_generators(PruningProfile* profile, const char* object,
                                     const char* stage, int pass, std::size_t map,
                                     std::size_t count) {
    if (profile && profile->trace_generators)
        profile->trace_generators(object, stage, pass, map, count);
}

// A null destination leaves ordinary runs uninstrumented (no clock reads).
class PruningTimer {
    using Clock = std::chrono::steady_clock;
    PruningTiming* timing_;
    Clock::time_point start_;
    PruningProfile* profile_ = nullptr;
    PruningTiming PruningProfile::* field_ = nullptr;
public:
    explicit PruningTimer(PruningTiming* timing) : timing_(timing) {
        if (timing_) start_ = Clock::now();
    }
    PruningTimer(PruningProfile* profile, PruningTiming PruningProfile::* field)
        : timing_(profile ? &(profile->*field) : nullptr), profile_(profile), field_(field) {
        if (profile_ && profile_->trace_step) profile_->trace_step(field_, false, 0);
        if (timing_) start_ = Clock::now();
    }
    ~PruningTimer() {
        if (timing_) {
            const auto seconds = std::chrono::duration<double>(Clock::now() - start_).count();
            timing_->seconds += seconds;
            ++timing_->calls;
            if (profile_ && profile_->trace_step) profile_->trace_step(field_, true, seconds);
        }
    }
    PruningTimer(const PruningTimer&) = delete;
    PruningTimer& operator=(const PruningTimer&) = delete;
};

// Run an operation and return its result; also record its time and call count
// in the selected field when profile is provided.
template <typename Function>
decltype(auto) measure(PruningProfile* profile, PruningTiming PruningProfile::* field, Function&& operation) {
    PruningTimer timer(profile, field);
    return operation();
}

void print_pruning_comparison(std::ostream& out, const PruningProfile& matrix, const PruningProfile& module);
void print_module_pruning_profile(std::ostream& out, const PruningProfile& profile,
                                  const PruningOptions& options);
} // namespace stable_decomposition
