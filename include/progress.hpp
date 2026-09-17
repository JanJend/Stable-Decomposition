#pragma once
#include <grlina/progress.hpp>
#include <chrono>

namespace stable_decomposition {
using graded_linalg::timed_with_progress;
void print_progress(int iteration, std::size_t current, std::size_t total);

struct PruningTiming {
    std::size_t calls = 0;
    double seconds = 0;
};

struct PruningProfile {
    int i_iterations = 0, k_iterations = 0;
    PruningTiming total, basis, i_phase, k_phase;
    PruningTiming i_preimage, k_preimage, i_image, k_image;
    PruningTiming i_intersection, k_sum, i_reduce, k_reduce, i_equal, k_equal;
    PruningTiming quotient, presentation, minimize;
};

// A null destination leaves ordinary runs uninstrumented (no clock reads).
class PruningTimer {
    using Clock = std::chrono::steady_clock;
    PruningTiming* timing_;
    Clock::time_point start_;
public:
    explicit PruningTimer(PruningTiming* timing) : timing_(timing) {
        if (timing_) start_ = Clock::now();
    }
    ~PruningTimer() {
        if (timing_) {
            timing_->seconds += std::chrono::duration<double>(Clock::now() - start_).count();
            ++timing_->calls;
        }
    }
    PruningTimer(const PruningTimer&) = delete;
    PruningTimer& operator=(const PruningTimer&) = delete;
};

template <typename Function>
decltype(auto) measure(PruningProfile* profile, PruningTiming PruningProfile::* field, Function&& operation) {
    PruningTimer timer(profile ? &(profile->*field) : nullptr);
    return operation();
}

void print_pruning_comparison(std::ostream& out, const PruningProfile& matrix, const PruningProfile& module);
} // namespace stable_decomposition
