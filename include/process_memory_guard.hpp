#pragma once
#include <atomic>
#include <chrono>
#include <condition_variable>
#include <cstdlib>
#include <iostream>
#include <mutex>
#include <sys/resource.h>
#include <thread>

namespace stable_decomposition {

// Optional sampled peak-RSS guard for reproducible resource-limited runs.
class ProcessMemoryGuard {
    std::atomic<bool> finished_{false};
    std::mutex mutex_;
    std::condition_variable changed_;
    std::thread thread_;
public:
    explicit ProcessMemoryGuard(double limit_mb) {
        if (limit_mb <= 0) return;
        thread_ = std::thread([this, limit_mb] {
            while (!finished_.load()) {
                rusage usage{};
                if (getrusage(RUSAGE_SELF, &usage)) std::_Exit(76);
#ifdef __APPLE__
                const double megabytes = usage.ru_maxrss / (1024.0 * 1024.0);
#else
                const double megabytes = usage.ru_maxrss / 1024.0;
#endif
                if (megabytes > limit_mb) {
                    std::cerr << "PRUNING_MEMORY_LIMIT " << megabytes << " MiB\n";
                    std::_Exit(75);
                }
                std::unique_lock<std::mutex> lock(mutex_);
                changed_.wait_for(lock, std::chrono::milliseconds(200),
                                  [this] { return finished_.load(); });
            }
        });
    }
    ~ProcessMemoryGuard() {
        {
            std::lock_guard<std::mutex> lock(mutex_);
            finished_ = true;
        }
        changed_.notify_all();
        if (thread_.joinable()) thread_.join();
    }
};

} // namespace stable_decomposition
