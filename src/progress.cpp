#include "progress.hpp"
#include <iomanip>

namespace stable_decomposition {

void print_pruning_comparison(std::ostream& out, const PruningProfile& matrix, const PruningProfile& module) {
  const auto flags = out.flags();
  const auto precision = out.precision();
  out << std::fixed << std::setprecision(6)
      << "\nPruning timing comparison (matrix first, module second)\n"
      << "Total pruning: matrix=" << matrix.total.seconds << " s, module=" << module.total.seconds << " s\n";
  if (matrix.total.seconds > 0)
    out << "Module / matrix total time: " << module.total.seconds / matrix.total.seconds << "\n";
  auto phase = [&](const char* name, int old_n, int new_n, double old_s, double new_s) {
    out << name << " iterations: matrix=" << old_n << ", module=" << new_n << '\n'
        << "  Phase time / average per iteration (s): matrix=" << old_s << " / " << (old_n ? old_s/old_n : 0)
        << ", module=" << new_s << " / " << (new_n ? new_s/new_n : 0) << '\n';
  };
  phase("I", matrix.i_iterations, module.i_iterations, matrix.i_phase.seconds, module.i_phase.seconds);
  phase("K", matrix.k_iterations, module.k_iterations, matrix.k_phase.seconds, module.k_phase.seconds);
  out << std::left << std::setw(25) << "Operation" << std::right
      << std::setw(12) << "Old calls" << std::setw(14) << "Old total s" << std::setw(14) << "Old avg ms"
      << std::setw(12) << "New calls" << std::setw(14) << "New total s" << std::setw(14) << "New avg ms" << '\n';
  auto row = [&](const char* name, PruningTiming PruningProfile::*field) {
    out << std::left << std::setw(25) << name << std::right;
    for (const auto* p : {&matrix, &module}) {
      const auto& t = p->*field;
      out << std::setw(12) << t.calls << std::setw(14) << t.seconds
          << std::setw(14) << (t.calls ? 1000*t.seconds/t.calls : 0);
    }
    out << '\n';
  };
  row("End_2eps/0", &PruningProfile::basis);
  row("I preimage", &PruningProfile::i_preimage);
  row("I intersection", &PruningProfile::i_intersection);
  row("I image / composition", &PruningProfile::i_image);
  row("I generator reduction", &PruningProfile::i_reduce);
  row("I convergence check", &PruningProfile::i_equal);
  row("K preimage", &PruningProfile::k_preimage);
  row("K image / composition", &PruningProfile::k_image);
  row("K sum / append", &PruningProfile::k_sum);
  row("K generator reduction", &PruningProfile::k_reduce);
  row("K convergence check", &PruningProfile::k_equal);
  row("Quotient construction", &PruningProfile::quotient);
  row("Presentation", &PruningProfile::presentation);
  row("Minimization", &PruningProfile::minimize);
  out << "Times are wall time. Total excludes input preparation, file I/O, images and AIDA.\n"
      << "Phase averages include convergence checks and setup; zero passes have average zero.\n"
      << "Operation times include work inside that API call; nested library operations are not counted separately.\n"
      << "Matrix I preimage uses f restricted to I; module I uses ambient preimage then intersection.\n"
      << "Matrix K takes preimages per map; module K takes one after summing per iteration.\n"
      << "Progress printing is disabled in both timed runs. Timings include instrumentation overhead.\n";
  out.flags(flags);
  out.precision(precision);
}

void print_progress(int iteration_I, size_t current, size_t total) {
  if(total != 0){
  int progress = (current * 100) / total;
  int bars = progress / 5; // 20 segments (100/5)

  std::cout << "\r Iteration " << iteration_I << " [";
  for (int i = 0; i < 20; ++i) {
    std::cout << (i < bars ? "█" : " ");
  }
  std::cout << "] " << progress << "%" << std::flush;
  }
}

} // namespace stable_decomposition
