#include "progress.hpp"
#include <iomanip>

namespace stable_decomposition {

void print_pruning_comparison(std::ostream& out, const PruningProfile& matrix, const PruningProfile& module) {
  const auto flags = out.flags();
  const auto precision = out.precision();
  out << std::defaultfloat << std::setprecision(3) << '\n'
      << std::left << std::setw(25) << "Operation" << std::right
      << std::setw(10) << "Old calls" << std::setw(10) << "New calls"
      << std::setw(11) << "Old s" << std::setw(11) << "New s"
      << std::setw(12) << "Old avg ms" << std::setw(12) << "New avg ms" << '\n';
  auto row = [&](const char* name, PruningTiming a, PruningTiming b, bool average = true) {
    out << std::left << std::setw(25) << name << std::right
        << std::setw(10) << a.calls << std::setw(10) << b.calls
        << std::setw(11) << a.seconds << std::setw(11) << b.seconds;
    for (const auto& t : {a, b}) {
      if (average && t.calls) out << std::setw(12) << 1000*t.seconds/t.calls;
      else out << std::setw(12) << "-";
    }
    out << '\n';
  };
  row("Total", matrix.total, module.total, false);
  row("I iterations", {static_cast<std::size_t>(matrix.i_iterations), matrix.i_phase.seconds},
                      {static_cast<std::size_t>(module.i_iterations), module.i_phase.seconds});
  row("K iterations", {static_cast<std::size_t>(matrix.k_iterations), matrix.k_phase.seconds},
                      {static_cast<std::size_t>(module.k_iterations), module.k_phase.seconds});
  row("End_2eps/0", matrix.basis, module.basis);
  row("I preimage", matrix.i_preimage, module.i_preimage);
  row("I intersection", matrix.i_intersection, module.i_intersection);
  row("I image / composition", matrix.i_image, module.i_image);
  row("I generator reduction", matrix.i_reduce, module.i_reduce);
  row("I convergence check", matrix.i_equal, module.i_equal);
  row("K preimage", matrix.k_preimage, module.k_preimage);
  row("K image / composition", matrix.k_image, module.k_image);
  row("K sum / append", matrix.k_sum, module.k_sum);
  row("K generator reduction", matrix.k_reduce, module.k_reduce);
  row("K convergence check", matrix.k_equal, module.k_equal);
  row("Quotient construction", matrix.quotient, module.quotient);
  row("Presentation", matrix.presentation, module.presentation);
  row("Minimization", matrix.minimize, module.minimize);
  out << "Module image/preimage/sum times include generator minimization.\n";
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
