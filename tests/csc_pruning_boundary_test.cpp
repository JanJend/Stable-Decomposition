// Compact storage is converted at the pruning boundary, never edited inside.
#define GRLINA_CSC_COMPILE_WARNINGS 0
#define GRLINA_CSC_RUNTIME_WARNINGS 1
#include <grlina/mat_grid.hpp>
#include "pruning.hpp"
#include <cassert>
#include <sstream>

using namespace graded_linalg;
namespace sd = stable_decomposition;

int main() {
    using Editable = R2GradedSparseMatrix<int>;
    using Compact = R2GradedSparseMatrix<int, CSCStorage<int>>;
    Editable matrix(sd::Mat(2, 1, {{0}, {0}}, {{0.25, 1.5}, {1.25, -2}}, {{0.25, -2}}));
    Module<Editable> source(matrix);
    source.compute_projective_resolution();
    assert(source.projective_resolution().size() > 1);
    source.set_injective_resolution(ChainComplex<Editable>(vec<Editable>{matrix}));
    const auto expected = sd::pruning(source, 0.5);
    auto compact = sd::convert_module_storage<CSCStorage<int>>(source);
    assert(compact.has_complete_projective_resolution());
    assert(compact.projective_resolution().size() == source.projective_resolution().size());
    assert(compact.has_injective_resolution());
    assert(compact.injective_resolution().size() == source.injective_resolution().size());
    assert(compact.presentation().row_degrees == source.presentation().row_degrees);
    auto parent = std::make_shared<const Module<Compact>>(compact);
    auto original = parent.get();
    std::weak_ptr<const Module<Compact>> lifetime = parent;
    const auto before = compact.presentation().data;
    std::ostringstream warnings;
    auto* old_buffer = std::cerr.rdbuf(warnings.rdbuf());
    auto actual = sd::pruning(compact, 0.5);
    auto pair = sd::pruning_pair(*parent, 0.5);
    std::cerr.rdbuf(old_buffer);
    assert(warnings.str().find("CSC performance warning") == std::string::npos);
    assert(actual.presentation().row_degrees == expected.presentation().row_degrees);
    assert(actual.presentation().col_degrees == expected.presentation().col_degrees);
    assert(compact.presentation().data.shares_storage_with(before));
    assert(pair.first.parent().get() == original && pair.second.parent().get() == original);
    assert(pair.first.generator_map().generator_lift().is_graded_matrix());
    assert(pair.second.generator_map().generator_lift().is_graded_matrix());
    parent.reset();
    assert(!lifetime.expired());
    assert(!pair.first.generator_map().generator_lift().get_col(0).empty());
    assert(pair.second.number_of_generators() == 0);

    const long double precise = 1.234567890123456789L;
    using Precise = GridGradedSparseMatrix<long double, 2, long>;
    Precise precision(0, 1, {}, {{precise, -precise}});
    auto precision_result = sd::pruning(Module<Precise>(precision), 0);
    assert(precision_result.presentation().real_row_degree(0)[0] == precise);
    assert(precision_result.presentation().grids == precision.grids);
}
