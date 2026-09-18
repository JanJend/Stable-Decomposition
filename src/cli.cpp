#include "cli.hpp"
#include "pruning.hpp"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include <grlina/draw_hf.hpp>
#include <grlina/isomorphism_test.hpp>
#ifdef PRUNING_WITH_AIDA
#include <aida_interface.hpp>
#endif
#include <filesystem>
#include <fstream>
#include <map>
#include <sstream>

namespace stable_decomposition {
namespace fs = std::filesystem;
namespace {

std::ofstream output_file(const fs::path& path) {
    if (!path.parent_path().empty()) fs::create_directories(path.parent_path());
    std::ofstream out(path);
    if (!out) throw std::runtime_error("Cannot write " + path.string());
    out.exceptions(std::ios::badbit | std::ios::failbit);
    return out;
}

} // namespace

void write_module(const Module& module, const fs::path& path) {
    auto out = output_file(path);
    module.to_stream(out);
}


void write_hilbert_images(Module& input, Module& output, const fs::path& prefix, int size) {
    auto degrees = input.support_degrees();
    const auto output_degrees = output.support_degrees();
    degrees.insert(degrees.end(), output_degrees.begin(), output_degrees.end());
    if (degrees.empty()) degrees.push_back({0,0});
    double xmin=degrees[0].first, xmax=xmin, ymin=degrees[0].second, ymax=ymin;
    for (const auto& d : degrees) {
        xmin=std::min(xmin,d.first); xmax=std::max(xmax,d.first);
        ymin=std::min(ymin,d.second); ymax=std::max(ymax,d.second);
    }
    const double dx = xmax > xmin ? (xmax-xmin)*0.05 : 0.1;
    const double dy = ymax > ymin ? (ymax-ymin)*0.05 : 0.1;
    xmin-=dx; xmax+=dx; ymin-=dy; ymax+=dy;
    vec<double> xs(size), ys(size);
    for (int i=0; i<size; ++i) {
        xs[i]=xmin+(xmax-xmin)*i/(size-1);
        ys[i]=ymin+(ymax-ymin)*i/(size-1);
    }
    std::cout << "Computing Hilbert images on a shared " << size << " x " << size << " grid...\n";
    // Explicitly complete the resolutions: counting presentation degrees alone
    // would omit syzygies and can give incorrect dimensions.
    input.compute_projective_resolution();
    output.compute_projective_resolution();
    const auto before = input.hilbert_function_on_grid(xs,ys);
    const auto after = output.hilbert_function_on_grid(xs,ys);
    int positive_min=std::numeric_limits<int>::max(), positive_max=0;
    for (const auto* grid : {&before,&after}) for (const auto& row : grid->values) for (int value : row)
        if (value > 0) { positive_min=std::min(positive_min,value); positive_max=std::max(positive_max,value); }
    if (positive_max == 0) positive_min=1;
    const auto before_path=prefix.string()+"_input_hilbert.png";
    const auto after_path=prefix.string()+"_output_hilbert.png";
    graded_linalg::save_hilbert_png(before,before_path,positive_min,positive_max,"Input Hilbert function");
    graded_linalg::save_hilbert_png(after,after_path,positive_min,positive_max,"Output Hilbert function");
    std::cout << "Saved Hilbert images: " << before_path << "\n                      " << after_path << '\n';
}

#ifdef PRUNING_WITH_AIDA
namespace {
struct Decomposition {
    aida::Block_list blocks;
    std::map<std::string,int> types, signatures;
};

Decomposition decompose(const Module& module, const fs::path& path = {}) {
    Decomposition result;
    Module working(module.presentation());
    working.minimize(); // AIDA requires a minimal presentation.
    if (working.number_of_generators() != 0) {
        working.sort_compatibly();
        working.mutable_presentation().compute_col_batches(); // Sets k_max before AIDA loads its tables.
        aida::AIDA_functor decomposer;
        decomposer.config.sort = true;
        decomposer(working,result.blocks);
    }
    std::ofstream file;
    if (!path.empty()) {
        file = output_file(path);
        file << "scc2020sum\n" << result.blocks.size() << '\n';
    }
    for (auto& block : result.blocks) {
        ++result.types[block.get_type()];
        if (!path.empty()) {
            file << '\n' << block.get_type() << '\n';
            block.to_stream_r2(file);
        }
        Module summand(static_cast<const Mat&>(block));
        summand.minimize();
        // Compare graded Betti signatures, not coefficient matrices or presumed
        // isomorphism classes. Minimal presentations can have different bases.
        std::ostringstream key;
        key << std::setprecision(17);
        for (const auto* grades : {&summand.presentation().row_degrees,&summand.presentation().col_degrees}) {
            auto sorted=*grades;
            std::sort(sorted.begin(),sorted.end(),graded_linalg::Degree_traits<graded_linalg::r2degree>::lex_lambda());
            key << '[';
            for (const auto& d : sorted) key << '(' << d.first << ',' << d.second << ')';
            key << ']';
        }
        ++result.signatures[key.str()];
    }
    return result;
}

} // namespace

void compare_decompositions(const Module& input, const Module& output, const fs::path& prefix) {
    const auto before_path=prefix.string()+"_input_decomposition.sccsum";
    const auto after_path=prefix.string()+"_output_decomposition.sccsum";
    std::cout << "Running AIDA on the input module...\n" << std::flush;
    auto before=decompose(input,before_path);
    std::cout << "Running AIDA on the pruned module...\n" << std::flush;
    auto after=decompose(output,after_path);
    int shared=0;
    for (const auto& entry : before.signatures) {
        auto found=after.signatures.find(entry.first);
        if (found != after.signatures.end()) shared+=std::min(entry.second,found->second);
    }
    std::ostringstream summary;
    summary << "AIDA decomposition comparison\n"
            << "Summands: input=" << before.blocks.size() << ", output=" << after.blocks.size() << '\n'
            << "Type            Input    Output\n";
    for (const std::string type : {"free","cyclic","interval","non-interval"})
        summary << std::left << std::setw(16) << type << std::setw(9) << before.types[type] << after.types[type] << '\n';
    summary << "Summands with matching graded Betti signatures: " << shared << '\n'
            << "Compared at their actual grades, including pruning's final shift.\n"
            << "Matching signatures do not prove isomorphism. Different decompositions are expected after pruning.\n";
    const auto report_path=prefix.string()+"_decomposition_comparison.txt";
    auto report=output_file(report_path);
    report << summary.str() << "\nSignature multiplicities (generator degrees then relation degrees)\n"
           << "Input  Output  Signature\n";
    auto signatures=before.signatures;
    signatures.insert(after.signatures.begin(),after.signatures.end());
    for (const auto& entry : signatures)
        report << before.signatures[entry.first] << "      " << after.signatures[entry.first] << "       " << entry.first << '\n';
    std::cout << summary.str() << "Saved decompositions: " << before_path << "\n                      " << after_path
              << "\nSaved comparison: " << report_path << '\n';
}
#endif

Module compare_pruning(const Module& input, double epsilon, const fs::path& prefix) {
    // Prepare independent copies of the same presentation before timing either run.
    Mat old_input = input.presentation();
    Module new_input(input.presentation());
    PruningProfile old_profile, new_profile;
    std::cout << "Running matrix pruning...\n" << std::flush;
    Module old_output(pruning_profiled(old_input, epsilon, false, &old_profile));
    std::cout << "Running module pruning...\n" << std::flush;
    Module new_output = pruning(std::move(new_input), epsilon, false, &new_profile);

    std::ostringstream report;
    print_pruning_comparison(report, old_profile, new_profile);
    // Normalize copies with the same minimizer; all checks are outside pruning timers.
    Module old_checked(old_output.presentation()), new_checked(new_output.presentation());
    old_checked.minimize();
    new_checked.minimize();
    const auto& old_p = old_checked.presentation();
    const auto& new_p = new_checked.presentation();
    auto compare_degrees = [&](const char* label, const auto& a, const auto& b) {
        std::map<graded_linalg::r2degree, std::pair<std::size_t, std::size_t>> counts;
        for (const auto& d : a) ++counts[d].first;
        for (const auto& d : b) ++counts[d].second;
        const bool same = std::all_of(counts.begin(), counts.end(), [](const auto& entry) {
            return entry.second.first == entry.second.second;
        });
        report << label << " degrees: " << (same ? "match" : "DIFFER") << '\n';
        for (const auto& entry : counts) if (entry.second.first != entry.second.second)
            report << std::setprecision(17) << "  (" << entry.first.first << ", " << entry.first.second
                   << "): matrix=" << entry.second.first << ", module=" << entry.second.second << '\n';
        return same;
    };
    bool match = compare_degrees("Generator", old_p.row_degrees, new_p.row_degrees);
    match = compare_degrees("Relation", old_p.col_degrees, new_p.col_degrees) && match;
    const bool degrees_match = match;
#ifdef PRUNING_WITH_AIDA
    std::cout << "Running AIDA on both pruning results...\n" << std::flush;
    const auto old_decomposition = decompose(old_checked);
    const auto new_decomposition = decompose(new_checked);
    const bool counts_match = old_decomposition.blocks.size() == new_decomposition.blocks.size();
    const bool types_match = old_decomposition.types == new_decomposition.types;
    const bool signatures_match = old_decomposition.signatures == new_decomposition.signatures;
    if (counts_match && types_match && signatures_match) {
        report << "AIDA: match; " << old_decomposition.blocks.size() << " indecomposables (";
        const char* separator = "";
        for (const auto& entry : old_decomposition.types) {
            report << separator << entry.second << ' ' << entry.first;
            separator = ", ";
        }
        report << "); degree signatures match.\n";
    } else {
        report << "AIDA: DIFFER; indecomposables matrix=" << old_decomposition.blocks.size()
               << ", module=" << new_decomposition.blocks.size() << '\n';
        auto types = old_decomposition.types;
        types.insert(new_decomposition.types.begin(), new_decomposition.types.end());
        for (const auto& entry : types) {
            auto count = [&](const auto& decomposition) {
                auto found = decomposition.types.find(entry.first);
                return found == decomposition.types.end() ? 0 : found->second;
            };
            if (count(old_decomposition) != count(new_decomposition))
                report << "  " << entry.first << ": matrix=" << count(old_decomposition)
                       << ", module=" << count(new_decomposition) << '\n';
        }
        report << "  Degree signatures: " << (signatures_match ? "match" : "DIFFER") << '\n';
    }
    match = match && counts_match && types_match && signatures_match;
#else
    report << "AIDA: unavailable.\n";
#endif
    auto multiplicity = [](const auto& degrees) {
        std::map<graded_linalg::r2degree, std::size_t> counts;
        std::size_t largest = 0;
        for (const auto& d : degrees) largest = std::max(largest, ++counts[d]);
        return largest;
    };
    const auto k = std::max(multiplicity(old_p.row_degrees), multiplicity(new_p.row_degrees));
    // Benchmarked up to k=4; larger blocks can require exponential enumeration.
    constexpr std::size_t exact_multiplicity_limit = 4;
    if (!degrees_match || k <= exact_multiplicity_limit) {
        const bool isomorphic = graded_linalg::is_isomorphic(old_p, new_p, true);
        report << "Isomorphism: " << (isomorphic ? "match" : "DIFFER") << ".\n";
        match = match && isomorphic;
    } else {
        report << "Isomorphism: skipped (generator multiplicity " << k << " > "
               << exact_multiplicity_limit << ").\n";
    }
    report << (match ? "Checks passed.\n" : "Checks FAILED.\n");
    const fs::path path = prefix.string() + "_pruning_comparison.txt";
    auto file = output_file(path);
    file << report.str();
    std::cout << report.str() << "Output directory: " << fs::absolute(path).parent_path().string()
              << "\nSaved comparison: " << path.filename().string() << '\n';
    return new_output;
}

} // namespace stable_decomposition
