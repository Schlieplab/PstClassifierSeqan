#include <benchmark/benchmark.h>

#include "../src/distances/cv.hpp"
#include "../src/kl_tree_map.hpp"

static void CV(benchmark::State &state) {
  std::filesystem::path first_path{"./../trees/CM008035.1.tree"};
  std::filesystem::path second_path{"./../trees/NC_009067.tree"};
  pst::KullbackLieblerTreeMap<seqan3::dna5> first{first_path};
  pst::KullbackLieblerTreeMap<seqan3::dna5> second{second_path};

  for (auto _ : state) {
    pst::cv(first, second);
  }
}
BENCHMARK(CV);

BENCHMARK_MAIN();