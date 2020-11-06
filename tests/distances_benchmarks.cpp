#include <benchmark/benchmark.h>

#include "../src/distances/cv.hpp"
#include "../src/kl_tree_map.hpp"

static void CV(benchmark::State &state) {
  std::filesystem::path first_path{"./../trees/CM008035.1.tree"};
  std::filesystem::path second_path{"./../trees/NC_009067.tree"};
  pst::KullbackLieblerTreeMap<seqan3::dna5> first{first_path};
  pst::KullbackLieblerTreeMap<seqan3::dna5> second{second_path};

  for (auto _ : state) {
    pst::distances::cv(first, second);
  }
}
BENCHMARK(CV);

static void SharedContexts(benchmark::State &state) {
  std::filesystem::path first_path{"./../trees/CM008035.1.tree"};
  std::filesystem::path second_path{"./../trees/NC_009067.tree"};
  pst::KullbackLieblerTreeMap<seqan3::dna5> first{first_path};
  pst::KullbackLieblerTreeMap<seqan3::dna5> second{second_path};

  for (auto _ : state) {
    pst::distances::details::get_shared_contexts(first, second);
  }
}
BENCHMARK(SharedContexts);

static void CompositionVector(benchmark::State &state) {
  std::filesystem::path first_path{"./../trees/CM008035.1.tree"};
  std::filesystem::path second_path{"./../trees/NC_009067.tree"};
  pst::KullbackLieblerTreeMap<seqan3::dna5> first{first_path};
  pst::KullbackLieblerTreeMap<seqan3::dna5> second{second_path};

  auto contexts = pst::distances::details::get_shared_contexts(first, second);

  for (auto _ : state) {
    pst::distances::details::composition_vector(first, contexts, 2);
  }
}
BENCHMARK(CompositionVector);

static void GetTerminalNodes(benchmark::State &state) {
  std::filesystem::path first_path{"./../trees/CM008035.1.tree"};
  pst::KullbackLieblerTreeMap<seqan3::dna5> first{first_path};

  for (auto _ : state) {
    first.get_terminal_nodes();
  }
}
BENCHMARK(GetTerminalNodes);

static void CountTerminalNodes(benchmark::State &state) {
  std::filesystem::path first_path{"./../trees/CM008035.1.tree"};
  pst::KullbackLieblerTreeMap<seqan3::dna5> first{first_path};

  for (auto _ : state) {
    first.count_terminal_nodes();
  }
}
BENCHMARK(CountTerminalNodes);

static void GetTransitionProbability(benchmark::State &state) {
  std::filesystem::path first_path{"./../trees/CM008035.1.tree"};
  pst::KullbackLieblerTreeMap<seqan3::dna5> first{first_path};

  std::string context{"GTGTGT"};
  int char_rank = 4;
  for (auto _ : state) {
    first.get_transition_probability(context, char_rank);
  }
}
BENCHMARK(GetTransitionProbability);

static void GetClosestState(benchmark::State &state) {
  std::filesystem::path first_path{"./../trees/CM008035.1.tree"};
  pst::KullbackLieblerTreeMap<seqan3::dna5> first{first_path};

  std::string context{"ACGGTGTGT"};
  for (auto _ : state) {
    first.get_closest_state(context);
  }
}
BENCHMARK(GetClosestState);

BENCHMARK_MAIN();