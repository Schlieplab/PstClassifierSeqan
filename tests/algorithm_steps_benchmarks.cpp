#include <benchmark/benchmark.h>
#include <random>
#include <vector>

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/alphabet/nucleotide/dna5.hpp>

#include "../../src/kl_tree.hpp"
#include "../../src/kl_tree_map.hpp"

#include "random_sequence.hpp"

using seqan3::operator""_dna5;

class KLTreeixture : public benchmark::Fixture {
public:
  void SetUp(const ::benchmark::State &state) override {
    sequence = random_sequence(1000000);
  }
  std::vector<seqan3::dna5> sequence;
};

BENCHMARK_F(KLTreeixture, MapSupportPruning)(benchmark::State &state) {
  for (auto _ : state) {
    state.PauseTiming();
    pst::KullbackLieblerTreeMap tree{"test", sequence, 15,    10, 1.2,
                                     0,      "cutoff", false, 0};
    state.ResumeTiming();

    tree.support_pruning();
  }
}

BENCHMARK_F(KLTreeixture, MapSimilarityPruning)(benchmark::State &state) {
  pst::KullbackLieblerTreeMap tree{"test", sequence, 15,    10, 1.2,
                                   0,      "cutoff", false, 0};
  tree.support_pruning();
  robin_hood::unordered_set<std::string> status{tree.status};

  for (auto _ : state) {
    state.PauseTiming();
    robin_hood::unordered_set<std::string> status_copy{status};
    tree.status = status_copy;

    state.ResumeTiming();

    tree.similarity_pruning();
  }
}

BENCHMARK_F(KLTreeixture, MapSupportPruningParallel)(benchmark::State &state) {
  for (auto _ : state) {
    state.PauseTiming();
    pst::KullbackLieblerTreeMap tree{"test", sequence, 15,   10, 1.2,
                                     0,      "cutoff", true, 2};
    state.ResumeTiming();

    tree.support_pruning();
  }
}

BENCHMARK_F(KLTreeixture, MapSimilarityPruningParallel)
(benchmark::State &state) {
  pst::KullbackLieblerTreeMap tree{"test", sequence, 15,   10, 1.2,
                                   0,      "cutoff", true, 2};
  tree.support_pruning();
  robin_hood::unordered_set<std::string> status{tree.status};

  for (auto _ : state) {
    state.PauseTiming();
    robin_hood::unordered_set<std::string> status_copy{status};
    tree.status = status_copy;
    state.ResumeTiming();

    tree.similarity_pruning();
  }
}

BENCHMARK_F(KLTreeixture, TreeGetPstLeaves)(benchmark::State &state) {
  pst::KullbackLieblerTreeMap tree{"test", sequence, 15,    10, 1.2,
                                   0,      "cutoff", false, 0};
  tree.support_pruning();
  for (auto _ : state) {
    benchmark::DoNotOptimize(tree.get_pst_leaves());
  }
}

BENCHMARK_F(KLTreeixture, TreeSupportPruning)(benchmark::State &state) {
  for (auto _ : state) {
    state.PauseTiming();
    pst::KullbackLieblerTree tree{"test", sequence, 15,    10, 1.2,
                                  0,      "cutoff", false, 0};
    state.ResumeTiming();
    tree.support_pruning();
  }
}

BENCHMARK_F(KLTreeixture, TreeSimilarityPruning)(benchmark::State &state) {
  for (auto _ : state) {
    state.PauseTiming();
    pst::KullbackLieblerTree tree{"test", sequence, 15,    10, 1.2,
                                  0,      "cutoff", false, 0};
    tree.support_pruning();
    tree.add_suffix_links();
    tree.add_reverse_suffix_links();
    state.ResumeTiming();

    tree.similarity_pruning();
  }
}

BENCHMARK_F(KLTreeixture, TreeSuffixLinks)(benchmark::State &state) {
  for (auto _ : state) {
    state.PauseTiming();
    pst::KullbackLieblerTree tree{"test", sequence, 15,    10, 1.2,
                                  0,      "cutoff", false, 0};
    tree.support_pruning();
    state.ResumeTiming();

    tree.add_suffix_links();
    tree.add_reverse_suffix_links();
  }
}

BENCHMARK_F(KLTreeixture, TreeSupportPruningParallel)(benchmark::State &state) {
  for (auto _ : state) {
    state.PauseTiming();
    pst::KullbackLieblerTree tree{"test", sequence, 15,   10, 1.2,
                                  0,      "cutoff", true, 2};
    state.ResumeTiming();
    tree.support_pruning();
  }
}

BENCHMARK_F(KLTreeixture, TreeSimilarityPruningParallel)
(benchmark::State &state) {
  for (auto _ : state) {
    state.PauseTiming();
    pst::KullbackLieblerTree tree{"test", sequence, 15,   10, 1.2,
                                  0,      "cutoff", true, 2};
    tree.support_pruning();
    tree.add_suffix_links();
    tree.add_reverse_suffix_links();
    state.ResumeTiming();

    tree.similarity_pruning();
  }
}

BENCHMARK_F(KLTreeixture, TreeSuffixLinksParallel)(benchmark::State &state) {
  for (auto _ : state) {
    state.PauseTiming();
    pst::KullbackLieblerTree tree{"test", sequence, 15,   10, 1.2,
                                  0,      "cutoff", true, 2};
    tree.support_pruning();
    state.ResumeTiming();

    tree.add_suffix_links();
    tree.add_reverse_suffix_links();
  }
}

BENCHMARK_MAIN();