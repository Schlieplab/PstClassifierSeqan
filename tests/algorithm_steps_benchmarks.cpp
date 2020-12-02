#include <benchmark/benchmark.h>
#include <random>
#include <vector>

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/alphabet/nucleotide/dna5.hpp>

#include "../../src/kl_tree.hpp"
#include "../../src/kl_tree_map.hpp"

using seqan3::operator""_dna5;

class KLTreeixture : public benchmark::Fixture {
public:
  void SetUp(const ::benchmark::State &state) override {
    std::random_device rd;
    gen = std::mt19937{rd()};
    sequence = random_sequence(1000000);
  }
  std::vector<seqan3::dna5> sequence;
  std::mt19937 gen;

  std::vector<seqan3::dna5> random_sequence(int length) {
    std::vector<seqan3::dna5> local_sequence{"ACGT"_dna5};
    std::uniform_int_distribution<> distrib(0, 3);
    for (int i = 0; i < length; i++) {
      auto rand = distrib(this->gen);
      auto c_dna4 = seqan3::assign_rank_to(rand, seqan3::dna4{});
      auto c_dna5 = seqan3::assign_char_to(c_dna4.to_char(), seqan3::dna5{});
      local_sequence.push_back(c_dna5);
    }

    return local_sequence;
  }
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