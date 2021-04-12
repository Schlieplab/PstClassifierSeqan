#include <benchmark/benchmark.h>
#include <random>
#include <vector>

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/alphabet/nucleotide/dna5.hpp>

#include "../../src/kl_tree.hpp"
#include "../../src/kl_tree_map.hpp"

#include "random_sequence.hpp"

using seqan3::operator""_dna5;

using counts_map = robin_hood::unordered_map<
    std::string,
    std::tuple<size_t, std::array<double, seqan3::alphabet_size<seqan3::dna5>>,
               bool>>;

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
  counts_map counts{tree.counts};

  for (auto _ : state) {
    state.PauseTiming();
    counts_map counts_copy{counts};
    tree.counts = counts;

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
  counts_map counts{tree.counts};

  for (auto _ : state) {
    state.PauseTiming();
    counts_map counts_copy{counts};
    tree.counts = counts;
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

BENCHMARK_F(KLTreeixture, TreeExplicitSuffixLinks)(benchmark::State &state) {
  pst::KullbackLieblerTree tree{"test", sequence, 15,    10, 1.2,
                                0,      "cutoff", false, 0};
  tree.support_pruning();

  for (auto _ : state) {
    state.PauseTiming();
    tree.suffix_links.resize(tree.table.size() / 2, (size_t)-1);
    state.ResumeTiming();

    lst::details::add_explicit_suffix_links(
        tree.sequence, tree.suffixes, tree.table, tree.suffix_links, false, 0);
  }
}

BENCHMARK_F(KLTreeixture, TreeLeafSuffixLinks)(benchmark::State &state) {
  pst::KullbackLieblerTree tree{"test", sequence, 15,    10, 1.2,
                                0,      "cutoff", false, 0};
  tree.support_pruning();

  for (auto _ : state) {
    state.PauseTiming();
    tree.suffix_links.resize(tree.table.size() / 2, (size_t)-1);
    lst::details::add_explicit_suffix_links(
        tree.sequence, tree.suffixes, tree.table, tree.suffix_links, false, 0);
    state.ResumeTiming();

    lst::details::add_leaf_suffix_links(
        tree.sequence, tree.suffixes, tree.table, tree.suffix_links, false, 0);
  }
}

BENCHMARK_F(KLTreeixture, TreeImplicitSuffixLinks)(benchmark::State &state) {
  pst::KullbackLieblerTree tree{"test", sequence, 15,    10, 1.2,
                                0,      "cutoff", false, 0};
  tree.support_pruning();
  for (auto _ : state) {
    state.PauseTiming();

    tree.suffix_links.resize(tree.table.size() / 2, (size_t)-1);

    lst::details::add_explicit_suffix_links(
        tree.sequence, tree.suffixes, tree.table, tree.suffix_links, false, 0);
    lst::details::add_leaf_suffix_links(
        tree.sequence, tree.suffixes, tree.table, tree.suffix_links, false, 0);

    state.ResumeTiming();

    lst::details::add_implicit_suffix_links(
        tree.sequence, tree.suffixes, tree.table, tree.suffix_links, false, 0);
  }
}

BENCHMARK_F(KLTreeixture, FindSuffixMatch)(benchmark::State &state) {
  pst::KullbackLieblerTree tree{"test", sequence, 15,    10, 1.2,
                                0,      "cutoff", false, 0};

  tree.suffix_links.resize(tree.table.size() / 2, (size_t)-1);

  lst::details::add_explicit_suffix_links(
      tree.sequence, tree.suffixes, tree.table, tree.suffix_links, false, 0);
  lst::details::add_leaf_suffix_links(tree.sequence, tree.suffixes, tree.table,
                                      tree.suffix_links, false, 0);

  for (auto _ : state) {
    find_suffix_match(2, 1, 0, tree.sequence, tree.suffixes, tree.table);
  }
}

BENCHMARK_F(KLTreeixture, TreeSupportPruningParallel)(benchmark::State &state) {
  for (auto _ : state) {
    state.PauseTiming();
    pst::KullbackLieblerTree tree{"test", sequence, 15,   10, 1.2,
                                  0,      "cutoff", true, 1};
    state.ResumeTiming();
    tree.support_pruning();
  }
}

BENCHMARK_F(KLTreeixture, TreeSimilarityPruningParallel)
(benchmark::State &state) {
  pst::KullbackLieblerTree tree{"test", sequence, 15,   10, 1.2,
                                0,      "cutoff", true, 1};
  tree.support_pruning();
  tree.add_suffix_links();
  tree.add_reverse_suffix_links();

  for (auto _ : state) {
    tree.similarity_pruning();
  }
}

BENCHMARK_F(KLTreeixture, TreeSuffixLinksParallel)(benchmark::State &state) {
  pst::KullbackLieblerTree tree{"test", sequence, 15,   10, 1.2,
                                0,      "cutoff", true, 1};
  tree.support_pruning();
  for (auto _ : state) {
    tree.add_suffix_links();
    tree.add_reverse_suffix_links();
  }
}

BENCHMARK_F(KLTreeixture, ToTree)(benchmark::State &state) {
  pst::KullbackLieblerTree tree{"test", sequence, 15,   10, 1.2,
                                0,      "cutoff", true, 1};
  tree.construct_tree();
  for (auto _ : state) {
    benchmark::DoNotOptimize(tree.to_tree());
  }
}

BENCHMARK_F(KLTreeixture, ToTreeMap)(benchmark::State &state) {
  pst::KullbackLieblerTreeMap tree{"test", sequence, 15,   10, 1.2,
                                   0,      "cutoff", true, 1};
  tree.construct_tree();
  for (auto _ : state) {
    benchmark::DoNotOptimize(tree.to_tree());
  }
}

BENCHMARK_F(KLTreeixture, BreadthFirstIterationP)(benchmark::State &state) {
  pst::KullbackLieblerTreeMap tree{"test", sequence, 15,   10, 1.2,
                                   0,      "cutoff", true, 1};
  tree.construct_tree();
  for (auto _ : state) {
    tree.breadth_first_iteration_p([&](const std::string &label, size_t level) {
      tree.assign_node_probabilities(label);
      return true;
    });
  }
}

BENCHMARK_F(KLTreeixture, BreadthFirstIterationPNoWork)
(benchmark::State &state) {
  pst::KullbackLieblerTreeMap tree{"test", sequence, 15,   10, 1.2,
                                   0,      "cutoff", true, 1};
  tree.construct_tree();
  for (auto _ : state) {
    tree.breadth_first_iteration_p(
        [&](const std::string &label, size_t level) { return true; });
  }
}

BENCHMARK_F(KLTreeixture, BreadthAssignNodeProbabilities)
(benchmark::State &state) {
  pst::KullbackLieblerTreeMap tree{"test", sequence, 15,   10, 1.2,
                                   0,      "cutoff", true, 1};
  tree.construct_tree();
  for (auto _ : state) {
    tree.assign_node_probabilities("");
  }
}

BENCHMARK_MAIN();
