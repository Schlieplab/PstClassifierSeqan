#include <benchmark/benchmark.h>
#include <random>
#include <vector>

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/alphabet/nucleotide/dna5.hpp>

#include "../../src/search/lazy_suffix_tree.hpp"

#include "../random_sequence.hpp"

using seqan3::operator""_dna5;

class LSTFixture : public benchmark::Fixture {
public:
  void SetUp(const ::benchmark::State &state) override {
    sequence = random_sequence(1000000);
  }
  std::vector<seqan3::dna5> sequence;
};

BENCHMARK_F(LSTFixture, ExpandAll)(benchmark::State &state) {
  for (auto _ : state) {
    lst::LazySuffixTree tree{sequence, false};
    tree.expand_all();
  }
}

BENCHMARK_F(LSTFixture, ExpandRoot)(benchmark::State &state) {

  for (auto _ : state) {
    state.PauseTiming();

    lst::details::Table<> table{{0, lst::details::Flag::RIGHT_MOST_CHILD},
                                {2, lst::details::Flag::NONE}};

    std::vector<size_t> suffixes(sequence.size() + 1);

    state.ResumeTiming();

    lst::details::expand_root(sequence, suffixes, table);
  }
}

BENCHMARK_F(LSTFixture, ExpandNodeLevel1)(benchmark::State &state) {
  lst::LazySuffixTree tree{sequence, false};
  auto first_level_child = tree.table[0].second;

  for (auto _ : state) {
    state.PauseTiming();

    std::vector<size_t> suffixes_copy(tree.suffixes.begin(),
                                      tree.suffixes.end());
    lst::details::Table<> table_copy = tree.table;

    state.ResumeTiming();

    lst::details::expand_node(first_level_child, tree.sequence, suffixes_copy,
                              table_copy);
  }
}

BENCHMARK_F(LSTFixture, ExpandNodeLevel2)(benchmark::State &state) {
  lst::LazySuffixTree tree{sequence, false};
  auto first_level_child = tree.table[0].second;

  lst::details::expand_node(first_level_child, tree.sequence, tree.suffixes,
                            tree.table);
  auto second_level_child = tree.table[first_level_child + 1].value;

  for (auto _ : state) {
    state.PauseTiming();

    std::vector<size_t> suffixes_copy(tree.suffixes.begin(),
                                      tree.suffixes.end());
    lst::details::Table<> table_copy = tree.table;

    state.ResumeTiming();

    lst::details::expand_node(second_level_child, tree.sequence, suffixes_copy,
                              table_copy);
  }
}

BENCHMARK_F(LSTFixture, ExpandNodeLevel3)(benchmark::State &state) {
  lst::LazySuffixTree tree{sequence, false};
  auto first_level_child = tree.table[0].second;

  lst::details::expand_node(first_level_child, tree.sequence, tree.suffixes,
                            tree.table);
  auto second_level_child = tree.table[first_level_child + 1].value;
  lst::details::expand_node(second_level_child, tree.sequence, tree.suffixes,
                            tree.table);
  auto third_level_child = tree.table[second_level_child + 1].value;

  for (auto _ : state) {
    state.PauseTiming();

    std::vector<size_t> suffixes_copy(tree.suffixes.begin(),
                                      tree.suffixes.end());
    lst::details::Table<> table_copy = tree.table;

    state.ResumeTiming();

    lst::details::expand_node(third_level_child, tree.sequence, suffixes_copy,
                              table_copy);
  }
}

BENCHMARK_F(LSTFixture, LongestCommonPrefix)(benchmark::State &state) {
  lst::LazySuffixTree tree{sequence, false};
  auto first_level_child = tree.table[0].second;
  auto lower_bound = tree.table[first_level_child].value;
  auto upper_bound = tree.table[first_level_child + 1].value;

  for (auto _ : state) {
    benchmark::DoNotOptimize(lst::details::longest_common_prefix(
        lower_bound, upper_bound, tree.sequence, tree.suffixes));
  }
}

BENCHMARK_F(LSTFixture, AddLCPToSuffixes)(benchmark::State &state) {
  lst::LazySuffixTree tree{sequence, false};
  auto first_level_child = tree.table[0].second;
  auto lower_bound = tree.table[first_level_child].value;
  auto upper_bound = tree.table[first_level_child + 1].value;

  for (auto _ : state) {
    state.PauseTiming();

    std::vector<size_t> suffixes_copy(tree.suffixes.begin(),
                                      tree.suffixes.end());

    state.ResumeTiming();

    lst::details::add_lcp_to_suffixes(lower_bound, upper_bound, 4,
                                      suffixes_copy);
  }
}

BENCHMARK_F(LSTFixture, CountSuffixes)(benchmark::State &state) {
  lst::LazySuffixTree tree{sequence, false};
  auto first_level_child = tree.table[0].second;

  auto lower_bound = tree.table[first_level_child].value;
  auto upper_bound = tree.table[first_level_child + 1].value;
  auto lcp = lst::details::longest_common_prefix(lower_bound, upper_bound,
                                                 tree.sequence, tree.suffixes);
  lst::details::add_lcp_to_suffixes(lower_bound, upper_bound, lcp,
                                    tree.suffixes);

  for (auto _ : state) {
    benchmark::DoNotOptimize(lst::details::count_suffixes(
        lower_bound, upper_bound, tree.sequence, tree.suffixes));
  }
}

BENCHMARK_F(LSTFixture, AddLcpAndCountSuffixes)(benchmark::State &state) {
  lst::LazySuffixTree tree{sequence, false};
  auto first_level_child = tree.table[0].second;

  auto lower_bound = tree.table[first_level_child].value;
  auto upper_bound = tree.table[first_level_child + 1].value;
  auto lcp = lst::details::longest_common_prefix(lower_bound, upper_bound,
                                                 tree.sequence, tree.suffixes);
  lst::details::add_lcp_to_suffixes(lower_bound, upper_bound, lcp,
                                    tree.suffixes);

  for (auto _ : state) {
    benchmark::DoNotOptimize(lst::details::add_lcp_and_count_suffixes(
        lower_bound, upper_bound, lcp, tree.sequence, tree.suffixes));
  }
}

BENCHMARK_F(LSTFixture, SortSuffixes)(benchmark::State &state) {
  lst::LazySuffixTree tree{sequence, false};
  auto first_level_child = tree.table[0].second;

  auto lower_bound = tree.table[first_level_child].value;
  auto upper_bound = tree.table[first_level_child + 1].value;
  auto lcp = lst::details::longest_common_prefix(lower_bound, upper_bound,
                                                 tree.sequence, tree.suffixes);

  auto [counts, count] = lst::details::add_lcp_and_count_suffixes(
      lower_bound, upper_bound, lcp, tree.sequence, tree.suffixes);

  for (auto _ : state) {
    state.PauseTiming();
    std::vector<size_t> suffixes_copy(tree.suffixes.begin(),
                                      tree.suffixes.end());
    state.ResumeTiming();

    lst::details::sort_suffixes(counts, lower_bound, upper_bound, tree.sequence,
                                suffixes_copy);
  }
}

BENCHMARK_F(LSTFixture, SuffixPointers)(benchmark::State &state) {
  lst::LazySuffixTree tree{sequence, false};
  auto first_level_child = tree.table[0].second;

  auto lower_bound = tree.table[first_level_child].value;
  auto upper_bound = tree.table[first_level_child + 1].value;
  auto lcp = lst::details::longest_common_prefix(lower_bound, upper_bound,
                                                 tree.sequence, tree.suffixes);
  lst::details::add_lcp_to_suffixes(lower_bound, upper_bound, lcp,
                                    tree.suffixes);

  auto [counts, count] = lst::details::count_suffixes(
      lower_bound, upper_bound, tree.sequence, tree.suffixes);

  for (auto _ : state) {
    benchmark::DoNotOptimize(
        lst::details::suffix_pointers<seqan3::dna5>(counts));
  }
}

BENCHMARK_F(LSTFixture, CopySuffixes)(benchmark::State &state) {
  lst::LazySuffixTree tree{sequence, false};
  auto first_level_child = tree.table[0].second;

  auto lower_bound = tree.table[first_level_child].value;
  auto upper_bound = tree.table[first_level_child + 1].value;

  for (auto _ : state) {
    benchmark::DoNotOptimize(
        std::vector<size_t>(tree.suffixes.begin() + lower_bound,
                            tree.suffixes.begin() + upper_bound));
  }
}

BENCHMARK_F(LSTFixture, AddChildren)(benchmark::State &state) {
  lst::LazySuffixTree tree{sequence, false};
  auto first_level_child = tree.table[0].second;

  auto lower_bound = tree.table[first_level_child].value;
  auto upper_bound = tree.table[first_level_child + 1].value;
  auto lcp = lst::details::longest_common_prefix(lower_bound, upper_bound,
                                                 tree.sequence, tree.suffixes);
  lst::details::add_lcp_to_suffixes(lower_bound, upper_bound, lcp,
                                    tree.suffixes);

  auto [counts, count] = lst::details::count_suffixes(
      lower_bound, upper_bound, tree.sequence, tree.suffixes);

  lst::details::sort_suffixes(counts, lower_bound, upper_bound, tree.sequence,
                              tree.suffixes);

  for (auto _ : state) {
    state.PauseTiming();
    lst::details::Table<> table_copy = tree.table;
    state.ResumeTiming();
    lst::details::add_children<seqan3::dna5>(counts, lower_bound, tree.suffixes,
                                             table_copy);
  }
}

BENCHMARK_F(LSTFixture, GetCharacterRank)(benchmark::State &state) {
  for (auto _ : state) {
    benchmark::DoNotOptimize(lst::details::get_character_rank(sequence, 5));
  }
}

BENCHMARK_F(LSTFixture, GetEdgeLcpExpanded)(benchmark::State &state) {
  lst::LazySuffixTree tree{sequence, false};
  auto first_level_child = tree.table[0].second;

  lst::details::expand_node(first_level_child, tree.sequence, tree.suffixes,
                            tree.table);

  for (auto _ : state) {
    benchmark::DoNotOptimize(lst::details::get_edge_lcp(
        first_level_child, tree.sequence, tree.suffixes, tree.table));
  }
}

BENCHMARK_F(LSTFixture, GetEdgeLcpNotExpanded)(benchmark::State &state) {
  lst::LazySuffixTree tree{sequence, false};
  auto first_level_child = tree.table[0].second;

  for (auto _ : state) {
    benchmark::DoNotOptimize(lst::details::get_edge_lcp(
        first_level_child, tree.sequence, tree.suffixes, tree.table));
  }
}

BENCHMARK_F(LSTFixture, IteratingOverSuffixes)(benchmark::State &state) {
  lst::LazySuffixTree tree{sequence, false};
  auto first_level_child = tree.table[0].second;
  auto lower_bound = tree.table[first_level_child].value;
  auto upper_bound = tree.table[first_level_child + 1].value;

  for (auto _ : state) {
    for (auto i = lower_bound; i < upper_bound; i++) {
      benchmark::DoNotOptimize(tree.suffixes[i]);
    }
  }
}

BENCHMARK_F(LSTFixture, NodeOccurrences)(benchmark::State &state) {
  lst::LazySuffixTree tree{sequence, false};
  auto first_level_child = tree.table[0].second;
  lst::details::expand_node(first_level_child, tree.sequence, tree.suffixes,
                            tree.table);

  for (auto _ : state) {
    benchmark::DoNotOptimize(node_occurrences(first_level_child, tree.table,
                                              tree.sequence, tree.suffixes));
  }
}

BENCHMARK_MAIN();