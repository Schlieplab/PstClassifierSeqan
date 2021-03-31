#include <gtest/gtest.h>

#include "../../../src/search/lazy_suffix_tree/construction.hpp"
#include "../../../src/search/lazy_suffix_tree/iteration.hpp"
#include <vector>

using namespace lst::details;

#include <seqan3/alphabet/nucleotide/dna5.hpp>
using seqan3::operator""_dna5;

class IterationTests : public ::testing::Test {
protected:
  void SetUp() override {}
  sequence_t<seqan3::dna5> sequence{"CACAC"_dna5};

  Table<> table{
      {0, 1, Flag::RIGHT_MOST_CHILD},
      {1, 4, Flag::NONE},
      {0, 6, Flag::NONE},
      {5, 0, Flag(Flag::LEAF | Flag::RIGHT_MOST_CHILD)},
      {3, 0, Flag::LEAF},
      {5, 0, Flag(Flag::LEAF | Flag::RIGHT_MOST_CHILD)},
      {1, 8, Flag::NONE},
      {5, 0, Flag(Flag::LEAF | Flag::RIGHT_MOST_CHILD)},
      {3, 0, Flag::LEAF},
      {5, 0, Flag(Flag::LEAF | Flag::RIGHT_MOST_CHILD)},
  };
  std::vector<size_t> suffixes{3, 5, 3, 5, 5, 5};
};

TEST_F(IterationTests, GetEdgeLCP) {
  size_t root_edge_lcp = get_edge_lcp(0, sequence, suffixes, table);
  EXPECT_EQ(root_edge_lcp, 0);

  size_t ac_edge_lcp = get_edge_lcp(2, sequence, suffixes, table);
  EXPECT_EQ(ac_edge_lcp, 1);
}

TEST_F(IterationTests, NodeOccurrences) {
  size_t root_occurrences = node_occurrences(0, table, sequence, suffixes);
  EXPECT_EQ(root_occurrences, 6);

  size_t ac_occurrences = node_occurrences(1, table, sequence, suffixes);
  EXPECT_EQ(ac_occurrences, 2);
}

TEST_F(IterationTests, IterateChildren) {
  std::vector<size_t> visited{};

  iterate_children(0, table, [&](size_t index) { visited.push_back(index); });

  std::vector<size_t> expected_visited{1, 2, 3};

  EXPECT_EQ(visited, expected_visited);
}

TEST_F(IterationTests, BreadthFirstIteration) {
  std::vector<size_t> visited{};
  breadth_first_iteration(
      sequence, suffixes, table, true,
      [&](size_t index, size_t lcp, size_t edge_lcp, size_t node_count,
          lst::details::alphabet_array<size_t, seqan3::dna5> &child_counts)
          -> bool {
        visited.push_back(index);
        return true;
      },
      [](size_t n, size_t l, size_t &e) {});

  std::vector<size_t> expected_visited{1, 2, 3, 4, 5, 6, 7, 8, 9};

  EXPECT_EQ(visited, expected_visited);
}

TEST_F(IterationTests, BreadthFirstIterationParallel) {
  std::set<size_t> expected_visited{1, 2, 3, 4, 5, 6, 7, 8, 9};
  for (size_t i = 0; i < 6; i++) {
    std::set<size_t> visited{};
    breadth_first_iteration_parallel(
        sequence, suffixes, table, true,
        [&](size_t index, size_t lcp, size_t edge_lcp, size_t node_count,
            lst::details::alphabet_array<size_t, seqan3::dna5> &child_counts)
            -> bool {
          visited.insert(index);
          return true;
        },
        []() {}, 0, [](size_t n, size_t l, size_t &e) {});
    EXPECT_EQ(visited, expected_visited);
  }
}
