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
      {0, Flag::RIGHT_MOST_CHILD},
      {2, Flag::NONE},
      {1, Flag::NONE},
      {8, Flag::NONE},
      {0, Flag::NONE},
      {12, Flag::NONE},
      {5, Flag(Flag::LEAF | Flag::RIGHT_MOST_CHILD)},
      {0, Flag::NONE}, // Added to allow for explicit labels in the tree
      {3, Flag::LEAF},
      {0, Flag::NONE}, // Added to allow for explicit labels in the tree
      {5, Flag(Flag::LEAF | Flag::RIGHT_MOST_CHILD)},
      {0, Flag::NONE}, // Added to allow for explicit labels in the tree
      {1, Flag::NONE},
      {16, Flag::NONE},
      {5, Flag(Flag::LEAF | Flag::RIGHT_MOST_CHILD)},
      {0, Flag::NONE}, // Added to allow for explicit labels in the tree
      {3, Flag::LEAF},
      {0, Flag::NONE}, // Added to allow for explicit labels in the tree
      {5, Flag(Flag::LEAF | Flag::RIGHT_MOST_CHILD)},
      {0, Flag::NONE}, // Added to allow for explicit labels in the tree
  };
  std::vector<size_t> suffixes{3, 5, 3, 5, 5, 5};
};

TEST_F(IterationTests, GetEdgeLCP) {
  size_t root_edge_lcp = get_edge_lcp(0, sequence, suffixes, table);
  EXPECT_EQ(root_edge_lcp, 0);

  size_t ac_edge_lcp = get_edge_lcp(2, sequence, suffixes, table);
  EXPECT_EQ(ac_edge_lcp, 2);
}

TEST_F(IterationTests, NodeOccurrences) {
  size_t root_occurrences = node_occurrences(0, table, sequence, suffixes);
  EXPECT_EQ(root_occurrences, 6);

  size_t ac_occurrences = node_occurrences(2, table, sequence, suffixes);
  EXPECT_EQ(ac_occurrences, 2);
}

TEST_F(IterationTests, IterateChildren) {
  std::vector<size_t> visited{};

  iterate_children(0, table, [&](size_t index) { visited.push_back(index); });

  std::vector<size_t> expected_visited{2, 4, 6};

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

  std::vector<size_t> expected_visited{2, 4, 6, 8, 10, 12, 14, 16, 18};

  EXPECT_EQ(visited, expected_visited);
}

TEST_F(IterationTests, BreadthFirstIterationParallel) {
  std::set<size_t> expected_visited{2, 4, 6, 8, 10, 12, 14, 16, 18};
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
