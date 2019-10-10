#include <gtest/gtest.h>

#include "../../src/search/lazy_suffix_tree.hpp"
#include "../../src/search/lazy_suffix_tree/construction.hpp"
#include "../../src/search/lazy_suffix_tree/suffix_links.hpp"

using namespace lst::details;

#include <seqan3/alphabet/all.hpp>
using seqan3::operator""_dna5;

class LazySuffixTreeTest : public ::testing::Test {
protected:
  void SetUp() override {
    std::vector<seqan3::dna5> sequence_{"CACAC"_dna5};
    sequence = sequence_t<seqan3::dna5>{
        sequence_ | seqan3::view::convert<seqan3::gapped<seqan3::dna5>>};
    sequence.push_back(seqan3::gap{});
  }

  sequence_t<seqan3::dna5> sequence{};

  std::vector<int> table{0, 2, 1, 8,  0, 12, 5, 0, 3, 0,
                         5, 0, 1, 16, 5, 0,  3, 0, 5, 0};
  std::vector<Flag> flags{
      Flag::RIGHT_MOST_CHILD,
      Flag::NONE,
      Flag::NONE,
      Flag::NONE,
      Flag::NONE,
      Flag::NONE,
      Flag(Flag::LEAF | Flag::RIGHT_MOST_CHILD),
      Flag::NONE, // Added to allow for explicit labels in the tree
      Flag::LEAF,
      Flag::NONE, // Added to allow for explicit labels in the tree
      Flag(Flag::LEAF | Flag::RIGHT_MOST_CHILD),
      Flag::NONE, // Added to allow for explicit labels in the tree
      Flag::NONE,
      Flag::NONE,
      Flag(Flag::LEAF | Flag::RIGHT_MOST_CHILD),
      Flag::NONE, // Added to allow for explicit labels in the tree
      Flag::LEAF,
      Flag::NONE, // Added to allow for explicit labels in the tree
      Flag(Flag::LEAF | Flag::RIGHT_MOST_CHILD),
      Flag::NONE, // Added to allow for explicit labels in the tree
  };
  std::vector<int> suffixes{3, 5, 3, 5, 5, 5};

  std::vector<int> unfinished_suffixes{1, 3, 0, 2, 4, 5};
  ;
  std::vector<int> unfinished_table{0, 2, 2, 5, 5, 0};
  std::vector<Flag> unfinished_flags{
      Flag::UNEVALUATED,
      Flag::NONE,
      Flag::UNEVALUATED,
      Flag::NONE,
      Flag(Flag::LEAF | Flag::RIGHT_MOST_CHILD),
      Flag::NONE, // Added to allow for explicit labels in the tree
  };
};

TEST_F(LazySuffixTreeTest, GetEdgeLCP) {
  int root_edge_lcp = get_edge_lcp(0, sequence, suffixes, table, flags);
  EXPECT_EQ(root_edge_lcp, 0);

  int ac_edge_lcp = get_edge_lcp(2, sequence, suffixes, table, flags);
  EXPECT_EQ(ac_edge_lcp, 2);
}

TEST_F(LazySuffixTreeTest, NodeOccurrences) {
  int root_occurrences = node_occurrences(0, table, flags);
  EXPECT_EQ(root_occurrences, 6);

  int ac_occurrences = node_occurrences(2, table, flags);
  EXPECT_EQ(ac_occurrences, 2);
}

TEST_F(LazySuffixTreeTest, IterateChildren) {
  std::vector<int> visited{};

  iterate_children(0, table, flags,
                   [&](int index) { visited.push_back(index); });

  std::vector<int> expected_visited{2, 4, 6};

  EXPECT_EQ(visited, expected_visited);
}

TEST_F(LazySuffixTreeTest, BreadthFirstIteration) {
  std::vector<int> visited{};
  breadth_first_iteration(sequence, suffixes, table, flags,
                          [&](int index, int lcp, int edge_lcp) -> bool {
                            visited.push_back(index);
                            return true;
                          });

  std::vector<int> expected_visited{2, 4, 6, 8, 10, 12, 14, 16, 18};

  EXPECT_EQ(visited, expected_visited);
}
