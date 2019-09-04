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
    sequence = sequence_ | seqan3::view::convert<seqan3::gapped<seqan3::dna5>>;
    sequence.push_back(seqan3::gap{});
  }

  std::vector<seqan3::gapped<seqan3::dna5>> sequence{};

  std::vector<int> table{0, 2, 1, 8,  0, 12, 5, 0, 3, 0,
                         5, 0, 1, 16, 5, 0,  3, 0, 5, 0};
  std::vector<Flag> flags{
      Flag::RightMostChild,
      Flag::None,
      Flag::None,
      Flag::None,
      Flag::None,
      Flag::None,
      Flag(Flag::Leaf | Flag::RightMostChild),
      Flag::None, // Added to allow for explicit labels in the tree
      Flag::Leaf,
      Flag::None, // Added to allow for explicit labels in the tree
      Flag(Flag::Leaf | Flag::RightMostChild),
      Flag::None, // Added to allow for explicit labels in the tree
      Flag::None,
      Flag::None,
      Flag(Flag::Leaf | Flag::RightMostChild),
      Flag::None, // Added to allow for explicit labels in the tree
      Flag::Leaf,
      Flag::None, // Added to allow for explicit labels in the tree
      Flag(Flag::Leaf | Flag::RightMostChild),
      Flag::None, // Added to allow for explicit labels in the tree
  };
  std::vector<int> suffixes{3, 5, 3, 5, 5, 5};

  std::vector<int> unfinished_suffixes{1, 3, 0, 2, 4, 5};
  ;
  std::vector<int> unfinished_table{0, 2, 2, 5, 5, 0};
  std::vector<Flag> unfinished_flags{
      Flag::Unevaluated,
      Flag::None,
      Flag::Unevaluated,
      Flag::None,
      Flag(Flag::Leaf | Flag::RightMostChild),
      Flag::None, // Added to allow for explicit labels in the tree
  };
};

TEST_F(LazySuffixTreeTest, TreeHeight) {
  std::vector<int> depths(table.size() / 2, -1);
  int height = tree_height(depths, sequence, suffixes, table, flags);

  EXPECT_EQ(6, height);

  std::vector<int> expected_depths{0, 2, 1, 1, 5, 3, 3, 2, 6, 4};

  EXPECT_EQ(depths, expected_depths);
}

TEST_F(LazySuffixTreeTest, LeafIndex) {
  int leaf_index = get_leaf_index(6, 0, suffixes, table, flags);
  EXPECT_EQ(leaf_index, 5);

  int root_leaf_index = get_leaf_index(0, 0, suffixes, table, flags);
  EXPECT_EQ(root_leaf_index, -1);

  int unfinished_leaf_index = get_leaf_index(
      0, 0, unfinished_suffixes, unfinished_table, unfinished_flags);
  EXPECT_EQ(unfinished_leaf_index, 1);
}
