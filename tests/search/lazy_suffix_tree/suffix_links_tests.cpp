#include <gtest/gtest.h>

#include <vector>

#include "../../../src/search/lazy_suffix_tree.hpp"
#include "../../../src/search/lazy_suffix_tree/construction.hpp"
#include "../../../src/search/lazy_suffix_tree/suffix_links.hpp"

using namespace lst::details;

#include <seqan3/alphabet/nucleotide/dna5.hpp>

using seqan3::operator""_dna5;

class SuffixLinksTests : public ::testing::Test {
protected:
  void SetUp() override {
  }

  sequence_t<seqan3::dna5> sequence{"CACAC"_dna5};

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

TEST_F(SuffixLinksTests, TreeHeight) {
  std::vector<int> depths(table.size() / 2, -1);
  int height =
      tree_height<seqan3::dna5>(depths, sequence, suffixes, table, flags);

  EXPECT_EQ(6, height);

  std::vector<int> expected_depths{0, 2, 1, 1, 5, 3, 3, 2, 6, 4};

  EXPECT_EQ(depths, expected_depths);
}

TEST_F(SuffixLinksTests, LeafIndex) {
  int leaf_index = get_leaf_index(6, 0, suffixes, table, flags);
  EXPECT_EQ(leaf_index, 5);

  int root_leaf_index = get_leaf_index(0, 0, suffixes, table, flags);
  EXPECT_EQ(root_leaf_index, -1);

  int unfinished_leaf_index = get_leaf_index(
      0, 0, unfinished_suffixes, unfinished_table, unfinished_flags);
  EXPECT_EQ(unfinished_leaf_index, 1);
}
