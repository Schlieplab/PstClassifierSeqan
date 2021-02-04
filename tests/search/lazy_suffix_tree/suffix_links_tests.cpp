#include <gtest/gtest.h>

#include <vector>

#include "../../../src/search/lazy_suffix_tree/construction.hpp"
#include "../../../src/search/lazy_suffix_tree/suffix_links.hpp"

using namespace lst::details;

#include <seqan3/alphabet/nucleotide/dna5.hpp>

using seqan3::operator""_dna5;

class SuffixLinksTests : public ::testing::Test {
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

  std::vector<size_t> unfinished_suffixes{1, 3, 0, 2, 4, 5};

  Table<> unfinished_table{
      {0, Flag::UNEVALUATED},
      {2, Flag::NONE},
      {2, Flag::UNEVALUATED},
      {5, Flag::NONE},
      {5, Flag(Flag::LEAF | Flag::RIGHT_MOST_CHILD)},
      {0, Flag::NONE}, // Added to allow for explicit labels in the tree
  };
};

TEST_F(SuffixLinksTests, TreeHeight) {
  size_t height = tree_height<seqan3::dna5>(sequence, suffixes, table);

  EXPECT_EQ(3, height);
}

TEST_F(SuffixLinksTests, TreeHeightParallel) {
  size_t height =
      tree_height_parallel<seqan3::dna5>(sequence, suffixes, table, 2);

  EXPECT_EQ(3, height);
}

TEST_F(SuffixLinksTests, LeafIndex) {
  size_t leaf_index = get_leaf_index(6, 0, suffixes, table);
  EXPECT_EQ(leaf_index, 5);

  size_t root_leaf_index = get_leaf_index(0, 0, suffixes, table);
  EXPECT_EQ(root_leaf_index, -1);

  size_t unfinished_leaf_index =
      get_leaf_index(0, 0, unfinished_suffixes, unfinished_table);
  EXPECT_EQ(unfinished_leaf_index, 1);
}
