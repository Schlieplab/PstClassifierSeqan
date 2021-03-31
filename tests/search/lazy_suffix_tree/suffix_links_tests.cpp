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

  std::vector<size_t> unfinished_suffixes{1, 3, 0, 2, 4, 5};

  Table<> unfinished_table{
      {0, 1, Flag::UNEVALUATED},
      {2, 5, Flag::UNEVALUATED},
      {5, 0, Flag(Flag::LEAF | Flag::RIGHT_MOST_CHILD)},
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
  size_t leaf_index = get_leaf_index(3, 0, suffixes, table);
  EXPECT_EQ(leaf_index, 5);

  size_t root_leaf_index = get_leaf_index(0, 0, suffixes, table);
  EXPECT_EQ(root_leaf_index, -1);

  size_t unfinished_leaf_index =
      get_leaf_index(0, 0, unfinished_suffixes, unfinished_table);
  EXPECT_EQ(unfinished_leaf_index, 1);
}
