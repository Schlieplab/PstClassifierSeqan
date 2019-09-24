#include <gtest/gtest.h>

#include "../../src/search/lazy_suffix_tree/construction.hpp"
using namespace lst::details;

#include <seqan3/alphabet/nucleotide/dna5.hpp>
using seqan3::operator""_dna5;
#include <seqan3/alphabet/nucleotide/dna4.hpp>
using seqan3::operator""_dna4;
#include <seqan3/range/view/convert.hpp>

#include <seqan3/core/debug_stream.hpp>

class LazySuffixTreeTest : public ::testing::Test {
protected:
  void SetUp() override {
    std::vector<seqan3::dna5> sequence_{"CACAC"_dna5};
    sequence = sequence_t<seqan3::dna5>{
        sequence_ | seqan3::view::convert<seqan3::gapped<seqan3::dna5>>};
    sequence.push_back(seqan3::gap{});
  }

  sequence_t<seqan3::dna5> sequence{};
  std::vector<int> suffixes{0, 1, 2, 3, 4, 5};
  std::vector<int> table{0, 2, 2, 5, 5};
  std::vector<Flag> flags{Flag::Unevaluated, Flag::None, Flag::Unevaluated,
                          Flag::None, Flag(Flag::Leaf | Flag::RightMostChild)};
};

TEST_F(LazySuffixTreeTest, CountSuffixes) {
  std::array<int, 6> expected_counts{2, 3, 0, 0, 0, 1};
  EXPECT_EQ(count_suffixes(0, 6, sequence, suffixes), expected_counts);

  EXPECT_THROW(count_suffixes(0, 16, sequence, suffixes),
               std::invalid_argument);
  EXPECT_THROW(count_suffixes(-1, 5, sequence, suffixes),
               std::invalid_argument);
}

TEST_F(LazySuffixTreeTest, LongestCommonPrefix) {
  table = std::vector<int>{0, 2, 2, 5, 5};
  suffixes = std::vector<int>{1, 3, 0, 2, 4, 5};
  EXPECT_EQ(longest_common_prefix(1, 2, sequence, suffixes), 2);
  EXPECT_EQ(longest_common_prefix(2, 4, sequence, suffixes), 3);
  EXPECT_EQ(longest_common_prefix(3, 5, sequence, suffixes), 1);
}

TEST_F(LazySuffixTreeTest, SuffixPointers) {
  std::array<int, 6> counts{2, 3, 0, 0, 0, 1};
  auto actual_pointers = suffix_pointers<seqan3::dna5>(counts);

  std::array<int, 6> expected_pointers{0, 2, 5, 5, 5, 5};
  EXPECT_EQ(actual_pointers, expected_pointers);
}

TEST_F(LazySuffixTreeTest, SortSuffixes) {
  int lower_bound = 0;
  int upper_bound = 6;
  std::array<int, 6> counts =
      count_suffixes(lower_bound, upper_bound, sequence, suffixes);
  sort_suffixes(counts, lower_bound, upper_bound, sequence, suffixes);

  std::vector<int> expected_suffixes{1, 3, 0, 2, 4, 5};
  EXPECT_EQ(suffixes, expected_suffixes);
}

TEST_F(LazySuffixTreeTest, AddLcpToSuffixes) {
  int lower_bound = 0;
  int upper_bound = 6;
  auto counts = count_suffixes(lower_bound, upper_bound, sequence, suffixes);
  sort_suffixes(counts, lower_bound, upper_bound, sequence, suffixes);

  int lcp = longest_common_prefix(lower_bound, upper_bound, sequence, suffixes);
  add_lcp_to_suffixes(lower_bound, upper_bound, lcp, suffixes);

  std::vector<int> expected_suffixes{1, 3, 0, 2, 4, 5};
  EXPECT_EQ(suffixes, expected_suffixes);

  int lcp_0_2 = longest_common_prefix(0, 2, sequence, suffixes);
  add_lcp_to_suffixes(0, 2, lcp_0_2, suffixes);
  expected_suffixes = {3, 5, 0, 2, 4, 5};
  EXPECT_EQ(suffixes, expected_suffixes);
}

TEST_F(LazySuffixTreeTest, AddChildren) {
  table = std::vector<int>{};
  flags = std::vector<Flag>{};
  std::array<int, 6> counts = count_suffixes(0, 6, sequence, suffixes);
  add_children<seqan3::dna5>(counts, 0, suffixes, table, flags);

  // Changed to allow for explicit labels in the tree
  std::vector<int> expected_table{0, 2, 2, 5, 5, 0};
  EXPECT_EQ(table, expected_table);

  std::vector<Flag> expected_flags{
      Flag::Unevaluated,
      Flag::None,
      Flag::Unevaluated,
      Flag::None,
      Flag(Flag::Leaf | Flag::RightMostChild),
      Flag::None, // Added to allow for explicit labels in the tree
  };
  EXPECT_EQ(flags, expected_flags);
}

TEST_F(LazySuffixTreeTest, ExpandRoot) {
  table = std::vector<int>{};
  flags = std::vector<Flag>{};
  expand_root(sequence, suffixes, table, flags);

  // Changed to allow for explicit labels in the tree
  std::vector<int> expected_table{0, 2, 2, 5, 5, 0};
  EXPECT_EQ(table, expected_table);

  std::vector<Flag> expected_flags{
      Flag::Unevaluated,
      Flag::None,
      Flag::Unevaluated,
      Flag::None,
      Flag(Flag::Leaf | Flag::RightMostChild),
      Flag::None, // Added to allow for explicit labels in the tree
  };
  EXPECT_EQ(flags, expected_flags);

  std::vector<int> expected_suffixes{1, 3, 0, 2, 4, 5};
  EXPECT_EQ(suffixes, expected_suffixes);
}

TEST_F(LazySuffixTreeTest, ExpandTest) {
  table = std::vector<int>{};
  flags = std::vector<Flag>{};
  expand_root(sequence, suffixes, table, flags);
  expand_node(0, sequence, suffixes, table, flags);

  // Changed to allow for explicit labels in the tree
  std::vector<int> expected_table{1, 6, 2, 5, 5, 0, 3, 0, 5, 0};

  std::vector<Flag> expected_flags{
      Flag::None,
      Flag::None,
      Flag::Unevaluated,
      Flag::None,
      Flag(Flag::Leaf | Flag::RightMostChild),
      Flag::None, // Added to allow for explicit labels in the tree
      Flag::Leaf,
      Flag::None, // Added to allow for explicit labels in the tree
      Flag(Flag::Leaf | Flag::RightMostChild),
      Flag::None, // Added to allow for explicit labels in the tree
  };

  EXPECT_EQ(table, expected_table);
  EXPECT_EQ(flags, expected_flags);
}
