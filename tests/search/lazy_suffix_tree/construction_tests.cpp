#include <gtest/gtest.h>

#include <iostream>

#include "../../../src/search/lazy_suffix_tree/construction.hpp"
using namespace lst::details;

#include <seqan3/alphabet/gap/all.hpp>
#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/core/debug_stream.hpp>

using seqan3::operator""_dna5;

class ConstructionTests : public ::testing::Test {
protected:
  void SetUp() override {}

  sequence_t<seqan3::dna5> sequence{"CACAC"_dna5};
  std::vector<size_t> suffixes{0, 1, 2, 3, 4, 5};
  Table<> table{{0, 2, Flag::UNEVALUATED},
                {2, 5, Flag::UNEVALUATED},
                {5, 0, Flag(Flag::LEAF | Flag::RIGHT_MOST_CHILD)}};
};

TEST_F(ConstructionTests, CountSuffixes) {
  std::array<size_t, 6> expected_counts{2, 3, 0, 0, 0, 1};
  auto [counts, _] = count_suffixes(0, 6, sequence, suffixes);
  EXPECT_EQ(counts, expected_counts);

  EXPECT_DEBUG_DEATH(count_suffixes(0, 16, sequence, suffixes), "");
}

TEST_F(ConstructionTests, LongestCommonPrefix) {
  suffixes = std::vector<size_t>{1, 3, 0, 2, 4, 5};
  EXPECT_EQ(longest_common_prefix(1, 2, sequence, suffixes), 2);
  EXPECT_EQ(longest_common_prefix(2, 4, sequence, suffixes), 3);
  EXPECT_EQ(longest_common_prefix(3, 5, sequence, suffixes), 1);
}

TEST_F(ConstructionTests, SuffixPointers) {
  std::array<size_t, 6> counts{2, 3, 0, 0, 0, 1};
  auto actual_pointers = suffix_pointers<seqan3::dna5>(counts);

  std::array<size_t, 6> expected_pointers{0, 2, 5, 5, 5, 5};
  EXPECT_EQ(actual_pointers, expected_pointers);
}

TEST_F(ConstructionTests, SortSuffixes) {
  size_t lower_bound = 0;
  size_t upper_bound = 6;
  auto [counts, _] =
      count_suffixes(lower_bound, upper_bound, sequence, suffixes);
  sort_suffixes(counts, lower_bound, upper_bound, sequence, suffixes);

  std::vector<size_t> expected_suffixes{1, 3, 0, 2, 4, 5};
  EXPECT_EQ(suffixes, expected_suffixes);
}

TEST_F(ConstructionTests, AddLcpToSuffixes) {
  size_t lower_bound = 0;
  size_t upper_bound = 6;
  auto [counts, _] =
      count_suffixes(lower_bound, upper_bound, sequence, suffixes);
  sort_suffixes(counts, lower_bound, upper_bound, sequence, suffixes);

  size_t lcp =
      longest_common_prefix(lower_bound, upper_bound, sequence, suffixes);
  add_lcp_to_suffixes(lower_bound, upper_bound, lcp, suffixes);

  std::vector<size_t> expected_suffixes{1, 3, 0, 2, 4, 5};
  EXPECT_EQ(suffixes, expected_suffixes);

  size_t lcp_0_2 = longest_common_prefix(0, 2, sequence, suffixes);
  add_lcp_to_suffixes(0, 2, lcp_0_2, suffixes);
  expected_suffixes = {3, 5, 0, 2, 4, 5};
  EXPECT_EQ(suffixes, expected_suffixes);
}

TEST_F(ConstructionTests, AddChildren) {
  table = Table{};
  auto [counts, _] = count_suffixes(0, 6, sequence, suffixes);
  add_children<seqan3::dna5>(counts, 0, suffixes, table);

  std::vector<size_t> tree_values{};
  std::vector<Flag> tree_flags{};
  for (auto &e : table.table) {
    tree_values.push_back(e.first);
    tree_values.push_back(e.second);
    tree_flags.push_back(e.flag);
  }

  // Changed to allow for explicit labels in the tree
  std::vector<size_t> expected_table{0, 2, 2, 5, 5, 0};

  EXPECT_EQ(tree_values, expected_table);

  std::vector<Flag> expected_flags{
      Flag::UNEVALUATED,
      Flag::UNEVALUATED,
      Flag(Flag::LEAF | Flag::RIGHT_MOST_CHILD),
  };
  EXPECT_EQ(tree_flags, expected_flags);
}

TEST_F(ConstructionTests, ExpandRoot) {
  table = Table{};
  expand_root(sequence, suffixes, table);

  std::vector<size_t> tree_values{};
  std::vector<Flag> tree_flags{};
  for (auto &e : table.table) {
    tree_values.push_back(e.first);
    tree_values.push_back(e.second);
    tree_flags.push_back(e.flag);
  }

  // Changed to allow for explicit labels in the tree
  std::vector<size_t> expected_table{0, 2, 2, 5, 5, 0};
  EXPECT_EQ(tree_values, expected_table);

  std::vector<Flag> expected_flags{
      Flag::UNEVALUATED,
      Flag::UNEVALUATED,
      Flag(Flag::LEAF | Flag::RIGHT_MOST_CHILD),
  };
  EXPECT_EQ(tree_flags, expected_flags);

  std::vector<size_t> expected_suffixes{1, 3, 0, 2, 4, 5};
  EXPECT_EQ(suffixes, expected_suffixes);
}

TEST_F(ConstructionTests, ExpandTest) {
  table = Table{};
  expand_root(sequence, suffixes, table);
  expand_node(0, sequence, suffixes, table);

  // Changed to allow for explicit labels in the tree
  std::vector<size_t> expected_table{1, 3, 2, 5, 5, 0, 3, 0, 5, 0};

  std::vector<Flag> expected_flags{
      Flag::NONE,
      Flag::UNEVALUATED,
      Flag(Flag::LEAF | Flag::RIGHT_MOST_CHILD),
      Flag::LEAF,
      Flag(Flag::LEAF | Flag::RIGHT_MOST_CHILD),
  };

  std::vector<size_t> tree_values{};
  std::vector<Flag> tree_flags{};
  for (auto &e : table.table) {
    tree_values.push_back(e.first);
    tree_values.push_back(e.second);
    tree_flags.push_back(e.flag);
  }

  EXPECT_EQ(tree_values, expected_table);
  EXPECT_EQ(tree_flags, expected_flags);
}

TEST_F(ConstructionTests, GetCharacter) {
  auto c_gapped = seqan3::gapped<seqan3::dna5>{'C'_dna5};
  auto c = get_character(sequence, 0);
  EXPECT_EQ(c, c_gapped);

  auto gap_gapped = seqan3::gapped<seqan3::dna5>{seqan3::gap{}};
  auto gap = get_character(sequence, sequence.size());
  EXPECT_EQ(gap, gap_gapped);
}