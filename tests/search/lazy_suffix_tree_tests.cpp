#include <gtest/gtest.h>

#include "../../src/search/lazy_suffix_tree.cpp"

#include <seqan3/alphabet/nucleotide/dna5.hpp>
using seqan3::operator""_dna5;
#include <seqan3/alphabet/nucleotide/dna4.hpp>
using seqan3::operator""_dna4;
#include <seqan3/range/view/convert.hpp>

#include <seqan3/core/debug_stream.hpp>

class LazySuffixTreeTest : public ::testing::Test {
protected:
  lst::LazySuffixTree<seqan3::dna5> tree{"CACAC"_dna5};

  template <seqan3::Alphabet alphabet_t = seqan3::dna5>
  std::vector<seqan3::gapped<alphabet_t>>
  make_gapped(std::vector<alphabet_t> sequence) {
    std::vector<seqan3::gapped<alphabet_t>> sequence_ =
        sequence | seqan3::view::convert<seqan3::gapped<alphabet_t>>;

    return sequence_;
  }

  template <seqan3::Alphabet alphabet_t = seqan3::dna5>
  std::vector<seqan3::gapped<alphabet_t>>
  make_gapped_with_gap(std::vector<alphabet_t> sequence) {
    std::vector<seqan3::gapped<alphabet_t>> sequence_ =
        sequence | seqan3::view::convert<seqan3::gapped<alphabet_t>>;

    sequence_.push_back(seqan3::gap{});
    return sequence_;
  }
};

TEST_F(LazySuffixTreeTest, SimpleTest) {
  tree.expand_all();
  std::vector<int> expected_table{1, 5, 0, 7, 5, 3, 5, 1, 10, 5, 3, 5};
  EXPECT_EQ(tree.table, expected_table);

  std::vector<int> expected_suffixes{3, 5, 3, 5, 5, 5};
  EXPECT_EQ(tree.suffixes, expected_suffixes);

  std::vector<lst::Flag> expected_flags{
      lst::Flag::None,
      lst::Flag::None,
      lst::Flag::None,
      lst::Flag::None,
      lst::Flag(lst::Flag::Leaf | lst::Flag::RightMostChild),
      lst::Flag::Leaf,
      lst::Flag(lst::Flag::Leaf | lst::Flag::RightMostChild),
      lst::Flag::None,
      lst::Flag::None,
      lst::Flag(lst::Flag::Leaf | lst::Flag::RightMostChild),
      lst::Flag::Leaf,
      lst::Flag(lst::Flag::Leaf | lst::Flag::RightMostChild)};
  EXPECT_EQ(tree.flags, expected_flags);
}

TEST_F(LazySuffixTreeTest, CountSuffixes) {
  std::array<int, 6> expected_counts{2, 3, 0, 0, 0, 1};
  EXPECT_EQ(tree.count_suffixes(0, 6), expected_counts);

  lst::LazySuffixTree<seqan3::dna5> tree2{"ACGACGACGTNNNN"_dna5};
  expected_counts = {3, 3, 3, 4, 1, 1};
  EXPECT_EQ(tree2.count_suffixes(0, 15), expected_counts);
  EXPECT_THROW(tree2.count_suffixes(0, 16), std::invalid_argument);
  EXPECT_THROW(tree2.count_suffixes(-1, 13), std::invalid_argument);
}

TEST_F(LazySuffixTreeTest, LongestCommonPrefix) {
  EXPECT_EQ(tree.longest_common_prefix({1, 3}), 2);
  EXPECT_EQ(tree.longest_common_prefix({0, 2}), 3);
  EXPECT_EQ(tree.longest_common_prefix({0, 2, 4}), 1);

  lst::LazySuffixTree<seqan3::dna5> tree2{"ACGTTACGACGTTATTT"_dna5};
  EXPECT_EQ(tree2.longest_common_prefix({0, 5, 8}), 3);
}

TEST_F(LazySuffixTreeTest, ConstructorTest) {
  std::vector<int> expected_suffixes{0, 1, 2, 3, 4, 5};

  std::vector<seqan3::gapped<seqan3::dna5>> seq{};
  seq.push_back('C'_dna5);
  seq.push_back('A'_dna5);
  seq.push_back('C'_dna5);
  seq.push_back('A'_dna5);
  seq.push_back('C'_dna5);
  seq.push_back(seqan3::gap{});

  EXPECT_EQ(tree.sequence, seq);
  EXPECT_EQ(tree.suffixes, expected_suffixes);
}

TEST_F(LazySuffixTreeTest, SuffixPointers) {
  std::array<int, 6> counts{2, 3, 0, 0, 0, 1};
  auto actual_pointers = tree.suffix_pointers(counts);

  std::array<int, 6> expected_pointers{0, 2, 5, 5, 5, 5};
  EXPECT_EQ(actual_pointers, expected_pointers);
}

TEST_F(LazySuffixTreeTest, SortSuffixes) {
  int lower_bound = 0;
  int upper_bound = 6;
  std::array<int, 6> counts = tree.count_suffixes(lower_bound, upper_bound);
  tree.sort_suffixes(counts, lower_bound, upper_bound);

  std::vector<int> expected_suffixes{1, 3, 0, 2, 4, 5};
  EXPECT_EQ(tree.suffixes, expected_suffixes);
}

TEST_F(LazySuffixTreeTest, AddLcpToSuffixes) {
  int lower_bound = 0;
  int upper_bound = 6;
  auto counts = tree.count_suffixes(lower_bound, upper_bound);
  tree.sort_suffixes(counts, lower_bound, upper_bound);
  tree.add_lcp_to_suffixes(lower_bound, upper_bound);

  std::vector<int> expected_suffixes{1, 3, 0, 2, 4, 5};
  EXPECT_EQ(tree.suffixes, expected_suffixes);

  tree.add_lcp_to_suffixes(0, 2);
  expected_suffixes = {3, 5, 0, 2, 4, 5};
  EXPECT_EQ(tree.suffixes, expected_suffixes);
}

TEST_F(LazySuffixTreeTest, AddChildren) {
  std::array<int, 6> counts = tree.count_suffixes(0, 6);
  tree.add_children(counts, 0);

  std::vector<int> expected_table{0, 2, 2, 5, 5};
  EXPECT_EQ(tree.table, expected_table);

  std::vector<lst::Flag> expected_flags{
      lst::Flag::Unevaluated, lst::Flag::None, lst::Flag::Unevaluated,
      lst::Flag::None, lst::Flag(lst::Flag::Leaf | lst::Flag::RightMostChild)};
  EXPECT_EQ(tree.flags, expected_flags);
}

TEST_F(LazySuffixTreeTest, ExpandRoot) {
  tree.expand_root();

  std::vector<int> expected_table{0, 2, 2, 5, 5};
  EXPECT_EQ(tree.table, expected_table);

  std::vector<lst::Flag> expected_flags{
      lst::Flag::Unevaluated, lst::Flag::None, lst::Flag::Unevaluated,
      lst::Flag::None, lst::Flag(lst::Flag::Leaf | lst::Flag::RightMostChild)};
  EXPECT_EQ(tree.flags, expected_flags);

  std::vector<int> expected_suffixes{1, 3, 0, 2, 4, 5};
  EXPECT_EQ(tree.suffixes, expected_suffixes);
}

TEST_F(LazySuffixTreeTest, ExpandTest) {
  tree.expand_root();
  tree.expand_node(0);

  std::vector<int> expected_table{1, 5, 2, 5, 5, 3, 5};

  std::vector<lst::Flag> expected_flags{
      lst::Flag::None,
      lst::Flag::None,
      lst::Flag::Unevaluated,
      lst::Flag::None,
      lst::Flag(lst::Flag::Leaf | lst::Flag::RightMostChild),
      lst::Flag::Leaf,
      lst::Flag(lst::Flag::Leaf | lst::Flag::RightMostChild)};

  EXPECT_EQ(tree.table, expected_table);
  EXPECT_EQ(tree.flags, expected_flags);
}

TEST_F(LazySuffixTreeTest, LongerTests) {
  lst::LazySuffixTree<seqan3::dna5> tree{"CACACACACACACACACACACACAGT"_dna5};
  tree.expand_all();
  EXPECT_TRUE(tree.table.size() != 0);

  lst::LazySuffixTree<seqan3::dna4> tree4{"ACGTACGTACGTACGTACGTACGT"_dna4};
  tree4.expand_all();
  EXPECT_TRUE(tree4.table.size() != 0);
}

TEST_F(LazySuffixTreeTest, LabelSets) {
  std::vector<std::vector<seqan3::gapped<seqan3::dna5>>> expected_labels{
      make_gapped("AC"_dna5),          make_gapped("C"_dna5),
      make_gapped_with_gap(""_dna5),   make_gapped_with_gap("ACAC"_dna5),
      make_gapped_with_gap("AC"_dna5), make_gapped("CAC"_dna5),
      make_gapped_with_gap("C"_dna5),  make_gapped_with_gap("CACAC"_dna5),
      make_gapped_with_gap("CAC"_dna5)};
  auto labels = tree.get_all_labels();
  EXPECT_EQ(labels, expected_labels);

  std::vector<std::vector<seqan3::gapped<seqan3::dna4>>> double_expected_labels{
      make_gapped<seqan3::dna4>("A"_dna4),
      make_gapped_with_gap<seqan3::dna4>("TAA"_dna4),
      make_gapped_with_gap<seqan3::dna4>(""_dna4),
      make_gapped_with_gap<seqan3::dna4>("AA"_dna4),
      make_gapped_with_gap<seqan3::dna4>("ATAA"_dna4),
      make_gapped_with_gap<seqan3::dna4>("A"_dna4)};

  lst::LazySuffixTree<seqan3::dna4> double_tree{"ATAA"_dna4};
  auto double_labels = double_tree.get_all_labels();

  EXPECT_EQ(double_labels, double_expected_labels);
}

TEST_F(LazySuffixTreeTest, Count) {
  tree.expand_root();
  int a_count = tree.count_occurrences(0);
  EXPECT_EQ(a_count, 2);
  int b_count = tree.count_occurrences(2);
  EXPECT_EQ(b_count, 3);

  tree.expand_all();
  a_count = tree.count_occurrences(0);
  EXPECT_EQ(a_count, 2);
  b_count = tree.count_occurrences(2);
  EXPECT_EQ(b_count, 3);

  lst::LazySuffixTree<seqan3::dna4> double_tree{"ATAA"_dna4};
  EXPECT_THROW(double_tree.count_occurrences(0), std::invalid_argument);

  double_tree.expand_root();
  a_count = double_tree.count_occurrences(0);
  EXPECT_EQ(a_count, 3);

  double_tree.expand_all();
  a_count = double_tree.count_occurrences(0);
  EXPECT_EQ(a_count, 3);
}
