#include <gtest/gtest.h>

#include <tuple>
#include <vector>

#include "../../src/search/lazy_suffix_tree.hpp"
#include "../../src/search/lazy_suffix_tree/construction.hpp"

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/range/container/bitcompressed_vector.hpp>
#include <seqan3/range/views/convert.hpp>
#include <seqan3/range/views/to.hpp>

#include <seqan3/core/debug_stream.hpp>

using namespace lst::details;

using seqan3::operator""_dna5;
using seqan3::operator""_dna4;

class LazySuffixTreeTest : public ::testing::Test {
protected:
  void SetUp() override {
    tree = lst::LazySuffixTree<seqan3::dna5>{sequence};

    dna_tree = lst::LazySuffixTree<seqan3::dna5>{dna_sequence};

    tree4 = lst::LazySuffixTree{sequence4};
  }

  sequence_t<seqan3::dna5> sequence{"CACAC"_dna5};
  lst::LazySuffixTree<seqan3::dna5> tree{};

  sequence_t<seqan3::dna5> dna_sequence{"ACGATCGCT"_dna5};
  lst::LazySuffixTree<seqan3::dna5> dna_tree{};

  sequence_t<seqan3::dna4> sequence4{"ACGTACGTACGTACGTACGTACGT"_dna4};
  lst::LazySuffixTree<seqan3::dna4> tree4{};
};

TEST_F(LazySuffixTreeTest, SimpleTest) {
  tree.expand_all();

  std::vector<int> expected_table{0, 2, 1, 8,  0, 12, 5, 0, 3, 0,
                                  5, 0, 1, 16, 5, 0,  3, 0, 5, 0};
  EXPECT_EQ(tree.table, expected_table);

  std::vector<int> expected_suffixes{3, 5, 3, 5, 5, 5};
  EXPECT_EQ(tree.suffixes, expected_suffixes);

  std::vector<Flag> expected_flags{
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

  EXPECT_EQ(tree.flags, expected_flags);
}

TEST_F(LazySuffixTreeTest, LongerTests) {
  sequence_t<seqan3::dna5> longer_sequence{"CACACACACACACACACACACACAGT"_dna5};
  lst::LazySuffixTree<seqan3::dna5> tree{longer_sequence};
  tree.expand_all();
  EXPECT_TRUE(tree.table.size() != 0);

  tree4.expand_all();
  EXPECT_TRUE(tree4.table.size() != 0);
}

TEST_F(LazySuffixTreeTest, LabelSets) {
  std::vector<std::tuple<sequence_t<seqan3::dna5>, int>> expected_labels{
      std::make_tuple(sequence_t<seqan3::dna5>{"AC"_dna5}, 2),
      std::make_tuple(sequence_t<seqan3::dna5>{"C"_dna5}, 3),
      std::make_tuple(sequence_t<seqan3::dna5>{""_dna5}, 1),
      std::make_tuple(sequence_t<seqan3::dna5>{"ACAC"_dna5}, 1),
      std::make_tuple(sequence_t<seqan3::dna5>{"AC"_dna5}, 1),
      std::make_tuple(sequence_t<seqan3::dna5>{"CAC"_dna5}, 2),
      std::make_tuple(sequence_t<seqan3::dna5>{"C"_dna5}, 1),
      std::make_tuple(sequence_t<seqan3::dna5>{"CACAC"_dna5}, 1),
      std::make_tuple(sequence_t<seqan3::dna5>{"CAC"_dna5}, 1)};
  auto labels = tree.get_all_labels();

  EXPECT_EQ(labels, expected_labels);

  std::vector<std::tuple<sequence_t<seqan3::dna4>, int>> double_expected_labels{
      std::make_tuple(sequence_t<seqan3::dna4>{"A"_dna4}, 3),
      std::make_tuple(sequence_t<seqan3::dna4>{"TAA"_dna4}, 1),
      std::make_tuple(sequence_t<seqan3::dna4>{""_dna4}, 1),
      std::make_tuple(sequence_t<seqan3::dna4>{"AA"_dna4}, 1),
      std::make_tuple(sequence_t<seqan3::dna4>{"ATAA"_dna4}, 1),
      std::make_tuple(sequence_t<seqan3::dna4>{"A"_dna4}, 1)};

  sequence_t<seqan3::dna4> double_sequence{"ATAA"_dna4};
  lst::LazySuffixTree<seqan3::dna4> double_tree{double_sequence};
  auto double_labels = double_tree.get_all_labels();
  EXPECT_EQ(double_labels, double_expected_labels);
}

TEST_F(LazySuffixTreeTest, Search) {
  auto ac_indicies = tree.search("A"_dna5);
  std::vector<int> expected_ac_indicies{1, 3};
  EXPECT_EQ(ac_indicies, expected_ac_indicies);

  auto c_indicies = tree.search("C"_dna5);
  std::vector<int> expected_c_indicies{0, 2, 4};
  EXPECT_EQ(c_indicies, expected_c_indicies);

  tree.expand_all();

  ac_indicies = tree.search("AC"_dna5);
  EXPECT_EQ(ac_indicies, expected_ac_indicies);

  c_indicies = tree.search("C"_dna5);
  EXPECT_EQ(c_indicies, (std::vector<int>{4, 0, 2}));

  auto indices = tree.search("ACGT"_dna5);
  EXPECT_EQ(indices, (std::vector<int>{}));

  EXPECT_EQ(tree.search(""_dna5).size(), 6);

  EXPECT_EQ(tree4.search("ACGT"_dna4), (std::vector<int>{0, 4, 8, 12, 16, 20}));

  EXPECT_EQ(tree4.search("ACGTACGTACGTACGTACGTACGT"_dna4),
            (std::vector<int>{0}));
  EXPECT_EQ(tree4.search("CGTACGTACGTACGTACGTACGT"_dna4),
            (std::vector<int>{1}));
  EXPECT_EQ(tree4.search("GTACGTACGTACGTACGTACGT"_dna4), (std::vector<int>{2}));
  EXPECT_EQ(tree4.search("CCCC"_dna4), (std::vector<int>{}));
  EXPECT_EQ(tree4.search("ACGTCT"_dna4), (std::vector<int>{}));

  tree4.expand_all();
  EXPECT_EQ(tree4.search("ACGT"_dna4), (std::vector<int>{20, 16, 12, 8, 0, 4}));

  seqan3::bitcompressed_vector<seqan3::dna5> long_sequence{
      "ACTAGCTAGCTACGCGCTATCATCATTTACGACTAGCAGCCTACTACATTATATAGCGCTAGCATCAGTCAGCACTACTACAGCAGCAGCATCACGACTAGCTACGATCAGCATCGATCGATCATTATCGACTAG"_dna5};
  lst::LazySuffixTree<seqan3::dna5> tree5{long_sequence};
  EXPECT_EQ(tree5.search("TAG"_dna5).size(), 7);
  EXPECT_EQ(tree5.search("GAC"_dna5).size(), 3);
  EXPECT_EQ(tree5.search("ATTAT"_dna5).size(), 2);
  EXPECT_EQ(tree5.search("AT"_dna5).size(), 14);
  EXPECT_EQ(tree5.search("ACCCC"_dna5).size(), 0);
  EXPECT_EQ(tree5.search("TCTCTCTCTCT"_dna5).size(), 0);
  EXPECT_EQ(tree5.search("ATCT"_dna5).size(), 0);
}

TEST_F(LazySuffixTreeTest, ExpandImplicitNodes) {
  tree.expand_all();
  tree.expand_implicit_nodes();

  std::vector<std::tuple<sequence_t<seqan3::dna5>, int>> expected_labels{
      std::make_tuple(sequence_t<seqan3::dna5>{"A"_dna5}, 2),
      std::make_tuple(sequence_t<seqan3::dna5>{"C"_dna5}, 3),
      std::make_tuple(sequence_t<seqan3::dna5>{""_dna5}, 1),
      std::make_tuple(sequence_t<seqan3::dna5>{"AC"_dna5}, 2),
      std::make_tuple(sequence_t<seqan3::dna5>{"CA"_dna5}, 2),
      std::make_tuple(sequence_t<seqan3::dna5>{"C"_dna5}, 1),
      std::make_tuple(sequence_t<seqan3::dna5>{"ACA"_dna5}, 1),
      std::make_tuple(sequence_t<seqan3::dna5>{"AC"_dna5}, 1),
      std::make_tuple(sequence_t<seqan3::dna5>{"CAC"_dna5}, 2),
      std::make_tuple(sequence_t<seqan3::dna5>{"ACAC"_dna5}, 1),
      std::make_tuple(sequence_t<seqan3::dna5>{"CACA"_dna5}, 1),
      std::make_tuple(sequence_t<seqan3::dna5>{"CAC"_dna5}, 1),
      std::make_tuple(sequence_t<seqan3::dna5>{"ACAC"_dna5}, 1),
      std::make_tuple(sequence_t<seqan3::dna5>{"CACAC"_dna5}, 1),
      std::make_tuple(sequence_t<seqan3::dna5>{"CACAC"_dna5}, 1)};

  auto labels = tree.get_all_labels();

  EXPECT_EQ(labels, expected_labels);
}

TEST_F(LazySuffixTreeTest, SuffixLinks) {
  tree.expand_all();

  tree.add_suffix_links();
  std::vector<int> expected_links{-1, 4, 0, 0, 18, 14, 2, 6, 8, 10};

  EXPECT_EQ(tree.suffix_links, expected_links);

  tree.expand_implicit_nodes();
  tree.add_suffix_links();

  std::vector<int> implicit_expected_links{-1, 0,  0, 0,  12, 14, 2,  6,
                                           8,  10, 4, 26, 18, 20, 22, 24};

  EXPECT_EQ(tree.suffix_links, implicit_expected_links);

  dna_tree.expand_all();
  dna_tree.expand_implicit_nodes();
  dna_tree.add_suffix_links();

  std::vector<int> implicit_expected_dna_links{
      -1, 0,  0,  0,  0,  0,  4,  8,  6,  8,  2,  4,  4,  10, 20, 22, 16,
      28, 84, 86, 88, 90, 92, 94, 24, 76, 78, 80, 82, 26, 14, 48, 50, 52,
      54, 56, 18, 58, 16, 30, 96, 98, 60, 62, 64, 66, 68, 70, 72, 74};

  EXPECT_EQ(dna_tree.suffix_links, implicit_expected_dna_links);
}

TEST_F(LazySuffixTreeTest, UnevaluatedSuffixLinks) {
  tree4.search("AC"_dna4);
  tree4.add_suffix_links();
  std::vector<int> expected_tree4_links{-1, 4, 6, 8, 0, 0, -1, -1};
  EXPECT_EQ(tree4.suffix_links, expected_tree4_links);
}

TEST_F(LazySuffixTreeTest, ReverseSuffixLinks) {
  dna_tree.expand_all();
  dna_tree.add_reverse_suffix_links();

  std::vector<
      std::array<int, seqan3::alphabet_size<seqan3::gapped<seqan3::dna5>>>>
      expected_reverse_links{
          {2, 4, 6, -1, 8, 10},     // root
          {-1, -1, -1, -1, -1, -1}, // A, 2
          {-1, -1, -1, -1, -1, -1}, // C, 4
          {-1, 16, -1, -1, -1, -1}, // G, 6
          {-1, -1, -1, -1, -1, -1}, // T, 8
          {-1, -1, -1, -1, 26, -1}, // -, 10
          {-1, -1, -1, -1, -1, -1}, // ACGATCGCT-, 12
          {-1, -1, 20, -1, -1, -1}, // ATCGCT-, 14
          {-1, -1, -1, -1, -1, -1}, // CG, 16
          {-1, -1, 22, -1, -1, -1}, // CT-, 18
          {-1, 28, -1, -1, -1, -1}, // GATCGCT-, 20
          {-1, 30, -1, -1, -1, -1}, // GCT-, 22
          {14, -1, -1, -1, -1, -1}, // TCGCT-, 24
          {-1, 18, -1, -1, -1, -1}, // T-, 26
          {12, -1, -1, -1, -1, -1}, // CGATCGCT-, 28
          {-1, -1, -1, -1, 24, -1}  // CGCT-, 30
      };

  EXPECT_EQ(dna_tree.reverse_suffix_links, expected_reverse_links);
}
