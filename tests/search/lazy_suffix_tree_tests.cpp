#include <gtest/gtest.h>

#include "../../src/search/lazy_suffix_tree.hpp"
#include "../../src/search/lazy_suffix_tree_construction.hpp"

#include <seqan3/alphabet/nucleotide/dna5.hpp>
using seqan3::operator""_dna5;
#include <seqan3/alphabet/nucleotide/dna4.hpp>
using seqan3::operator""_dna4;
#include <seqan3/range/view/convert.hpp>

#include <seqan3/core/debug_stream.hpp>

using namespace lst::details;

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

  std::vector<Flag> expected_flags{Flag::None,
                                   Flag::None,
                                   Flag::None,
                                   Flag::None,
                                   Flag(Flag::Leaf | Flag::RightMostChild),
                                   Flag::Leaf,
                                   Flag(Flag::Leaf | Flag::RightMostChild),
                                   Flag::None,
                                   Flag::None,
                                   Flag(Flag::Leaf | Flag::RightMostChild),
                                   Flag::Leaf,
                                   Flag(Flag::Leaf | Flag::RightMostChild)};
  EXPECT_EQ(tree.flags, expected_flags);
}

TEST_F(LazySuffixTreeTest, ConstructorTest) {
  std::vector<int> expected_suffixes{1, 3, 0, 2, 4, 5};

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

TEST_F(LazySuffixTreeTest, LongerTests) {
  lst::LazySuffixTree<seqan3::dna5> tree{"CACACACACACACACACACACACAGT"_dna5};
  tree.expand_all();
  EXPECT_TRUE(tree.table.size() != 0);

  lst::LazySuffixTree<seqan3::dna4> tree4{"ACGTACGTACGTACGTACGTACGT"_dna4};
  tree4.expand_all();
  EXPECT_TRUE(tree4.table.size() != 0);
}

TEST_F(LazySuffixTreeTest, LabelSets) {
  std::vector<std::tuple<std::vector<seqan3::gapped<seqan3::dna5>>, int>>
      expected_labels{std::make_tuple(make_gapped("AC"_dna5), 2),
                      std::make_tuple(make_gapped("C"_dna5), 3),
                      std::make_tuple(make_gapped_with_gap(""_dna5), 1),
                      std::make_tuple(make_gapped_with_gap("ACAC"_dna5), 1),
                      std::make_tuple(make_gapped_with_gap("AC"_dna5), 1),
                      std::make_tuple(make_gapped("CAC"_dna5), 2),
                      std::make_tuple(make_gapped_with_gap("C"_dna5), 1),
                      std::make_tuple(make_gapped_with_gap("CACAC"_dna5), 1),
                      std::make_tuple(make_gapped_with_gap("CAC"_dna5), 1)};
  auto labels = tree.get_all_labels();
  EXPECT_EQ(labels, expected_labels);

  std::vector<std::tuple<std::vector<seqan3::gapped<seqan3::dna4>>, int>>
      double_expected_labels{
          std::make_tuple(make_gapped<seqan3::dna4>("A"_dna4), 3),
          std::make_tuple(make_gapped_with_gap<seqan3::dna4>("TAA"_dna4), 1),
          std::make_tuple(make_gapped_with_gap<seqan3::dna4>(""_dna4), 1),
          std::make_tuple(make_gapped_with_gap<seqan3::dna4>("AA"_dna4), 1),
          std::make_tuple(make_gapped_with_gap<seqan3::dna4>("ATAA"_dna4), 1),
          std::make_tuple(make_gapped_with_gap<seqan3::dna4>("A"_dna4), 1)};

  lst::LazySuffixTree<seqan3::dna4> double_tree{"ATAA"_dna4};
  auto double_labels = double_tree.get_all_labels();
  EXPECT_EQ(double_labels, double_expected_labels);
}

TEST_F(LazySuffixTreeTest, Search) {
  auto ac_indicies = tree.search("A"_dna5);
  std::vector<int> expected_ac_indicies{1, 3};
  EXPECT_EQ(ac_indicies, expected_ac_indicies);

  auto c_indicies = tree.search("C"_dna5);
  std::vector<int> expected_c_indicies{4, 0, 2};
  EXPECT_EQ(c_indicies, expected_c_indicies);

  tree.expand_all();

  ac_indicies = tree.search("AC"_dna5);
  EXPECT_EQ(ac_indicies, expected_ac_indicies);

  c_indicies = tree.search("C"_dna5);
  EXPECT_EQ(c_indicies, (std::vector<int>{4, 0, 2}));

  EXPECT_EQ(tree.search("ACGT"_dna5), (std::vector<int>{}));

  EXPECT_EQ(tree.search(""_dna5).size(), 6);

  lst::LazySuffixTree<seqan3::dna4> tree4{"ACGTACGTACGTACGTACGTACGT"_dna4};
  EXPECT_EQ(tree4.search("ACGT"_dna4), (std::vector<int>{20, 16, 12, 8, 0, 4}));

  EXPECT_EQ(tree4.search("ACGTACGTACGTACGTACGTACGT"_dna4),
            (std::vector<int>{0}));
  EXPECT_EQ(tree4.search("CGTACGTACGTACGTACGTACGT"_dna4),
            (std::vector<int>{1}));
  EXPECT_EQ(tree4.search("GTACGTACGTACGTACGTACGT"_dna4), (std::vector<int>{2}));
  EXPECT_EQ(tree4.search("CCCC"_dna4), (std::vector<int>{}));
  EXPECT_EQ(tree4.search("ACGTCT"_dna4), (std::vector<int>{}));

  tree4.expand_all();
  EXPECT_EQ(tree4.search("ACGT"_dna4), (std::vector<int>{20, 16, 12, 8, 0, 4}));

  lst::LazySuffixTree<seqan3::dna5> tree5{
      "ACTAGCTAGCTACGCGCTATCATCATTTACGACTAGCAGCCTACTACATTATATAGCGCTAGCATCAGTCAGCACTACTACAGCAGCAGCATCACGACTAGCTACGATCAGCATCGATCGATCATTATCGACTAG"_dna5};
  EXPECT_EQ(tree5.search("TAG"_dna5).size(), 7);
  EXPECT_EQ(tree5.search("GAC"_dna5).size(), 3);
  EXPECT_EQ(tree5.search("ATTAT"_dna5).size(), 2);
  EXPECT_EQ(tree5.search("AT"_dna5).size(), 14);
  EXPECT_EQ(tree5.search("ACCCC"_dna5).size(), 0);
  EXPECT_EQ(tree5.search("TCTCTCTCTCT"_dna5).size(), 0);
  EXPECT_EQ(tree5.search("ATCT"_dna5).size(), 0);
}

TEST_F(LazySuffixTreeTest, Find) {
  auto [a_index, a_lcp] = tree.find("A"_dna5);
  EXPECT_EQ(a_index, 0);
  EXPECT_EQ(a_lcp, 0);

  auto [c_index, c_lcp] = tree.find("C"_dna5);
  EXPECT_EQ(c_index, 2);
  EXPECT_EQ(c_lcp, 0);

  auto [tt_index, tt_lcp] = tree.find("TT"_dna5);
  EXPECT_EQ(tt_index, -1);
  EXPECT_EQ(tt_lcp, -1);

  lst::LazySuffixTree<seqan3::dna4> tree4{"ACGTACGTACGTACGTACGTACGT"_dna4};
  auto [acgt_index, acgt_lcp] = tree4.find("ACGT"_dna4);
  EXPECT_EQ(acgt_index, 0);
  EXPECT_EQ(acgt_lcp, 0);
}

TEST_F(LazySuffixTreeTest, NodeOccurrences) {
  EXPECT_EQ(tree.node_occurrences(0), 2);

  EXPECT_EQ(tree.node_occurrences(2), 3);

  tree.expand_all();

  EXPECT_EQ(tree.node_occurrences(0), 2);

  EXPECT_EQ(tree.node_occurrences(2), 3);
}
