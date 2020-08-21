#include <gtest/gtest.h>

#include <tuple>
#include <vector>

#include "../../src/search/lazy_suffix_tree.hpp"

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/range/container/bitcompressed_vector.hpp>


using namespace lst::details;

using seqan3::operator""_dna5;
using seqan3::operator""_dna4;

class LazySuffixTreeTest : public ::testing::Test {
protected:
  void SetUp() override {
    tree = lst::LazySuffixTree<seqan3::dna5>{sequence, false};
    tree_parallel = lst::LazySuffixTree<seqan3::dna5>{sequence, true};

    dna_tree = lst::LazySuffixTree<seqan3::dna5>{dna_sequence, false};
    dna_tree_parallel = lst::LazySuffixTree<seqan3::dna5>{dna_sequence, true};

    tree4 = lst::LazySuffixTree{sequence4, false};
    tree4_parallel = lst::LazySuffixTree{sequence4, true};
  }

  sequence_t<seqan3::dna5> sequence{"CACAC"_dna5};
  lst::LazySuffixTree<seqan3::dna5> tree{};
  lst::LazySuffixTree<seqan3::dna5> tree_parallel{};

  sequence_t<seqan3::dna5> dna_sequence{"ACGATCGCT"_dna5};
  lst::LazySuffixTree<seqan3::dna5> dna_tree{};
  lst::LazySuffixTree<seqan3::dna5> dna_tree_parallel{};

  sequence_t<seqan3::dna4> sequence4{"ACGTACGTACGTACGTACGTACGT"_dna4};
  lst::LazySuffixTree<seqan3::dna4> tree4{};
  lst::LazySuffixTree<seqan3::dna4> tree4_parallel{};
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
  std::set<std::tuple<std::string, int>> expected_labels{
      std::make_tuple("AC", 2), std::make_tuple("C", 3),
      std::make_tuple("", 1),   std::make_tuple("ACAC", 1),
      std::make_tuple("AC", 1), std::make_tuple("CAC", 2),
      std::make_tuple("C", 1),  std::make_tuple("CACAC", 1),
      std::make_tuple("CAC", 1)};

  std::vector<lst::LazySuffixTree<seqan3::dna5>> trees{tree, tree_parallel};
  for (auto &tree : trees) {
    auto labels = tree.get_all_labels();
    std::set<std::tuple<std::string, int>> labels_set{labels.begin(),
                                                      labels.end()};

    EXPECT_EQ(labels_set, expected_labels);
  }

  std::set<std::tuple<std::string, int>> double_expected_labels{
      std::make_tuple("A", 3),    std::make_tuple("TAA", 1),
      std::make_tuple("", 1),     std::make_tuple("AA", 1),
      std::make_tuple("ATAA", 1), std::make_tuple("A", 1)};

  sequence_t<seqan3::dna4> double_sequence{"ATAA"_dna4};
  bool multicores[2]{true, false};
  for (bool multi_core : multicores) {
    lst::LazySuffixTree<seqan3::dna4> double_tree{double_sequence, multi_core};
    auto double_labels = double_tree.get_all_labels();
    std::set<std::tuple<std::string, int>> double_labels_set{
        double_labels.begin(), double_labels.end()};
    EXPECT_EQ(double_labels_set, double_expected_labels);
  }
}

TEST_F(LazySuffixTreeTest, Search) {
  std::vector<lst::LazySuffixTree<seqan3::dna5>> trees{tree, tree_parallel};
  for (auto &t : trees) {
    auto ac_indicies = t.search("A"_dna5);
    std::vector<int> expected_ac_indicies{1, 3};
    EXPECT_EQ(ac_indicies, expected_ac_indicies);

    auto c_indicies = t.search("C"_dna5);
    std::vector<int> expected_c_indicies{0, 2, 4};
    EXPECT_EQ(c_indicies, expected_c_indicies);

    t.expand_all();

    ac_indicies = t.search("AC"_dna5);
    EXPECT_EQ(ac_indicies, expected_ac_indicies);

    c_indicies = t.search("C"_dna5);
    EXPECT_EQ(c_indicies, (std::vector<int>{4, 0, 2}));

    auto indices = t.search("ACGT"_dna5);
    EXPECT_EQ(indices, (std::vector<int>{}));

    EXPECT_EQ(t.search(""_dna5).size(), 6);
  }

  std::vector<lst::LazySuffixTree<seqan3::dna4>> dna4_trees{tree4,
                                                            tree4_parallel};
  for (auto &t : dna4_trees) {
    EXPECT_EQ(t.search("ACGT"_dna4), (std::vector<int>{0, 4, 8, 12, 16, 20}));

    EXPECT_EQ(t.search("ACGTACGTACGTACGTACGTACGT"_dna4), (std::vector<int>{0}));
    EXPECT_EQ(t.search("CGTACGTACGTACGTACGTACGT"_dna4), (std::vector<int>{1}));
    EXPECT_EQ(t.search("GTACGTACGTACGTACGTACGT"_dna4), (std::vector<int>{2}));
    EXPECT_EQ(t.search("CCCC"_dna4), (std::vector<int>{}));
    EXPECT_EQ(t.search("ACGTCT"_dna4), (std::vector<int>{}));

    t.expand_all();
    EXPECT_EQ(t.search("ACGT"_dna4), (std::vector<int>{20, 16, 12, 8, 0, 4}));
  }

  seqan3::bitcompressed_vector<seqan3::dna5> long_sequence{
      "ACTAGCTAGCTACGCGCTATCATCATTTACGACTAGCAGCCTACTACATTATATAGCGCTAGCATCAGTCAGCACTACTACAGCAGCAGCATCACGACTAGCTACGATCAGCATCGATCGATCATTATCGACTAG"_dna5};
  ;
  std::vector<lst::LazySuffixTree<seqan3::dna5>> dna5_trees{
      lst::LazySuffixTree<seqan3::dna5>{long_sequence},
      lst::LazySuffixTree<seqan3::dna5>{long_sequence, false}};
  for (auto &t : dna5_trees) {
    EXPECT_EQ(t.search("TAG"_dna5).size(), 7);
    EXPECT_EQ(t.search("GAC"_dna5).size(), 3);
    EXPECT_EQ(t.search("ATTAT"_dna5).size(), 2);
    EXPECT_EQ(t.search("AT"_dna5).size(), 14);
    EXPECT_EQ(t.search("ACCCC"_dna5).size(), 0);
    EXPECT_EQ(t.search("TCTCTCTCTCT"_dna5).size(), 0);
    EXPECT_EQ(t.search("ATCT"_dna5).size(), 0);
  }
}

TEST_F(LazySuffixTreeTest, ExpandImplicitNodes) {
  std::vector<lst::LazySuffixTree<seqan3::dna5>> trees{tree, tree_parallel};
  for (auto &t : trees) {
    t.expand_all();
    t.expand_implicit_nodes();

    std::set<std::tuple<std::string, int>> expected_labels{
        std::make_tuple("A", 2),    std::make_tuple("C", 3),
        std::make_tuple("", 1),     std::make_tuple("AC", 2),
        std::make_tuple("CA", 2),   std::make_tuple("C", 1),
        std::make_tuple("ACA", 1),  std::make_tuple("AC", 1),
        std::make_tuple("CAC", 2),  std::make_tuple("ACAC", 1),
        std::make_tuple("CACA", 1), std::make_tuple("CAC", 1),
        std::make_tuple("ACAC", 1), std::make_tuple("CACAC", 1),
        std::make_tuple("CACAC", 1)};

    auto labels = t.get_all_labels();
    std::set<std::tuple<std::string, int>> labels_set{labels.begin(),
                                                      labels.end()};

    EXPECT_EQ(labels_set, expected_labels);
  }
}

TEST_F(LazySuffixTreeTest, SuffixLinks) {
  std::vector<lst::LazySuffixTree<seqan3::dna5>> trees{tree, tree_parallel};
  for (auto &t : trees) {
    t.expand_all();

    t.add_suffix_links();
    std::vector<int> expected_links{
        -1, // root
        4,  // AC
        0,  // C
        0,  // -
        18, // ACAC
        14, // AC
        2,  // CAC
        6,  // C
        8,  // CACAC
        10  // CAC
    };

    EXPECT_EQ(t.suffix_links, expected_links);

    t.expand_implicit_nodes();
    t.add_suffix_links();

    std::vector<int> implicit_expected_links{
        -1, // 0: root
        0,  // 2: A
        0,  // 4: C
        0,  // 6: -
        12, // 8: ACA
        14, // 10: AC-
        2,  // 12: CA
        6,  // 14: C-
        8,  // 16: CACA
        10, // 18: CAC-
        4,  // 20: AC
        26, // 22: ACAC
        18, // 24: ACAC-
        20, // 26: CAC
        22, // 28: CACAC
        24  // 30: CACAC-
    };

    EXPECT_EQ(t.suffix_links, implicit_expected_links);
  }

  dna_tree.expand_all();
  dna_tree.expand_implicit_nodes();
  dna_tree.add_suffix_links();

  std::vector<int> implicit_expected_dna_links{
      -1, 0,  0,  0,  0,  0,  4,  8,  6,  8,  2,  4,  4,  10, 20, 22, 16,
      28, 84, 86, 88, 90, 92, 94, 24, 76, 78, 80, 82, 26, 14, 48, 50, 52,
      54, 56, 18, 58, 16, 30, 96, 98, 60, 62, 64, 66, 68, 70, 72, 74};

  EXPECT_EQ(dna_tree.suffix_links, implicit_expected_dna_links);

  dna_tree_parallel.expand_all();
  dna_tree_parallel.expand_implicit_nodes();
  dna_tree_parallel.add_suffix_links();

  std::set<int> implicit_expected_dna_links_set{
      -1, 0,  0,  0,  0,  0,  4,  8,  6,  8,  2,  4,  4,  10, 20, 22, 16,
      28, 84, 86, 88, 90, 92, 94, 24, 76, 78, 80, 82, 26, 14, 48, 50, 52,
      54, 56, 18, 58, 16, 30, 96, 98, 60, 62, 64, 66, 68, 70, 72, 74};

  std::set<int> suffix_links_set{dna_tree_parallel.suffix_links.begin(),
                                 dna_tree_parallel.suffix_links.end()};

  EXPECT_EQ(suffix_links_set, implicit_expected_dna_links_set);
}

TEST_F(LazySuffixTreeTest, ReverseSuffixLinks) {
  std::vector<lst::LazySuffixTree<seqan3::dna5>> dna_trees{dna_tree,
                                                           dna_tree_parallel};
  for (auto &t : dna_trees) {
    t.expand_all();
    t.add_reverse_suffix_links();

    std::set<
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

    std::set<
        std::array<int, seqan3::alphabet_size<seqan3::gapped<seqan3::dna5>>>>
        reverse_links_set{t.reverse_suffix_links.begin(),
                          t.reverse_suffix_links.end()};

    EXPECT_EQ(reverse_links_set, expected_reverse_links);
  }
}
