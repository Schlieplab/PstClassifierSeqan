#include <gtest/gtest.h>

#include <tuple>
#include <vector>

#include "../../src/search/lazy_suffix_tree.hpp"

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/alphabet/nucleotide/dna5.hpp>

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

  Table<> expected_table{{0, 1, Flag::RIGHT_MOST_CHILD},
                         {1, 4, Flag::NONE},
                         {0, 6, Flag::NONE},
                         {5, 0, Flag(Flag::LEAF | Flag::RIGHT_MOST_CHILD)},
                         {3, 0, Flag::LEAF},
                         {5, 0, Flag(Flag::LEAF | Flag::RIGHT_MOST_CHILD)},
                         {1, 8, Flag::NONE},
                         {5, 0, Flag(Flag::LEAF | Flag::RIGHT_MOST_CHILD)},
                         {3, 0, Flag::LEAF},
                         {5, 0, Flag(Flag::LEAF | Flag::RIGHT_MOST_CHILD)}};

  std::vector<size_t> expected_values{};
  for (auto &e : expected_table.table) {
    expected_values.push_back(e.first);
    expected_values.push_back(e.second);
  }
  std::vector<size_t> tree_values{};
  for (auto &e : tree.table.table) {
    tree_values.push_back(e.first);
    tree_values.push_back(e.second);
  }
  std::vector<Flag> expected_flags{};
  for (auto &e : expected_table.table) {
    expected_flags.push_back(e.flag);
  }
  std::vector<Flag> tree_flags{};
  for (auto &e : tree.table.table) {
    tree_flags.push_back(e.flag);
  }
  EXPECT_EQ(tree_values, expected_values);
  EXPECT_EQ(tree_flags, expected_flags);

  std::vector<size_t> expected_suffixes{3, 5, 3, 5, 5, 5};
  EXPECT_EQ(tree.suffixes, expected_suffixes);
}

TEST_F(LazySuffixTreeTest, LongerTests) {
  sequence_t<seqan3::dna5> longer_sequence{"CACACACACACACACACACACACAGT"_dna5};
  lst::LazySuffixTree<seqan3::dna5> tree{longer_sequence};
  tree.expand_all();
  EXPECT_FALSE(tree.table.table.empty());

  tree4.expand_all();
  EXPECT_FALSE(tree4.table.table.empty());
}

TEST_F(LazySuffixTreeTest, LabelSets) {
  std::set<std::tuple<std::string, size_t>> expected_labels{
      std::make_tuple("AC", 2), std::make_tuple("C", 3),
      std::make_tuple("", 1),   std::make_tuple("ACAC", 1),
      std::make_tuple("AC", 1), std::make_tuple("CAC", 2),
      std::make_tuple("C", 1),  std::make_tuple("CACAC", 1),
      std::make_tuple("CAC", 1)};

  std::vector<lst::LazySuffixTree<seqan3::dna5>> trees{tree, tree_parallel};
  for (auto &tree : trees) {
    auto labels = tree.get_all_labels();
    std::set<std::tuple<std::string, size_t>> labels_set{labels.begin(),
                                                         labels.end()};

    EXPECT_EQ(labels_set, expected_labels);
  }

  std::set<std::tuple<std::string, size_t>> double_expected_labels{
      std::make_tuple("A", 3),    std::make_tuple("TAA", 1),
      std::make_tuple("", 1),     std::make_tuple("AA", 1),
      std::make_tuple("ATAA", 1), std::make_tuple("A", 1)};

  sequence_t<seqan3::dna4> double_sequence{"ATAA"_dna4};
  bool multicores[2]{true, false};
  for (bool multi_core : multicores) {
    lst::LazySuffixTree<seqan3::dna4> double_tree{double_sequence, multi_core};
    auto double_labels = double_tree.get_all_labels();
    std::set<std::tuple<std::string, size_t>> double_labels_set{
        double_labels.begin(), double_labels.end()};
    EXPECT_EQ(double_labels_set, double_expected_labels);
  }
}

TEST_F(LazySuffixTreeTest, Search) {
  std::vector<lst::LazySuffixTree<seqan3::dna5>> trees{tree, tree_parallel};
  for (auto &t : trees) {
    auto ac_indicies = t.search("A"_dna5);
    std::vector<size_t> expected_ac_indicies{1, 3};
    EXPECT_EQ(ac_indicies, expected_ac_indicies);

    auto c_indicies = t.search("C"_dna5);
    std::vector<size_t> expected_c_indicies{0, 2, 4};
    EXPECT_EQ(c_indicies, expected_c_indicies);

    t.expand_all();

    ac_indicies = t.search("AC"_dna5);
    EXPECT_EQ(ac_indicies, expected_ac_indicies);

    c_indicies = t.search("C"_dna5);
    EXPECT_EQ(c_indicies, (std::vector<size_t>{4, 0, 2}));

    auto indices = t.search("ACGT"_dna5);
    EXPECT_EQ(indices, (std::vector<size_t>{}));

    EXPECT_EQ(t.search(""_dna5).size(), 6);
  }

  std::vector<lst::LazySuffixTree<seqan3::dna4>> dna4_trees{tree4,
                                                            tree4_parallel};
  for (auto &t : dna4_trees) {
    EXPECT_EQ(t.search("ACGT"_dna4),
              (std::vector<size_t>{0, 4, 8, 12, 16, 20}));

    EXPECT_EQ(t.search("ACGTACGTACGTACGTACGTACGT"_dna4),
              (std::vector<size_t>{0}));
    EXPECT_EQ(t.search("CGTACGTACGTACGTACGTACGT"_dna4),
              (std::vector<size_t>{1}));
    EXPECT_EQ(t.search("GTACGTACGTACGTACGTACGT"_dna4),
              (std::vector<size_t>{2}));
    EXPECT_EQ(t.search("CCCC"_dna4), (std::vector<size_t>{}));
    EXPECT_EQ(t.search("ACGTCT"_dna4), (std::vector<size_t>{}));

    t.expand_all();
    EXPECT_EQ(t.search("ACGT"_dna4),
              (std::vector<size_t>{20, 16, 12, 8, 0, 4}));
  }

  sequence_t<seqan3::dna5> long_sequence{
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

    std::set<std::tuple<std::string, size_t>> expected_labels{
        std::make_tuple("A", 2),    std::make_tuple("C", 3),
        std::make_tuple("", 1),     std::make_tuple("AC", 2),
        std::make_tuple("CA", 2),   std::make_tuple("C", 1),
        std::make_tuple("ACA", 1),  std::make_tuple("AC", 1),
        std::make_tuple("CAC", 2),  std::make_tuple("ACAC", 1),
        std::make_tuple("CACA", 1), std::make_tuple("CAC", 1),
        std::make_tuple("ACAC", 1), std::make_tuple("CACAC", 1),
        std::make_tuple("CACAC", 1)};

    auto labels = t.get_all_labels();
    std::set<std::tuple<std::string, size_t>> labels_set{labels.begin(),
                                                         labels.end()};

    EXPECT_EQ(labels_set, expected_labels);
  }
}

TEST_F(LazySuffixTreeTest, SuffixLinks) {
  std::vector<lst::LazySuffixTree<seqan3::dna5>> trees{tree, tree_parallel};
  for (auto &t : trees) {
    t.expand_all();

    t.add_suffix_links();
    std::vector<size_t> expected_links{
        max_size, // root
        2,        // AC
        0,        // C
        0,        // -
        9,        // ACAC
        7,        // AC
        1,        // CAC
        3,        // C
        4,        // CACAC
        5         // CAC
    };

    EXPECT_EQ(t.suffix_links, expected_links);

    t.expand_implicit_nodes();
    t.add_suffix_links();

    std::vector<size_t> implicit_expected_links{
        max_size, // 0: root
        0,        // 2: A
        0,        // 4: C
        0,        // 6: -
        6,        // 8: ACA
        7,        // 10: AC-
        1,        // 12: CA
        3,        // 14: C-
        4,        // 16: CACA
        5,        // 18: CAC-
        2,        // 20: AC
        13,       // 22: ACAC
        9,        // 24: ACAC-
        10,       // 26: CAC
        11,       // 28: CACAC
        12        // 30: CACAC-
    };

    EXPECT_EQ(t.suffix_links, implicit_expected_links);
  }

  dna_tree.expand_all();
  dna_tree.expand_implicit_nodes();
  dna_tree.add_suffix_links();

  std::vector<size_t> implicit_expected_dna_links{
      max_size, 0,  0,  0,  0,  0,  2,  4,  3,  4,  1,  2,  2,  5,  10, 11, 8,
      14,       42, 43, 44, 45, 46, 47, 12, 38, 39, 40, 41, 13, 7,  24, 25, 26,
      27,       28, 9,  29, 8,  15, 48, 49, 30, 31, 32, 33, 34, 35, 36, 37};

  EXPECT_EQ(dna_tree.suffix_links, implicit_expected_dna_links);

  dna_tree_parallel.expand_all();
  dna_tree_parallel.expand_implicit_nodes();
  dna_tree_parallel.add_suffix_links();

  std::set<size_t> implicit_expected_dna_links_set{
      implicit_expected_dna_links.begin(), implicit_expected_dna_links.end()};

  std::set<size_t> suffix_links_set{dna_tree_parallel.suffix_links.begin(),
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
        std::array<size_t, seqan3::alphabet_size<seqan3::gapped<seqan3::dna5>>>>
        expected_reverse_links{
            {1, 2, 3, max_size, 4, 5}, // root
            {max_size, max_size, max_size, max_size, max_size,
             max_size}, // A, 2
            {max_size, max_size, max_size, max_size, max_size,
             max_size},                                            // C, 4
            {max_size, 8, max_size, max_size, max_size, max_size}, // G, 6
            {max_size, max_size, max_size, max_size, max_size,
             max_size},                                             // T, 8
            {max_size, max_size, max_size, max_size, 13, max_size}, // -, 10
            {max_size, max_size, max_size, max_size, max_size,
             max_size}, // ACGATCGCT-, 12
            {max_size, max_size, 10, max_size, max_size,
             max_size}, // ATCGCT-, 14
            {max_size, max_size, max_size, max_size, max_size,
             max_size},                                             // CG, 16
            {max_size, max_size, 11, max_size, max_size, max_size}, // CT-, 18
            {max_size, 14, max_size, max_size, max_size,
             max_size}, // GATCGCT-, 20
            {max_size, 15, max_size, max_size, max_size, max_size}, // GCT-, 22
            {7, max_size, max_size, max_size, max_size, max_size}, // TCGCT-, 24
            {max_size, 9, max_size, max_size, max_size, max_size}, // T-, 26
            {6, max_size, max_size, max_size, max_size,
             max_size}, // CGATCGCT-, 28
            {max_size, max_size, max_size, max_size, 12, max_size} // CGCT-, 30
        };

    std::set<
        std::array<size_t, seqan3::alphabet_size<seqan3::gapped<seqan3::dna5>>>>
        reverse_links_set{t.reverse_suffix_links.begin(),
                          t.reverse_suffix_links.end()};

    EXPECT_EQ(reverse_links_set, expected_reverse_links);
  }
}
