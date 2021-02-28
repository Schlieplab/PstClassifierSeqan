#include <gtest/gtest.h>

#include "../src/kl_tree.hpp"
#include "../src/probabilistic_suffix_tree.hpp"

#include <array>
#include <map>
#include <string>
#include <tuple>
#include <vector>

#include <robin_hood.h>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/io/sequence_file/input.hpp>
#include <seqan3/range/container/bitcompressed_vector.hpp>
#include <seqan3/range/views/to.hpp>
#include <seqan3/range/views/to_char.hpp>

#include "random_sequence.hpp"

class ProbabilisticSuffixTreeTest : public ::testing::Test {
protected:
  void SetUp() override {
    using seqan3::operator""_dna4;
    using seqan3::operator""_dna5;
    sequence = lst::details::sequence_t<seqan3::dna4>{"GATTATA"_dna4};
    probabilisticSuffixTree = pst::ProbabilisticSuffixTree<seqan3::dna4>{
        "TEST", sequence, 2, 1, 192, "parameters", false, 2};
    probabilisticSuffixTree.construct_tree();

    probabilisticSuffixTreeParallel =
        pst::ProbabilisticSuffixTree<seqan3::dna4>{
            "TEST", sequence, 2, 1, 192, "parameters", true, 2};
    probabilisticSuffixTreeParallel.construct_tree();

    long_sequence = lst::details::sequence_t<seqan3::dna5>{
        "CGCGCGGCGCGTCGTCGAGCTCGCCGCCCAGGGCGCGCCCCTCGGCGCGATCCGCGCCGCCCTCAACGACGTGACCCCCGCCGCCGACAAGGGCGAAGCCTTCGTCGGACGCCCCGACTGGCTGGGCGAGCTGTGGACCGCGCGCCGCACGGACCGGCCCCTGATCGACTCGATCACCAAAAAGGCCCTGCCCCGTGCGACCAAGGTCAAGGGCTGGCGCTGGAAGAAGCGCCCCGAGGTCGCCGACTACACGGGCAACAAAACCGAAATCTCCTCCAACGAGATCGAGACCGAGCCCGTCGAGGCCGCCGTCAAACGCATCGCCGCAGGGTGGGACACCGACCGCATTTTCGTCGACCTCGGCGACGGCGACATGATCGAGAGCCTGTGGGAGGGCGCCCGCGAGGACTACGCGATCAAGACCGAGGCCGCCGTCACCACCGGCCTCAAGACCGCAGCGACGAAGCTCACCGGCGCGCCCACCGAACTGGACAAGGCCCTCGTGGTCCTCGGCTCCAAGGCCGCCGCCATCGGCTCCCGCCTGAACTTCGTCGCCTTCGGCGCCGACGTGTGGAGCAAGTTCACCGCGCTGACCCGCGACCAGGTGCCGTGGTGGATCACCAACGGCGACCGCCTCAACCTCTCGACCGCGACCGGCGAGGTCAACGGCCTGCGCCTGTTCGTGGACCCGACCCTCGCGGCGGGCGACATCCTCGCGGGCGACACCCGGTCCGCGACGTTCTACGAGGAGCCGACCCCCATCAGGGTCAACGCCATCGACCTGCCCAAGGGCGGCGTGGACCTCGGCCTGTTCGGGTACCACGCCCTTCTCGTGAACGACCCGAACTCGTTGTTCATCATCACGGCGGGCTGACCCCATGACCCCCGACGACCTCGCCACGCGGGCCGCCGCGTGGGCGAAGCTCCCGGGCGGCGTGGACGACGCCATGAGGGCGTGCGCAGCCGCAGTGCACGCCCTCGTGGCCGCCCTGCCCGTCACGCAGGGCCGCCCCGCCTGGCGCGAGGACACGGCCCTCGGAGCGGTCATGCTCACCGCCCGCCTGCACCGCCGCCGCAACAGCCCGGCGGGCATCGAGTCCCTGACCGAGATGGGCGCGACCTACGTGAGCCGCTACGACAGCGACATCGCGCGCCTGCTGCGCATCGACGCCTTCGTCGGGCCCGTCGCCATCTGAGGGGGGCCACGAGATGAACCCGCTCTACGCGGCCGCCCAGGATGTGGCCGACATGCTCGCCGCGGCCGGAGTCCACACGGTCACCGACCCCAGGGACATCGAGCCGCCGTGCGCGTGGGTCAGCCCCAGCCGCATCGCCTACCCCACGCTCGCTGGCCGCCCCCGCACCGTCGAGTGGGAGGTGTACCTCATCGCACCCGACAGCGGCGCGCCCCTCTTCCCCCTCGGCGACCTCATCGACCGGGCCGCCACCGTCTTCCCCGGCATCGAGGCCCGCACCCTCGGCCTGACAATCCCCAACCTCAGCCCGGACCCCCTGCCCGCGATCACGTTCACCATCGAAACAGAAACGGACTAAACCCATGGCAGTGAAAACCCTCACCCTCGGCCCCGGCAAACTCAGCTTCGGCGCCCCCGAGTCCCTGACCCACGCCGCCGCCCAGGTCACCAAATGCGCCGTCAAGCCCACCGCGAAGCAGGGCGACTCCGTGGCCGTCCTGTCCGGCGACCGCGTGCCCGGCGACCGCACCGAAGCCGCGACCCTGGAGTTCACGATCTACCAGGACTTCGGCGAGGCCGAATCCTTCGTCGAATGGACCTGGGCCAACGCAGGGAAGGAACTCCCCTTCGAGTTCATCCCCGCCGACAAGCACGACAAGGCCGTGCGCGGCCGCGTCACGATCGAGCGGTCCGACATCGGCGGCGAGGTCGGCGTCAAGGTCACCGCCGACCTGGAGTTCACCTGCACCACCATGCCCACCATCGAGCCCAAAACCAAGATCGGGCACTGAGGTGGCCGACTACTCCGGCGTCAAGATCGACGGCGCGCGCCGCCTCCGCTCGACCCTGCGCAAAGCGGGCGCGGACATGCGCGACATGCGGGAGGTGAACCGCGTCGTCGCCGGCATCGTCGTCGGCGCGGCCACCGCCCGCGTCCCCCGACGCACCGGGGCCCTGGCCGCCACCGTGCGCGCAGGGGCCACCCAGGCCGCAGCCATCGGCCGCGCCGGGAACAACCGCCGCACCGGCGTCCCCTACGCCAACCCCATCCACTGGGGATGGCACCGCCACCGAATCCGCCCTAACCCGTTCCTCAGCCTCGCCGCCCAGGACACCGAACCCCAGTGGTTCGGCGTCTACGCCGACCGCATCGAACGCCTCATCAACAGCATCGAAGGAGCCTGACCCATGTCGAGCATCAAAGCCATCAACGTCGAGGTAGTCACCTCCGCCGTGACCGGCGACCTCGCCGCCGTCACCGTCCGCACCGACAACCGCGACCGCATCGCCTGGGACCTCGCGAGGGGCCGCAACAAATGGCCCCAAGCACAGGAGGCCCCCAGCCTGTGGGCCACCCACATCGCCTACACCGCGCTCCGCCGCACCGGCGAAGTCAGCTGCTCGTTCGAGGAGTTCTCCGAGGCAACTGTGAGCGCCGAACCCGAGGTCATCGACGTGGACCCTACCCGGACGGCGACCGCCGGGGCCTGATCGTCGCCCTGGCCCTCGCCACCCGCATCCCCATGAGCGAGTGGGAGACCCGCCCCGACGAGGACATCGCCACCGCACTGCAACTGCTAGAAGAGAGGAGGAGCTGACTTGGCGTCGAAAACCGCCATCCTGAGCGTCCGCGTCGTCTCCGACGTGAAGGACGCCACCAAGGGACTGGACGACGTGGCCGACAAGACCGGCCGCCTGGAGGACGGCCTCAAACGGGCCGCCGCCCCCGCCGGGATCGCCGTCGCCGCCCTCGCCGGGATCGGCAAGGCCGCCACCGACTCCGCCAGCGAGTTGCAGCAGAGCGCGGGCGCCGTCGAATCCGTATTCGGCGGGCACGCCGCCGCCGTCCAGGACGCCGCCAAGACCGCCGCCTCCAGCGTCGGCCTGGCAGCAAGCGAGTACCAGAACATGAGCGCGGTCCTGGGCGCCCAGCTCAAGAACATGGGCACCCCCATGGAGGACCTGGCCGGATCGACCCAGAACCTCATAGGCCTGGGCTCCGACCTCGCCGCCACCTTCGGGGGAACCACCGCCGACGCCGTGAGCGCCATCTCAGCCCTCCTCCGGGGCGAGCGCGACCCCATCGAGCGCTACGGCGTCTCGATCAAACAGTCGGATATCAACGCGCGTCTGGCCGCCGAGGGCATGGACAAGCTGGAAGGCGCGGCCAAGACCCAGGCCGAAGCCCAGGCCGCCCTCGCCCTGCTCACCGAGCAGACCGCATCTGCGCAAGGCCAGTTCGCGCGCGAGACCGACACGATGGCCGGGAGCCAGCAGATCGCCGCCGCCCAGTTCGAGAACGCAAAAGCCGCCCTCGGGGAGAAGCTGCTGCCCGTCGTCACGCAGTTCATGGAGGCCATGAGCGGGGCGGCTCAATGGGTCGCCCAGAACAGCGATGCGCTGCTCGTCCTCGGCGGCGTCGTCGGAACCATCGCGGGCGTGATCCTCGCCGCCAACGCCGCCATGGGCGTGTGGACCGCAGTCCAGACGACCGCCAGAGTCGCGACGGCCGCCTGGACCGGCGTCCAGGCCGCGTTCAACGCGGTCATGGCCCTGAACCCGATCACACTGGTGGTCATCGCCATCGGGGCCCTGGTCGCCGCCGTCGTCGTCGCCTACAACAAGTCCGAGGCGTTCCGCAACGCTGTATCCGCGCTGTGGGACGCCATCAAAGCGGGGGGCGGCTGGATAGTCGATCACGTCATCAAACCCATCGGAGACGCTTTCAACGCCGTTGTTGATGCCGTGAAATCCGTCTACGACTGGGTGAAGAACCTGTTCAGCGGATTCCAGCTCCCAGGCTGGCTATCGAGCGTACTCAGCTGGTTCGGGCTCGAAGCCCCCTCCGGCCCCGAATCCGGGGCCATCCTCGCCGCCACCGGGACAACGGACGCCCCCCTCGCGCGCCTCGCATCGTGGGCGCTCGCCCCCCGCACCGGCTCCAGCCCCACCCCCGCCGGCAGCGTCGTGAACATCACTGTGAACGGCGCCCTGGACCCTGACGCGGTGGCCCGCCAGATCGGCCGCATCCTGTCGCGCCGCGACCTCATCAACGGCACCGAGCAGATCGTGGGGGCGACCCTATGAGCGTCAGCGCAAGCCTCCGCGTCGCCGCAGGCGGCCTCGGCGGTGTCATCAACGCCGCCGCCGACAAGTACCCGACGACGGTGACCGTCCTCGACGACCTCACCGTCACGTGGGGCCGCGATAGCGTGGTCTCACACCCGGACCCGTCGTCAATGACCGCGACAATCGCCCTCGTGGACACTGTGCCCGACTGGCTGCGCGTCGGCGCGCTGGCGACCGTCAACGCGGTCGCCCGGACCGAGGAGTCCCAGCGGTCCTATATGCGGCTCCTGCCCTGGCGCGCCATCGAACCGGGCACCGGGTGGCGCCAGCAGGTCACCCCCGATCCGCCCGGCGCGTGGGTCGGGAGCCTCCCGGTGTTCGCCGCCGCCGGAGCAGACTCCGGGATCGGATGGTTCATAGCGCCAGGGGTGCAGCCCCCCTCCGATACCCCCGAGACCACCCAGTGGGCCGCGAACGCAAAGACGACAGCGGGCAAGCCAGTCACGTTCACAATCACCGTCCCCGAGCTCACCGGCGCGACCGTGCGCGCCTTCCCGCTCACTTACCGGCGCCCCGGCGGGCTCTACACCCGCGCCCCCGGCATGGCGATCGAGCTCTCCCCGGAGAAGTACACGCCGGGCACTGTCGAGTACTCGGGCACCTGGACCCCCGAGGCGACGGGCCTTTACGTCGGCGCCTATCTGCATATCCAGCTGCACAAGGCCCCGGCCTGGACGACGATCCCCCGCGAGCGGACCTGGCGGGCCGCACCCGGCACGTGGGCGGACGCCGGCGGCCGCGCAACCGTCACTGACGTTCACATCGCCGGGACCTCCGGCCACGTCGCTGAGCACGCCGTCGAGGTCTTCACCGGCCGGGTCCAGTCCCTCCGCGTCGAATGGTCCGAGCGCCTGTCCCGACCCATCGCGCGGATCACCGCAGTGGACAAGCTCGCCGATCTGAATGGTACCTACATCGGCGACACGCCGTGGGGCGAGGAATCCTGGAAGCTGCGCGCAGAACGCATCCTGAAACAGGCCCTCGGACCAGCAGACACACTGGAGGGCGAACCCGGCAATTGGCTGGGGACGATCCGCCCCCGAGACGTGGACCACCGCAGCGCCGGTGAGCTCATGAGGAACACACTCGCCTCGTGCGCCGCCACCGCATTCCCCGTCAGCTGGCGCAAATGGCGGGTAATCCCATTCATCTACAAGGGGAGCGATCAATCAATCACGATCCCCGGGCGCGCCATCCGCCGCGACGGCGTACAGGTCAGCACTGACGAATCCGCGAACATCTCGACAATTCAAGCCACGTATTTCGATGTGACCTACGACGGGAAAACCGGGAGGGTGAAAGACGTTATAGAGCGCACGACTACGCGGAAGAATACACCGGCAAATGAGGGACCCCCCAGGTCCATCAAAATGAAAACCGAACTATCCCGCAGTAACGAAGCGAGCGAACTCACCCGGATCATCGGGAAATACGTGAACGTGAACCAGTGGATCATCAGCGCCCTATCAGTGAAGCACGACCGGATCAGCGAGGACGCCCTTGTGCGCCTACTGTCCGCCACCGAGCGCATCGCCCAGCAGGTCGTCCTCACCGGCCTCCCACGATGGTTCCCAGCAGCGACGATGCGCGGCATCGTCATCGGCGGGTCCCTCACCATGCACCGCGGCCACTGGACCCCCACCCTCCGCATCGCAAACACACCCGACTAGAAAGAGTACCCATGCCATCGACAACCCCCCGGGGTCTTCCCTACGCAATCCCCACCGACGCCCAAGCCGCATTCCCCGACGCCGTGTCAAAACCCATCGCCGAATGGATCGAGGCGAACCTCCCGGTCATGCAAGCCGGGACCATCGCCTACCCCGCCCTCGGCTCCCAAGACCAGACGGGAGAGTACACAGTCACGTTTCCCAAACCCTTCCCAGTCACGCCCCGGATATTCATGCAAGCCGATAACCAGCGCCTCACAATCGCCGTATGGAATATCAGCCGCACCGGGTTCAAATGGATGGCCCGCAACAACAGCAACGGCAATTCGTCCTCTGGAGCGGCCTCGTGGTTCGCCGTTAGCGGCGCCACCGGACAGTAACGAAAGGAAACAAAGGAAATGACCACAGCGGTCGACGTGTTCACCGCCCGCCTCGCCTGGATGATGACCCAAGCCGACGGCGGCTACTCCCAGCCCAACCGCCTCGACGTGCGCCGCACGCGCGGCGTGTGGGACCCCGGCTTCCAGTTCGAGGGGGACTGCTCCTCCTGCGTCCTGGAGGCCGCCCACCAGGCGGGCCTGCCCACAGGCTCTGCGTCCTACACGGGCGACATGCGCGCGGGTCTGGAGGCCGTGGGATGGGCCGTCATCCCCTACGCCGCGACCGGCGGGGACCTTGACAACCTCGCCGACGGCGACGTGCTCCTATCCGAGGCCGCGAGCGGCGGCGTCGGCCATACCGGCGGCCTCATCCCCGGCGGCCTCGTCGCCGAGGCGTGGATCGACGGTCACGGAGACATCATGGGCTCCGCAGGCGGGGACGGGCCCGGCGACGACACCGGCGGGGAAACCCGCGCAGTGCCGTTCTATTCCCACCCCTACACAGTGCGGGGGCTCTGGACGCACGTCCTGCGCCCCCCAGCCCTCGACGCCGCAGACTCGCCCGCCGAACCCACCCCCACAACGAAAGGAATCCCCAATATGTTCGGAATCACCTACACGGCAAACGCCTTCGGCGGTATCACCGCCTACGTCCTCATCCACGAGTCCGCCGGTGCCGACGCCCTTGACCGCGTTCAGGCCCAGGTGTACAACAGCGTCCTTCCCAACGGCTTCACCGAGGTCCCTGAGCACCACGCCGAAATGCTCATCCGCGAGTCGTGGGTGCGCCACAACCGCATCGCCAACGCCGTCGCCGCGACCACTCGCGTAGACATCAACGAGGCCACCGCCCGCGTTCTCGCCGCCGTCAAGGAAGGAGCTGCCAAGTGAACGCCATCACGTCCCAGACCCCCGATGATCCCACGCCGCAGCCCATCTCCTGGCTTACACCCGCAGTGCGGCGTTACATCTACAACGTCACTATCGCCGCCCTCGGCGTCGCTCTCGTCTACGGCGTCGTCGATGGCCAGCACGCCGCCGCCTGGGAGGCCCTGGCCCTCGCTGTCGTGGGCCTCGCCCGCGCGCACGTCCCCGGAGACCCCCAATGAGCGACGCCTCGGCGGCCGTCGAGGTCATCGCCGCCATAGGCGGCCTCGGCGGCCTCGGGGCCGCGCTCTCAGCCGTCGCCTCCCTCATGGAGGCCCGCAGGGTCCGCGCCAGTATCCCCGCCGCCGCCGACCGCACCGAGGAGGCCATCGACGCCCTGCGCTCCGATGTCCGCGCCATCGACCGCCGCATCGGACACGAGCTCGGCGAAATCCGCCGAGCCGCCGACCGGGAACACGCCGACTATGACGCACGTCTCAGACGATTGGAGGGGTCATGAGTTGACACTCCCACTCAGGTGGTTATAAACTATAACCATCAGGAAGCCATAAGGGGCAAGCCTGAAACCTGAAGGGAGCGCGAAAATGCGCAAGTCAATCGAACTCACCTGGACCCCCGAAACCCGCGTTTGGGGCACAAGTGGGAACACCAGCGTGGCCGTCGGCACCGGAACCCTGGACGGCCGCCGCCTCGCCGTCTACGCCTTCCCCCAGTCCGATCACTGGTCGTTCTGGTCGCAAATCGAGCGTCCCGGCGGTGGTTCGACCTCGATAGAAATCAGGTCCACATTGCCCGCAGGCACCGTTCCATCGGTCCTCGGCCCGAACGGCGCCATCCGCGAAACGACGACAATCGAACTCTGACCCAACGCCCCGGCGCACACGCGCCGGGGCACTGTCACGGAAGGGCAAAGGATGAACGGCATCGAACTGCGCGCCCGCCGCGAAGCACTCGGGCTCTCACAAACCAAGTTCGCGAAAATGTGCGAAACGACTCAAGTGACCGTTTCCCGCTGGGAGAACGGCACCCGCGAACCAAGGAACGACATTGCAATACACCTGCTGATGGCAAATATCGAAGACGCCGCTATCGACCTCATCGAGGACCTGCTAGAGCTCGCCGAAGACGAAGAACTCCTAACAGCAACGCCCGACCTCCAGCTCACGGTGTACAACGACGAAGCCCGCTACGCCGCCGGGGAGCCCGTCTGGTCAAAACGCCTCCCCATGGAGACGCACCGCGTGTGCGCCGCCCGAGCCGCCGCCCTCCTCGGCGCCGAGGACGGCACACACGTCACACTCATCGAGGGCTGAGCGCCCTCACGCGATAGCCTCAACGACCTCGCGCACATCATCGTCGGGAACGAGGATGTAGCGAAGCGTCGTCGAGGGTGAAGCGTGCCCGAGCGCGCGCTGCACCGCGACGAGGTTCCTCGTGCGCGCGAACCCCGTCGAGGCGAAGGCGTGCCGGAGGGCGTGCATCGTCACGCCCTCCGGCAGCGCCCGGCCGACCAGCTTCCCCACCCACGCCGGGGAGAGGTGTCCATGGTCCGCCCCCGGGAAGAAGAACCCCGGGTCATGATCGAGCAGCTCATCGGCGAGAGAATGCGGGAGGGGGATCACGCGGGTTTTCCCGCCCTTTCCGTGGACCACGAGGGACCAGCCCGCCAGATCGCGCACAAGATCGCGCGTGTGCGCGCGGGCGACCTCACCCCGCCGCATCCCCAACTCCGCCGCCATGCGCACCATGAGACGCACTCGTGGATCCGTCGCCCGGCGCCCCACCGCGATGGCGCCCGGCGTCGCCGGCCTCGGCGCAGGGTCCGACTGTCTCACCGACGGCACCGGCGGCGCTACCTCGATGTATCCGACCCCCTGGGCCCACCGGTAGAACTGGTCGACGCTTTGATGCGCGCTCCGACGCGTATCCCGCGCCCAATCATGCGCCCCGGACCACTCGATCACCGTGAGCGGCCCAACCTCCCACGGGCCCGCCCGCAGATCGCGGGCGAACCTACTCACCCACTCGATCCGCAGTCGGATAGTCTCGGCCCGCCGGCCGGCCGCCGCCAGCGCTGTAGTCCACTCGCCTATAGGACCGGCCCATCCGGCGGGTACCGGTCGCGGTTTCATGCTCATGATGATTACCCTGCATCTATTCGGCCCCAGTGTCGCGTCATCACGCCGCCGGGACCCACGATGTAGGATCGAGCCGTGGGCCACCAATCCGTAGGTTGGGGGTTCGAGTCCCCCTGGGCCTACTCGCCCGCCCCAGCCCCGCCGGCTGGGGCGGTTGCCGCACTGCGCCGGTACCGACGCGCTGCCGACGGATGGAGGTCGCCATGTCCGACGACGCGGACAGGGCCCACGGCGCCCTGGCGGGACTCGCGCTGGGGGACGCCCTTGGCATGCCGACCCAGGCGATGACCGCCGATCAGATCAGGCTGACCTACGGGTGGGTGGACGCCCTGGTGCCAGCCGACGCCTCGCAGCCCTACGCGCCCGGCATGCCCGCCGGCAGCGTCACGGACGACACGGAGCAGGCGCTGCTCGTCGCCGGCCTGCTGGTATCGGGCGGGGGCGGCATCGACCCCCACGCCTTCTCCCGCGCCCTGCTGGACTGGGAGGACTCGATGGCGGCCCGCGGTTCCCTCGACCTGCTGGGCCCATCGACGAAGGCCGCCCTGGAACGCGTGCGGGCCGGGGAGGACCCCCTCCGCGTGGGCGGCGCGGGCACCACCAACGGCTCAGCGATGCGAGTCGCGCCCGTCGGGATCGCCTCCTCCACCCGGGATCCGCGTTTCGCCGACACCGTGTGGGAGTCGTGCCGCGTCACCCACGCCACCGAACAGGGCTTCCACGCCGCCGCGCTCGTGGCGGCAGCGGTCTCCCTCGGCATCGACGGAGCAGGGGCGGACAGCCCTTCGGACTCCGCCCGCGCCTCCTTGGAACGCGCCCTGGCCCTCGTGGAGGCACTCGGGCGCCGGGGGGCGCGGACGCCCCAGCCGGACGTGTGCGAGCGGACCCGCTACGCGCTGCGGTTCGCGCGCGCCCGCGACCCCGCCCCCGGTACTGCCGACGACGACCGGGCATTCGCCGGGGCACTGCGGGCACGCGTCGGCGCCTCCGTGGAGGCCGCCCAGTCCGTCCCCGCAGCATTCGCCATCGCCTGGCGCTACGCCGCCGATCCGTGGCGGGGCCTGTGCGTCGCCGCCAACCTCGGCGGTGACACCGACACGATCGGCGCTATCACCGGCGCCGTGCTCGGCGCCGCCCTGGGGGCCCGGTGCTGGCCCGCCCAGGAGCTGGAACGAGTGGAGGCCGTCTCCGGGCTGCGGCTGCGCGAGACCGCCGACGGTTTGCTCCGCCTGCGCGCCCACGGATCCCGACTGCCCGCCCACGGGGAGCCGGTCGCAGCACCGCAGGAGGGCAGGGTCGTCCTGCTCGGGCAGGTGGTCGTCGACCTCGCACTGCTGGCGCCGCGCGTGCCCGCTCCCGGCGGCGACGTGTTCGCAGAGGACGCGGGCATGCACGCGGGCGGGGGCTTCAACGTGCTGGCCGCTGCGCGCCGGATGGGAGCGGAGGCAGTGAGCCTGTCCGGCGTCGGGGACGGCGGATTCGCCTCGATCATCACCGCTGCGTTGGAGCGCATCGGCGCCTCCTGCGAGGGACCGCGCGTCGCGGGAACGGACTCGGGGTACTGCGTGGCCATCACGGACGGCGACGGCGAGCGCACCTTCGTCTCGACCAGGGGCGCGGAGGCCCGCCTGCCGCGCGGGTCGTGGTCCGCCCACGCGGCCCGCTTGCGCAGCGGGGACGTGGTGCACGTGGACGGCTACGCGCTGGCCCATCCGGCCAACACCGCAGCGCTGCGGGAGTTCCTCTCGGCGCACCTGCCCGCAGGGCTCCGCGCGATCGTCGACGTGTCGCCCGTCGTCGGCGATGTGGACCTCGACGACCTGCTTGCCCTGCGGGCCCTGGCCCCCCTGTGGTCCATGAACGAGCGCGAGGCGGGGATCCTCGCGGGCCGCCTCGCGCGGGCGTCCGCCGCTCCCCCGCACGGAGGCGCTCCCCCGGGGGAGGCGACACCACCGGCCGGAGCGGCCCCCGGGA"_dna5};
  }

  lst::details::sequence_t<seqan3::dna4> sequence;
  lst::details::sequence_t<seqan3::dna5> long_sequence;
  pst::ProbabilisticSuffixTree<seqan3::dna4> probabilisticSuffixTree;
  pst::ProbabilisticSuffixTree<seqan3::dna4> probabilisticSuffixTreeParallel;
};

TEST_F(ProbabilisticSuffixTreeTest, ConstructorStatus) {
  std::vector<bool> expected_status{
      true,  // root
      true,  // A
      true,  // G
      true,  // T
      false, // -
      true,  // AT
      false, // A-
      true,  // GA
      true,  // TA
      true,  // TT
      false, // ATA
      false, // ATTATA
      false, // GATTATA-
      false, // TATA
      false, // TA
      false  // TTATA
  };
  probabilisticSuffixTree.print();

  std::vector<bool> tree_status{};
  for (auto &v : probabilisticSuffixTree.entries) {
    tree_status.push_back(v.included);
  }

  EXPECT_EQ(tree_status, expected_status);
}

TEST_F(ProbabilisticSuffixTreeTest, ConstructorSuffixLinks) {
  std::vector<size_t> expected_suffix_links{
      max_size, // root
      0,        // A
      0,        // G
      0,        // T
      0,        // -
      6,        // AT
      8,        // A-
      2,        // GA
      2,        // TA
      6,        // TT
      28,       // ATA
      30,       // ATTATA
      22,       // GATTATA-
      20,       // TATA
      12,       // TA
      26        // TTATA
  };

  EXPECT_EQ(probabilisticSuffixTree.suffix_links, expected_suffix_links);
}

TEST_F(ProbabilisticSuffixTreeTest, ConstructorProbabilities) {
  std::vector<std::array<float, seqan3::alphabet_size<seqan3::dna4>>>
      expected_probabilities{
          {4.0 / 11.0, 1.0 / 11.0, 2.0 / 11.0, 4.0 / 11.0}, // root
          {1.0 / 6.0, 1.0 / 6.0, 1.0 / 6.0, 3.0 / 6.0},     // A
          {2.0 / 5.0, 1.0 / 5.0, 1.0 / 5.0, 1.0 / 5.0},     // G
          {3.0 / 7.0, 1.0 / 7.0, 1.0 / 7.0, 2.0 / 7.0},     // T
          {0, 0, 0, 0},                                     // -
          {2.0 / 6.0, 1.0 / 6.0, 1.0 / 6.0, 2.0 / 6.0},     // AT
          {0, 0, 0, 0},                                     // A-
          {1.0 / 5.0, 1.0 / 5.0, 1.0 / 5.0, 2.0 / 5.0},     // GA
          {1.0 / 5.0, 1.0 / 5.0, 1.0 / 5.0, 2.0 / 5.0},     // TA
          {2.0 / 5.0, 1.0 / 5.0, 1.0 / 5.0, 1.0 / 5.0},     // TT
          {0, 0, 0, 0},                                     // ATA
          {0, 0, 0, 0},                                     // ATTATA
          {0, 0, 0, 0},                                     // GATTATA-
          {0, 0, 0, 0},                                     // TATA
          {0, 0, 0, 0},                                     // TA
          {0, 0, 0, 0}                                      // TTATA
      };

  for (size_t i = 0; i < probabilisticSuffixTree.entries.size() &&
                     i < expected_probabilities.size();
       i++) {
    for (size_t j = 0; j < 4; j++) {
      EXPECT_FLOAT_EQ(probabilisticSuffixTree.entries[i].probabilities[j],
                      expected_probabilities[i][j]);
    }
  }
}

TEST_F(ProbabilisticSuffixTreeTest, PrunedKL) {
  auto kl_tree = pst::KullbackLieblerTree<seqan3::dna4>{
      "TEST", sequence, 3, 2, 0.3, 192, "cutoff", false, 2};
  kl_tree.construct_tree();

  std::vector<bool> expected_status{
      true,  // root
      true,  // A
      false, // G
      false, // T
      false, // -
      false, // AT
      false, // A-
      false, // GA
      false, // GATTATA-
      false, // TA
      false, // TT
      false, // ATA
      false, // ATT
      false, // TAT
      false, // TA
      false, // TTATA
      false, //  ATA
      false, // ATTATA
      false, // TATA
  };

  std::vector<bool> kl_status{};
  for (auto &v : kl_tree.entries) {
    kl_status.push_back(v.included);
  }

  EXPECT_EQ(kl_status, expected_status);
}

TEST_F(ProbabilisticSuffixTreeTest, PrunedParameters) {
  auto kl_tree = pst::KullbackLieblerTree<seqan3::dna4>{
      "TEST", sequence, 3, 2, 0.0, 6, "parameters", false, 2};
  kl_tree.construct_tree();
  std::vector<bool> expected_status{
      true,  // root
      true,  // A
      false, // G
      false, // T
      false, // -
      false, // AT
      false, // A-
      false, // GA
      false, // GATTATA-
      false, // TA
      false, // TT
      false, // ATA
      false, // ATT
      false, // TAT
      false, // TA
      false, // TTATA
      false, //  ATA
      false, // ATTATA
      false, // TATA
  };

  std::vector<bool> kl_status{};
  for (auto &v : kl_tree.entries) {
    kl_status.push_back(v.included);
  }

  EXPECT_EQ(kl_status, expected_status);
}

TEST_F(ProbabilisticSuffixTreeTest, Print) {
  auto pst_unpruned = pst::ProbabilisticSuffixTree<seqan3::dna4>{
      "TEST", sequence, 10, 1, 192, "parameters", false};
  pst_unpruned.construct_tree();
  pst_unpruned.print();
  seqan3::debug_stream << std::endl;

  probabilisticSuffixTree.print();
  seqan3::debug_stream << std::endl;

  auto pst_pruned = pst::KullbackLieblerTree<seqan3::dna4>{
      "TEST", sequence, 3, 2, 1.2, 0, "threshold", false, 2};
  pst_pruned.construct_tree();

  pst_pruned.print();
  seqan3::debug_stream << std::endl;
}

TEST_F(ProbabilisticSuffixTreeTest, MemoryAllocationException) {
  // Crashes when allocating memory in expand suffix links
  EXPECT_NO_FATAL_FAILURE({
    auto pst = pst::ProbabilisticSuffixTree<seqan3::dna5>(
        "TEST", long_sequence, 15, 4, 48, "parameters", false);
    pst.construct_tree();
  });
}

TEST_F(ProbabilisticSuffixTreeTest, CorrectNumberOfParameters) {
  size_t sought_n_parameters = 192;

  pst::ProbabilisticSuffixTree<seqan3::dna5> pst{
      "TEST", long_sequence, 15, 4, sought_n_parameters, "parameters", false};
  pst.construct_tree();

  size_t n_terminal = pst.count_terminal_nodes();
  EXPECT_EQ(n_terminal * 3, sought_n_parameters);
}

TEST_F(ProbabilisticSuffixTreeTest, PSTBreadthFirstIteration) {
  std::vector<size_t> visited{};

  probabilisticSuffixTree.pst_breadth_first_iteration(
      0, 0, [&](size_t index, size_t level) {
        visited.push_back(index);
        return true;
      });

  std::vector<size_t> expected_visited{0, 2, 4, 6, 14, 16, 10, 18};

  EXPECT_EQ(visited, expected_visited);
}

TEST_F(ProbabilisticSuffixTreeTest, PSTBreadthFirstIterationSubtree) {
  std::vector<size_t> visited{};

  probabilisticSuffixTree.pst_breadth_first_iteration(
      2, 1, [&](size_t index, size_t level) {
        visited.push_back(index);
        return true;
      });

  std::vector<size_t> expected_visited{2, 14, 16};

  EXPECT_EQ(visited, expected_visited);
}

void test_suffix_links(pst::ProbabilisticSuffixTree<seqan3::dna5> tree) {
  tree.construct_tree();

  std::map<size_t, std::string> labels{};

  tree.breadth_first_iteration_sequential([&](size_t node_index, size_t lcp,
                                              size_t edge_lcp,
                                              size_t node_count) -> bool {
    labels[node_index] = tree.node_label(node_index, lcp, edge_lcp);
    return true;
  });

  tree.pst_breadth_first_iteration(
      [&](size_t node_index, size_t level) -> bool {
        if (node_index == 0) {
          return true;
        }
        auto label = labels[node_index];

        auto parent = tree.get_pst_parent(node_index);

        EXPECT_NE(parent, -1);

        auto parent_label = labels[parent];
        auto expected_parent_label = label.substr(1);

        EXPECT_EQ(expected_parent_label, parent_label);

        return true;
      });
}

TEST_F(ProbabilisticSuffixTreeTest, SuffixLinksCorrect) {
  size_t sought_n_parameters{30300};

  pst::ProbabilisticSuffixTree<seqan3::dna5> tree{
      "TEST", long_sequence, 15, 4, sought_n_parameters, "parameters", false};
  test_suffix_links(tree);
}

TEST_F(ProbabilisticSuffixTreeTest, SuffixLinksCorrectParallel) {
  size_t sought_n_parameters{30300};

  // If it succeeds 1000 times, we have no race conditions?
  for (int i = 0; i < 10; i++) {
    pst::ProbabilisticSuffixTree<seqan3::dna5> tree{
        "TEST",       long_sequence, 15, 4, sought_n_parameters,
        "parameters", true,          1};
    test_suffix_links(tree);
  }
}

struct my_traits : seqan3::sequence_file_input_default_traits_dna {
  template <typename alph>
  using sequence_container =
      std::vector<alph>; // must be defined as a template!
};
void test_benchmark(std::string filename, bool parallel, int parallel_depth) {
  seqan3::sequence_file_input<my_traits> file_in{filename};
  std::vector<lst::details::sequence_t<seqan3::dna5>> sequences{};
  std::vector<std::string> ids{};

  for (auto &[seq, id, qual] : file_in) {
    sequences.push_back(lst::details::sequence_t<seqan3::dna5>{std::move(seq)});
    ids.push_back(id);
  }

  size_t max_depth = 15;
  size_t min_count = 100;
  float threshold = 1.2;

  pst::KullbackLieblerTree<seqan3::dna5> pst{
      ids[0],    sequences[0], max_depth,     min_count,
      threshold, parallel,     parallel_depth};
  EXPECT_NO_FATAL_FAILURE(pst.construct_tree());
}

// TEST(ProbabilisticSuffixTreeLiveTest, HumanHerpesvirus5) {
//  std::string filename{"../../fasta/NC_006273.2.fa"};
//  test_benchmark(filename, true, 1);
//}
//
// TEST(ProbabilisticSuffixTreeLiveTest, SaccharomycesCerevisiae) {
//  std::string filename{"../../fasta/CM010781.1.fa"};
//  test_benchmark(filename, true, 1);
//}
//
// TEST(ProbabilisticSuffixTreeLiveTest, EColi) {
//  std::string filename{"../../fasta/CP007136.1.fa"};
//  test_benchmark(filename, true, 1);
//}

std::unordered_map<std::string, size_t>
get_label_index_map(pst::ProbabilisticSuffixTree<seqan3::dna5> &tree) {
  std::unordered_map<std::string, size_t> map{};

  static std::mutex labels_mutex{};

  tree.breadth_first_iteration([&](size_t node_index, size_t lcp,
                                   size_t edge_lcp, size_t node_count) -> bool {
    std::lock_guard labels_lock{labels_mutex};
    auto label = tree.node_label(node_index, lcp, edge_lcp);

    map[label] = node_index;

    return true;
  });

  return map;
}

void compare_trees(pst::ProbabilisticSuffixTree<seqan3::dna5> &left,
                   pst::ProbabilisticSuffixTree<seqan3::dna5> &right) {
  // Can't compare node indices since the trees can be created in different
  // orders.

  auto left_index_map = get_label_index_map(left);
  auto right_index_map = get_label_index_map(right);

  // Check that the labels are the same
  for (auto &[str, v] : left_index_map) {
    EXPECT_TRUE(right_index_map.find(str) != right_index_map.end())
        << "right does not contain " << str << std::endl;
  }

  // Check that the labels are the same
  for (auto &[str, v] : right_index_map) {
    EXPECT_TRUE(left_index_map.find(str) != left_index_map.end())
        << "left does not contain " << str << std::endl;
  }

  // Check that the counts are the same
  for (auto &[str, left_idx] : left_index_map) {
    auto right_idx = right_index_map[str];

    auto left_count = left.get_counts(left_idx);
    auto right_count = right.get_counts(right_idx);

    EXPECT_EQ(left_count, right_count) << str << " didn't match!" << std::endl;
  }

  // Check that the probabilities are the same
  for (auto &[str, left_idx] : left_index_map) {
    int right_idx = right_index_map[str];

    auto left_probabilities = left.get_probabilities(left_idx);
    auto right_probabilities = right.get_probabilities(right_idx);

    for (size_t i = 0; i < right_probabilities.size(); i++) {
      EXPECT_FLOAT_EQ(left_probabilities[i], right_probabilities[i]);
    }
  }
}

TEST_F(ProbabilisticSuffixTreeTest, ParallelAndSequentialSame) {
  std::vector<seqan3::dna5> sequence = random_sequence(100000);
  // If parallel successfully runs 1000 times, is it correct?
  for (int i = 0; i < 10; i++) {
    pst::KullbackLieblerTree<seqan3::dna5> tree{
        "TEST", sequence, 15, 10, 1.2, 0, "cutoff", false, 1};
    tree.construct_tree();

    pst::KullbackLieblerTree<seqan3::dna5> parallel_tree{
        "TEST", sequence, 15, 10, 1.2, 0, "cutoff", true, 1};
    parallel_tree.construct_tree();

    compare_trees(tree, parallel_tree);
  }
}

TEST_F(ProbabilisticSuffixTreeTest, ResourceTemporarilyUnavailableError) {
  std::vector<seqan3::dna5> sequence = random_sequence(1000000);
  pst::KullbackLieblerTree tree{"test", sequence, 15,    10, 1.2,
                                0,      "cutoff", false, 0};
  tree.support_pruning();
  tree.add_suffix_links();
  tree.add_reverse_suffix_links();
  tree.similarity_pruning();
}

std::unordered_map<std::string, size_t>
get_label_count_map(pst::ProbabilisticSuffixTree<seqan3::dna5> &tree) {
  std::unordered_map<std::string, size_t> map{};

  static std::mutex labels_mutex{};

  tree.breadth_first_iteration([&](size_t node_index, size_t lcp,
                                   size_t edge_lcp, size_t node_count) -> bool {
    std::lock_guard labels_lock{labels_mutex};
    auto label = tree.node_label(node_index, lcp, edge_lcp);

    map[label] = tree.get_counts(node_index);

    return true;
  });

  return map;
}

void correct_counts(pst::ProbabilisticSuffixTree<seqan3::dna5> tree,
                    lst::details::sequence_t<seqan3::dna5> seq) {
  auto tree_counts = get_label_count_map(tree);
  robin_hood::unordered_map<std::string, int> counts{};

  std::string sequence =
      seq | seqan3::views::to_char | seqan3::views::to<std::string>;

  for (size_t i = 0; i < sequence.size(); i++) {
    for (int j = 0; j < 16 && j + i < sequence.size(); j++) {
      auto kmer = sequence.substr(i, j);
      if (counts.find(kmer) != counts.end()) {
        counts[kmer] += 1;
      } else {
        counts[kmer] = 1;
      }
    }
  }

  for (auto &[k, v] : tree_counts) {
    EXPECT_EQ(counts[k], v) << k;
  }
}

TEST_F(ProbabilisticSuffixTreeTest, CorrectCounts) {
  pst::ProbabilisticSuffixTree<seqan3::dna5> tree{
      "TEST", long_sequence, 15, 100, 24601, "parameters", false, 1};
  tree.construct_tree();

  pst::ProbabilisticSuffixTree<seqan3::dna5> parallel_tree{
      "TEST", long_sequence, 15, 100, 24601, "parameters", true, 1};
  parallel_tree.construct_tree();

  correct_counts(tree, long_sequence);
  correct_counts(parallel_tree, long_sequence);
}
