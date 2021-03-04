#include <gtest/gtest.h>

#include "../src/kl_tree_map.hpp"
#include "../src/probabilistic_suffix_tree_map.hpp"

#include <robin_hood.h>

#include <array>
#include <string>
#include <tuple>

#include <vector>

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/io/sequence_file/input.hpp>
#include <seqan3/range/views/to.hpp>
#include <seqan3/range/views/to_char.hpp>
#include <seqan3/std/filesystem>

class ProbabilisticSuffixTreeTestMap : public ::testing::Test {
protected:
  void SetUp() override {
    using seqan3::operator""_dna4;
    using seqan3::operator""_dna5;

    sequence = lst::details::sequence_t<seqan3::dna4>{"GATTATA"_dna4};
    long_sequence = lst::details::sequence_t<seqan3::dna5>{
        "CGCGCGGCGCGTCGTCGAGCTCGCCGCCCAGGGCGCGCCCCTCGGCGCGATCCGCGCCGCCCTCAACGACGTGACCCCCGCCGCCGACAAGGGCGAAGCCTTCGTCGGACGCCCCGACTGGCTGGGCGAGCTGTGGACCGCGCGCCGCACGGACCGGCCCCTGATCGACTCGATCACCAAAAAGGCCCTGCCCCGTGCGACCAAGGTCAAGGGCTGGCGCTGGAAGAAGCGCCCCGAGGTCGCCGACTACACGGGCAACAAAACCGAAATCTCCTCCAACGAGATCGAGACCGAGCCCGTCGAGGCCGCCGTCAAACGCATCGCCGCAGGGTGGGACACCGACCGCATTTTCGTCGACCTCGGCGACGGCGACATGATCGAGAGCCTGTGGGAGGGCGCCCGCGAGGACTACGCGATCAAGACCGAGGCCGCCGTCACCACCGGCCTCAAGACCGCAGCGACGAAGCTCACCGGCGCGCCCACCGAACTGGACAAGGCCCTCGTGGTCCTCGGCTCCAAGGCCGCCGCCATCGGCTCCCGCCTGAACTTCGTCGCCTTCGGCGCCGACGTGTGGAGCAAGTTCACCGCGCTGACCCGCGACCAGGTGCCGTGGTGGATCACCAACGGCGACCGCCTCAACCTCTCGACCGCGACCGGCGAGGTCAACGGCCTGCGCCTGTTCGTGGACCCGACCCTCGCGGCGGGCGACATCCTCGCGGGCGACACCCGGTCCGCGACGTTCTACGAGGAGCCGACCCCCATCAGGGTCAACGCCATCGACCTGCCCAAGGGCGGCGTGGACCTCGGCCTGTTCGGGTACCACGCCCTTCTCGTGAACGACCCGAACTCGTTGTTCATCATCACGGCGGGCTGACCCCATGACCCCCGACGACCTCGCCACGCGGGCCGCCGCGTGGGCGAAGCTCCCGGGCGGCGTGGACGACGCCATGAGGGCGTGCGCAGCCGCAGTGCACGCCCTCGTGGCCGCCCTGCCCGTCACGCAGGGCCGCCCCGCCTGGCGCGAGGACACGGCCCTCGGAGCGGTCATGCTCACCGCCCGCCTGCACCGCCGCCGCAACAGCCCGGCGGGCATCGAGTCCCTGACCGAGATGGGCGCGACCTACGTGAGCCGCTACGACAGCGACATCGCGCGCCTGCTGCGCATCGACGCCTTCGTCGGGCCCGTCGCCATCTGAGGGGGGCCACGAGATGAACCCGCTCTACGCGGCCGCCCAGGATGTGGCCGACATGCTCGCCGCGGCCGGAGTCCACACGGTCACCGACCCCAGGGACATCGAGCCGCCGTGCGCGTGGGTCAGCCCCAGCCGCATCGCCTACCCCACGCTCGCTGGCCGCCCCCGCACCGTCGAGTGGGAGGTGTACCTCATCGCACCCGACAGCGGCGCGCCCCTCTTCCCCCTCGGCGACCTCATCGACCGGGCCGCCACCGTCTTCCCCGGCATCGAGGCCCGCACCCTCGGCCTGACAATCCCCAACCTCAGCCCGGACCCCCTGCCCGCGATCACGTTCACCATCGAAACAGAAACGGACTAAACCCATGGCAGTGAAAACCCTCACCCTCGGCCCCGGCAAACTCAGCTTCGGCGCCCCCGAGTCCCTGACCCACGCCGCCGCCCAGGTCACCAAATGCGCCGTCAAGCCCACCGCGAAGCAGGGCGACTCCGTGGCCGTCCTGTCCGGCGACCGCGTGCCCGGCGACCGCACCGAAGCCGCGACCCTGGAGTTCACGATCTACCAGGACTTCGGCGAGGCCGAATCCTTCGTCGAATGGACCTGGGCCAACGCAGGGAAGGAACTCCCCTTCGAGTTCATCCCCGCCGACAAGCACGACAAGGCCGTGCGCGGCCGCGTCACGATCGAGCGGTCCGACATCGGCGGCGAGGTCGGCGTCAAGGTCACCGCCGACCTGGAGTTCACCTGCACCACCATGCCCACCATCGAGCCCAAAACCAAGATCGGGCACTGAGGTGGCCGACTACTCCGGCGTCAAGATCGACGGCGCGCGCCGCCTCCGCTCGACCCTGCGCAAAGCGGGCGCGGACATGCGCGACATGCGGGAGGTGAACCGCGTCGTCGCCGGCATCGTCGTCGGCGCGGCCACCGCCCGCGTCCCCCGACGCACCGGGGCCCTGGCCGCCACCGTGCGCGCAGGGGCCACCCAGGCCGCAGCCATCGGCCGCGCCGGGAACAACCGCCGCACCGGCGTCCCCTACGCCAACCCCATCCACTGGGGATGGCACCGCCACCGAATCCGCCCTAACCCGTTCCTCAGCCTCGCCGCCCAGGACACCGAACCCCAGTGGTTCGGCGTCTACGCCGACCGCATCGAACGCCTCATCAACAGCATCGAAGGAGCCTGACCCATGTCGAGCATCAAAGCCATCAACGTCGAGGTAGTCACCTCCGCCGTGACCGGCGACCTCGCCGCCGTCACCGTCCGCACCGACAACCGCGACCGCATCGCCTGGGACCTCGCGAGGGGCCGCAACAAATGGCCCCAAGCACAGGAGGCCCCCAGCCTGTGGGCCACCCACATCGCCTACACCGCGCTCCGCCGCACCGGCGAAGTCAGCTGCTCGTTCGAGGAGTTCTCCGAGGCAACTGTGAGCGCCGAACCCGAGGTCATCGACGTGGACCCTACCCGGACGGCGACCGCCGGGGCCTGATCGTCGCCCTGGCCCTCGCCACCCGCATCCCCATGAGCGAGTGGGAGACCCGCCCCGACGAGGACATCGCCACCGCACTGCAACTGCTAGAAGAGAGGAGGAGCTGACTTGGCGTCGAAAACCGCCATCCTGAGCGTCCGCGTCGTCTCCGACGTGAAGGACGCCACCAAGGGACTGGACGACGTGGCCGACAAGACCGGCCGCCTGGAGGACGGCCTCAAACGGGCCGCCGCCCCCGCCGGGATCGCCGTCGCCGCCCTCGCCGGGATCGGCAAGGCCGCCACCGACTCCGCCAGCGAGTTGCAGCAGAGCGCGGGCGCCGTCGAATCCGTATTCGGCGGGCACGCCGCCGCCGTCCAGGACGCCGCCAAGACCGCCGCCTCCAGCGTCGGCCTGGCAGCAAGCGAGTACCAGAACATGAGCGCGGTCCTGGGCGCCCAGCTCAAGAACATGGGCACCCCCATGGAGGACCTGGCCGGATCGACCCAGAACCTCATAGGCCTGGGCTCCGACCTCGCCGCCACCTTCGGGGGAACCACCGCCGACGCCGTGAGCGCCATCTCAGCCCTCCTCCGGGGCGAGCGCGACCCCATCGAGCGCTACGGCGTCTCGATCAAACAGTCGGATATCAACGCGCGTCTGGCCGCCGAGGGCATGGACAAGCTGGAAGGCGCGGCCAAGACCCAGGCCGAAGCCCAGGCCGCCCTCGCCCTGCTCACCGAGCAGACCGCATCTGCGCAAGGCCAGTTCGCGCGCGAGACCGACACGATGGCCGGGAGCCAGCAGATCGCCGCCGCCCAGTTCGAGAACGCAAAAGCCGCCCTCGGGGAGAAGCTGCTGCCCGTCGTCACGCAGTTCATGGAGGCCATGAGCGGGGCGGCTCAATGGGTCGCCCAGAACAGCGATGCGCTGCTCGTCCTCGGCGGCGTCGTCGGAACCATCGCGGGCGTGATCCTCGCCGCCAACGCCGCCATGGGCGTGTGGACCGCAGTCCAGACGACCGCCAGAGTCGCGACGGCCGCCTGGACCGGCGTCCAGGCCGCGTTCAACGCGGTCATGGCCCTGAACCCGATCACACTGGTGGTCATCGCCATCGGGGCCCTGGTCGCCGCCGTCGTCGTCGCCTACAACAAGTCCGAGGCGTTCCGCAACGCTGTATCCGCGCTGTGGGACGCCATCAAAGCGGGGGGCGGCTGGATAGTCGATCACGTCATCAAACCCATCGGAGACGCTTTCAACGCCGTTGTTGATGCCGTGAAATCCGTCTACGACTGGGTGAAGAACCTGTTCAGCGGATTCCAGCTCCCAGGCTGGCTATCGAGCGTACTCAGCTGGTTCGGGCTCGAAGCCCCCTCCGGCCCCGAATCCGGGGCCATCCTCGCCGCCACCGGGACAACGGACGCCCCCCTCGCGCGCCTCGCATCGTGGGCGCTCGCCCCCCGCACCGGCTCCAGCCCCACCCCCGCCGGCAGCGTCGTGAACATCACTGTGAACGGCGCCCTGGACCCTGACGCGGTGGCCCGCCAGATCGGCCGCATCCTGTCGCGCCGCGACCTCATCAACGGCACCGAGCAGATCGTGGGGGCGACCCTATGAGCGTCAGCGCAAGCCTCCGCGTCGCCGCAGGCGGCCTCGGCGGTGTCATCAACGCCGCCGCCGACAAGTACCCGACGACGGTGACCGTCCTCGACGACCTCACCGTCACGTGGGGCCGCGATAGCGTGGTCTCACACCCGGACCCGTCGTCAATGACCGCGACAATCGCCCTCGTGGACACTGTGCCCGACTGGCTGCGCGTCGGCGCGCTGGCGACCGTCAACGCGGTCGCCCGGACCGAGGAGTCCCAGCGGTCCTATATGCGGCTCCTGCCCTGGCGCGCCATCGAACCGGGCACCGGGTGGCGCCAGCAGGTCACCCCCGATCCGCCCGGCGCGTGGGTCGGGAGCCTCCCGGTGTTCGCCGCCGCCGGAGCAGACTCCGGGATCGGATGGTTCATAGCGCCAGGGGTGCAGCCCCCCTCCGATACCCCCGAGACCACCCAGTGGGCCGCGAACGCAAAGACGACAGCGGGCAAGCCAGTCACGTTCACAATCACCGTCCCCGAGCTCACCGGCGCGACCGTGCGCGCCTTCCCGCTCACTTACCGGCGCCCCGGCGGGCTCTACACCCGCGCCCCCGGCATGGCGATCGAGCTCTCCCCGGAGAAGTACACGCCGGGCACTGTCGAGTACTCGGGCACCTGGACCCCCGAGGCGACGGGCCTTTACGTCGGCGCCTATCTGCATATCCAGCTGCACAAGGCCCCGGCCTGGACGACGATCCCCCGCGAGCGGACCTGGCGGGCCGCACCCGGCACGTGGGCGGACGCCGGCGGCCGCGCAACCGTCACTGACGTTCACATCGCCGGGACCTCCGGCCACGTCGCTGAGCACGCCGTCGAGGTCTTCACCGGCCGGGTCCAGTCCCTCCGCGTCGAATGGTCCGAGCGCCTGTCCCGACCCATCGCGCGGATCACCGCAGTGGACAAGCTCGCCGATCTGAATGGTACCTACATCGGCGACACGCCGTGGGGCGAGGAATCCTGGAAGCTGCGCGCAGAACGCATCCTGAAACAGGCCCTCGGACCAGCAGACACACTGGAGGGCGAACCCGGCAATTGGCTGGGGACGATCCGCCCCCGAGACGTGGACCACCGCAGCGCCGGTGAGCTCATGAGGAACACACTCGCCTCGTGCGCCGCCACCGCATTCCCCGTCAGCTGGCGCAAATGGCGGGTAATCCCATTCATCTACAAGGGGAGCGATCAATCAATCACGATCCCCGGGCGCGCCATCCGCCGCGACGGCGTACAGGTCAGCACTGACGAATCCGCGAACATCTCGACAATTCAAGCCACGTATTTCGATGTGACCTACGACGGGAAAACCGGGAGGGTGAAAGACGTTATAGAGCGCACGACTACGCGGAAGAATACACCGGCAAATGAGGGACCCCCCAGGTCCATCAAAATGAAAACCGAACTATCCCGCAGTAACGAAGCGAGCGAACTCACCCGGATCATCGGGAAATACGTGAACGTGAACCAGTGGATCATCAGCGCCCTATCAGTGAAGCACGACCGGATCAGCGAGGACGCCCTTGTGCGCCTACTGTCCGCCACCGAGCGCATCGCCCAGCAGGTCGTCCTCACCGGCCTCCCACGATGGTTCCCAGCAGCGACGATGCGCGGCATCGTCATCGGCGGGTCCCTCACCATGCACCGCGGCCACTGGACCCCCACCCTCCGCATCGCAAACACACCCGACTAGAAAGAGTACCCATGCCATCGACAACCCCCCGGGGTCTTCCCTACGCAATCCCCACCGACGCCCAAGCCGCATTCCCCGACGCCGTGTCAAAACCCATCGCCGAATGGATCGAGGCGAACCTCCCGGTCATGCAAGCCGGGACCATCGCCTACCCCGCCCTCGGCTCCCAAGACCAGACGGGAGAGTACACAGTCACGTTTCCCAAACCCTTCCCAGTCACGCCCCGGATATTCATGCAAGCCGATAACCAGCGCCTCACAATCGCCGTATGGAATATCAGCCGCACCGGGTTCAAATGGATGGCCCGCAACAACAGCAACGGCAATTCGTCCTCTGGAGCGGCCTCGTGGTTCGCCGTTAGCGGCGCCACCGGACAGTAACGAAAGGAAACAAAGGAAATGACCACAGCGGTCGACGTGTTCACCGCCCGCCTCGCCTGGATGATGACCCAAGCCGACGGCGGCTACTCCCAGCCCAACCGCCTCGACGTGCGCCGCACGCGCGGCGTGTGGGACCCCGGCTTCCAGTTCGAGGGGGACTGCTCCTCCTGCGTCCTGGAGGCCGCCCACCAGGCGGGCCTGCCCACAGGCTCTGCGTCCTACACGGGCGACATGCGCGCGGGTCTGGAGGCCGTGGGATGGGCCGTCATCCCCTACGCCGCGACCGGCGGGGACCTTGACAACCTCGCCGACGGCGACGTGCTCCTATCCGAGGCCGCGAGCGGCGGCGTCGGCCATACCGGCGGCCTCATCCCCGGCGGCCTCGTCGCCGAGGCGTGGATCGACGGTCACGGAGACATCATGGGCTCCGCAGGCGGGGACGGGCCCGGCGACGACACCGGCGGGGAAACCCGCGCAGTGCCGTTCTATTCCCACCCCTACACAGTGCGGGGGCTCTGGACGCACGTCCTGCGCCCCCCAGCCCTCGACGCCGCAGACTCGCCCGCCGAACCCACCCCCACAACGAAAGGAATCCCCAATATGTTCGGAATCACCTACACGGCAAACGCCTTCGGCGGTATCACCGCCTACGTCCTCATCCACGAGTCCGCCGGTGCCGACGCCCTTGACCGCGTTCAGGCCCAGGTGTACAACAGCGTCCTTCCCAACGGCTTCACCGAGGTCCCTGAGCACCACGCCGAAATGCTCATCCGCGAGTCGTGGGTGCGCCACAACCGCATCGCCAACGCCGTCGCCGCGACCACTCGCGTAGACATCAACGAGGCCACCGCCCGCGTTCTCGCCGCCGTCAAGGAAGGAGCTGCCAAGTGAACGCCATCACGTCCCAGACCCCCGATGATCCCACGCCGCAGCCCATCTCCTGGCTTACACCCGCAGTGCGGCGTTACATCTACAACGTCACTATCGCCGCCCTCGGCGTCGCTCTCGTCTACGGCGTCGTCGATGGCCAGCACGCCGCCGCCTGGGAGGCCCTGGCCCTCGCTGTCGTGGGCCTCGCCCGCGCGCACGTCCCCGGAGACCCCCAATGAGCGACGCCTCGGCGGCCGTCGAGGTCATCGCCGCCATAGGCGGCCTCGGCGGCCTCGGGGCCGCGCTCTCAGCCGTCGCCTCCCTCATGGAGGCCCGCAGGGTCCGCGCCAGTATCCCCGCCGCCGCCGACCGCACCGAGGAGGCCATCGACGCCCTGCGCTCCGATGTCCGCGCCATCGACCGCCGCATCGGACACGAGCTCGGCGAAATCCGCCGAGCCGCCGACCGGGAACACGCCGACTATGACGCACGTCTCAGACGATTGGAGGGGTCATGAGTTGACACTCCCACTCAGGTGGTTATAAACTATAACCATCAGGAAGCCATAAGGGGCAAGCCTGAAACCTGAAGGGAGCGCGAAAATGCGCAAGTCAATCGAACTCACCTGGACCCCCGAAACCCGCGTTTGGGGCACAAGTGGGAACACCAGCGTGGCCGTCGGCACCGGAACCCTGGACGGCCGCCGCCTCGCCGTCTACGCCTTCCCCCAGTCCGATCACTGGTCGTTCTGGTCGCAAATCGAGCGTCCCGGCGGTGGTTCGACCTCGATAGAAATCAGGTCCACATTGCCCGCAGGCACCGTTCCATCGGTCCTCGGCCCGAACGGCGCCATCCGCGAAACGACGACAATCGAACTCTGACCCAACGCCCCGGCGCACACGCGCCGGGGCACTGTCACGGAAGGGCAAAGGATGAACGGCATCGAACTGCGCGCCCGCCGCGAAGCACTCGGGCTCTCACAAACCAAGTTCGCGAAAATGTGCGAAACGACTCAAGTGACCGTTTCCCGCTGGGAGAACGGCACCCGCGAACCAAGGAACGACATTGCAATACACCTGCTGATGGCAAATATCGAAGACGCCGCTATCGACCTCATCGAGGACCTGCTAGAGCTCGCCGAAGACGAAGAACTCCTAACAGCAACGCCCGACCTCCAGCTCACGGTGTACAACGACGAAGCCCGCTACGCCGCCGGGGAGCCCGTCTGGTCAAAACGCCTCCCCATGGAGACGCACCGCGTGTGCGCCGCCCGAGCCGCCGCCCTCCTCGGCGCCGAGGACGGCACACACGTCACACTCATCGAGGGCTGAGCGCCCTCACGCGATAGCCTCAACGACCTCGCGCACATCATCGTCGGGAACGAGGATGTAGCGAAGCGTCGTCGAGGGTGAAGCGTGCCCGAGCGCGCGCTGCACCGCGACGAGGTTCCTCGTGCGCGCGAACCCCGTCGAGGCGAAGGCGTGCCGGAGGGCGTGCATCGTCACGCCCTCCGGCAGCGCCCGGCCGACCAGCTTCCCCACCCACGCCGGGGAGAGGTGTCCATGGTCCGCCCCCGGGAAGAAGAACCCCGGGTCATGATCGAGCAGCTCATCGGCGAGAGAATGCGGGAGGGGGATCACGCGGGTTTTCCCGCCCTTTCCGTGGACCACGAGGGACCAGCCCGCCAGATCGCGCACAAGATCGCGCGTGTGCGCGCGGGCGACCTCACCCCGCCGCATCCCCAACTCCGCCGCCATGCGCACCATGAGACGCACTCGTGGATCCGTCGCCCGGCGCCCCACCGCGATGGCGCCCGGCGTCGCCGGCCTCGGCGCAGGGTCCGACTGTCTCACCGACGGCACCGGCGGCGCTACCTCGATGTATCCGACCCCCTGGGCCCACCGGTAGAACTGGTCGACGCTTTGATGCGCGCTCCGACGCGTATCCCGCGCCCAATCATGCGCCCCGGACCACTCGATCACCGTGAGCGGCCCAACCTCCCACGGGCCCGCCCGCAGATCGCGGGCGAACCTACTCACCCACTCGATCCGCAGTCGGATAGTCTCGGCCCGCCGGCCGGCCGCCGCCAGCGCTGTAGTCCACTCGCCTATAGGACCGGCCCATCCGGCGGGTACCGGTCGCGGTTTCATGCTCATGATGATTACCCTGCATCTATTCGGCCCCAGTGTCGCGTCATCACGCCGCCGGGACCCACGATGTAGGATCGAGCCGTGGGCCACCAATCCGTAGGTTGGGGGTTCGAGTCCCCCTGGGCCTACTCGCCCGCCCCAGCCCCGCCGGCTGGGGCGGTTGCCGCACTGCGCCGGTACCGACGCGCTGCCGACGGATGGAGGTCGCCATGTCCGACGACGCGGACAGGGCCCACGGCGCCCTGGCGGGACTCGCGCTGGGGGACGCCCTTGGCATGCCGACCCAGGCGATGACCGCCGATCAGATCAGGCTGACCTACGGGTGGGTGGACGCCCTGGTGCCAGCCGACGCCTCGCAGCCCTACGCGCCCGGCATGCCCGCCGGCAGCGTCACGGACGACACGGAGCAGGCGCTGCTCGTCGCCGGCCTGCTGGTATCGGGCGGGGGCGGCATCGACCCCCACGCCTTCTCCCGCGCCCTGCTGGACTGGGAGGACTCGATGGCGGCCCGCGGTTCCCTCGACCTGCTGGGCCCATCGACGAAGGCCGCCCTGGAACGCGTGCGGGCCGGGGAGGACCCCCTCCGCGTGGGCGGCGCGGGCACCACCAACGGCTCAGCGATGCGAGTCGCGCCCGTCGGGATCGCCTCCTCCACCCGGGATCCGCGTTTCGCCGACACCGTGTGGGAGTCGTGCCGCGTCACCCACGCCACCGAACAGGGCTTCCACGCCGCCGCGCTCGTGGCGGCAGCGGTCTCCCTCGGCATCGACGGAGCAGGGGCGGACAGCCCTTCGGACTCCGCCCGCGCCTCCTTGGAACGCGCCCTGGCCCTCGTGGAGGCACTCGGGCGCCGGGGGGCGCGGACGCCCCAGCCGGACGTGTGCGAGCGGACCCGCTACGCGCTGCGGTTCGCGCGCGCCCGCGACCCCGCCCCCGGTACTGCCGACGACGACCGGGCATTCGCCGGGGCACTGCGGGCACGCGTCGGCGCCTCCGTGGAGGCCGCCCAGTCCGTCCCCGCAGCATTCGCCATCGCCTGGCGCTACGCCGCCGATCCGTGGCGGGGCCTGTGCGTCGCCGCCAACCTCGGCGGTGACACCGACACGATCGGCGCTATCACCGGCGCCGTGCTCGGCGCCGCCCTGGGGGCCCGGTGCTGGCCCGCCCAGGAGCTGGAACGAGTGGAGGCCGTCTCCGGGCTGCGGCTGCGCGAGACCGCCGACGGTTTGCTCCGCCTGCGCGCCCACGGATCCCGACTGCCCGCCCACGGGGAGCCGGTCGCAGCACCGCAGGAGGGCAGGGTCGTCCTGCTCGGGCAGGTGGTCGTCGACCTCGCACTGCTGGCGCCGCGCGTGCCCGCTCCCGGCGGCGACGTGTTCGCAGAGGACGCGGGCATGCACGCGGGCGGGGGCTTCAACGTGCTGGCCGCTGCGCGCCGGATGGGAGCGGAGGCAGTGAGCCTGTCCGGCGTCGGGGACGGCGGATTCGCCTCGATCATCACCGCTGCGTTGGAGCGCATCGGCGCCTCCTGCGAGGGACCGCGCGTCGCGGGAACGGACTCGGGGTACTGCGTGGCCATCACGGACGGCGACGGCGAGCGCACCTTCGTCTCGACCAGGGGCGCGGAGGCCCGCCTGCCGCGCGGGTCGTGGTCCGCCCACGCGGCCCGCTTGCGCAGCGGGGACGTGGTGCACGTGGACGGCTACGCGCTGGCCCATCCGGCCAACACCGCAGCGCTGCGGGAGTTCCTCTCGGCGCACCTGCCCGCAGGGCTCCGCGCGATCGTCGACGTGTCGCCCGTCGTCGGCGATGTGGACCTCGACGACCTGCTTGCCCTGCGGGCCCTGGCCCCCCTGTGGTCCATGAACGAGCGCGAGGCGGGGATCCTCGCGGGCCGCCTCGCGCGGGCGTCCGCCGCTCCCCCGCACGGAGGCGCTCCCCCGGGGGAGGCGACACCACCGGCCGGAGCGGCCCCCGGGA"_dna5};

    probabilisticSuffixTree = pst::ProbabilisticSuffixTreeMap<seqan3::dna4>{
        "TEST", sequence, 2, 1, 192, "parameters", false, 2};
  }

  lst::details::sequence_t<seqan3::dna4> sequence;
  lst::details::sequence_t<seqan3::dna5> long_sequence;
  pst::ProbabilisticSuffixTreeMap<seqan3::dna4> probabilisticSuffixTree;
};

TEST_F(ProbabilisticSuffixTreeTestMap, ConstructorStatus) {
  robin_hood::unordered_set<std::string> expected_status{
      "", "A", "G", "T", "AT", "GA", "TA", "TT"};

  probabilisticSuffixTree.construct_tree();

  probabilisticSuffixTree.print();

  EXPECT_EQ(probabilisticSuffixTree.status, expected_status);
}

TEST_F(ProbabilisticSuffixTreeTestMap, ConstructorProbabilities) {
  robin_hood::unordered_map<
      std::string, std::array<float, seqan3::alphabet_size<seqan3::dna4>>>
      expected_probabilities{
          {"", {4.0 / 11.0, 1.0 / 11.0, 2.0 / 11.0, 4.0 / 11.0}}, // root
          {"A", {1.0 / 6.0, 1.0 / 6.0, 1.0 / 6.0, 3.0 / 6.0}},    // A
          {"G", {2.0 / 5.0, 1.0 / 5.0, 1.0 / 5.0, 1.0 / 5.0}},    // G
          {"T", {3.0 / 7.0, 1.0 / 7.0, 1.0 / 7.0, 2.0 / 7.0}},    // T
          {"AT", {2.0 / 6.0, 1.0 / 6.0, 1.0 / 6.0, 2.0 / 6.0}},   // AT
          {"GA", {1.0 / 5.0, 1.0 / 5.0, 1.0 / 5.0, 2.0 / 5.0}},   // GA
          {"TA", {1.0 / 5.0, 1.0 / 5.0, 1.0 / 5.0, 2.0 / 5.0}},   // TA
          {"TT", {2.0 / 5.0, 1.0 / 5.0, 1.0 / 5.0, 1.0 / 5.0}},   // TT
      };

  probabilisticSuffixTree.construct_tree();

  std::vector<std::string> labels{"", "A", "G", "T", "AT", "GA", "TA", "TT"};

  for (auto &label : labels) {
    for (size_t j = 0; j < 4; j++) {
      EXPECT_FLOAT_EQ(probabilisticSuffixTree.probabilities[label][j],
                      expected_probabilities[label][j])
          << label << " " << j;
    }
  }
}

TEST_F(ProbabilisticSuffixTreeTestMap, PrunedKL) {
  auto kl_tree = pst::KullbackLieblerTreeMap<seqan3::dna4>{
      "TEST", sequence, 3, 2, 0.2, 192, "cutoff", false, 2};
  kl_tree.construct_tree();
  kl_tree.print();

  robin_hood::unordered_set<std::string> expected_status{"", "A"};

  EXPECT_EQ(kl_tree.status, expected_status);
}

TEST_F(ProbabilisticSuffixTreeTestMap, PrunedParameters) {
  auto kl_tree = pst::KullbackLieblerTreeMap<seqan3::dna4>{
      "TEST", sequence, 3, 2, 0.0, 6, "parameters", false, 2};
  kl_tree.construct_tree();
  robin_hood::unordered_set<std::string> expected_status{"", "A"};

  EXPECT_EQ(kl_tree.status, expected_status);
}

TEST_F(ProbabilisticSuffixTreeTestMap, Print) {
  auto pst_unpruned = pst::ProbabilisticSuffixTreeMap<seqan3::dna4>{
      "TEST", sequence, 10, 1, 192, "parameters", false};
  pst_unpruned.construct_tree();
  pst_unpruned.print();
  seqan3::debug_stream << std::endl;
  probabilisticSuffixTree.construct_tree();

  probabilisticSuffixTree.print();
  seqan3::debug_stream << std::endl;

  auto pst_pruned = pst::KullbackLieblerTreeMap<seqan3::dna4>{
      "TEST", sequence, 3, 2, 1.2, 0, "threshold", false, 2};
  pst_pruned.construct_tree();
  pst_pruned.print();
  seqan3::debug_stream << std::endl;
}

TEST_F(ProbabilisticSuffixTreeTestMap, MemoryAllocationException) {
  // Crashes when allocating memory in expand suffix links
  EXPECT_NO_FATAL_FAILURE({
    auto pst = pst::ProbabilisticSuffixTreeMap<seqan3::dna5>(
        "TEST", long_sequence, 15, 4, 48, "parameters", false);
    pst.construct_tree();
  });
}

TEST_F(ProbabilisticSuffixTreeTestMap, CorrectNumberOfParameters) {
  size_t sought_n_parameters = 192;

  pst::ProbabilisticSuffixTreeMap<seqan3::dna5> pst{
      "TEST", long_sequence, 15, 4, sought_n_parameters, "parameters", false};
  pst.construct_tree();

  size_t n_terminal = pst.count_terminal_nodes();
  EXPECT_EQ(n_terminal * 3, sought_n_parameters);
}

TEST_F(ProbabilisticSuffixTreeTestMap, PSTBreadthFirstIteration) {
  pst::ProbabilisticSuffixTreeMap<seqan3::dna4> tree{
      "TEST", sequence, 2, 1, 192, "parameters", false, 2};
  tree.construct_tree();
  std::set<std::string> visited{};

  tree.pst_breadth_first_iteration([&](const std::string &label, size_t level) {
    visited.insert(label);
    return true;
  });

  std::set<std::string> expected_visited{"",   "A",  "G",  "T",
                                         "AT", "GA", "TA", "TT"};

  EXPECT_EQ(visited, expected_visited);
}

TEST_F(ProbabilisticSuffixTreeTestMap, PSTBreadthFirstIterationSubtree) {
  std::set<std::string> visited{};
  probabilisticSuffixTree.construct_tree();

  probabilisticSuffixTree.pst_breadth_first_iteration(
      "A", 1, [&](const std::string &label, size_t level) {
        visited.insert(label);
        return true;
      });

  std::set<std::string> expected_visited{"A", "TA", "GA"};

  EXPECT_EQ(visited, expected_visited);
}

void test_suffix_links(pst::ProbabilisticSuffixTreeMap<seqan3::dna5> tree) {
  tree.construct_tree();

  tree.pst_breadth_first_iteration(
      [&](const std::string &label, size_t level) -> bool {
        if (label.empty()) {
          return true;
        }
        auto parent_label = tree.get_pst_parent(label);

        auto expected_parent_label = label.substr(1);

        EXPECT_EQ(expected_parent_label, parent_label);

        return true;
      });
}

TEST_F(ProbabilisticSuffixTreeTestMap, CorrectTree) {
  size_t sought_n_parameters{30300};

  pst::ProbabilisticSuffixTreeMap<seqan3::dna5> tree{
      "TEST", long_sequence, 15, 4, sought_n_parameters, "parameters", false};
  test_suffix_links(tree);
}

TEST_F(ProbabilisticSuffixTreeTestMap, CorrectTreeParallel) {
  size_t sought_n_parameters{30300};

  // If parallel code succeeds 1000 times, it is correct?
  for (size_t i = 0; i < 10; i++) {
    pst::ProbabilisticSuffixTreeMap<seqan3::dna5> tree{
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
    sequences.emplace_back(std::move(seq));
    ids.push_back(id);
  }

  size_t max_depth = 15;
  size_t min_count = 100;
  float threshold = 1.2;

  pst::KullbackLieblerTreeMap<seqan3::dna5> pst{
      ids[0],    sequences[0], max_depth,     min_count,
      threshold, parallel,     parallel_depth};
  EXPECT_NO_FATAL_FAILURE(pst.construct_tree());
}

void compare_trees(pst::ProbabilisticSuffixTreeMap<seqan3::dna5> left,
                   pst::ProbabilisticSuffixTreeMap<seqan3::dna5> right) {
  // Left can contain other, non-essential counts.
  for (auto &[k, v] : right.counts) {
    EXPECT_EQ(left.counts[k], v);
  }

  for (auto &v : left.status) {
    EXPECT_TRUE(right.status.find(v) != right.status.end());
  }

  for (auto &[k, v] : right.probabilities) {
    for (size_t i = 0; i < v.size(); i++) {
      EXPECT_FLOAT_EQ(left.probabilities[k][i], v[i]) << k << " " << i;
    }
  }
}

TEST_F(ProbabilisticSuffixTreeTestMap, FromTreeToStreamToTree) {
  std::filesystem::path filename{"./test.tree"};

  pst::ProbabilisticSuffixTreeMap<seqan3::dna5> tree{
      "TEST", long_sequence, 15, 4, 192, "parameters", false};
  tree.construct_tree();
  auto tree_out = tree.to_tree();
  std::ofstream out{filename};
  out << tree_out;
  out.close();

  pst::ProbabilisticSuffixTreeMap<seqan3::dna5> tree_stream{filename};

  compare_trees(tree, tree_stream);

  std::filesystem::remove(filename);
}

TEST_F(ProbabilisticSuffixTreeTestMap, FromTreeToStringToTree) {
  pst::ProbabilisticSuffixTreeMap<seqan3::dna5> tree{
      "TEST", long_sequence, 15, 4, 192, "parameters", false};
  tree.construct_tree();
  auto tree_out = tree.to_tree();

  pst::ProbabilisticSuffixTreeMap<seqan3::dna5> tree_string{tree_out};

  compare_trees(tree, tree_string);
}

TEST_F(ProbabilisticSuffixTreeTestMap, ParallelAndSequentialSame) {
  pst::ProbabilisticSuffixTreeMap<seqan3::dna5> tree{
      "TEST", long_sequence, 15, 4, 24601, "parameters", false, 1};
  tree.construct_tree();

  pst::ProbabilisticSuffixTreeMap<seqan3::dna5> parallel_tree{
      "TEST", long_sequence, 15, 4, 24601, "parameters", true, 1};
  parallel_tree.construct_tree();

  compare_trees(tree, parallel_tree);
}

void correct_counts(pst::ProbabilisticSuffixTreeMap<seqan3::dna5> &tree,
                    lst::details::sequence_t<seqan3::dna5> &seq) {
  robin_hood::unordered_map<std::string, int> counts{};

  std::string sequence =
      seq | seqan3::views::to_char | seqan3::views::to<std::string>;

  for (size_t i = 0; i < sequence.size(); i++) {
    for (int j = 0; j < 16 && j + i <= sequence.size(); j++) {
      auto kmer = sequence.substr(i, j);
      if (counts.find(kmer) != counts.end()) {
        counts[kmer] += 1;
      } else {
        counts[kmer] = 1;
      }
    }
  }

  for (auto &[k, v] : tree.counts) {
    EXPECT_EQ(counts[k], v) << k;
  }
}

TEST_F(ProbabilisticSuffixTreeTestMap, CorrectCounts) {
  pst::ProbabilisticSuffixTreeMap<seqan3::dna5> tree{
      "TEST", long_sequence, 15, 4, 24601, "parameters", false, 1};
  tree.construct_tree();

  pst::ProbabilisticSuffixTreeMap<seqan3::dna5> parallel_tree{
      "TEST", long_sequence, 15, 4, 24601, "parameters", true, 1};
  parallel_tree.construct_tree();

  correct_counts(tree, long_sequence);
  correct_counts(parallel_tree, long_sequence);
}
