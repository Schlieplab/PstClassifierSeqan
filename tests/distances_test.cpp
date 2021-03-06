#include <gtest/gtest.h>

#include <array>
#include <string>
#include <tuple>

#include <vector>

#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/std/filesystem>

#include "../src/distances/cv.hpp"
#include "../src/distances/score.hpp"
#include "../src/kl_tree_map.hpp"
#include "../src/probabilistic_suffix_tree_map.hpp"

#include "../src/kl_tree.hpp"
#include "random_sequence.hpp"

class DistancesTest : public ::testing::Test {
protected:
  void SetUp() override {
    first = pst::KullbackLieblerTreeMap<seqan3::dna5>{first_path};
    second = pst::KullbackLieblerTreeMap<seqan3::dna5>{second_path};
    third = pst::KullbackLieblerTreeMap<seqan3::dna5>{third_path};
    using seqan3::operator""_dna5;
    sequence = lst::details::sequence_t<seqan3::dna5>{
        "AAAAATTTTTTAAAAAATTTTTTAAAAAATTTTT"_dna5};
    seq = std::string{"AAAATAAAAATAAAAATAAAAAAATAAAAAAATAAAACAAAAGAAAACAAAAAAAA"
                      "AAAAAAATAAAAAAAC"
                      "AAAATAAAAATAAAAATAAAAAAATAAAAAAATAAAACAAAAGAAAACAAAAAAAA"
                      "AAAAAAATAAAAAAAC"
                      "AAAATAAAAATAAAAATAAAAAAATAAAAAAATAAAACAAAAGAAAACAAAAAAAA"
                      "AAAAAAATAAAAAAAC"
                      "AAAATAAAAATAAAAATAAAAAAATAAAAAAATAAAACAAAAGAAAACAAAAAAAA"
                      "AAAAAAATAAAAAAAC"
                      "AAAATAAAAATAAAAATAAAAAAATAAAAAAATAAAACAAAAGAAAACAAAAAAAA"
                      "AAAAAAATAAAAAAAC"
                      "AAAATAAAAATAAAAATAAAAAAATAAAAAAATAAAACAAAAGAAAACAAAAAAAA"
                      "AAAAAAATAAAAAAAC"
                      "AAAATAAAAATAAAAATAAAAAAATAAAAAAATAAAACAAAAGAAAACAAAAAAAA"
                      "AAAAAAATAAAAAAA"
                      "C"};

    binary_seq =
        std::string{"AAAATAAAAATAAAAATAAAAAAATAAAAAAATAAAAAAAAAAAAAAAAAAAA"
                    "AAAAAAATAAAAAAA"
                    "AAAATAAAAATAAAAATAAAAAAATAAAAAAATAAAAAAAAAAAAAAAAAAAA"
                    "AAAAAAATAAAAAAA"
                    "AAAATAAAAATAAAAATAAAAAAATAAAAAAATAAAAAAAAAAAAAAAAAAAA"
                    "AAAAAAATAAAAAAA"
                    "AAAATAAAAATAAAAATAAAAAAATAAAAAAATAAAAAAAAAAAAAAAAAAAA"
                    "AAAAAAATAAAAAAA"
                    "AAAATAAAAATAAAAATAAAAAAATAAAAAAATAAAAAAAAAAAAAAAAAAAA"
                    "AAAAAAATAAAAAAA"
                    "AAAATAAAAATAAAAATAAAAAAATAAAAAAATAAAAAAAAAAAAAAAAAAAA"
                    "AAAAAAATAAAAAAA"
                    "AAAATAAAAATAAAAATAAAAAAATAAAAAAATAAAAAAAAAAAAAAAAAAAA"
                    "AAAAAAATAAAAAAA"};
  }

  pst::KullbackLieblerTreeMap<seqan3::dna5> first;
  pst::KullbackLieblerTreeMap<seqan3::dna5> second;
  pst::KullbackLieblerTreeMap<seqan3::dna5> third;
  std::filesystem::path first_path{"./../trees/CM008035.1.tree"};
  std::filesystem::path second_path{"./../trees/CM008036.1.tree"};
  std::filesystem::path third_path{"./../trees/NC_009067.tree"};

  size_t max_depth = 1;
  size_t min_count = 1;
  float threshold = 3.9075;

  lst::details::sequence_t<seqan3::dna5> sequence;
  std::string seq;
  std::string binary_seq;
};

TEST_F(DistancesTest, CVSymmetry) {
  auto result_forward = pst::distances::cv(first, second);
  auto result_backward = pst::distances::cv(second, first);
  EXPECT_FLOAT_EQ(result_forward, result_backward);
}

TEST_F(DistancesTest, CVIdentity) {
  auto result_zero = pst::distances::cv(first, first);
  EXPECT_EQ(result_zero, 0);
}

TEST_F(DistancesTest, CVSnapshots) {
  float expected_dissimilar = 0.49773079;
  EXPECT_FLOAT_EQ(expected_dissimilar, pst::distances::cv(first, third));

  float expected_similar = 0.063189894;
  EXPECT_FLOAT_EQ(expected_similar, pst::distances::cv(first, second));
}

TEST_F(DistancesTest, AllContexts) {
  auto contexts = pst::distances::details::get_all_contexts<seqan3::dna5>(
      2, first.valid_characters);

  std::vector<std::string> expected_contexts{"AA", "AC", "AG", "AT", "CA", "CC",
                                             "CG", "CT", "GA", "GC", "GG", "GT",
                                             "TA", "TC", "TG", "TT"};

  EXPECT_EQ(contexts, expected_contexts);
}

TEST_F(DistancesTest, CVEstimationSymmetry) {
  auto result_forward = pst::distances::cv_estimation(first, second);
  auto result_backward = pst::distances::cv_estimation(second, first);
  EXPECT_FLOAT_EQ(result_forward, result_backward);
}

TEST_F(DistancesTest, CVEstimationIdentity) {
  auto result_zero = pst::distances::cv_estimation(first, first);
  EXPECT_EQ(result_zero, 0);
}

TEST_F(DistancesTest, BackgroudState) {
  std::string c1{"AACG"};
  std::string c2{"AACG"};
  std::string c3{"ACG"};
  std::string c4{"CG"};
  std::string c5{"G"};
  EXPECT_EQ(pst::distances::get_background_state(c1, 2), "CG");
  EXPECT_EQ(pst::distances::get_background_state(c2, 3), "ACG");
  EXPECT_EQ(pst::distances::get_background_state(c2, 1), "G");
  EXPECT_EQ(pst::distances::get_background_state(c3, 2), "CG");
  EXPECT_EQ(pst::distances::get_background_state(c4, 2), "CG");
  EXPECT_EQ(pst::distances::get_background_state(c5, 2), "G");
}

TEST_F(DistancesTest, ScoringHangs) {
  using seqan3::operator""_dna5;
  std::filesystem::path first_path{"./../trees/CM008035.1.tree"};
  pst::ProbabilisticSuffixTreeMap<seqan3::dna5> first{first_path};

  std::string seq{"ACGATCGATCGATCGATCGACACTACCAGCACATAGTAGCTAGCATGATCGACTACTAGC"
                  "ATCTACGGCTACGATCATCGATCGATCATATCAGCACTAGCACG"};
  std::vector<std::string> sequences(2000, seq);

  std::vector<pst::ProbabilisticSuffixTreeMap<seqan3::dna5>> trees(50, first);

  EXPECT_NO_FATAL_FAILURE({ pst::score_sequences(trees, sequences, 0); });
}

TEST_F(DistancesTest, LogLikelighoodGrowth) {
  using seqan3::operator""_dna5;
  lst::details::sequence_t<seqan3::dna5> sequence{
      "AAAATAAAAATAAAAATAAAAAAATAAAAAAATAAAACAAAAGAAAACAAAAAAAAAAAAAAATAAAAAAACAAAATAAAAATAAAAATAAAAAAATAAAAAAATAAAACAAAAGAAAACAAAAAAAAAAAAAAATAAAAAAACAAAATAAAAATAAAAATAAAAAAATAAAAAAATAAAACAAAAGAAAACAAAAAAAAAAAAAAATAAAAAAACAAAATAAAAATAAAAATAAAAAAATAAAAAAATAAAACAAAAGAAAACAAAAAAAAAAAAAAATAAAAAAACAAAATAAAAATAAAAATAAAAAAATAAAAAAATAAAACAAAAGAAAACAAAAAAAAAAAAAAATAAAAAAACAAAATAAAAATAAAAATAAAAAAATAAAAAAATAAAACAAAAGAAAACAAAAAAAAAAAAAAATAAAAAAACAAAATAAAAATAAAAATAAAAAAATAAAAAAATAAAACAAAAGAAAACAAAAAAAAAAAAAAATAAAAAAAC"_dna5};

  size_t max_depth_4 = 4;
  size_t max_depth_7 = 7;

  pst::KullbackLieblerTreeMap<seqan3::dna5> tree_simple{
      "TEST", sequence, max_depth, min_count, threshold, false, 1};
  tree_simple.construct_tree();

  pst::KullbackLieblerTreeMap<seqan3::dna5> tree_medium{
      "TEST", sequence, max_depth_4, min_count, threshold, false, 1};
  tree_medium.construct_tree();

  pst::KullbackLieblerTreeMap<seqan3::dna5> tree_complex{
      "TEST", sequence, max_depth_7, min_count, threshold, false, 1};
  tree_complex.construct_tree();

  auto log_likelihood_simple = pst::distances::log_likelihood(tree_simple, seq);
  auto log_likelihood_medium = pst::distances::log_likelihood(tree_medium, seq);
  auto log_likelihood_complex =
      pst::distances::log_likelihood(tree_complex, seq);

  EXPECT_GT(log_likelihood_medium, log_likelihood_simple);
  EXPECT_GT(log_likelihood_complex, log_likelihood_medium);
}

TEST_F(DistancesTest, LogLikelighoodHandCrafted0Order) {
  pst::KullbackLieblerTreeMap<seqan3::dna5> tree{
      "0-order", sequence, max_depth, min_count, threshold, false, 1};
  tree.counts[""] = {8, {0.8, 0.0, 0.0, 0.0, 0.2}, true};

  double log_likelihood = pst::distances::log_likelihood_part(
      tree, binary_seq, 0, binary_seq.size());
  double log_likelihood_manual = std::log(0.8) * 434 + std::log(0.2) * 42;

  EXPECT_FLOAT_EQ(log_likelihood, log_likelihood_manual);
}
TEST_F(DistancesTest, LogLikelighoodHandCrafted1Order) {
  pst::KullbackLieblerTreeMap<seqan3::dna5> tree{
      "1-order", sequence, max_depth, min_count, threshold, false, 1};
  tree.counts[""] = {8, {0.8, 0.0, 0.0, 0.0, 0.2}, true};
  tree.counts["A"] = {3, {0.9, 0.0, 0.0, 0.0, 0.1}, true};
  tree.counts["T"] = {5, {1.0, 0.0, 0.0, 0.0, 0.0}, true};

  double log_likelihood = pst::distances::log_likelihood(tree, binary_seq);

  double log_likelihood_manual =
      391 * std::log(0.9) + 42 * std::log(0.1) + 42 * std::log(1);

  EXPECT_FLOAT_EQ(log_likelihood, log_likelihood_manual);
}

std::unordered_map<std::string, size_t> count_2_mers(std::string &sequence) {
  std::unordered_map<std::string, size_t> counts{};

  for (size_t i = 0; i < sequence.size(); i++) {
    auto kmer = sequence.substr(i, 2);
    if (counts.find(kmer) != counts.end()) {
      counts[kmer] += 1;
    } else {
      counts[kmer] = 1;
    }
  }
  return counts;
}

TEST_F(DistancesTest, LogLikelighoodHandCrafted1OrderLong) {
  pst::KullbackLieblerTreeMap<seqan3::dna5> tree{
      "1-order", sequence, max_depth, min_count, threshold, false, 1};
  tree.counts[""] = {8, {0.6, 0.1, 0.1, 0.0, 0.1}, true};
  tree.counts["A"] = {3, {0.7, 0.1, 0.1, 0.0, 0.1}, true};
  tree.counts["T"] = {5, {0.4, 0.2, 0.1, 0.0, 0.3}, true};

  auto seq = random_sequence(100000000);
  std::string sequence =
      seq | seqan3::views::to_char | seqan3::views::to<std::string>;

  auto counted_2_mers = count_2_mers(sequence);

  auto log_likelihood = pst::distances::log_likelihood(tree, sequence);

  double log_likelihood_manual =
      counted_2_mers["AA"] * std::log(std::get<1>(tree.counts["A"])[0]) +
      counted_2_mers["AC"] * std::log(std::get<1>(tree.counts["A"])[1]) +
      counted_2_mers["AG"] * std::log(std::get<1>(tree.counts["A"])[2]) +
      counted_2_mers["AT"] * std::log(std::get<1>(tree.counts["A"])[4]) +
      counted_2_mers["TA"] * std::log(std::get<1>(tree.counts["T"])[0]) +
      counted_2_mers["TC"] * std::log(std::get<1>(tree.counts["T"])[1]) +
      counted_2_mers["TG"] * std::log(std::get<1>(tree.counts["T"])[2]) +
      counted_2_mers["TT"] * std::log(std::get<1>(tree.counts["T"])[4]) +
      (counted_2_mers["CA"] + counted_2_mers["GA"]) *
          std::log(std::get<1>(tree.counts[""])[0]) +
      (counted_2_mers["CC"] + counted_2_mers["GC"]) *
          std::log(std::get<1>(tree.counts[""])[1]) +
      (counted_2_mers["CG"] + counted_2_mers["GG"]) *
          std::log(std::get<1>(tree.counts[""])[2]) +
      (counted_2_mers["CT"] + counted_2_mers["GT"]) *
          std::log(std::get<1>(tree.counts[""])[4]);

  EXPECT_FLOAT_EQ(log_likelihood, log_likelihood_manual);
}

TEST_F(DistancesTest, SlidingWindows) {
  std::vector<int> window_sizes{10, 30};
  std::vector<std::vector<double>> log_likelihoods =
      pst::distances::sliding_windows(first, seq, window_sizes);

  for (auto &likelihoods : log_likelihoods) {
    for (auto val : likelihoods) {
    }
  }
}

TEST_F(DistancesTest, SlidingWindowsBackground) {
  std::vector<int> window_sizes{10, 30};
  std::vector<std::vector<double>> log_likelihoods =
      pst::distances::sliding_windows_background(first, seq, window_sizes, 0);

  for (auto &likelihoods : log_likelihoods) {
    for (auto val : likelihoods) {
      std::cout << val << " ";
    }
    std::cout << std::endl;
  }
}

TEST_F(DistancesTest, ReverseComplementWindows) {
  std::cout << pst::distances::details::get_reverse_complement("CCGGAATT");
}

TEST_F(DistancesTest, LogLikelihoodSame) {
  using seqan3::operator""_dna5;
  std::vector<seqan3::dna5> long_sequence{
      "CGCGCGGCGCGTCGTCGAGCTCGCCGCCCAGGGCGCGCCCCTCGGCGCGATCCGCGCCGCCCTCAACGACGTGACCCCCGCCGCCGACAAGGGCGAAGCCTTCGTCGGACGCCCCGACTGGCTGGGCGAGCTGTGGACCGCGCGCCGCACGGACCGGCCCCTGATCGACTCGATCACCAAAAAGGCCCTGCCCCGTGCGACCAAGGTCAAGGGCTGGCGCTGGAAGAAGCGCCCCGAGGTCGCCGACTACACGGGCAACAAAACCGAAATCTCCTCCAACGAGATCGAGACCGAGCCCGTCGAGGCCGCCGTCAAACGCATCGCCGCAGGGTGGGACACCGACCGCATTTTCGTCGACCTCGGCGACGGCGACATGATCGAGAGCCTGTGGGAGGGCGCCCGCGAGGACTACGCGATCAAGACCGAGGCCGCCGTCACCACCGGCCTCAAGACCGCAGCGACGAAGCTCACCGGCGCGCCCACCGAACTGGACAAGGCCCTCGTGGTCCTCGGCTCCAAGGCCGCCGCCATCGGCTCCCGCCTGAACTTCGTCGCCTTCGGCGCCGACGTGTGGAGCAAGTTCACCGCGCTGACCCGCGACCAGGTGCCGTGGTGGATCACCAACGGCGACCGCCTCAACCTCTCGACCGCGACCGGCGAGGTCAACGGCCTGCGCCTGTTCGTGGACCCGACCCTCGCGGCGGGCGACATCCTCGCGGGCGACACCCGGTCCGCGACGTTCTACGAGGAGCCGACCCCCATCAGGGTCAACGCCATCGACCTGCCCAAGGGCGGCGTGGACCTCGGCCTGTTCGGGTACCACGCCCTTCTCGTGAACGACCCGAACTCGTTGTTCATCATCACGGCGGGCTGACCCCATGACCCCCGACGACCTCGCCACGCGGGCCGCCGCGTGGGCGAAGCTCCCGGGCGGCGTGGACGACGCCATGAGGGCGTGCGCAGCCGCAGTGCACGCCCTCGTGGCCGCCCTGCCCGTCACGCAGGGCCGCCCCGCCTGGCGCGAGGACACGGCCCTCGGAGCGGTCATGCTCACCGCCCGCCTGCACCGCCGCCGCAACAGCCCGGCGGGCATCGAGTCCCTGACCGAGATGGGCGCGACCTACGTGAGCCGCTACGACAGCGACATCGCGCGCCTGCTGCGCATCGACGCCTTCGTCGGGCCCGTCGCCATCTGAGGGGGGCCACGAGATGAACCCGCTCTACGCGGCCGCCCAGGATGTGGCCGACATGCTCGCCGCGGCCGGAGTCCACACGGTCACCGACCCCAGGGACATCGAGCCGCCGTGCGCGTGGGTCAGCCCCAGCCGCATCGCCTACCCCACGCTCGCTGGCCGCCCCCGCACCGTCGAGTGGGAGGTGTACCTCATCGCACCCGACAGCGGCGCGCCCCTCTTCCCCCTCGGCGACCTCATCGACCGGGCCGCCACCGTCTTCCCCGGCATCGAGGCCCGCACCCTCGGCCTGACAATCCCCAACCTCAGCCCGGACCCCCTGCCCGCGATCACGTTCACCATCGAAACAGAAACGGACTAAACCCATGGCAGTGAAAACCCTCACCCTCGGCCCCGGCAAACTCAGCTTCGGCGCCCCCGAGTCCCTGACCCACGCCGCCGCCCAGGTCACCAAATGCGCCGTCAAGCCCACCGCGAAGCAGGGCGACTCCGTGGCCGTCCTGTCCGGCGACCGCGTGCCCGGCGACCGCACCGAAGCCGCGACCCTGGAGTTCACGATCTACCAGGACTTCGGCGAGGCCGAATCCTTCGTCGAATGGACCTGGGCCAACGCAGGGAAGGAACTCCCCTTCGAGTTCATCCCCGCCGACAAGCACGACAAGGCCGTGCGCGGCCGCGTCACGATCGAGCGGTCCGACATCGGCGGCGAGGTCGGCGTCAAGGTCACCGCCGACCTGGAGTTCACCTGCACCACCATGCCCACCATCGAGCCCAAAACCAAGATCGGGCACTGAGGTGGCCGACTACTCCGGCGTCAAGATCGACGGCGCGCGCCGCCTCCGCTCGACCCTGCGCAAAGCGGGCGCGGACATGCGCGACATGCGGGAGGTGAACCGCGTCGTCGCCGGCATCGTCGTCGGCGCGGCCACCGCCCGCGTCCCCCGACGCACCGGGGCCCTGGCCGCCACCGTGCGCGCAGGGGCCACCCAGGCCGCAGCCATCGGCCGCGCCGGGAACAACCGCCGCACCGGCGTCCCCTACGCCAACCCCATCCACTGGGGATGGCACCGCCACCGAATCCGCCCTAACCCGTTCCTCAGCCTCGCCGCCCAGGACACCGAACCCCAGTGGTTCGGCGTCTACGCCGACCGCATCGAACGCCTCATCAACAGCATCGAAGGAGCCTGACCCATGTCGAGCATCAAAGCCATCAACGTCGAGGTAGTCACCTCCGCCGTGACCGGCGACCTCGCCGCCGTCACCGTCCGCACCGACAACCGCGACCGCATCGCCTGGGACCTCGCGAGGGGCCGCAACAAATGGCCCCAAGCACAGGAGGCCCCCAGCCTGTGGGCCACCCACATCGCCTACACCGCGCTCCGCCGCACCGGCGAAGTCAGCTGCTCGTTCGAGGAGTTCTCCGAGGCAACTGTGAGCGCCGAACCCGAGGTCATCGACGTGGACCCTACCCGGACGGCGACCGCCGGGGCCTGATCGTCGCCCTGGCCCTCGCCACCCGCATCCCCATGAGCGAGTGGGAGACCCGCCCCGACGAGGACATCGCCACCGCACTGCAACTGCTAGAAGAGAGGAGGAGCTGACTTGGCGTCGAAAACCGCCATCCTGAGCGTCCGCGTCGTCTCCGACGTGAAGGACGCCACCAAGGGACTGGACGACGTGGCCGACAAGACCGGCCGCCTGGAGGACGGCCTCAAACGGGCCGCCGCCCCCGCCGGGATCGCCGTCGCCGCCCTCGCCGGGATCGGCAAGGCCGCCACCGACTCCGCCAGCGAGTTGCAGCAGAGCGCGGGCGCCGTCGAATCCGTATTCGGCGGGCACGCCGCCGCCGTCCAGGACGCCGCCAAGACCGCCGCCTCCAGCGTCGGCCTGGCAGCAAGCGAGTACCAGAACATGAGCGCGGTCCTGGGCGCCCAGCTCAAGAACATGGGCACCCCCATGGAGGACCTGGCCGGATCGACCCAGAACCTCATAGGCCTGGGCTCCGACCTCGCCGCCACCTTCGGGGGAACCACCGCCGACGCCGTGAGCGCCATCTCAGCCCTCCTCCGGGGCGAGCGCGACCCCATCGAGCGCTACGGCGTCTCGATCAAACAGTCGGATATCAACGCGCGTCTGGCCGCCGAGGGCATGGACAAGCTGGAAGGCGCGGCCAAGACCCAGGCCGAAGCCCAGGCCGCCCTCGCCCTGCTCACCGAGCAGACCGCATCTGCGCAAGGCCAGTTCGCGCGCGAGACCGACACGATGGCCGGGAGCCAGCAGATCGCCGCCGCCCAGTTCGAGAACGCAAAAGCCGCCCTCGGGGAGAAGCTGCTGCCCGTCGTCACGCAGTTCATGGAGGCCATGAGCGGGGCGGCTCAATGGGTCGCCCAGAACAGCGATGCGCTGCTCGTCCTCGGCGGCGTCGTCGGAACCATCGCGGGCGTGATCCTCGCCGCCAACGCCGCCATGGGCGTGTGGACCGCAGTCCAGACGACCGCCAGAGTCGCGACGGCCGCCTGGACCGGCGTCCAGGCCGCGTTCAACGCGGTCATGGCCCTGAACCCGATCACACTGGTGGTCATCGCCATCGGGGCCCTGGTCGCCGCCGTCGTCGTCGCCTACAACAAGTCCGAGGCGTTCCGCAACGCTGTATCCGCGCTGTGGGACGCCATCAAAGCGGGGGGCGGCTGGATAGTCGATCACGTCATCAAACCCATCGGAGACGCTTTCAACGCCGTTGTTGATGCCGTGAAATCCGTCTACGACTGGGTGAAGAACCTGTTCAGCGGATTCCAGCTCCCAGGCTGGCTATCGAGCGTACTCAGCTGGTTCGGGCTCGAAGCCCCCTCCGGCCCCGAATCCGGGGCCATCCTCGCCGCCACCGGGACAACGGACGCCCCCCTCGCGCGCCTCGCATCGTGGGCGCTCGCCCCCCGCACCGGCTCCAGCCCCACCCCCGCCGGCAGCGTCGTGAACATCACTGTGAACGGCGCCCTGGACCCTGACGCGGTGGCCCGCCAGATCGGCCGCATCCTGTCGCGCCGCGACCTCATCAACGGCACCGAGCAGATCGTGGGGGCGACCCTATGAGCGTCAGCGCAAGCCTCCGCGTCGCCGCAGGCGGCCTCGGCGGTGTCATCAACGCCGCCGCCGACAAGTACCCGACGACGGTGACCGTCCTCGACGACCTCACCGTCACGTGGGGCCGCGATAGCGTGGTCTCACACCCGGACCCGTCGTCAATGACCGCGACAATCGCCCTCGTGGACACTGTGCCCGACTGGCTGCGCGTCGGCGCGCTGGCGACCGTCAACGCGGTCGCCCGGACCGAGGAGTCCCAGCGGTCCTATATGCGGCTCCTGCCCTGGCGCGCCATCGAACCGGGCACCGGGTGGCGCCAGCAGGTCACCCCCGATCCGCCCGGCGCGTGGGTCGGGAGCCTCCCGGTGTTCGCCGCCGCCGGAGCAGACTCCGGGATCGGATGGTTCATAGCGCCAGGGGTGCAGCCCCCCTCCGATACCCCCGAGACCACCCAGTGGGCCGCGAACGCAAAGACGACAGCGGGCAAGCCAGTCACGTTCACAATCACCGTCCCCGAGCTCACCGGCGCGACCGTGCGCGCCTTCCCGCTCACTTACCGGCGCCCCGGCGGGCTCTACACCCGCGCCCCCGGCATGGCGATCGAGCTCTCCCCGGAGAAGTACACGCCGGGCACTGTCGAGTACTCGGGCACCTGGACCCCCGAGGCGACGGGCCTTTACGTCGGCGCCTATCTGCATATCCAGCTGCACAAGGCCCCGGCCTGGACGACGATCCCCCGCGAGCGGACCTGGCGGGCCGCACCCGGCACGTGGGCGGACGCCGGCGGCCGCGCAACCGTCACTGACGTTCACATCGCCGGGACCTCCGGCCACGTCGCTGAGCACGCCGTCGAGGTCTTCACCGGCCGGGTCCAGTCCCTCCGCGTCGAATGGTCCGAGCGCCTGTCCCGACCCATCGCGCGGATCACCGCAGTGGACAAGCTCGCCGATCTGAATGGTACCTACATCGGCGACACGCCGTGGGGCGAGGAATCCTGGAAGCTGCGCGCAGAACGCATCCTGAAACAGGCCCTCGGACCAGCAGACACACTGGAGGGCGAACCCGGCAATTGGCTGGGGACGATCCGCCCCCGAGACGTGGACCACCGCAGCGCCGGTGAGCTCATGAGGAACACACTCGCCTCGTGCGCCGCCACCGCATTCCCCGTCAGCTGGCGCAAATGGCGGGTAATCCCATTCATCTACAAGGGGAGCGATCAATCAATCACGATCCCCGGGCGCGCCATCCGCCGCGACGGCGTACAGGTCAGCACTGACGAATCCGCGAACATCTCGACAATTCAAGCCACGTATTTCGATGTGACCTACGACGGGAAAACCGGGAGGGTGAAAGACGTTATAGAGCGCACGACTACGCGGAAGAATACACCGGCAAATGAGGGACCCCCCAGGTCCATCAAAATGAAAACCGAACTATCCCGCAGTAACGAAGCGAGCGAACTCACCCGGATCATCGGGAAATACGTGAACGTGAACCAGTGGATCATCAGCGCCCTATCAGTGAAGCACGACCGGATCAGCGAGGACGCCCTTGTGCGCCTACTGTCCGCCACCGAGCGCATCGCCCAGCAGGTCGTCCTCACCGGCCTCCCACGATGGTTCCCAGCAGCGACGATGCGCGGCATCGTCATCGGCGGGTCCCTCACCATGCACCGCGGCCACTGGACCCCCACCCTCCGCATCGCAAACACACCCGACTAGAAAGAGTACCCATGCCATCGACAACCCCCCGGGGTCTTCCCTACGCAATCCCCACCGACGCCCAAGCCGCATTCCCCGACGCCGTGTCAAAACCCATCGCCGAATGGATCGAGGCGAACCTCCCGGTCATGCAAGCCGGGACCATCGCCTACCCCGCCCTCGGCTCCCAAGACCAGACGGGAGAGTACACAGTCACGTTTCCCAAACCCTTCCCAGTCACGCCCCGGATATTCATGCAAGCCGATAACCAGCGCCTCACAATCGCCGTATGGAATATCAGCCGCACCGGGTTCAAATGGATGGCCCGCAACAACAGCAACGGCAATTCGTCCTCTGGAGCGGCCTCGTGGTTCGCCGTTAGCGGCGCCACCGGACAGTAACGAAAGGAAACAAAGGAAATGACCACAGCGGTCGACGTGTTCACCGCCCGCCTCGCCTGGATGATGACCCAAGCCGACGGCGGCTACTCCCAGCCCAACCGCCTCGACGTGCGCCGCACGCGCGGCGTGTGGGACCCCGGCTTCCAGTTCGAGGGGGACTGCTCCTCCTGCGTCCTGGAGGCCGCCCACCAGGCGGGCCTGCCCACAGGCTCTGCGTCCTACACGGGCGACATGCGCGCGGGTCTGGAGGCCGTGGGATGGGCCGTCATCCCCTACGCCGCGACCGGCGGGGACCTTGACAACCTCGCCGACGGCGACGTGCTCCTATCCGAGGCCGCGAGCGGCGGCGTCGGCCATACCGGCGGCCTCATCCCCGGCGGCCTCGTCGCCGAGGCGTGGATCGACGGTCACGGAGACATCATGGGCTCCGCAGGCGGGGACGGGCCCGGCGACGACACCGGCGGGGAAACCCGCGCAGTGCCGTTCTATTCCCACCCCTACACAGTGCGGGGGCTCTGGACGCACGTCCTGCGCCCCCCAGCCCTCGACGCCGCAGACTCGCCCGCCGAACCCACCCCCACAACGAAAGGAATCCCCAATATGTTCGGAATCACCTACACGGCAAACGCCTTCGGCGGTATCACCGCCTACGTCCTCATCCACGAGTCCGCCGGTGCCGACGCCCTTGACCGCGTTCAGGCCCAGGTGTACAACAGCGTCCTTCCCAACGGCTTCACCGAGGTCCCTGAGCACCACGCCGAAATGCTCATCCGCGAGTCGTGGGTGCGCCACAACCGCATCGCCAACGCCGTCGCCGCGACCACTCGCGTAGACATCAACGAGGCCACCGCCCGCGTTCTCGCCGCCGTCAAGGAAGGAGCTGCCAAGTGAACGCCATCACGTCCCAGACCCCCGATGATCCCACGCCGCAGCCCATCTCCTGGCTTACACCCGCAGTGCGGCGTTACATCTACAACGTCACTATCGCCGCCCTCGGCGTCGCTCTCGTCTACGGCGTCGTCGATGGCCAGCACGCCGCCGCCTGGGAGGCCCTGGCCCTCGCTGTCGTGGGCCTCGCCCGCGCGCACGTCCCCGGAGACCCCCAATGAGCGACGCCTCGGCGGCCGTCGAGGTCATCGCCGCCATAGGCGGCCTCGGCGGCCTCGGGGCCGCGCTCTCAGCCGTCGCCTCCCTCATGGAGGCCCGCAGGGTCCGCGCCAGTATCCCCGCCGCCGCCGACCGCACCGAGGAGGCCATCGACGCCCTGCGCTCCGATGTCCGCGCCATCGACCGCCGCATCGGACACGAGCTCGGCGAAATCCGCCGAGCCGCCGACCGGGAACACGCCGACTATGACGCACGTCTCAGACGATTGGAGGGGTCATGAGTTGACACTCCCACTCAGGTGGTTATAAACTATAACCATCAGGAAGCCATAAGGGGCAAGCCTGAAACCTGAAGGGAGCGCGAAAATGCGCAAGTCAATCGAACTCACCTGGACCCCCGAAACCCGCGTTTGGGGCACAAGTGGGAACACCAGCGTGGCCGTCGGCACCGGAACCCTGGACGGCCGCCGCCTCGCCGTCTACGCCTTCCCCCAGTCCGATCACTGGTCGTTCTGGTCGCAAATCGAGCGTCCCGGCGGTGGTTCGACCTCGATAGAAATCAGGTCCACATTGCCCGCAGGCACCGTTCCATCGGTCCTCGGCCCGAACGGCGCCATCCGCGAAACGACGACAATCGAACTCTGACCCAACGCCCCGGCGCACACGCGCCGGGGCACTGTCACGGAAGGGCAAAGGATGAACGGCATCGAACTGCGCGCCCGCCGCGAAGCACTCGGGCTCTCACAAACCAAGTTCGCGAAAATGTGCGAAACGACTCAAGTGACCGTTTCCCGCTGGGAGAACGGCACCCGCGAACCAAGGAACGACATTGCAATACACCTGCTGATGGCAAATATCGAAGACGCCGCTATCGACCTCATCGAGGACCTGCTAGAGCTCGCCGAAGACGAAGAACTCCTAACAGCAACGCCCGACCTCCAGCTCACGGTGTACAACGACGAAGCCCGCTACGCCGCCGGGGAGCCCGTCTGGTCAAAACGCCTCCCCATGGAGACGCACCGCGTGTGCGCCGCCCGAGCCGCCGCCCTCCTCGGCGCCGAGGACGGCACACACGTCACACTCATCGAGGGCTGAGCGCCCTCACGCGATAGCCTCAACGACCTCGCGCACATCATCGTCGGGAACGAGGATGTAGCGAAGCGTCGTCGAGGGTGAAGCGTGCCCGAGCGCGCGCTGCACCGCGACGAGGTTCCTCGTGCGCGCGAACCCCGTCGAGGCGAAGGCGTGCCGGAGGGCGTGCATCGTCACGCCCTCCGGCAGCGCCCGGCCGACCAGCTTCCCCACCCACGCCGGGGAGAGGTGTCCATGGTCCGCCCCCGGGAAGAAGAACCCCGGGTCATGATCGAGCAGCTCATCGGCGAGAGAATGCGGGAGGGGGATCACGCGGGTTTTCCCGCCCTTTCCGTGGACCACGAGGGACCAGCCCGCCAGATCGCGCACAAGATCGCGCGTGTGCGCGCGGGCGACCTCACCCCGCCGCATCCCCAACTCCGCCGCCATGCGCACCATGAGACGCACTCGTGGATCCGTCGCCCGGCGCCCCACCGCGATGGCGCCCGGCGTCGCCGGCCTCGGCGCAGGGTCCGACTGTCTCACCGACGGCACCGGCGGCGCTACCTCGATGTATCCGACCCCCTGGGCCCACCGGTAGAACTGGTCGACGCTTTGATGCGCGCTCCGACGCGTATCCCGCGCCCAATCATGCGCCCCGGACCACTCGATCACCGTGAGCGGCCCAACCTCCCACGGGCCCGCCCGCAGATCGCGGGCGAACCTACTCACCCACTCGATCCGCAGTCGGATAGTCTCGGCCCGCCGGCCGGCCGCCGCCAGCGCTGTAGTCCACTCGCCTATAGGACCGGCCCATCCGGCGGGTACCGGTCGCGGTTTCATGCTCATGATGATTACCCTGCATCTATTCGGCCCCAGTGTCGCGTCATCACGCCGCCGGGACCCACGATGTAGGATCGAGCCGTGGGCCACCAATCCGTAGGTTGGGGGTTCGAGTCCCCCTGGGCCTACTCGCCCGCCCCAGCCCCGCCGGCTGGGGCGGTTGCCGCACTGCGCCGGTACCGACGCGCTGCCGACGGATGGAGGTCGCCATGTCCGACGACGCGGACAGGGCCCACGGCGCCCTGGCGGGACTCGCGCTGGGGGACGCCCTTGGCATGCCGACCCAGGCGATGACCGCCGATCAGATCAGGCTGACCTACGGGTGGGTGGACGCCCTGGTGCCAGCCGACGCCTCGCAGCCCTACGCGCCCGGCATGCCCGCCGGCAGCGTCACGGACGACACGGAGCAGGCGCTGCTCGTCGCCGGCCTGCTGGTATCGGGCGGGGGCGGCATCGACCCCCACGCCTTCTCCCGCGCCCTGCTGGACTGGGAGGACTCGATGGCGGCCCGCGGTTCCCTCGACCTGCTGGGCCCATCGACGAAGGCCGCCCTGGAACGCGTGCGGGCCGGGGAGGACCCCCTCCGCGTGGGCGGCGCGGGCACCACCAACGGCTCAGCGATGCGAGTCGCGCCCGTCGGGATCGCCTCCTCCACCCGGGATCCGCGTTTCGCCGACACCGTGTGGGAGTCGTGCCGCGTCACCCACGCCACCGAACAGGGCTTCCACGCCGCCGCGCTCGTGGCGGCAGCGGTCTCCCTCGGCATCGACGGAGCAGGGGCGGACAGCCCTTCGGACTCCGCCCGCGCCTCCTTGGAACGCGCCCTGGCCCTCGTGGAGGCACTCGGGCGCCGGGGGGCGCGGACGCCCCAGCCGGACGTGTGCGAGCGGACCCGCTACGCGCTGCGGTTCGCGCGCGCCCGCGACCCCGCCCCCGGTACTGCCGACGACGACCGGGCATTCGCCGGGGCACTGCGGGCACGCGTCGGCGCCTCCGTGGAGGCCGCCCAGTCCGTCCCCGCAGCATTCGCCATCGCCTGGCGCTACGCCGCCGATCCGTGGCGGGGCCTGTGCGTCGCCGCCAACCTCGGCGGTGACACCGACACGATCGGCGCTATCACCGGCGCCGTGCTCGGCGCCGCCCTGGGGGCCCGGTGCTGGCCCGCCCAGGAGCTGGAACGAGTGGAGGCCGTCTCCGGGCTGCGGCTGCGCGAGACCGCCGACGGTTTGCTCCGCCTGCGCGCCCACGGATCCCGACTGCCCGCCCACGGGGAGCCGGTCGCAGCACCGCAGGAGGGCAGGGTCGTCCTGCTCGGGCAGGTGGTCGTCGACCTCGCACTGCTGGCGCCGCGCGTGCCCGCTCCCGGCGGCGACGTGTTCGCAGAGGACGCGGGCATGCACGCGGGCGGGGGCTTCAACGTGCTGGCCGCTGCGCGCCGGATGGGAGCGGAGGCAGTGAGCCTGTCCGGCGTCGGGGACGGCGGATTCGCCTCGATCATCACCGCTGCGTTGGAGCGCATCGGCGCCTCCTGCGAGGGACCGCGCGTCGCGGGAACGGACTCGGGGTACTGCGTGGCCATCACGGACGGCGACGGCGAGCGCACCTTCGTCTCGACCAGGGGCGCGGAGGCCCGCCTGCCGCGCGGGTCGTGGTCCGCCCACGCGGCCCGCTTGCGCAGCGGGGACGTGGTGCACGTGGACGGCTACGCGCTGGCCCATCCGGCCAACACCGCAGCGCTGCGGGAGTTCCTCTCGGCGCACCTGCCCGCAGGGCTCCGCGCGATCGTCGACGTGTCGCCCGTCGTCGGCGATGTGGACCTCGACGACCTGCTTGCCCTGCGGGCCCTGGCCCCCCTGTGGTCCATGAACGAGCGCGAGGCGGGGATCCTCGCGGGCCGCCTCGCGCGGGCGTCCGCCGCTCCCCCGCACGGAGGCGCTCCCCCGGGGGAGGCGACACCACCGGCCGGAGCGGCCCCCGGGA"_dna5};
  std::vector<seqan3::dna5> short_sequence{
      "CGCGCGGCGCGTCGTCGAGCTCGCCGCCCAGGGCGC"_dna5};

  pst::KullbackLieblerTree<seqan3::dna5> tree{"Tree", long_sequence, 7, 2,
                                              3.9075, false,         2};
  tree.construct_tree();

  pst::KullbackLieblerTreeMap<seqan3::dna5> hash_map{
      "HashMap", long_sequence, 7, 2, 3.9075, true, 2};
  hash_map.construct_tree();

  auto tree_nll = pst::distances::negative_log_likelihood<seqan3::dna5>(
      tree, short_sequence);
  auto hash_map_nll = pst::distances::negative_log_likelihood<seqan3::dna5>(
      hash_map, short_sequence);

  EXPECT_FLOAT_EQ(tree_nll, hash_map_nll);
}
