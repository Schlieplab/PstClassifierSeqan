#include <gtest/gtest.h>

#include <array>
#include <filesystem>
#include <string>
#include <tuple>

#include <vector>

#include <seqan3/alphabet/nucleotide/dna5.hpp>

#include "../src/distances/cv.hpp"
#include "../src/kl_tree_map.hpp"

class DistancesTest : public ::testing::Test {
protected:
  void SetUp() override {
    first = pst::KullbackLieblerTreeMap<seqan3::dna5>{first_path};
    second = pst::KullbackLieblerTreeMap<seqan3::dna5>{second_path};
    third = pst::KullbackLieblerTreeMap<seqan3::dna5>{third_path};
  }

  pst::KullbackLieblerTreeMap<seqan3::dna5> first;
  pst::KullbackLieblerTreeMap<seqan3::dna5> second;
  pst::KullbackLieblerTreeMap<seqan3::dna5> third;
  std::filesystem::path first_path{"./../trees/CM008035.1.tree"};
  std::filesystem::path second_path{"./../trees/CM008036.1.tree"};
  std::filesystem::path third_path{"./../trees/NC_009067.tree"};
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
  float expected_dissimilar = 0.46000117222458986;
  EXPECT_FLOAT_EQ(expected_dissimilar, pst::distances::cv(first, third));

  float expected_similar = 0.10424618759761345;
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
  EXPECT_EQ(pst::distances::details::get_background_state(c1, 2), "CG");
  EXPECT_EQ(pst::distances::details::get_background_state(c2, 3), "ACG");
  EXPECT_EQ(pst::distances::details::get_background_state(c2, 1), "G");
  EXPECT_EQ(pst::distances::details::get_background_state(c3, 2), "CG");
  EXPECT_EQ(pst::distances::details::get_background_state(c4, 2), "CG");
  EXPECT_EQ(pst::distances::details::get_background_state(c5, 2), "G");
}