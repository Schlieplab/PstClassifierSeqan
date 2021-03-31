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

#include "../src/kl_tree.hpp"
#include "../src/kl_tree_map.hpp"
#include "random_sequence.hpp"
#include "test_utils.hpp"

class ProbabilisticSuffixTreeRegressionTest : public ::testing::Test {
protected:
  void SetUp() override {}
};

void same_contexts(pst::KullbackLieblerTree<seqan3::dna5> &tree,
                   pst::KullbackLieblerTreeMap<seqan3::dna5> &map) {
  auto tree_label_counts = get_label_count_map(tree);

  for (auto &[context, count] : tree_label_counts) {
    EXPECT_TRUE(map.counts.find(context) != map.counts.end())
        << "context " << context << " not in map";
    auto [map_count, map_prob, included] = map.counts[context];
    if (!included) {
      continue;
    }
    EXPECT_EQ(count, map_count)
        << "map count for " << context << " is incorrect. Should be " << count
        << " is " << map_count;
  }

  for (auto &[context, values] : map.counts) {
    if (context.size() == 0) {
      continue;
    }
    auto &[map_count, _p, included] = values;
    if (!included) {
      continue;
    }

    EXPECT_TRUE(tree_label_counts.find(context) != tree_label_counts.end())
        << "context " << context << " not in tree";
    auto count = tree_label_counts[context];
    EXPECT_EQ(count, map_count)
        << "tree count for " << context << " is incorrect. Should be " << count
        << " is " << map_count;
  }
}

TEST_F(ProbabilisticSuffixTreeRegressionTest, BothAlgorithmsSame) {
  std::vector<seqan3::dna5> sequence = random_sequence(100000);

  size_t max_depth = 15;
  size_t min_count = 10;
  float threshold = 1.2;

  pst::KullbackLieblerTree<seqan3::dna5> tree{
      "TEST", sequence, max_depth, min_count, threshold, 0, "cutoff", false, 1};
  tree.construct_tree();

  pst::KullbackLieblerTreeMap<seqan3::dna5> tree_map{
      "TEST", sequence, max_depth, min_count, threshold, 0, "cutoff", false, 1};
  tree_map.construct_tree();

  same_contexts(tree, tree_map);
}
