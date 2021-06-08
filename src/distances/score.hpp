#include <algorithm>
#include <functional>
#include <iostream>
#include <string>
#include <thread>

#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/argument_parser/argument_parser.hpp>
#include <seqan3/io/sequence_file/input.hpp>

#include "../probabilistic_suffix_tree_map.hpp"
#include "negative_log_likelihood.hpp"
#include "parallelize.hpp"

namespace pst {

using tree_t = pst::ProbabilisticSuffixTreeMap<seqan3::dna5>;

void score_trees_slice(
    size_t start_index, size_t stop_index,
    std::vector<std::vector<double>> &scores, std::vector<tree_t> &trees,
    std::vector<seqan3::dna5_vector> &sequences,
    const std::function<float(tree_t &, std::vector<seqan3::dna5> &)> &fun) {

  std::cout << "Sizes: " << trees.size() << " " << sequences.size()
            << std::endl;
  for (size_t j = 0; j < sequences.size(); j++) {
    for (size_t i = start_index; i < stop_index; i++) {
      std::cout << i << " " << j << std::endl;
      scores[j][i] = fun(trees[i], sequences[j]);
    }
  }
}

void score_sequences_slice(
    size_t start_index, size_t stop_index,
    std::vector<std::vector<double>> &scores, std::vector<tree_t> &trees,
    std::vector<seqan3::dna5_vector> &sequences,
    const std::function<double(tree_t &, std::vector<seqan3::dna5> &)> &fun) {

  for (size_t i = 0; i < trees.size(); i++) {
    for (size_t j = start_index; j < stop_index; j++) {
      scores[j][i] = fun(trees[i], sequences[j]);
    }
  }
}

std::vector<std::vector<double>>
score_sequences(std::vector<tree_t> &trees, std::vector<std::string> &sequences,
                size_t background_order) {
  std::vector<std::vector<double>> scores(sequences.size(),
                                          std::vector<double>(trees.size()));

  std::vector<std::vector<seqan3::dna5>> dna_sequences(sequences.size());
  for (size_t sequence_idx = 0; sequence_idx < sequences.size();
       sequence_idx++) {
    std::vector<seqan3::dna5> seq = sequences[sequence_idx] |
                                    seqan3::views::char_to<seqan3::dna5> |
                                    seqan3::views::to<std::vector>;
    dna_sequences[sequence_idx] = seq;
  }

  auto fun = [&](size_t start_index, size_t stop_index) {
    score_trees_slice(
        start_index, stop_index, std::ref(scores), std::ref(trees),
        std::ref(dna_sequences),
        pst::distances::negative_log_likelihood_symmetric<seqan3::dna5>);
  };

  pst::parallelize::parallelize(trees.size(), fun);

  return scores;
}

static std::vector<std::vector<double>>
score_cpp(std::vector<std::string> tree_strings,
          std::vector<std::string> sequences) {
  std::vector<tree_t> trees{};

  std::transform(tree_strings.begin(), tree_strings.end(),
                 std::back_inserter(trees),
                 [](std::string &tree) -> tree_t { return tree_t{tree}; });

  auto scores = score_sequences(trees, sequences, 0);

  return scores;
}

static std::vector<std::vector<double>>
sliding_windows_cpp(std::string tree_string, std::string sequence,
                    std::vector<int> window_sizes) {
  tree_t tree{tree_string};

  auto scores = pst::distances::sliding_windows(tree, sequence, window_sizes);

  return scores;
}

static std::vector<std::vector<double>>
sliding_windows_background_cpp(std::string tree_string, std::string sequence,
                               std::vector<int> window_sizes,
                               int background_order) {
  tree_t tree{tree_string};

  auto scores = pst::distances::sliding_windows_background(
      tree, sequence, window_sizes, background_order);

  return scores;
}

} // namespace pst
