#include <algorithm>
#include <filesystem>
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

void score_slice(
    int start_index, int stop_index, std::vector<double> &scores,
    std::vector<tree_t> &trees, std::vector<seqan3::dna5> &sequence,
    const std::function<float(tree_t &, std::vector<seqan3::dna5> &)> &fun) {

  for (int i = start_index; i < stop_index; i++) {
    scores[i] = fun(trees[i], sequence);
  }
}

std::vector<std::vector<double>>
score_sequences(std::vector<tree_t> &trees, std::vector<std::string> &sequences,
                int background_order) {
  std::vector<std::vector<double>> scores{};

  int sequence_idx = 0;

  for (auto &sequence : sequences) {
    std::vector<seqan3::dna5> seq = sequence |
                                    seqan3::views::char_to<seqan3::dna5> |
                                    seqan3::views::to<std::vector>;

    std::vector<double> scores_row(trees.size());

    auto fun = [&](int start_index, int stop_index) {
      score_slice(
          start_index, stop_index, std::ref(scores_row), std::ref(trees),
          std::ref(seq),
          pst::distances::negative_log_likelihood_symmetric<seqan3::dna5>);
    };

    pst::parallelize::parallelize(trees.size(), fun);

    scores.push_back(scores_row);
  }

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

} // namespace pst
