#include <algorithm>
#include <chrono>
#include <functional>
#include <iostream>
#include <mutex>
#include <string>
#include <thread>

#include "indicators/block_progress_bar.hpp"
#include "indicators/cursor_control.hpp"
#include "indicators/dynamic_progress.hpp"

#include "seqan3/alphabet/nucleotide/dna5.hpp"
#include "seqan3/argument_parser/argument_parser.hpp"
#include "seqan3/io/sequence_file/input.hpp"

#include "../probabilistic_suffix_tree_map.hpp"
#include "negative_log_likelihood.hpp"
#include "parallelize.hpp"

namespace pst {

using tree_t = pst::ProbabilisticSuffixTreeMap<seqan3::dna5>;

void score_trees_slice_with_progress(
    size_t start_index, size_t stop_index,
    std::vector<std::vector<double>> &scores, std::vector<tree_t> &trees,
    const std::vector<std::string> &sequences,
    const std::function<double(tree_t &, const std::string &)> &fun,
    indicators::DynamicProgress<indicators::BlockProgressBar> &bars) {

  std::stringstream ss;
  ss << "Computing negative log likelihood of trees from " << start_index
     << " to " << stop_index << "...";

  indicators::BlockProgressBar bar{indicators::option::BarWidth{50},
                                   indicators::option::PrefixText{ss.str()},
                                   indicators::option::ShowElapsedTime{true},
                                   indicators::option::ShowRemainingTime{true}};

  //  auto bars_i = bars.push_back(bar);

  for (size_t j = 0; j < sequences.size(); j++) {
    for (size_t i = start_index; i < stop_index; i++) {
      scores[j][i] = fun(trees[i], sequences[j]);
    }
    //    double progress = double(j) / double(sequences.size());
    //    bars[bars_i].set_progress(progress * 100);
  }

  //  bars[bars_i].mark_as_completed();
}

void score_trees_slice(
    size_t start_index, size_t stop_index,
    std::vector<std::vector<double>> &scores, std::vector<tree_t> &trees,
    std::vector<seqan3::dna5_vector> &sequences,
    const std::function<double(tree_t &, std::vector<seqan3::dna5> &)> &fun) {

  for (size_t j = 0; j < sequences.size(); j++) {
    for (size_t i = start_index; i < stop_index; i++) {
      scores[j][i] = fun(trees[i], sequences[j]);
    }
  }
}

void score_sequences_slice_with_progress(
    size_t start_index, size_t stop_index,
    std::vector<std::vector<double>> &scores, std::vector<tree_t> &trees,
    std::vector<seqan3::dna5_vector> &sequences,
    const std::function<double(tree_t &, std::vector<seqan3::dna5> &)> &fun,
    indicators::DynamicProgress<indicators::ProgressBar> &bars) {
  std::stringstream ss;
  ss << "\"Computing negative log likelihood of sequences from " << start_index
     << " to " << stop_index << " sequences for all trees... ";

  indicators::ProgressBar bar{indicators::option::BarWidth{50},
                              indicators::option::PrefixText{ss.str()},
                              indicators::option::ShowElapsedTime{true},
                              indicators::option::ShowRemainingTime{true}};

  // auto bars_i =  bars.push_back(bar);;

  for (size_t i = 0; i < trees.size(); i++) {
    for (size_t j = start_index; j < stop_index; j++) {
      scores[j][i] = fun(trees[i], sequences[j]);
    }
    double progress = double(i) / double(trees.size());
    // bars[bars_i].set_progress(progress * 100);
  }

  // bars[bars_i].mark_as_completed();
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
score_sequences(std::vector<tree_t> &trees,
                const std::vector<std::string> &sequences,
                size_t background_order) {
  std::vector<std::vector<double>> scores(sequences.size(),
                                          std::vector<double>(trees.size()));

  auto fun = [&](size_t start_index, size_t stop_index,
                 indicators::DynamicProgress<indicators::BlockProgressBar>
                     &bars) {
    score_trees_slice_with_progress(
        start_index, stop_index, scores, trees, sequences,
        pst::distances::negative_log_likelihood_symmetric_string<seqan3::dna5>,
        bars);
  };

  pst::parallelize::parallelize_with_progress(trees.size(), fun);

  return scores;
}

static std::vector<std::vector<double>>
score_cpp(const std::vector<std::string> &tree_strings,
          const std::vector<std::string> &sequences) {
  std::vector<tree_t> trees{};

  std::transform(tree_strings.begin(), tree_strings.end(),
                 std::back_inserter(trees),
                 [](const std::string &tree) -> tree_t {
                   return tree_t{tree, 1.0};
                 });

  auto scores = score_sequences(trees, sequences, 0);

  return scores;
}

static std::vector<std::vector<double>>
score_background_cpp(const std::vector<std::string> &tree_strings,
                     const std::vector<std::string> &sequences) {
  std::vector<tree_t> trees{};

  std::transform(tree_strings.begin(), tree_strings.end(),
                 std::back_inserter(trees),
                 [](const std::string &tree) -> tree_t {
                   return tree_t{tree, 1.0};
                 });

  std::vector<std::vector<double>> scores(sequences.size(),
                                          std::vector<double>(trees.size()));

  auto fun =
      [&](size_t start_index, size_t stop_index,
          indicators::DynamicProgress<indicators::BlockProgressBar> &bars) {
        score_trees_slice_with_progress(
            start_index, stop_index, scores, trees, sequences,
            pst::distances::negative_log_likelihood_symmetric_background_string<
                seqan3::dna5>,
            bars);
      };

  pst::parallelize::parallelize_with_progress(trees.size(), fun);

  return scores;

  return scores;
}

} // namespace pst
