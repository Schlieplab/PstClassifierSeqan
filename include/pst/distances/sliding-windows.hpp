#pragma once

#include <algorithm>
#include <future>
#include <map>
#include <string>

#include "seqan3/alphabet/concept.hpp"
#include "seqan3/alphabet/views/all.hpp"
#include "seqan3/std/ranges"

#include "../probabilistic_suffix_tree_map.hpp"
#include "negative_log_likelihood.hpp"
#include "parallelize.hpp"

namespace pst::distances::details {

template <seqan3::alphabet alphabet_t>
std::vector<std::vector<double>>
sliding_window_part(ProbabilisticSuffixTreeMap<alphabet_t> &tree,
                    std::string &sequence, const std::vector<int> &window_sizes,
                    const scoring::score_signature<alphabet_t> &score_fun,
                    size_t start, size_t end) {

  int max_window_size =
      *std::max_element(window_sizes.begin(), window_sizes.end());

  std::deque<double> scores{};

  std::vector<std::vector<double>> log_likelihoods(window_sizes.size());

  size_t order_max = tree.get_max_order();

  size_t range_start = std::max(start + order_max - max_window_size, (size_t)0);

  for (size_t i = range_start; i < end; i++) {
    size_t start_index = i - order_max;
    if (order_max > i) {
      start_index = 0;
    }
    std::string subsequence =
        sequence.substr(start_index, std::min(i, order_max));

    char char_ = sequence[i];
    auto [context, val] = tree.get_closest_state(subsequence);

    double score = score_fun(tree, context, val, char_);
    scores.push_back(score);
    if (scores.size() > max_window_size) {
      scores.pop_front();
    }

    if (i < start) {
      continue;
    }

    for (int j = 0; j < window_sizes.size(); j++) {
      long window_size = window_sizes[j];

      auto start_iterator = scores.end() - window_size;
      if (scores.size() < window_size) {
        start_iterator = scores.begin();
      }

      double log_likelihood =
          std::accumulate(start_iterator, scores.end(), 0.0);
      double length = std::count_if(start_iterator, scores.end(),
                                 [](double i) { return i != 0.0; });
      if (length == 0.0) {
        log_likelihoods[j].push_back(-1.0);

      } else {
        double nll = -log_likelihood / length;
        log_likelihoods[j].push_back(nll);
      }
    }
  }

  return log_likelihoods;
}

template <seqan3::alphabet alphabet_t>
std::vector<std::vector<double>>
sliding_windows_(ProbabilisticSuffixTreeMap<alphabet_t> &tree,
                 std::string &sequence, const std::vector<int> &window_sizes,
                 const scoring::score_signature<alphabet_t> &score_fun) {

  std::vector<std::vector<double>> log_likelihoods(window_sizes.size());

  auto bounds = pst::parallelize::get_bounds(sequence.size());
  std::vector<std::future<std::vector<std::vector<double>>>> part_futures{};

  for (auto &[start_index, stop_index] : bounds) {
    part_futures.push_back(
        std::async(std::launch::async, sliding_window_part<alphabet_t>,
                   std::ref(tree), std::ref(sequence), std::ref(window_sizes),
                   std::ref(score_fun), start_index, stop_index));
  }

  double log_likelihood = 0.0;
  for (auto &part_future : part_futures) {
    auto part_likelihoods = part_future.get();
    for (int j = 0; j < window_sizes.size(); j++) {
      log_likelihoods[j].insert(log_likelihoods[j].end(),
                                part_likelihoods[j].begin(),
                                part_likelihoods[j].end());
    }
  }

  return log_likelihoods;
}
} // namespace pst::distances::details

namespace pst::distances {

template <seqan3::alphabet alphabet_t>
std::vector<std::vector<double>>
sliding_windows(ProbabilisticSuffixTreeMap<alphabet_t> &tree,
                std::string &sequence, std::vector<int> &window_sizes) {
  return details::sliding_windows_<alphabet_t>(
      tree, sequence, window_sizes,
      details::scoring::log_transition_prob<alphabet_t>);
}

template <seqan3::alphabet alphabet_t>
std::vector<std::vector<double>> sliding_windows_background(
    ProbabilisticSuffixTreeMap<alphabet_t> &tree, std::string &sequence,
    const std::vector<int> &window_sizes, int background_order) {

  return details::sliding_windows_<alphabet_t>(
      tree, sequence, window_sizes,
      details::scoring::specialise_background_log_transition_prob<alphabet_t>(
          background_order));
}

template <seqan3::alphabet alphabet_t>
std::vector<double>
sliding_windows(ProbabilisticSuffixTreeMap<alphabet_t> &tree,
                std::string &sequence, int window_size) {
  std::vector<int> window_sizes{window_size};
  return sliding_windows(tree, sequence, window_sizes)[0];
}

template <seqan3::alphabet alphabet_t>
std::vector<double>
sliding_windows_background(ProbabilisticSuffixTreeMap<alphabet_t> &tree,
                           std::string &sequence, int window_size,
                           int background_order) {
  std::vector<int> window_sizes{window_size};
  return sliding_windows_background(tree, sequence, window_sizes,
                                    background_order)[0];
}
} // namespace pst::distances

namespace pst {
static std::vector<std::vector<double>>
sliding_windows_cpp(std::string tree_string, std::string sequence,
                    std::vector<int> window_sizes) {
  tree_t tree{tree_string};

  auto scores = pst::distances::sliding_windows(tree, sequence, window_sizes);

  return scores;
}

static std::vector<std::vector<double>>
sliding_windows_background_cpp(std::string tree_string, std::string sequence,
                               const std::vector<int> &window_sizes,
                               int background_order) {
  tree_t tree{tree_string};

  auto scores = pst::distances::sliding_windows_background(
      tree, sequence, window_sizes, background_order);

  return scores;
}
} // namespace pst
