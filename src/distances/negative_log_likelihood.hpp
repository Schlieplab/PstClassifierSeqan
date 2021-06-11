#pragma once

#include <algorithm>
#include <future>
#include <map>
#include <string>

#include <seqan3/alphabet/concept.hpp>
#include <seqan3/alphabet/views/all.hpp>
#include <seqan3/std/ranges>

#include "../probabilistic_suffix_tree.hpp"
#include "../probabilistic_suffix_tree_map.hpp"
#include "parallelize.hpp"

namespace pst::distances::details {

static std::map<char, char> convert{
    {'A', 'T'}, {'T', 'A'}, {'G', 'C'}, {'C', 'G'}};

std::string get_reverse_complement(const std::string &forward) {
  std::string reverse_string(forward);
  std::transform(reverse_string.begin(), reverse_string.end(),
                 reverse_string.begin(),
                 [](char c) -> char { return convert[c]; });

  std::reverse(reverse_string.begin(), reverse_string.end());
  return reverse_string;
}

template <seqan3::alphabet alphabet_t>
double log_transition_prob(ProbabilisticSuffixTreeMap<alphabet_t> &tree,
                           const std::string &context, char char_) {
  double probability = tree.get_transition_probability(context, char_);
  if (probability == 0.0) {
    return 0.0;
  } else {
    return std::log(probability);
  }
}

template <seqan3::alphabet alphabet_t>
double
background_log_transition_prob(ProbabilisticSuffixTreeMap<alphabet_t> &tree,
                               const std::string &context, char char_,
                               int background_order) {

  std::string background_context{""};
  if (background_order < context.size()) {
    background_context = tree.get_closest_state(
        context.substr(context.size() - background_order));
  }

  double probability = tree.get_transition_probability(context, char_);
  double background_probability =
      tree.get_transition_probability(background_context, char_);

  if (probability == 0.0 || background_probability == 0.0 ||
      (probability - background_probability) == 0.0) {
    return 0.0;
  } else {
    return std::log(std::abs(probability - background_probability)) -
           std::log(background_probability);
  }
}

template <seqan3::alphabet alphabet_t>
double log_likelihood_part(ProbabilisticSuffixTreeMap<alphabet_t> &tree,
                           std::vector<alphabet_t> &sequence_dna, size_t start,
                           size_t end) {
  double log_likelihood = 0.0;

  size_t order_max = tree.get_max_order();

  std::string subsequence{};
  subsequence.resize(order_max);

  for (size_t i = start + order_max; i < end; i++) {
    size_t start_index = i - order_max;
    size_t end_index = i;
    if (order_max > i) {
      start_index = 0;
      end_index = i;
    }

    size_t context_len = end_index - start_index;
    size_t len_diff = order_max - context_len;
    for (int j = 0; j < len_diff; j++) {
      subsequence[j] = ' ';
    }
    for (size_t j = start_index; j < end_index; j++) {
      subsequence[j - start_index + len_diff] = sequence_dna[j].to_char();
    }

    size_t char_ = sequence_dna[i].to_rank();
    auto context = tree.get_closest_state(subsequence);

    double probability = tree.get_transition_probability(context, char_);

    if (probability != 0.0) {
      log_likelihood += std::log(probability);
    }
  }

  return log_likelihood;
}

template <seqan3::alphabet alphabet_t>
double log_likelihood_part_dna(ProbabilisticSuffixTreeMap<alphabet_t> &tree,
                               std::vector<alphabet_t> &sequence_dna,
                               size_t start, size_t end) {
  return log_likelihood_part(tree, sequence_dna, start, end);
}

template <seqan3::alphabet alphabet_t>
double log_likelihood_part(ProbabilisticSuffixTreeMap<alphabet_t> &tree,
                           std::string &sequence, size_t start, size_t end) {

  double log_likelihood = 0.0;

  size_t order_max = tree.get_max_order();

  for (size_t i = start + order_max; i < end; i++) {
    size_t start_index = i - order_max;
    if (order_max > i) {
      start_index = 0;
    }
    std::string subsequence =
        sequence.substr(start_index, std::min(i, order_max));

    char char_ = sequence[i];
    auto context = tree.get_closest_state(subsequence);

    double probability = tree.get_transition_probability(context, char_);

    if (probability != 0.0) {
      log_likelihood += std::log(probability);
    }
  }

  return log_likelihood;
}

template <seqan3::alphabet alphabet_t>
double likelihood_context(
    ProbabilisticSuffixTreeMap<alphabet_t> &tree, std::string &sequence,
    const std::function<double(ProbabilisticSuffixTreeMap<alphabet_t> &,
                               const std::string &, char)> &score_fun) {

  double log_likelihood = 0.0;

  size_t order_max = tree.get_max_order();

  for (size_t i = 0; i < sequence.size(); i++) {
    size_t start_index = i - order_max;
    if (order_max > i) {
      start_index = 0;
    }
    std::string subsequence =
        sequence.substr(start_index, std::min(i, order_max));

    char char_ = sequence[i];
    auto context = tree.get_closest_state(subsequence);

    double score = score_fun(tree, context, char_);
    log_likelihood += score;
  }

  return std::exp(log_likelihood);
}

template <seqan3::alphabet alphabet_t>
std::vector<std::vector<double>> sliding_window_part(
    ProbabilisticSuffixTreeMap<alphabet_t> &tree, std::string &sequence,
    std::vector<int> &window_sizes,
    const std::function<double(ProbabilisticSuffixTreeMap<alphabet_t> &,
                               const std::string &, char)> &score_fun,
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
    auto context = tree.get_closest_state(subsequence);

    double score = score_fun(tree, context, char_);
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
      double nll =
          -log_likelihood / double(std::distance(start_iterator, scores.end()));

      log_likelihoods[j].push_back(nll);
    }
  }

  return log_likelihoods;
}

template <seqan3::alphabet alphabet_t>
std::vector<std::vector<double>> sliding_windows_(
    ProbabilisticSuffixTreeMap<alphabet_t> &tree, std::string &sequence,
    std::vector<int> &window_sizes,
    const std::function<double(ProbabilisticSuffixTreeMap<alphabet_t> &,
                               const std::string &, char)> &score_fun) {

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
double log_likelihood(ProbabilisticSuffixTreeMap<alphabet_t> &tree,
                      std::string &sequence) {
  auto bounds = pst::parallelize::get_bounds(sequence.size());
  std::vector<std::future<double>> part_futures{};

  for (auto &[start_index, stop_index] : bounds) {
    part_futures.push_back(std::async(
        std::launch::async, details::log_likelihood_part<alphabet_t>,
        std::ref(tree), std::ref(sequence), start_index, stop_index));
  }

  double log_likelihood = 0.0;
  for (auto &f : part_futures) {
    log_likelihood += f.get();
  }
  return log_likelihood;
}

template <seqan3::alphabet alphabet_t>
double log_likelihood(ProbabilisticSuffixTreeMap<alphabet_t> &tree,
                      std::vector<alphabet_t> &sequence_dna) {
  auto bounds = pst::parallelize::get_bounds(sequence_dna.size());
  std::vector<std::future<double>> part_futures{};

  for (auto &[start_index, stop_index] : bounds) {
    part_futures.push_back(std::async(
        std::launch::async, log_likelihood_part_dna<alphabet_t>,
        std::ref(tree), std::ref(sequence_dna), start_index, stop_index));
  }

  double log_likelihood = 0.0;
  for (auto &f : part_futures) {
    log_likelihood += f.get();
  }
  return log_likelihood;
}

template <seqan3::alphabet alphabet_t>
double negative_log_likelihood_p(ProbabilisticSuffixTreeMap<alphabet_t> &tree,
                                 std::vector<alphabet_t> &sequence_dna) {
  return -log_likelihood(tree, sequence_dna) / sequence_dna.size();
}

template <seqan3::alphabet alphabet_t>
double negative_log_likelihood(ProbabilisticSuffixTreeMap<alphabet_t> &tree,
                               std::vector<alphabet_t> &sequence_dna) {
  double log_likelihood = 0.0;

  size_t order_max = tree.get_max_order();

  std::string subsequence{};
  subsequence.resize(order_max);

  size_t length = 0;
  for (size_t i = 0; i < sequence_dna.size(); i++) {
    size_t start_index = i - order_max;
    size_t end_index = i;
    if (order_max > i) {
      start_index = 0;
      end_index = i;
    }

    size_t context_len = end_index - start_index;
    size_t len_diff = order_max - context_len;
    for (int j = 0; j < len_diff; j++) {
      subsequence[j] = ' ';
    }
    for (size_t j = start_index; j < end_index; j++) {
      subsequence[j - start_index + len_diff] = sequence_dna[j].to_char();
    }

    size_t char_ = sequence_dna[i].to_rank();

    auto context = tree.get_closest_state(subsequence);

    auto probability = tree.get_transition_probability(context, char_);
    if (probability != 0.0) {
      length += 1;
      log_likelihood += std::log(probability);
    }
  }

  return -log_likelihood / double(length);
}

template <seqan3::alphabet alphabet_t>
double negative_log_likelihood(ProbabilisticSuffixTree<alphabet_t> &tree,
                               std::vector<alphabet_t> &sequence_dna) {
  // Example:
  // Tree edges:  ACGT -> ACGTT
  // Reverse edges:  ACGT -> CGT
  // A step forward if the current context is ACGT and the next char is T:
  // 1. ACGT -> find prob of T.
  // 2. If ACGTT in tree, follow forward edge to ACGTT
  // 3. Else, go to reverse suffix link CGT check for CGTT.
  // 4. Repeat 3 until match.
  // 5. Try to extend to the left to find GCGTT

  double log_likelihood = 0.0;

  size_t current_node = 0;
  int node_length = 0;

  size_t length = 0;
  for (size_t i = 0; i < sequence_dna.size(); i++) {
    size_t char_ = sequence_dna[i].to_rank();

    auto probability = tree.get_transition_probability(current_node, char_);

    if (probability != 0.0) {
      length += 1;
      log_likelihood += std::log(probability);
    }

    size_t prev_index = current_node;
    size_t child_index = tree.go_forward(current_node, char_);
    node_length++;
    while (child_index == prev_index) {
      node_length--;
      prev_index = tree.get_pst_parent(prev_index);
      child_index = tree.go_forward(prev_index, char_);
    }

    if (i > node_length) {
      do {
        prev_index = child_index;
        child_index = tree.go_backward(child_index,
                                       sequence_dna[i - node_length].to_rank());
        if (child_index != prev_index) {
          node_length++;
        }
      } while (child_index != prev_index && i > node_length);
    }

    current_node = child_index;
  }

  return -log_likelihood / double(length);
}

template <seqan3::alphabet alphabet_t>
double
negative_log_likelihood_symmetric(ProbabilisticSuffixTreeMap<alphabet_t> &tree,
                                  std::vector<alphabet_t> &sequence) {

  std::vector<alphabet_t> reverse_sequence =
      sequence | std::views::reverse | seqan3::views::complement |
      seqan3::views::to<std::vector<alphabet_t>>;

  return std::min(negative_log_likelihood<alphabet_t>(tree, sequence),
                  negative_log_likelihood<alphabet_t>(tree, reverse_sequence));
}

template <seqan3::alphabet alphabet_t>
double
negative_log_likelihood_symmetric_p(ProbabilisticSuffixTreeMap<alphabet_t> &tree,
                                  std::vector<alphabet_t> &sequence) {

  std::vector<alphabet_t> reverse_sequence =
      sequence | std::views::reverse | seqan3::views::complement |
      seqan3::views::to<std::vector<alphabet_t>>;

  return std::min(
      negative_log_likelihood_p<alphabet_t>(tree, sequence),
      negative_log_likelihood_p<alphabet_t>(tree, reverse_sequence));
}

std::vector<std::vector<double>>
sliding_windows(ProbabilisticSuffixTreeMap<alphabet_t> &tree,
                std::string &sequence, std::vector<int> &window_sizes) {
  return details::sliding_windows_<alphabet_t>(
      tree, sequence, window_sizes,
      [](ProbabilisticSuffixTreeMap<alphabet_t> &tree,
         const std::string &context, char char_) -> double {
        return details::log_transition_prob<alphabet_t>(tree, context, char_);
      });
}

template <seqan3::alphabet alphabet_t>
std::vector<std::vector<double>> sliding_windows_background(
    ProbabilisticSuffixTreeMap<alphabet_t> &tree, std::string &sequence,
    std::vector<int> &window_sizes, int background_order) {
  return details::sliding_windows_<alphabet_t>(
      tree, sequence, window_sizes,
      [&](ProbabilisticSuffixTreeMap<alphabet_t> &tree,
          const std::string &context, char char_) -> double {
        return details::background_log_transition_prob<alphabet_t>(
            tree, context, char_, background_order);
      });
}

template <seqan3::alphabet alphabet_t>
std::vector<std::vector<double>>
sliding_windows(ProbabilisticSuffixTreeMap<alphabet_t> &tree,
                std::string &sequence, int window_sizes) {
  return sliding_windows(tree, sequence, {window_sizes});
}

template <seqan3::alphabet alphabet_t>
std::vector<std::vector<double>>
sliding_windows_background(ProbabilisticSuffixTreeMap<alphabet_t> &tree,
                           std::string &sequence, int window_sizes,
                           int background_order) {
  return sliding_windows_background(tree, sequence, {window_sizes},
                                    background_order);
}

} // namespace pst::distances
