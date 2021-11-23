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

namespace pst::distances::details::scoring {
template <seqan3::alphabet alphabet_t>
using score_signature = std::function<double(
    ProbabilisticSuffixTreeMap<alphabet_t> &, const std::string &, char)>;

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
  if (background_order != 0 && background_order < context.size()) {
    background_context = tree.get_closest_state(
        context.substr(context.size() - background_order));
  }

  double probability = tree.get_transition_probability(context, char_);
  double background_probability =
      tree.get_transition_probability(background_context, char_);

  if (probability == 0.0 || background_probability == 0.0) {
    return 0.0;
  } else {
    return std::log(probability) - std::log(std::sqrt(background_probability));
  }
}
} // namespace pst::distances::details::scoring

namespace pst::distances::details {

std::string get_reverse_complement(const std::string &forward) {
  static std::map<char, char> convert{
      {'A', 'T'}, {'T', 'A'}, {'G', 'C'}, {'C', 'G'}};

  std::string reverse_string(forward);
  std::transform(reverse_string.begin(), reverse_string.end(),
                 reverse_string.begin(),
                 [](char c) -> char { return convert[c]; });

  std::reverse(reverse_string.begin(), reverse_string.end());
  return reverse_string;
}

template <seqan3::alphabet alphabet_t>
std::tuple<double, size_t>
log_likelihood_part(ProbabilisticSuffixTreeMap<alphabet_t> &tree,
                    std::vector<alphabet_t> &sequence_dna, size_t start,
                    size_t end,
                    const scoring::score_signature<alphabet_t> &score_fun,
                    const size_t max_depth) {
  double log_likelihood = 0.0;
  size_t length = 0.0;

  std::string subsequence{};
  subsequence.resize(max_depth);

  for (size_t i = start; i < end; i++) {
    size_t start_index = i - max_depth;
    size_t end_index = i;
    if (max_depth > i) {
      start_index = 0;
      end_index = i;
    }

    size_t context_len = end_index - start_index;
    size_t len_diff = max_depth - context_len;
    for (int j = 0; j < len_diff; j++) {
      subsequence[j] = ' ';
    }
    for (size_t j = start_index; j < end_index; j++) {
      subsequence[j - start_index + len_diff] = sequence_dna[j].to_char();
    }

    char char_ = sequence_dna[i].to_char();
    auto context = tree.get_closest_state(subsequence);

    double score = score_fun(tree, context, char_);

    if (score != 0.0) {
      log_likelihood += score;
      length += 1;
    }
  }

  return {log_likelihood, length};
}

template <seqan3::alphabet alphabet_t>
std::tuple<double, size_t>
log_likelihood_part(ProbabilisticSuffixTreeMap<alphabet_t> &tree,
                    std::vector<alphabet_t> &sequence_dna, size_t start,
                    size_t end,
                    const scoring::score_signature<alphabet_t> &score_fun) {
  auto max_depth = tree.get_max_order();
  return log_likelihood_part(tree, sequence_dna, start, end, score_fun,
                             max_depth);
}

template <seqan3::alphabet alphabet_t>
std::tuple<double, size_t>
log_likelihood_part(ProbabilisticSuffixTreeMap<alphabet_t> &tree,
                    std::vector<alphabet_t> &sequence_dna, size_t start,
                    size_t end) {
  return log_likelihood_part(tree, sequence_dna, start, end,
                             scoring::log_transition_prob<alphabet_t>);
}

template <seqan3::alphabet alphabet_t>
double
log_likelihood_part_dna(ProbabilisticSuffixTreeMap<alphabet_t> &tree,
                        std::vector<alphabet_t> &sequence_dna, size_t start,
                        size_t end,
                        const scoring::score_signature<alphabet_t> &score_fun) {
  return std::get<0>(
      log_likelihood_part(tree, sequence_dna, start, end, score_fun));
}

template <seqan3::alphabet alphabet_t>
std::tuple<double, size_t>
log_likelihood_part(ProbabilisticSuffixTreeMap<alphabet_t> &tree,
                    std::string &sequence, size_t start, size_t end,
                    const scoring::score_signature<alphabet_t> &score_fun,
                    const int max_depth) {

  double log_likelihood = 0.0;
  size_t length = 0;

  std::string subsequence{};
  subsequence.resize(max_depth);

  for (size_t i = start + max_depth; i < end; i++) {
    size_t start_index = i - max_depth;
    size_t end_index = i;
    if (max_depth > i) {
      start_index = 0;
      end_index = i;
    }

    size_t context_len = end_index - start_index;
    size_t len_diff = max_depth - context_len;
    for (int j = 0; j < len_diff; j++) {
      subsequence[j] = ' ';
    }
    for (size_t j = start_index; j < end_index; j++) {
      subsequence[j - start_index + len_diff] = sequence[j];
    }

    char char_ = sequence[i];
    auto context = tree.get_closest_state(subsequence);

    double score = score_fun(tree, context, char_);

    if (score != 0.0) {
      log_likelihood += score;
      length += 1;
    }
  }

  return {log_likelihood, length};
}

template <seqan3::alphabet alphabet_t>
std::tuple<double, size_t>
log_likelihood_part(ProbabilisticSuffixTreeMap<alphabet_t> &tree,
                    std::string &sequence, size_t start, size_t end,
                    const scoring::score_signature<alphabet_t> &score_fun) {
  size_t order_max = tree.get_max_order();
  return log_likelihood_part(tree, sequence, start, end, score_fun, order_max);
}

template <seqan3::alphabet alphabet_t>
std::tuple<double, size_t>
log_likelihood_part(ProbabilisticSuffixTreeMap<alphabet_t> &tree,
                    std::string &sequence, size_t start, size_t end) {
  return log_likelihood_part(tree, sequence, start, end,
                             scoring::log_transition_prob<alphabet_t>);
}

template <seqan3::alphabet alphabet_t>
double log_likelihood_part_string(ProbabilisticSuffixTreeMap<alphabet_t> &tree,
                                  std::string &sequence, size_t start,
                                  size_t end) {
  return std::get<0>(log_likelihood_part(tree, sequence, start, end));
}

template <seqan3::alphabet alphabet_t>
double
likelihood_context(ProbabilisticSuffixTreeMap<alphabet_t> &tree,
                   std::string &sequence,
                   const scoring::score_signature<alphabet_t> &score_fun) {

  auto [score, _len] =
      log_likelihood_part(tree, sequence, 0, sequence.size(), score_fun);
  return std::exp(score);
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
        std::launch::async, details::log_likelihood_part_string<alphabet_t>,
        std::ref(tree), std::ref(sequence), start_index, stop_index));
  }

  double log_likelihood = 0.0;
  for (auto &f : part_futures) {
    log_likelihood += f.get();
  }
  return log_likelihood;
}

template <seqan3::alphabet alphabet_t>
double
log_likelihood(ProbabilisticSuffixTreeMap<alphabet_t> &tree,
               std::vector<alphabet_t> &sequence_dna,
               const details::scoring::score_signature<alphabet_t> &score_fun) {
  auto bounds = pst::parallelize::get_bounds(sequence_dna.size());
  std::vector<std::future<double>> part_futures{};

  for (auto &[start_index, stop_index] : bounds) {
    part_futures.push_back(
        std::async(std::launch::async,
                   details::log_likelihood_part_dna<alphabet_t>, std::ref(tree),
                   std::ref(sequence_dna), start_index, stop_index, score_fun));
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
  return log_likelihood(tree, sequence_dna,
                        details::scoring::log_transition_prob<alphabet_t>);
}

template <seqan3::alphabet alphabet_t>
double negative_log_likelihood_p(
    ProbabilisticSuffixTreeMap<alphabet_t> &tree,
    std::vector<alphabet_t> &sequence_dna,
    const details::scoring::score_signature<alphabet_t> &score_fun) {
  return -log_likelihood(tree, sequence_dna, score_fun) / sequence_dna.size();
}

template <seqan3::alphabet alphabet_t>
double negative_log_likelihood(
    ProbabilisticSuffixTreeMap<alphabet_t> &tree,
    std::vector<alphabet_t> &sequence_dna,
    const details::scoring::score_signature<alphabet_t> &score_fun) {

  auto [score, length] = details::log_likelihood_part(
      tree, sequence_dna, 0, sequence_dna.size(), score_fun);

  return -score / double(length);
}

template <seqan3::alphabet alphabet_t>
double negative_log_likelihood(ProbabilisticSuffixTreeMap<alphabet_t> &tree,
                               std::vector<alphabet_t> &sequence_dna) {
  return negative_log_likelihood(
      tree, sequence_dna, details::scoring::log_transition_prob<alphabet_t>);
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
double negative_log_likelihood_symmetric_(
    ProbabilisticSuffixTreeMap<alphabet_t> &tree,
    std::vector<alphabet_t> &sequence,
    const details::scoring::score_signature<alphabet_t> &score_fun) {
  std::vector<alphabet_t> reverse_sequence =
      sequence | std::views::reverse | seqan3::views::complement |
      seqan3::views::to<std::vector<alphabet_t>>;

  return std::min(
      negative_log_likelihood<alphabet_t>(tree, sequence, score_fun),
      negative_log_likelihood<alphabet_t>(tree, reverse_sequence, score_fun));
}

template <seqan3::alphabet alphabet_t>
double
negative_log_likelihood_symmetric(ProbabilisticSuffixTreeMap<alphabet_t> &tree,
                                  std::vector<alphabet_t> &sequence) {

  return negative_log_likelihood_symmetric_<alphabet_t>(
      tree, sequence,
      [&](ProbabilisticSuffixTreeMap<alphabet_t> &tree_,
          const std::string &context, char char_) -> double {
        return details::scoring::log_transition_prob<alphabet_t>(tree_, context,
                                                                 char_);
      });
}

template <seqan3::alphabet alphabet_t>
double negative_log_likelihood_symmetric_p(
    ProbabilisticSuffixTreeMap<alphabet_t> &tree,
    std::vector<alphabet_t> &sequence,
    const details::scoring::score_signature<alphabet_t> &score_fun) {

  std::vector<alphabet_t> reverse_sequence =
      sequence | std::views::reverse | seqan3::views::complement |
      seqan3::views::to<std::vector<alphabet_t>>;

  return std::min(
      negative_log_likelihood_p<alphabet_t>(tree, sequence, score_fun),
      negative_log_likelihood_p<alphabet_t>(tree, reverse_sequence, score_fun));
}

template <seqan3::alphabet alphabet_t>
double negative_log_likelihood_symmetric_p(
    ProbabilisticSuffixTreeMap<alphabet_t> &tree,
    std::vector<alphabet_t> &sequence,
    std::function<double(int &, const std::string &, char)> function) {

  return negative_log_likelihood_symmetric_p(
      tree, sequence, details::scoring::log_transition_prob);
}
} // namespace pst::distances
