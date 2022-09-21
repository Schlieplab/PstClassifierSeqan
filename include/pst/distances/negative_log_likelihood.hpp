#pragma once

#include <algorithm>
#include <future>
#include <map>
#include <string>

#include "seqan3/alphabet/concept.hpp"
#include "seqan3/alphabet/views/all.hpp"
#include "seqan3/std/ranges"

#include "../probabilistic_suffix_tree.hpp"
#include "../probabilistic_suffix_tree_map.hpp"
#include "parallelize.hpp"

namespace pst::distances::details::scoring {
template <seqan3::alphabet alphabet_t>
using score_signature =
    std::function<double(ProbabilisticSuffixTreeMap<alphabet_t> &,
                         const std::string &, const hashmap_value<alphabet_t> &,
                         char)>;

template <seqan3::alphabet alphabet_t>
double log_transition_prob(ProbabilisticSuffixTreeMap<alphabet_t> &tree,
                           const std::string &context,
                           const hashmap_value<alphabet_t> &val, char char_) {
  double probability = tree.get_transition_probability(val, char_);
  if (probability == 0.0) {
    return 0.0;
  } else {
    return std::log(probability);
  }
}

template <seqan3::alphabet alphabet_t>
double
transition_prob_max_adjusted(ProbabilisticSuffixTreeMap<alphabet_t> &tree,
                             const std::string &context,
                             const hashmap_value<alphabet_t> &val, char char_) {
  double probability = tree.get_transition_probability(val, char_);

  double max_prob = 0.0;
  double min_prob = 1.0;

  for (auto &c : tree.valid_characters) {
    double p = tree.get_transition_probability(val, c);
    max_prob = std::max(p, max_prob);
    min_prob = std::min(p, min_prob);
  }
  if (probability == 0.0) {
    return 0.0;
  } else {
    return probability;
  }
}

template <seqan3::alphabet alphabet_t>
double background_log_transition_prob(
    ProbabilisticSuffixTreeMap<alphabet_t> &tree, const std::string &context,
    const hashmap_value<alphabet_t> &val, char char_, int background_order) {

  std::string_view context_in;
  if (background_order != 0 && background_order < context.size()) {
    context_in = context.substr(context.size() - background_order);
  }

  auto [background_context, background_val] =
      tree.get_closest_state(context_in);

  double probability = tree.get_transition_probability(val, char_);
  double background_probability =
      tree.get_transition_probability(background_val, char_);

  if (probability == 0.0 || background_probability == 0.0) {
    return 0.0;
  } else {
    return std::log(probability) - std::log(std::sqrt(background_probability));
  }
}

template <seqan3::alphabet alphabet_t>
score_signature<alphabet_t>
specialise_background_log_transition_prob(int background_order) {
  return [&](ProbabilisticSuffixTreeMap<alphabet_t> &tree,
             const std::string &context, const hashmap_value<alphabet_t> &val,
             char char_) -> double {
    return background_log_transition_prob<alphabet_t>(tree, context, val, char_,
                                                      background_order);
  };
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
                    const std::vector<alphabet_t> &sequence_dna, size_t start,
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
    auto [context, val] = tree.get_closest_state(subsequence);

    double score = score_fun(tree, context, val, char_);

    log_likelihood += score;
    length += 1;
  }

  return {log_likelihood, length};
}

template <seqan3::alphabet alphabet_t>
std::tuple<double, size_t>
log_likelihood_part(ProbabilisticSuffixTreeMap<alphabet_t> &tree,
                    const std::vector<alphabet_t> &sequence_dna, size_t start,
                    size_t end,
                    const scoring::score_signature<alphabet_t> &score_fun) {
  auto max_depth = tree.get_max_order();
  return log_likelihood_part(tree, sequence_dna, start, end, score_fun,
                             max_depth);
}

template <seqan3::alphabet alphabet_t>
std::tuple<double, size_t>
log_likelihood_part(ProbabilisticSuffixTreeMap<alphabet_t> &tree,
                    const std::vector<alphabet_t> &sequence_dna, size_t start,
                    size_t end) {
  return log_likelihood_part(tree, sequence_dna, start, end,
                             scoring::log_transition_prob<alphabet_t>);
}

template <seqan3::alphabet alphabet_t>
double
log_likelihood_part_dna(ProbabilisticSuffixTreeMap<alphabet_t> &tree,
                        const std::vector<alphabet_t> &sequence_dna, size_t start,
                        size_t end,
                        const scoring::score_signature<alphabet_t> &score_fun) {
  return std::get<0>(
      log_likelihood_part(tree, sequence_dna, start, end, score_fun));
}

template <seqan3::alphabet alphabet_t>
std::tuple<double, size_t>
log_likelihood_part(ProbabilisticSuffixTreeMap<alphabet_t> &tree,
                    const std::string &sequence, size_t start, size_t end,
                    const scoring::score_signature<alphabet_t> &score_fun,
                    const int max_depth) {
  double log_likelihood = 0.0;
  size_t length = 0;

  std::string_view sequence_view{sequence};

  for (size_t i = start; i < end; i++) {
    size_t start_index = i - max_depth;
    size_t end_index = i;
    if (max_depth > i) {
      start_index = 0;
      end_index = i;
    }

    auto subsequence =
        sequence_view.substr(start_index, end_index - start_index);

    char char_ = sequence[i];
    auto [context, val] = tree.get_closest_state(subsequence);

    double score = score_fun(tree, context, val, char_);

    log_likelihood += score;
    length += 1;
  }

  return {log_likelihood, length};
}

template <seqan3::alphabet alphabet_t>
std::tuple<double, size_t>
log_likelihood_part(ProbabilisticSuffixTreeMap<alphabet_t> &tree,
                    const std::string &sequence, size_t start, size_t end,
                    const scoring::score_signature<alphabet_t> &score_fun) {
  size_t order_max = tree.get_max_order();
  return log_likelihood_part(tree, sequence, start, end, score_fun, order_max);
}

template <seqan3::alphabet alphabet_t>
std::tuple<double, size_t>
log_likelihood_part(ProbabilisticSuffixTreeMap<alphabet_t> &tree,
                    const std::string &sequence, size_t start, size_t end) {
  return log_likelihood_part(tree, sequence, start, end,
                             scoring::log_transition_prob<alphabet_t>);
}

template <seqan3::alphabet alphabet_t>
std::tuple<double, size_t>
log_likelihood_part(ProbabilisticSuffixTreeMap<alphabet_t> &tree,
                    const std::string &sequence) {
  return log_likelihood_part(tree, sequence, 0, sequence.size(),
                             scoring::log_transition_prob<alphabet_t>);
}

template <seqan3::alphabet alphabet_t>
double log_likelihood_part_string(ProbabilisticSuffixTreeMap<alphabet_t> &tree,
                                  const std::string &sequence, size_t start,
                                  size_t end) {
  return std::get<0>(log_likelihood_part(tree, sequence, start, end));
}

template <seqan3::alphabet alphabet_t>
double
likelihood_context(ProbabilisticSuffixTreeMap<alphabet_t> &tree,
                   const std::string &sequence,
                   const scoring::score_signature<alphabet_t> &score_fun) {

  auto [score, _len] =
      log_likelihood_part(tree, sequence, 0, sequence.size(), score_fun);
  return std::exp(score);
}

template <seqan3::alphabet alphabet_t>
double likelihood_context(ProbabilisticSuffixTreeMap<alphabet_t> &tree,
                          const std::string &sequence) {

  auto [score, _len] =
      log_likelihood_part(tree, sequence, 0, sequence.size(),
                          scoring::log_transition_prob<alphabet_t>);
  return std::exp(score);
}

} // namespace pst::distances::details

namespace pst::distances {

template <seqan3::alphabet alphabet_t>
double log_likelihood_s(ProbabilisticSuffixTreeMap<alphabet_t> &tree,
                        const std::string &sequence) {
  size_t order_max = tree.get_max_order();

  auto [log_likelihood, _l] =
      details::log_likelihood_part<alphabet_t>(tree, sequence);
  return log_likelihood;
}

template <seqan3::alphabet alphabet_t>
double max_adjusted_likelihood_s(ProbabilisticSuffixTreeMap<alphabet_t> &tree,
                                 const std::string &sequence) {
  size_t order_max = tree.get_max_order();

  auto [log_likelihood, _l] = details::log_likelihood_part<alphabet_t>(
      tree, sequence, 0, sequence.size(),
      details::scoring::transition_prob_max_adjusted<alphabet_t>);
  return log_likelihood;
}

template <seqan3::alphabet alphabet_t>
double log_likelihood(ProbabilisticSuffixTreeMap<alphabet_t> &tree,
                      const std::string &sequence) {
  size_t order_max = tree.get_max_order();

  auto bounds = pst::parallelize::get_bounds(sequence.size());
  std::vector<std::future<std::tuple<double, size_t>>> part_futures{};

  auto fun = [&](size_t start_index, size_t stop_index) {
    return details::log_likelihood_part<alphabet_t>(
        std::ref(tree), std::ref(sequence), start_index, stop_index,
        details::scoring::log_transition_prob<alphabet_t>, order_max);
  };

  for (auto [start_index, stop_index] : bounds) {
    part_futures.push_back(
        std::async(std::launch::async, fun, start_index, stop_index));
  }

  double log_likelihood = 0.0;
  for (auto &f : part_futures) {
    auto [log_lik, len] = f.get();
    log_likelihood += log_lik;
  }
  return log_likelihood;
}

template <seqan3::alphabet alphabet_t>
double
log_likelihood(ProbabilisticSuffixTreeMap<alphabet_t> &tree,
               const std::vector<alphabet_t> &sequence_dna,
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
    const std::vector<alphabet_t> &sequence_dna,
    const details::scoring::score_signature<alphabet_t> &score_fun) {

  auto [score, length] = details::log_likelihood_part(
      tree, sequence_dna, 0, sequence_dna.size(), score_fun);

  return -score / double(length);
}

template <seqan3::alphabet alphabet_t>
double negative_log_likelihood(
    ProbabilisticSuffixTreeMap<alphabet_t> &tree, const std::string &sequence,
    const details::scoring::score_signature<alphabet_t> &score_fun) {
  auto [score, length] = details::log_likelihood_part(
      tree, sequence, 0, sequence.size(), score_fun);
  return -score / double(length);
}

template <seqan3::alphabet alphabet_t>
double negative_log_likelihood(ProbabilisticSuffixTreeMap<alphabet_t> &tree,
                               const std::string &sequence) {
  auto [score, length] = details::log_likelihood_part(tree, sequence);
  return -score / double(length);
}

template <seqan3::alphabet alphabet_t>
double negative_log_likelihood(ProbabilisticSuffixTreeMap<alphabet_t> &tree,
                               const std::vector<alphabet_t> &sequence_dna) {
  return negative_log_likelihood(
      tree, sequence_dna, details::scoring::log_transition_prob<alphabet_t>);
}

template <seqan3::alphabet alphabet_t>
double negative_log_likelihood(ProbabilisticSuffixTree<alphabet_t> &tree,
                               const std::vector<alphabet_t> &sequence_dna) {
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

  return -log_likelihood / double(length - 1);
}

template <seqan3::alphabet alphabet_t>
double negative_log_likelihood_symmetric_(
    ProbabilisticSuffixTreeMap<alphabet_t> &tree,
    const std::vector<alphabet_t> &sequence,
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
                                  const std::vector<alphabet_t> &sequence) {

  return negative_log_likelihood_symmetric_<alphabet_t>(
      tree, sequence, details::scoring::log_transition_prob<alphabet_t>);
}

template <seqan3::alphabet alphabet_t>
double
negative_log_likelihood_symmetric_string(ProbabilisticSuffixTreeMap<alphabet_t> &tree,
                                  const std::string &sequence) {
  auto reverse_sequence = details::get_reverse_complement(sequence);

  return std::min(
      negative_log_likelihood<alphabet_t>(tree, sequence),
      negative_log_likelihood<alphabet_t>(tree, reverse_sequence));
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
    std::vector<alphabet_t> &sequence) {

  return negative_log_likelihood_symmetric_p(
      tree, sequence, details::scoring::log_transition_prob);
}

template <seqan3::alphabet alphabet_t>
inline double
negative_log_likelihood(ProbabilisticSuffixTreeMap<alphabet_t> &left,
                        ProbabilisticSuffixTreeMap<alphabet_t> &right,
                        size_t order) {
  auto left_sequence = left.generate_sequence(order);
  auto right_sequence = right.generate_sequence(order);

  auto right_nll = negative_log_likelihood(right, left_sequence);
  auto left_nll = negative_log_likelihood(left, right_sequence);

  return (right_nll + left_nll) / 2;
}

template <seqan3::alphabet alphabet_t>
inline double negative_log_likelihood_background(
    ProbabilisticSuffixTreeMap<alphabet_t> &left,
    ProbabilisticSuffixTreeMap<alphabet_t> &right, size_t order,
    size_t background_order) {
  auto left_sequence = left.generate_sequence(order);
  auto right_sequence = right.generate_sequence(order);

  auto transition_fun =
      details::scoring::specialise_background_log_transition_prob<alphabet_t>(
          background_order);

  auto right_nll =
      negative_log_likelihood(right, left_sequence, transition_fun);
  auto left_nll = negative_log_likelihood(left, right_sequence, transition_fun);

  return (right_nll + left_nll) / 2;
}

} // namespace pst::distances
