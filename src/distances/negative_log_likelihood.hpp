#pragma once

#include <algorithm>
#include <future>
#include <string>

#include <seqan3/alphabet/concept.hpp>
#include <seqan3/range/views/complement.hpp>
#include <seqan3/range/views/to.hpp>
#include <seqan3/range/views/to_char.hpp>
#include <seqan3/std/ranges>

#include "../probabilistic_suffix_tree.hpp"
#include "../probabilistic_suffix_tree_map.hpp"
#include "parallelize.hpp"

namespace pst::distances {

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
double log_likelihood(ProbabilisticSuffixTreeMap<alphabet_t> &tree,
                      std::string &sequence) {

  auto bounds = pst::parallelize::get_bounds(sequence.size());
  std::vector<std::future<double>> part_futures{};

  for (auto &[start_index, stop_index] : bounds) {
    part_futures.push_back(std::async(
        std::launch::async, log_likelihood_part<alphabet_t>, std::ref(tree),
        std::ref(sequence), start_index, stop_index));
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
        std::launch::async, details::log_likelihood_part_dna<alphabet_t>,
        std::ref(tree), std::ref(sequence_dna), start_index, stop_index));
  }

  double log_likelihood = 0.0;
  for (auto &f : part_futures) {
    log_likelihood += f.get();
  }
  return log_likelihood;
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
double
negative_log_likelihood_symmetric(ProbabilisticSuffixTreeMap<alphabet_t> &tree,
                                  std::vector<alphabet_t> &sequence) {

  std::vector<alphabet_t> reverse_sequence =
      sequence | std::views::reverse | seqan3::views::complement |
      seqan3::views::to<std::vector<alphabet_t>>;

  return std::min(negative_log_likelihood<alphabet_t>(tree, sequence),
                  negative_log_likelihood<alphabet_t>(tree, reverse_sequence));
}
} // namespace pst::distances
