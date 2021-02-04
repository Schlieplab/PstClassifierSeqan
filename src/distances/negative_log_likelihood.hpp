#pragma once

#include <algorithm>
#include <string>

#include <seqan3/alphabet/concept.hpp>
#include <seqan3/range/views/complement.hpp>
#include <seqan3/range/views/to.hpp>
#include <seqan3/range/views/to_char.hpp>
#include <seqan3/std/ranges>

#include "../probabilistic_suffix_tree_map.hpp"
#include "composition_vectors.hpp"

namespace pst::distances {
template <seqan3::alphabet alphabet_t>
float negative_log_likelihood(ProbabilisticSuffixTreeMap<alphabet_t> &tree,
                              std::vector<alphabet_t> &sequence_dna) {
  float log_likelihood = 0.0;

  size_t order_max = tree.get_max_order();

  std::string sequence =
      sequence_dna | seqan3::views::to_char | seqan3::views::to<std::string>;

  size_t length = 0;
  for (size_t i = 0; i < sequence.size(); i++) {
    std::string subsequence = sequence.substr(
        std::max(size_t(0), i - order_max), std::min(i, order_max));

    char char_ = sequence[i];

    auto context = tree.get_closest_state(subsequence);

    auto probability = tree.get_transition_probability(context, char_);
    if (probability != 0) {
      length += 1;
      log_likelihood += std::log(probability);
    }
  }

  return -log_likelihood / float(length);
}

template <seqan3::alphabet alphabet_t>
float negative_log_likelihood_symmetric(
    ProbabilisticSuffixTreeMap<alphabet_t> &tree,
    std::vector<alphabet_t> &sequence) {

  std::vector<alphabet_t> reverse_sequence =
      sequence | std::views::reverse | seqan3::views::complement |
      seqan3::views::to<std::vector<alphabet_t>>;

  return std::min(negative_log_likelihood<alphabet_t>(tree, sequence),
                  negative_log_likelihood<alphabet_t>(tree, reverse_sequence));
}
} // namespace pst::distances
