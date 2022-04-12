#pragma once

#include <robin_hood.h>

#include <Eigen/Dense>
#include <algorithm>
#include <functional>
#include <numeric>

#include <seqan3/alphabet/concept.hpp>
#include <seqan3/alphabet/nucleotide/dna5.hpp>

#include "../probabilistic_suffix_tree_map.hpp"
#include "composition_vectors.hpp"
#include "negative_log_likelihood.hpp"

namespace pst::distances::details::d2star {
std::array<char, 5> rank_to_char{'A', 'C', 'G', 'N', 'T'};

template <seqan3::alphabet alphabet_t>
double get_component(ProbabilisticSuffixTreeMap<alphabet_t> &tree,
                     const hashmap_value<alphabet_t> &v,
                     const std::string &context, const size_t char_rank,
                     const size_t root_count) {
  char c = rank_to_char[char_rank];
  auto full_context = context + c;

  const double len = root_count - full_context.size() + 1;

  const auto [background_log_likelihood, _s] = log_likelihood_part(
      tree, full_context, 0, full_context.length(),
      pst::distances::details::scoring::log_transition_prob<alphabet_t>, 0);
  const double background_likelihood = std::exp(background_log_likelihood);

  if (background_likelihood == 0.0) {
    return 0.0;
  } else {
    const double prob = v.next_symbol_probabilities[char_rank];
    const double count = v.count;

    const double frequency = count * prob;

    const double numerator = frequency - len * background_likelihood;
    const double denominator = std::sqrt(background_likelihood * len);

    return numerator / denominator;
  }
}

template <seqan3::alphabet alphabet_t>
inline double core_d2star(ProbabilisticSuffixTreeMap<alphabet_t> &left,
                          ProbabilisticSuffixTreeMap<alphabet_t> &right,
                          size_t background_order) {
  double dot_product = 0.0;

  double left_norm = 0.0;
  double right_norm = 0.0;

  double left_root_count = left.counts[""].count;
  double right_root_count = right.counts[""].count;

  pst::distances::details::iterate_included_in_both<alphabet_t>(
      left, right, [&](auto &context, auto &left_v, auto &right_v) {
        const auto background_context =
            pst::distances::details::get_background_context(context, 0);

        auto left_background_v = left.counts[background_context];
        auto right_background_v = right.counts[background_context];

        for (auto &char_rank : left.valid_characters) {
          double left_component_value =
              get_component(left, left_v, context, char_rank, left_root_count);

          double right_component_value = get_component(
              right, right_v, context, char_rank, right_root_count);

          dot_product += left_component_value * right_component_value;
          left_norm += std::pow(left_component_value, 2.0);
          right_norm += std::pow(right_component_value, 2.0);
        }
      });

  left_norm = std::sqrt(left_norm);
  right_norm = std::sqrt(right_norm);

  if (left_norm == 0 || right_norm == 0) {
    return 1.0;
  } else {
    double D2star = dot_product / (left_norm * right_norm);

    double d2star = 0.5 * (1 - D2star);
    return d2star;
  }
}

} // namespace pst::distances::details::d2star

namespace pst::distances {

template <seqan3::alphabet alphabet_t>
inline double d2star(ProbabilisticSuffixTreeMap<alphabet_t> &left,
                     ProbabilisticSuffixTreeMap<alphabet_t> &right,
                     size_t background_order = 0) {
  return details::d2star::core_d2star<alphabet_t>(left, right,
                                                  background_order);
}

double d2star_cpp(std::string left_tree_string, std::string right_tree_string,
                  int background_order) {
  pst::ProbabilisticSuffixTreeMap<seqan3::dna5> left_tree{left_tree_string};
  pst::ProbabilisticSuffixTreeMap<seqan3::dna5> right_tree{right_tree_string};

  return d2star<seqan3::dna5>(left_tree, right_tree, background_order);
}

} // namespace pst::distances
