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

namespace pst::distances::details::dvstar {
template <seqan3::alphabet alphabet_t>
double get_component(ProbabilisticSuffixTreeMap<alphabet_t> &tree,
                     const hashmap_value<alphabet_t> &v,
                     const std::string &context,
                     const hashmap_value<alphabet_t> &background_v,
                     const size_t char_rank) {
  const double background_prob = std::get<1>(background_v)[char_rank];
  if (background_prob == 0.0) {
    return 0.0;
  } else {
    const double prob = std::get<1>(v)[char_rank];
    return prob / std::sqrt(background_prob);
  }
}

template <seqan3::alphabet alphabet_t>
double
get_component(ProbabilisticSuffixTreeMap<alphabet_t> &tree,
              const hashmap_value<alphabet_t> &v, const std::string &context,
              const std::string &background_context, const size_t char_rank) {
  const double background_prob = std::get<1>(tree.counts[context])[char_rank];
  if (background_prob == 0.0) {
    return 0.0;
  } else {
    const double prob = std::get<1>(v)[char_rank];
    return prob / std::sqrt(background_prob);
  }
}

template <seqan3::alphabet alphabet_t>
inline double core_dvstar(ProbabilisticSuffixTreeMap<alphabet_t> &left,
                          ProbabilisticSuffixTreeMap<alphabet_t> &right,
                          size_t background_order) {

  double dot_product = 0.0;

  double left_norm = 0.0;
  double right_norm = 0.0;

  pst::distances::details::iterate_included_in_both(
      left, right, [&](auto &context, auto &left_v, auto &right_v) {
        const auto background_context =
            pst::distances::details::get_background_context(context, 0);

        auto left_background_v = left.counts[background_context];
        auto right_background_v = right.counts[background_context];

        for (auto &char_rank : left.valid_characters) {
          double left_component_value = get_component(
              left, left_v, context, left_background_v, char_rank);

          double right_component_value = get_component(
              right, right_v, context, right_background_v, char_rank);

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
    double Dvstar = dot_product / (left_norm * right_norm);

    double dvstar = 0.5 * (1 - Dvstar);
    return dvstar;
  }
}

} // namespace pst::distances::details::dvstar

namespace pst::distances {

template <seqan3::alphabet alphabet_t>
inline double dvstar(ProbabilisticSuffixTreeMap<alphabet_t> &left,
                     ProbabilisticSuffixTreeMap<alphabet_t> &right,
                     size_t background_order = 0) {
  return details::dvstar::core_dvstar<alphabet_t>(left, right,
                                                  background_order);
}

double dvstar_cpp(const std::string &left_tree_string,
                  const std::string &right_tree_string, int background_order) {
  pst::ProbabilisticSuffixTreeMap<seqan3::dna5> left_tree{left_tree_string,
                                                          1.0};
  pst::ProbabilisticSuffixTreeMap<seqan3::dna5> right_tree{right_tree_string,
                                                           1.0};

  return dvstar<seqan3::dna5>(left_tree, right_tree, background_order);
}

} // namespace pst::distances
