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

namespace pst::distances::details {
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
std::tuple<Eigen::VectorXd, Eigen::VectorXd>
get_d2star_vectors(ProbabilisticSuffixTreeMap<alphabet_t> &left,
                   ProbabilisticSuffixTreeMap<alphabet_t> &right,
                   size_t background_order) {
  size_t number_of_entries = left.counts.size() * left.valid_characters.size();
  Eigen::VectorXd left_vector(number_of_entries);
  Eigen::VectorXd right_vector(number_of_entries);

  Eigen::Index i = 0;
  iterate_included_in_both(
      left, right, [&](auto &context, auto &left_v, auto &right_v) {
        const auto background_context =
            get_background_context(context, background_order);
        auto left_background_v = left.counts[background_context];
        auto right_background_v = right.counts[background_context];

        for (auto &char_rank : left.valid_characters) {
          left_vector(i) = get_component(left, left_v, context,
                                         left_background_v, char_rank);

          right_vector(i) = get_component(right, right_v, context,
                                          right_background_v, char_rank);
          i++;
        }
      });

  return {left_vector.head(i), right_vector.head(i)};
}

template <seqan3::alphabet alphabet_t>
inline double core_d2star(ProbabilisticSuffixTreeMap<alphabet_t> &left,
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
          double left_component_value = pst::distances::details::get_component(
              left, left_v, context, left_background_v, char_rank);

          double right_component_value = pst::distances::details::get_component(
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
    double D2star = dot_product / (left_norm * right_norm);

    double d2star = 0.5 * (1 - D2star);
    return d2star;
  }
}

} // namespace pst::distances::details

namespace pst::distances {

template <seqan3::alphabet alphabet_t>
inline double d2star(ProbabilisticSuffixTreeMap<alphabet_t> &left,
                     ProbabilisticSuffixTreeMap<alphabet_t> &right,
                     size_t background_order = 0) {
  return details::core_d2star<alphabet_t>(left, right, background_order);
}

double d2star_cpp(std::string left_tree_string, std::string right_tree_string,
                  int background_order) {
  pst::ProbabilisticSuffixTreeMap<seqan3::dna5> left_tree{left_tree_string};
  pst::ProbabilisticSuffixTreeMap<seqan3::dna5> right_tree{right_tree_string};

  return d2star<seqan3::dna5>(left_tree, right_tree, background_order);
}

} // namespace pst::distances
