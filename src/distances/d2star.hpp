#pragma once

#include <robin_hood.h>

#include <Eigen/Dense>
#include <algorithm>
#include <numeric>

#include <seqan3/alphabet/concept.hpp>
#include <seqan3/alphabet/nucleotide/dna5.hpp>

#include "../probabilistic_suffix_tree_map.hpp"
#include "composition_vectors.hpp"

namespace pst::distances::details {
template <seqan3::alphabet alphabet_t>
inline double core_d2star(ProbabilisticSuffixTreeMap<alphabet_t> &left,
                          ProbabilisticSuffixTreeMap<alphabet_t> &right,
                          std::vector<std::string> &contexts,
                          size_t background_order) {
  Eigen::VectorXd left_vector =
      adjusted_transition_frequency_vector(left, contexts, background_order);
  Eigen::VectorXd right_vector =
      adjusted_transition_frequency_vector(right, contexts, background_order);

  double D2star = left_vector.dot(right_vector) /
                  (left_vector.norm() * right_vector.norm());

  double d2star = 0.5 * (1 - D2star);
  return d2star;
}

} // namespace pst::distances::details

namespace pst::distances {

template <seqan3::alphabet alphabet_t>
inline double d2star(ProbabilisticSuffixTreeMap<alphabet_t> &left,
                     ProbabilisticSuffixTreeMap<alphabet_t> &right,
                     size_t background_order = 0) {
  auto contexts = details::get_shared_contexts(left, right);
  return details::core_d2star<alphabet_t>(left, right, contexts,
                                          background_order);
}

double d2star_cpp(std::string left_tree_string, std::string right_tree_string,
                  int background_order) {
  pst::ProbabilisticSuffixTreeMap<seqan3::dna5> left_tree{left_tree_string};
  pst::ProbabilisticSuffixTreeMap<seqan3::dna5> right_tree{right_tree_string};

  return d2star<seqan3::dna5>(left_tree, right_tree, background_order);
}

} // namespace pst::distances
