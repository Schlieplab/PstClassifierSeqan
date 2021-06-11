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
inline double core_d2(ProbabilisticSuffixTreeMap<alphabet_t> &left,
                      ProbabilisticSuffixTreeMap<alphabet_t> &right,
                      std::vector<std::string> &contexts) {
  Eigen::VectorXd left_vector = word_frequency_vector(left, contexts);
  Eigen::VectorXd right_vector = word_frequency_vector(right, contexts);

  double D2 = left_vector.dot(right_vector) /
              (left_vector.norm() * right_vector.norm());

  double d2 = 0.5 * (1 - D2);
  return d2;
}

} // namespace pst::distances::details

namespace pst::distances {

template <seqan3::alphabet alphabet_t>
inline double d2(ProbabilisticSuffixTreeMap<alphabet_t> &left,
                 ProbabilisticSuffixTreeMap<alphabet_t> &right) {
  auto contexts = details::get_shared_contexts(left, right);
  return details::core_d2<alphabet_t>(left, right, contexts);
}

} // namespace pst::distances
