#pragma once

#include "robin_hood.h"

#include "Eigen/Dense"
#include <algorithm>
#include <numeric>

#include "seqan3/alphabet/concept.hpp"
#include "seqan3/alphabet/nucleotide/dna5.hpp"

#include "../probabilistic_suffix_tree_map.hpp"
#include "composition_vectors.hpp"

namespace pst::distances::details {

template <seqan3::alphabet alphabet_t>
std::tuple<Eigen::VectorXd, Eigen::VectorXd>
get_d2_vectors(ProbabilisticSuffixTreeMap<alphabet_t> &left,
               ProbabilisticSuffixTreeMap<alphabet_t> &right) {
  size_t number_of_entries = left.counts.size() * left.valid_characters.size();
  Eigen::VectorXd left_vector(number_of_entries);
  Eigen::VectorXd right_vector(number_of_entries);

  Eigen::Index i = 0;
  iterate_included_in_both<alphabet_t>(
      left, right, [&](auto context, auto left_v, auto right_v) {
        for (auto &char_rank : left.valid_characters) {
          left_vector(i) = left_v.next_symbol_probabilities[char_rank];
          right_vector(i) = right_v.next_symbol_probabilities[char_rank];
          i++;
        }
      });

  return {left_vector.head(i), right_vector.head(i)};
}

template <seqan3::alphabet alphabet_t>
inline double core_d2(ProbabilisticSuffixTreeMap<alphabet_t> &left,
                      ProbabilisticSuffixTreeMap<alphabet_t> &right) {
  auto [left_vector, right_vector] = get_d2_vectors(left, right);

  auto left_norm = left_vector.norm();
  auto right_norm = right_vector.norm();

  if (left_norm == 0 || right_norm == 0) {
    return 1.0;
  } else {
    double D2 = left_vector.dot(right_vector) / (left_norm * right_norm);

    double d2 = 0.5 * (1 - D2);
    return d2;
  }
}

} // namespace pst::distances::details

namespace pst::distances {

template <seqan3::alphabet alphabet_t>
inline double d2(ProbabilisticSuffixTreeMap<alphabet_t> &left,
                 ProbabilisticSuffixTreeMap<alphabet_t> &right) {
  return details::core_d2<alphabet_t>(left, right);
}

double d2_cpp(std::string left_tree_string, std::string right_tree_string) {
  pst::ProbabilisticSuffixTreeMap<seqan3::dna5> left_tree{left_tree_string};
  pst::ProbabilisticSuffixTreeMap<seqan3::dna5> right_tree{right_tree_string};

  return d2<seqan3::dna5>(left_tree, right_tree);
}

} // namespace pst::distances
