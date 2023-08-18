#pragma once

#include "robin_hood.h"

#include "Eigen/Dense"
#include <Eigen/Core>
#include <algorithm>
#include <cmath>
#include <functional>
#include <numeric>

#include "seqan3/alphabet/concept.hpp"
#include "seqan3/alphabet/nucleotide/dna5.hpp"

#include "../probabilistic_suffix_tree_map.hpp"
#include "composition_vectors.hpp"

namespace pst::distances::details::dv {

double normalise_dv(double dot_product, double left_norm, double right_norm) {
  left_norm = std::sqrt(left_norm);
  right_norm = std::sqrt(right_norm);

  if (left_norm == 0 || right_norm == 0) {
    return 1.0;
  } else {
    double Dvstar = dot_product / (left_norm * right_norm);

    double dv = 0.5 * (1 - Dvstar);

    double angular_distance = 2 * std::acos(Dvstar) / M_PI;
    if (std::isnan(angular_distance)) {
      return 0.0;
    } else {
      return angular_distance;
    }
  }
}

template <seqan3::alphabet alphabet_t>
double get_root_count(ProbabilisticSuffixTreeMap<alphabet_t> &vlmc) {
  if (vlmc.root_state.count == 0) {
    vlmc.set_root_state();
  }
  return vlmc.root_state.count;
}

template <seqan3::alphabet alphabet_t>
std::array<std::vector<double>, 2>
core_dv_f(ProbabilisticSuffixTreeMap<alphabet_t> &left,
          ProbabilisticSuffixTreeMap<alphabet_t> &right,
          const std::string &context, const hashmap_value<alphabet_t> &left_v,
          const hashmap_value<alphabet_t> &right_v) {
  double dot_product = 0.0;

  double left_norm = 0.0;
  double right_norm = 0.0;

  std::array<std::vector<double>, 2> components{};

  for (auto &char_rank : left.valid_characters) {
    const double left_component_value =
        left_v.next_symbol_probabilities[char_rank];

    const double right_component_value =
        right_v.next_symbol_probabilities[char_rank];

    components[0].push_back(left_component_value);
    components[1].push_back(right_component_value);
  }

  return components;
}

template <seqan3::alphabet alphabet_t>
inline double core_dv(ProbabilisticSuffixTreeMap<alphabet_t> &left,
                      ProbabilisticSuffixTreeMap<alphabet_t> &right) {
  double dot_product = 0.0;

  double left_norm = 0.0;
  double right_norm = 0.0;

  pst::distances::details::iterate_included_in_both<alphabet_t>(
      left, right, [&](auto &context, auto &left_v, auto &right_v) {
        auto [left_components, right_components] =
            core_dv_f<alphabet_t>(left, right, context, left_v, right_v);

        for (int i = 0; i < left_components.size(); i++) {
          dot_product += left_components[i] * right_components[i];
          left_norm += std::pow(left_components[i], 2.0);
          right_norm += std::pow(right_components[i], 2.0);
        }
      });

  return normalise_dv(dot_product, left_norm, right_norm);
}
} // namespace pst::distances::details::dv

namespace pst::distances {

template <seqan3::alphabet alphabet_t>
inline double dv(ProbabilisticSuffixTreeMap<alphabet_t> &left,
                 ProbabilisticSuffixTreeMap<alphabet_t> &right) {
  return details::dv::core_dv<alphabet_t>(left, right);
}

double dv_cpp(const std::string &left_tree_string,
              const std::string &right_tree_string) {
  pst::ProbabilisticSuffixTreeMap<seqan3::dna5> left_tree{left_tree_string,
                                                          1.0};
  pst::ProbabilisticSuffixTreeMap<seqan3::dna5> right_tree{right_tree_string,
                                                           1.0};

  return dv<seqan3::dna5>(left_tree, right_tree);
}

} // namespace pst::distances
