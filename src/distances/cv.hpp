#pragma once

#include <robin_hood.h>

#include <Eigen/Dense>
#include <algorithm>
#include <numeric>

#include <seqan3/alphabet/concept.hpp>

#include "../probabilistic_suffix_tree_map.hpp"

namespace pst::distances::details {
template <seqan3::alphabet alphabet_t>
inline std::vector<std::string>
get_shared_contexts(ProbabilisticSuffixTreeMap<alphabet_t> &left,
                    ProbabilisticSuffixTreeMap<alphabet_t> &right) {

  auto left_terminal = left.get_terminal_nodes();
  auto right_terminal = right.get_terminal_nodes();

  std::vector<std::string> contexts{};
  // terminal nodes have to be sorted for this to work!
  std::set_union(left_terminal.begin(), left_terminal.end(),
                 right_terminal.begin(), right_terminal.end(),
                 std::back_inserter(contexts));
  return contexts;
}

template <seqan3::alphabet alphabet_t>
inline Eigen::VectorXd
composition_vector(ProbabilisticSuffixTreeMap<alphabet_t> &tree,
                   std::vector<std::string> &contexts, int background_order) {
  int number_of_entries = contexts.size() * tree.valid_characters.size();
  Eigen::VectorXd components(number_of_entries);

  int i = 0;
  for (auto &context : contexts) {
    auto state = tree.get_closest_state(context);

    int background = std::min(background_order, int(state.size()));
    const auto background_state = state.substr(background);

    for (auto &char_rank : tree.valid_characters) {
      const float background_prob =
          tree.get_transition_probability(background_state, char_rank);
      if (background_prob == 0.0) {
        components(i) = 0;
      } else {
        const float prob = tree.get_transition_probability(state, char_rank);
        components(i) = (prob - background_prob) / background_prob;
      }
      i++;
    }
  }
  return components;
}

inline float square_and_sum(Eigen::VectorXd &composition_vector) {
  return composition_vector.squaredNorm();
}

inline float multiply_and_sum(Eigen::VectorXd &left, Eigen::VectorXd &right) {
  return left.dot(right);
}

} // namespace pst::distances::details
namespace pst::distances {

template <seqan3::alphabet alphabet_t>
inline float cv(ProbabilisticSuffixTreeMap<alphabet_t> &left,
                ProbabilisticSuffixTreeMap<alphabet_t> &right,
                int background_order = 2) {
  auto contexts = details::get_shared_contexts(left, right);
  auto left_composition_vector =
      details::composition_vector(left, contexts, background_order);
  auto right_composition_vector =
      details::composition_vector(right, contexts, background_order);

  auto left_normalisation = details::square_and_sum(left_composition_vector);
  auto right_normalisation = details::square_and_sum(right_composition_vector);

  float normalisation = std::sqrt(left_normalisation * right_normalisation);

  if (normalisation == 0.0) {
    return 0.0;
  }

  float numerator = details::multiply_and_sum(left_composition_vector,
                                              right_composition_vector);
  float correlation = numerator / normalisation;

  return (1 - correlation) / 2;
}

} // namespace pst::distances

