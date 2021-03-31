#pragma once

#include <Eigen/Dense>
#include <string>

#include "../probabilistic_suffix_tree_map.hpp"

namespace pst::distances {

std::string get_background_state(std::string &state, size_t background_order) {
  if (state.size() <= background_order) {
    return state;
  } else {
    size_t background = state.size() - background_order;
    return state.substr(background);
  }
}

template <seqan3::alphabet alphabet_t>
inline Eigen::VectorXd
composition_vector(ProbabilisticSuffixTreeMap<alphabet_t> &tree,
                   std::vector<std::string> &contexts,
                   size_t background_order) {

  size_t number_of_entries = contexts.size() * tree.valid_characters.size();
  Eigen::VectorXd components(number_of_entries);

  size_t i = 0;
  for (auto &context : contexts) {
    auto state = tree.get_closest_state(context);

    const auto background_state = get_background_state(state, background_order);

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

template <seqan3::alphabet alphabet_t>
inline Eigen::VectorXd composition_vector_state_probability_scaled(
    ProbabilisticSuffixTreeMap<alphabet_t> &tree,
    std::vector<std::string> &contexts, size_t background_order) {

  size_t number_of_entries = contexts.size() * tree.valid_characters.size();
  Eigen::VectorXd components(number_of_entries);

  size_t i = 0;
  for (auto &context : contexts) {
    auto state = tree.get_closest_state(context);

    float state_probability = float(std::get<0>(tree.counts[state])) /
                              float(std::get<0>(tree.counts[""]));

    const auto background_state = get_background_state(state, background_order);

    for (auto &char_rank : tree.valid_characters) {
      const float background_prob =
          tree.get_transition_probability(background_state, char_rank);
      if (background_prob == 0.0) {
        components(i) = 0;
      } else {
        const float prob = tree.get_transition_probability(state, char_rank);
        components(i) =
            state_probability * (prob - background_prob) / background_prob;
      }
      i++;
    }
  }
  return components;
}

} // namespace pst::distances
