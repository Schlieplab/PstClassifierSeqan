#pragma once

#include <robin_hood.h>

#include <algorithm>
#include <numeric>

#include <seqan3/alphabet/concept.hpp>

#include "../probabilistic_suffix_tree_map.hpp"

namespace pst::details {
template <seqan3::alphabet alphabet_t>
std::vector<std::string>
get_shared_contexts(ProbabilisticSuffixTreeMap<alphabet_t> left,
                    ProbabilisticSuffixTreeMap<alphabet_t> right) {

  auto left_terminal = left.get_terminal_nodes();
  auto right_terminal = right.get_terminal_nodes();

  robin_hood::unordered_set<std::string> context_set{};
  for (auto &v : right_terminal) {
    context_set.insert(std::move(v));
  }

  for (auto &v : left_terminal) {
    context_set.insert(std::move(v));
  }

  std::vector<std::string> contexts{context_set.begin(), context_set.end()};
  return contexts;
}

template <seqan3::alphabet alphabet_t>
std::vector<float>
composition_vector(ProbabilisticSuffixTreeMap<alphabet_t> tree,
                   std::vector<std::string> &contexts, int background_order) {
  std::vector<float> components(contexts.size() * tree.valid_characters.size());

  int i = 0;
  for (auto &context : contexts) {
    auto state = tree.get_closest_state(context);

    int background = std::min(background_order, int(state.size()));
    const auto background_state = state.substr(background);

    for (auto &char_rank : tree.valid_characters) {
      const float background_prob =
          tree.get_transition_probability(background_state, char_rank);
      const float prob = tree.get_transition_probability(state, char_rank);
      ;

      if (background_prob == 0.0) {
        components[i] = 0;
      } else {
        components[i] = (prob - background_prob) / background_prob;
      }
      i++;
    }
  }
  return components;
}

float square_and_sum(std::vector<float> &composition_vector) {
  return std::transform_reduce(composition_vector.begin(),
                               composition_vector.end(), 0.0, std::plus<>(),
                               [](float x) { return x * x; });
}

float multiply_and_sum(std::vector<float> &left, std::vector<float> &right) {
  return std::transform_reduce(right.begin(), right.end(), left.begin(), 0.0,
                               std::plus<>(), std::multiplies<>());
}

} // namespace pst::details
namespace pst {
template <seqan3::alphabet alphabet_t>
float cv(ProbabilisticSuffixTreeMap<alphabet_t> left,
         ProbabilisticSuffixTreeMap<alphabet_t> right,
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

} // namespace pst
