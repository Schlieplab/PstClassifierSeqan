#pragma once

#include <Eigen/Dense>
#include <functional>
#include <string>

#include <seqan3/alphabet/concept.hpp>

#include "../probabilistic_suffix_tree_map.hpp"
#include "negative_log_likelihood.hpp"

namespace pst::distances::details {

std::string get_background_context(const std::string &state,
                                   const size_t background_order) {
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
    auto [state, val] = tree.get_closest_state(context);

    const auto background_context =
        get_background_context(state, background_order);

    for (auto &char_rank : tree.valid_characters) {
      const double background_prob =
          tree.get_transition_probability(background_context, char_rank);
      if (background_prob == 0.0) {
        components(i) = 0;
      } else {
        const double prob = tree.get_transition_probability(val, char_rank);
        components(i) = (prob - background_prob) / background_prob;
      }
      i++;
    }
  }
  return components;
}

template <seqan3::alphabet alphabet_t>
inline Eigen::VectorXd adjusted_transition_frequency_vector(
    ProbabilisticSuffixTreeMap<alphabet_t> &tree,
    const std::vector<std::string> &contexts, const size_t background_order) {

  size_t number_of_entries = contexts.size() * tree.valid_characters.size();
  Eigen::VectorXd components(number_of_entries);

  size_t i = 0;
  for (auto &context : contexts) {
    auto state = tree.get_closest_state(context);

    const auto background_context =
        get_background_context(context, background_order);

    for (auto &char_rank : tree.valid_characters) {
      const double background_prob =
          tree.get_transition_probability(background_context, char_rank);
      if (background_prob == 0.0) {
        components(i) = 0;
      } else {
        const double prob = tree.get_transition_probability(context, char_rank);
        components(i) = prob / std::sqrt(background_prob);
      }
      i++;
    }
  }
  return components;
}

template <seqan3::alphabet alphabet_t>
inline Eigen::VectorXd
adjusted_word_frequency_vector(ProbabilisticSuffixTreeMap<alphabet_t> &tree,
                               std::vector<std::string> &contexts,
                               size_t background_order) {

  size_t number_of_entries = contexts.size() * tree.valid_characters.size();
  Eigen::VectorXd components(number_of_entries);

  size_t i = 0;
  for (auto &context : contexts) {
    auto state = tree.get_closest_state(context);

    const auto background_context =
        get_background_context(state, background_order);

    for (auto &char_ : tree.valid_character_chars) {
      std::string extended_context{context + char_};

      const double background_prob = likelihood_context<alphabet_t>(
          tree, extended_context,
          [&](ProbabilisticSuffixTreeMap<alphabet_t> &tree,
              const std::string &context, char char_) -> double {
            return scoring::background_log_transition_prob<alphabet_t>(
                tree, context, char_, background_order);
          });

      if (background_prob == 0.0) {
        components(i) = 0;
      } else {
        const double prob = likelihood_context<alphabet_t>(
            tree, extended_context,
            [](ProbabilisticSuffixTreeMap<alphabet_t> &tree,
               const std::string &context, char char_) -> double {
              return scoring::log_transition_prob<alphabet_t>(tree, context,
                                                              char_);
            });

        const double count = tree.get_count("") - context.size() + 1;
        components(i) = count * prob / std::sqrt(count * background_prob);
      }
      i++;
    }
  }
  return components;
}

template <seqan3::alphabet alphabet_t>
inline Eigen::VectorXd
word_frequency_vector(ProbabilisticSuffixTreeMap<alphabet_t> &tree,
                      std::vector<std::string> &contexts) {

  size_t number_of_entries = contexts.size();
  Eigen::VectorXd components(number_of_entries);

  size_t i = 0;
  for (auto &context : contexts) {
    const double prob = likelihood_context<alphabet_t>(
        tree, context,
        [](ProbabilisticSuffixTreeMap<alphabet_t> &tree,
           const std::string &context, char char_) -> double {
          return scoring::log_transition_prob<alphabet_t>(tree, context, char_);
        });

    const double count = tree.get_count("") - context.size() + 1;

    components(i) = count * prob;

    i++;
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

    double state_probability = double(tree.counts[state].count) /
                               double(tree.counts[""].count);

    const auto background_context =
        get_background_context(state, background_order);

    for (auto &char_rank : tree.valid_characters) {
      const double background_prob =
          tree.get_transition_probability(background_context, char_rank);
      if (background_prob == 0.0) {
        components(i) = 0;
      } else {
        const double prob = tree.get_transition_probability(state, char_rank);
        components(i) =
            state_probability * (prob - background_prob) / background_prob;
      }
      i++;
    }
  }
  return components;
}

template <seqan3::alphabet alphabet_t>
inline std::vector<std::string>
get_shared_contexts(ProbabilisticSuffixTreeMap<alphabet_t> &left,
                    ProbabilisticSuffixTreeMap<alphabet_t> &right) {

  auto left_contexts = left.get_sorted_contexts();
  auto right_contexts = right.get_sorted_contexts();

  std::vector<std::string> contexts{};
  // terminal nodes have to be sorted for this to work!
  std::set_intersection(left_contexts.begin(), left_contexts.end(),
                        right_contexts.begin(), right_contexts.end(),
                        std::back_inserter(contexts));

  if (contexts.empty()) {
    contexts.emplace_back("");
  }
  return contexts;
}

template <seqan3::alphabet alphabet_t>
std::tuple<bool, hashmap_value<alphabet_t> &>
is_included_in_both(ProbabilisticSuffixTreeMap<alphabet_t> &left,
                    ProbabilisticSuffixTreeMap<alphabet_t> &right,
                    const std::string &context,
                    const hashmap_value<alphabet_t> &left_v) {
  bool left_included = left_v.is_included;
  if (left_included) {
    auto right_iter = right.counts.find(context);
    if (right_iter != right.counts.end()) {
      bool right_included = right_iter->second.is_included;
      return {right_included, right_iter->second};
    }
  }
  hashmap_value<alphabet_t> v{};
  return {false, v};
}

template <seqan3::alphabet alphabet_t>
void iterate_included_in_both(
    ProbabilisticSuffixTreeMap<alphabet_t> &left,
    ProbabilisticSuffixTreeMap<alphabet_t> &right,
    const std::function<void(const std::string &,
                             const hashmap_value<alphabet_t> &,
                             const hashmap_value<alphabet_t> &)> &f) {
  for (auto &[context, left_v] : left.counts) {
    auto [included_in_both, right_v] =
        is_included_in_both(left, right, context, left_v);

    if (included_in_both) {
      f(context, left_v, right_v);
    }
  }
}

} // namespace pst::distances::details
