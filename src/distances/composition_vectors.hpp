#pragma once

#include <Eigen/Dense>
#include <string>

#include "../probabilistic_suffix_tree_map.hpp"
#include "negative_log_likelihood.hpp"

namespace pst::distances::details {

std::string get_background_state(std::string &state, size_t background_order) {
  if (background_order == 0) {
    return "";
  }

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
      const double background_prob =
          tree.get_transition_probability(background_state, char_rank);
      if (background_prob == 0.0) {
        components(i) = 0;
      } else {
        const double prob = tree.get_transition_probability(state, char_rank);
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
    std::vector<std::string> &contexts, size_t background_order) {

  size_t number_of_entries = contexts.size() * tree.valid_characters.size();
  Eigen::VectorXd components(number_of_entries);

  size_t i = 0;
  for (auto &context : contexts) {
    auto state = tree.get_closest_state(context);

    const auto background_state = get_background_state(state, background_order);

    for (auto &char_rank : tree.valid_characters) {
      const double background_prob =
          tree.get_transition_probability(background_state, char_rank);
      if (background_prob == 0.0) {
        components(i) = 0;
      } else {
        const double prob = tree.get_transition_probability(state, char_rank);
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

    const auto background_state = get_background_state(state, background_order);

    for (auto &char_ : tree.valid_character_chars) {
      std::string extended_context{context + char_};

      const double background_prob = likelihood_context<alphabet_t>(
          tree, extended_context,
          [&](ProbabilisticSuffixTreeMap<alphabet_t> &tree,
              const std::string &context, char char_) -> double {
            return background_log_transition_prob<alphabet_t>(
                tree, context, char_, background_order);
          });

      if (background_prob == 0.0) {
        components(i) = 0;
      } else {
        const double prob = likelihood_context<alphabet_t>(
            tree, extended_context,
            [](ProbabilisticSuffixTreeMap<alphabet_t> &tree,
               const std::string &context, char char_) -> double {
              return log_transition_prob<alphabet_t>(tree, context, char_);
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

  size_t number_of_entries = contexts.size() * tree.valid_characters.size();
  Eigen::VectorXd components(number_of_entries);

  size_t i = 0;
  for (auto &context : contexts) {

    const double prob = likelihood_context<alphabet_t>(
        tree, context,
        [](ProbabilisticSuffixTreeMap<alphabet_t> &tree,
           const std::string &context, char char_) -> double {
          return log_transition_prob<alphabet_t>(tree, context, char_);
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

    double state_probability = double(std::get<0>(tree.counts[state])) /
                               double(std::get<0>(tree.counts[""]));

    const auto background_state = get_background_state(state, background_order);

    for (auto &char_rank : tree.valid_characters) {
      const double background_prob =
          tree.get_transition_probability(background_state, char_rank);
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

  auto left_terminal = left.get_terminal_nodes();
  auto right_terminal = right.get_terminal_nodes();

  std::vector<std::string> contexts{};
  // terminal nodes have to be sorted for this to work!
  std::set_union(left_terminal.begin(), left_terminal.end(),
                 right_terminal.begin(), right_terminal.end(),
                 std::back_inserter(contexts));

  if (contexts.empty()) {
    contexts.emplace_back("");
  }
  return contexts;
}

} // namespace pst::distances::details
