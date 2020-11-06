#pragma once

#include <robin_hood.h>

#include <Eigen/Dense>
#include <algorithm>
#include <numeric>

#include <seqan3/alphabet/concept.hpp>

#include "../probabilistic_suffix_tree_map.hpp"
#include "composition_vectors.hpp"

namespace pst::distances::details {

static bool use_cache = false;

std::vector<std::string> all_contexts{};

template <seqan3::alphabet alphabet_t>
std::vector<std::string>
get_all_contexts(int order, robin_hood::unordered_set<int> &valid_characters) {
  int n_valid_characters = valid_characters.size();
  int n_kmers = std::pow(n_valid_characters, order);

  if (all_contexts.size() == n_kmers) {
    return all_contexts;
  }

  std::vector<std::string> contexts(n_kmers, std::string(order, ' '));

  for (int n = 0; n < order; n++) {
    int n_char = 0;
    for (auto &char_rank : valid_characters) {
      int step_size = n_kmers / int(std::pow(n_valid_characters, (n + 1)));

      int start = step_size * n_char;
      int stop = step_size * (n_char + 1);

      while (stop <= n_kmers) {
        for (int i = start; i < stop; i++) {
          alphabet_t c = seqan3::assign_rank_to(char_rank, alphabet_t{});
          contexts[i][n] = c.to_char();
        }
        start = start + step_size * (n_valid_characters);
        stop = stop + step_size * (n_valid_characters);
      }

      n_char++;
    }
  }

  all_contexts = contexts;

  return contexts;
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

inline float square_and_sum(Eigen::VectorXd &composition_vector) {
  return composition_vector.squaredNorm();
}

inline float multiply_and_sum(Eigen::VectorXd &left, Eigen::VectorXd &right) {
  return left.dot(right);
}

float cosine_dissimilarity(Eigen::VectorXd &left, Eigen::VectorXd &right) {
  auto left_normalisation = details::square_and_sum(left);
  auto right_normalisation = details::square_and_sum(right);

  float normalisation = std::sqrt(left_normalisation * right_normalisation);

  if (normalisation == 0.0) {
    return 0.0;
  }

  float numerator = details::multiply_and_sum(left, right);
  float correlation = numerator / normalisation;

  return (1 - correlation) / 2;
}

template <seqan3::alphabet alphabet_t>
inline float core_cv(ProbabilisticSuffixTreeMap<alphabet_t> &left,
                     ProbabilisticSuffixTreeMap<alphabet_t> &right,
                     std::vector<std::string> &contexts, int background_order) {
  static std::unordered_map<std::string, Eigen::VectorXd> cache{};

  Eigen::VectorXd left_composition_vector =
      composition_vector(left, contexts, background_order);
  Eigen::VectorXd right_composition_vector =
      composition_vector(right, contexts, background_order);

  return cosine_dissimilarity(left_composition_vector,
                              right_composition_vector);
}

} // namespace pst::distances::details
namespace pst::distances {

template <seqan3::alphabet alphabet_t>
inline float cv(ProbabilisticSuffixTreeMap<alphabet_t> &left,
                ProbabilisticSuffixTreeMap<alphabet_t> &right,
                int background_order = 2) {
  auto contexts = details::get_shared_contexts(left, right);
  return details::core_cv<alphabet_t>(left, right, contexts, background_order);
}

template <seqan3::alphabet alphabet_t>
inline float cv_estimation(ProbabilisticSuffixTreeMap<alphabet_t> &left,
                           ProbabilisticSuffixTreeMap<alphabet_t> &right,
                           int order = 6, int background_order = 2) {
  if (details::all_contexts.empty()) {
    details::get_all_contexts<alphabet_t>(order, left.valid_characters);
  }
  details::use_cache = true;
  return details::core_cv<alphabet_t>(left, right, details::all_contexts,
                                      background_order);
}

} // namespace pst::distances

namespace pst::distances::properties {
constexpr bool symmetric_cv = true;
}
