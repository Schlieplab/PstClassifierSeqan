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
get_all_contexts(size_t order,
                 robin_hood::unordered_set<size_t> &valid_characters) {
  size_t n_valid_characters = valid_characters.size();
  size_t n_kmers = std::pow(n_valid_characters, order);

  if (all_contexts.size() == n_kmers) {
    return all_contexts;
  }

  std::vector<std::string> contexts(n_kmers, std::string(order, ' '));

  for (size_t n = 0; n < order; n++) {
    size_t n_char = 0;
    for (auto &char_rank : valid_characters) {
      size_t step_size =
          n_kmers / size_t(std::pow(n_valid_characters, (n + 1)));

      size_t start = step_size * n_char;
      size_t stop = step_size * (n_char + 1);

      while (stop <= n_kmers) {
        for (size_t i = start; i < stop; i++) {
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

inline double square_and_sum(Eigen::VectorXd &composition_vector) {
  return composition_vector.squaredNorm();
}

inline double multiply_and_sum(Eigen::VectorXd &left, Eigen::VectorXd &right) {
  return left.dot(right);
}

double cosine_dissimilarity(Eigen::VectorXd &left, Eigen::VectorXd &right) {
  auto left_normalisation = details::square_and_sum(left);
  auto right_normalisation = details::square_and_sum(right);

  double normalisation = std::sqrt(left_normalisation * right_normalisation);

  if (normalisation == 0.0) {
    return 0.0;
  }

  double numerator = details::multiply_and_sum(left, right);
  double correlation = numerator / normalisation;

  return (1 - correlation) / 2;
}

template <seqan3::alphabet alphabet_t>
inline double core_cv(ProbabilisticSuffixTreeMap<alphabet_t> &left,
                      ProbabilisticSuffixTreeMap<alphabet_t> &right,
                      std::vector<std::string> &contexts,
                      size_t background_order) {
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
inline double cv(ProbabilisticSuffixTreeMap<alphabet_t> &left,
                 ProbabilisticSuffixTreeMap<alphabet_t> &right,
                 size_t background_order = 2) {
  auto contexts = details::get_shared_contexts(left, right);
  return details::core_cv<alphabet_t>(left, right, contexts, background_order);
}

template <seqan3::alphabet alphabet_t>
inline double cv_estimation(ProbabilisticSuffixTreeMap<alphabet_t> &left,
                            ProbabilisticSuffixTreeMap<alphabet_t> &right,
                            size_t order = 6, size_t background_order = 2) {
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
