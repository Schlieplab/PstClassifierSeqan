#pragma once

#include <seqan3/alphabet/concept.hpp>
#include <seqan3/alphabet/nucleotide/dna5.hpp>

#include "../probabilistic_suffix_tree_map.hpp"
#include "cv.hpp"
#include "negative_log_likelihood.hpp"

// Contains an implementation of the KL-divergence given by Ron et al. in 1996.

namespace pst::distances::details {
template <seqan3::alphabet alphabet_t>
double kl_divergence(ProbabilisticSuffixTreeMap<alphabet_t> &left,
                     ProbabilisticSuffixTreeMap<alphabet_t> &right,
                     size_t order) {

  if (details::all_contexts.empty()) {
    details::get_all_contexts<alphabet_t>(order, left.valid_characters);
  }
  details::use_cache = true;

  double kl_sum = 0.0;
  for (auto &context : all_contexts) {
    double left_prob = likelihood_context(left, context);
    double right_prob = likelihood_context(right, context);

    kl_sum += left_prob * std::log(left_prob / right_prob);
  }

  return kl_sum / double(order);
}

template <seqan3::alphabet alphabet_t>
double kl_divergence_both(ProbabilisticSuffixTreeMap<alphabet_t> &left,
                          ProbabilisticSuffixTreeMap<alphabet_t> &right) {

  double left_root_count = left.counts[""].count;
  double right_root_count = right.counts[""].count;

  double kl_sum = 0.0;
  double context_sum = 0;
  double n_contexts = 0;
  pst::distances::details::iterate_included_in_both<alphabet_t>(
      left, right, [&](auto &context, auto &left_v, auto &right_v) {
        double left_prob = left_v.count / left_root_count;
        double right_prob = right_v.count / right_root_count;

        if (right_prob != 0.0 && left_prob != 0.0) {
          kl_sum += left_prob * std::log(left_prob / right_prob);
          n_contexts += 1;
          context_sum += context.size();
        }
      });

    double v = kl_sum / double(context_sum / n_contexts);
   return v;
}
} // namespace pst::distances::details

namespace pst::distances {

template <seqan3::alphabet alphabet_t>
double kl_divergence(ProbabilisticSuffixTreeMap<alphabet_t> &left,
                     ProbabilisticSuffixTreeMap<alphabet_t> &right,
                     size_t order) {

  return details::kl_divergence(left, right, order);
}

template <seqan3::alphabet alphabet_t>
inline double
symmetric_kl_divergence(ProbabilisticSuffixTreeMap<alphabet_t> &left,
                        ProbabilisticSuffixTreeMap<alphabet_t> &right,
                        size_t order) {

  return (details::kl_divergence(left, right, order) +
          details::kl_divergence(right, left, order)) /
         2;
}

template <seqan3::alphabet alphabet_t>
inline double
symmetric_kl_divergence_both(ProbabilisticSuffixTreeMap<alphabet_t> &left,
                             ProbabilisticSuffixTreeMap<alphabet_t> &right) {
  double v = (details::kl_divergence_both(left, right) +
              details::kl_divergence_both(right, left)) /
             2;
  return v;
}

} // namespace pst::distances
