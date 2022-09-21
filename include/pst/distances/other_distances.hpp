#pragma once

#include "../probabilistic_suffix_tree_map.hpp"
#include "composition_vectors.hpp"

namespace pst::distances::details {
using seqan3::operator""_dna5;

std::random_device rd;
std::mt19937 gen = std::mt19937{rd()};

std::string random_sequence(size_t length) {
  std::stringstream local_sequence{};
  std::uniform_int_distribution<> distrib(0, 3);
  std::array<char, 4> map{'A', 'C', 'G', 'T'};
  for (size_t i = 0; i < length; i++) {
    auto rand = distrib(gen);
    local_sequence << map[rand];
  }

  return local_sequence.str();
}
} // namespace pst::distances::details

namespace pst::distances {

template <seqan3::alphabet alphabet_t>
inline double jensenshannon(ProbabilisticSuffixTreeMap<alphabet_t> &left,
                            ProbabilisticSuffixTreeMap<alphabet_t> &right,
                            size_t background_order = 0) {

  auto left_root = left.counts[""];
  auto right_root = right.counts[""];

  double left_kl_sum = 0.0;
  double right_kl_sum = 0.0;

  pst::distances::details::iterate_included_in_both<alphabet_t>(
      left, right, [&](auto &context, auto &left_v, auto &right_v) {
        for (auto &char_rank : left.valid_characters) {
          double left_prob =
              left_v.next_symbol_probabilities[char_rank] /
              std::sqrt(left_root.next_symbol_probabilities[char_rank]);
          double right_prob =
              right_v.next_symbol_probabilities[char_rank] /
              std::sqrt(right_root.next_symbol_probabilities[char_rank]);

          double avg_next_prob =
              (right_v.next_symbol_probabilities[char_rank] +
               right_v.next_symbol_probabilities[char_rank]) /
              2;

          if (right_prob != 0.0 && left_prob != 0.0 && avg_next_prob != 0.0) {
            left_kl_sum += left_prob * std::log(left_prob / avg_next_prob);
            right_kl_sum += right_prob * std::log(right_prob / avg_next_prob);
          }
        }
      });

  return std::sqrt((left_kl_sum + right_kl_sum) / 2);
}

template <seqan3::alphabet alphabet_t>
inline double canberra(ProbabilisticSuffixTreeMap<alphabet_t> &left,
                       ProbabilisticSuffixTreeMap<alphabet_t> &right,
                       size_t background_order = 0) {

  auto left_root = left.counts[""];
  auto right_root = right.counts[""];
  double sum = 0.0;

  pst::distances::details::iterate_included_in_both<alphabet_t>(
      left, right, [&](auto &context, auto &left_v, auto &right_v) {
        for (auto &char_rank : left.valid_characters) {
          double left_prob =
              left_v.next_symbol_probabilities[char_rank] /
              std::sqrt(left_root.next_symbol_probabilities[char_rank]);
          double right_prob =
              right_v.next_symbol_probabilities[char_rank] /
              std::sqrt(right_root.next_symbol_probabilities[char_rank]);

          double denominator = std::abs(left_prob - right_prob);
          double numerator = left_prob + right_prob;

          sum += denominator / numerator;
        }
      });

  return sum;
}

template <seqan3::alphabet alphabet_t>
inline double correlation(ProbabilisticSuffixTreeMap<alphabet_t> &left,
                          ProbabilisticSuffixTreeMap<alphabet_t> &right) {

  auto left_root = left.counts[""];
  auto right_root = right.counts[""];
  double sum = 0.0;

  Eigen::VectorXd left_vector{};
  Eigen::VectorXd right_vector{};

  pst::distances::details::iterate_included_in_both<alphabet_t>(
      left, right, [&](auto &context, auto &left_v, auto &right_v) {
        for (auto &char_rank : left.valid_characters) {
          double left_prob =
              left_v.next_symbol_probabilities[char_rank] /
              std::sqrt(left_root.next_symbol_probabilities[char_rank]);
          double right_prob =
              right_v.next_symbol_probabilities[char_rank] /
              std::sqrt(right_root.next_symbol_probabilities[char_rank]);

          auto left_i = left_vector.size();
          left_vector.conservativeResize(left_i + 1);
          left_vector(left_i) = left_prob;

          auto right_i = right_vector.size();
          right_vector.conservativeResize(right_i + 1);
          right_vector(right_i) = right_prob;
        }
      });

  auto left_norm = (left_vector.array() - left_vector.mean());
  auto right_norm = (right_vector.array() - right_vector.mean());

  return 1 - left_norm.cwiseProduct(right_norm).sum() /
                 (std::sqrt(left_norm.square().sum()) *
                  std::sqrt(left_norm.square().sum()));
}

template <seqan3::alphabet alphabet_t>
double jaccard_estimation(ProbabilisticSuffixTreeMap<alphabet_t> &left,
                          ProbabilisticSuffixTreeMap<alphabet_t> &right) {
  int k = 21;

  auto left_root = left.counts[""];
  auto right_root = right.counts[""];
  double left_root_count = double(left_root.count);
  double right_root_count = double(right_root.count);

  double left_expected_frequency = left_root_count / std::pow(4.0, double(k));
  double right_expected_frequency = right_root_count / std::pow(4.0, double(k));

  double intersection = 0;
  double union_ = 0;

  for (int i = 0; i < 10000; i++) {
    auto seq = details::random_sequence(k);
    double left_prob = std::exp(log_likelihood_s(left, seq));
    double right_prob = std::exp(log_likelihood_s(right, seq));

    double est_left_count = left_prob * left_root_count;
    double est_right_count = right_prob * right_root_count;

    if ((est_left_count >= left_expected_frequency) ||
        (est_right_count >= right_expected_frequency)) {
      union_ += 1;
    }
    if ((est_left_count >= left_expected_frequency) &&
        (est_right_count >= right_expected_frequency)) {
      intersection += 1;
    }
  }
  double jaccard_v = intersection / union_;

  return 10 *
         (1 - std::pow(2.0 * jaccard_v / (jaccard_v + 1.0), 1.0 / double(k)));
}

template <seqan3::alphabet alphabet_t>
double jaccard_estimation_prob(ProbabilisticSuffixTreeMap<alphabet_t> &left,
                               ProbabilisticSuffixTreeMap<alphabet_t> &right) {
  int k = 11;

  auto left_root = left.counts[""];
  auto right_root = right.counts[""];

  double left_root_count = double(left_root.count);
  double right_root_count = double(right_root.count);

  double left_expected_frequency = left_root_count / std::pow(4.0, double(k));
  double right_expected_frequency = right_root_count / std::pow(4.0, double(k));

//  Eigen::VectorXd left_vector(1000);
//  Eigen::VectorXd right_vector(1000);
  double estimated_n_mut = 0.0;

  double n = 20000;
  for (int i = 0; i < n; i++) {
    auto seq = details::random_sequence(k);
    double left_prob = log_likelihood_s(left, seq);
    double right_prob = log_likelihood_s(right, seq);
    double est_left_count = left_prob * left_root_count;
    double est_right_count = right_prob * right_root_count;

//    left_vector(i) = left_prob;
//    right_vector(i) = right_prob;
    estimated_n_mut += std::abs(-left_prob + right_prob);
  }

  //  right_vector.normalize();
  //  left_vector.normalize();

//  double jaccard_sum = left_vector.cwiseMin(right_vector).sum() /
//                       left_vector.cwiseMax(right_vector).sum();

  // jaccard_sum *= 10.0;

  // jaccard_sum *= ((double(left_root.count) + double(right_root.count))
  // / 2.0);

  return estimated_n_mut * std::log((left_expected_frequency + right_expected_frequency) / 2.0)  / n; // 1 - std::pow(2.0 * jaccard_sum / (jaccard_sum
            // + 1.0), 1.0 / double(k));
}

} // namespace pst::distances