#pragma once

#include <algorithm>
#include <cmath>
#include <iterator>
#include <queue>
#include <string>

#include <seqan3/alphabet/concept.hpp>

#include "probabilistic_suffix_tree.hpp"

namespace pst {
template <seqan3::alphabet alphabet_t>
class PeresShieldsTree : public ProbabilisticSuffixTree<alphabet_t> {
  /**! \brief PST pruned using the Peres-Shields estimator.
   *
   */
public:
  PeresShieldsTree() = default;
  PeresShieldsTree(PeresShieldsTree const &) = default;
  ~PeresShieldsTree() = default;

  /*!\brief Constructor which assumes default values for all parameters.
   * \param[in] id The id of the model.
   * \param[in] sequence The text to construct from.
   */
  PeresShieldsTree(std::string id,
                   seqan3::bitcompressed_vector<alphabet_t> &sequence)
      : ProbabilisticSuffixTree<alphabet_t>(id, sequence, 15, 100, 192,
                                            "cutoff") {}

  /*!\brief Constructor.
   * \param[in] id_ The id of the model.
   * \param[in] sequence_ The text to construct from.
   * \param[in] max_depth_ Max length of a branch/context in the tree.
   * \param[in] freq_ Min frequency of each context/node in the tree.
   * \param[in] number_of_parameters_ Number of parameters to keep in the model,
   * for "parameters" pruning.
   * \param[in] pruning_method_ Pruning method, either "cutoff" or "parameters"
   * (which depend on the `number_of_parameters_`)
   */
  PeresShieldsTree(std::string id,
                   seqan3::bitcompressed_vector<alphabet_t> &sequence,
                   size_t max_depth, size_t freq)
      : ProbabilisticSuffixTree<alphabet_t>(id, sequence, max_depth, freq, 192,
                                            "cutoff") {}

  /*!\brief Constructor.
   * \param[in] id_ The id of the model.
   * \param[in] sequence_ The text to construct from.
   * \param[in] max_depth_ Max length of a branch/context in the tree.
   * \param[in] freq_ Min frequency of each context/node in the tree.
   * \param[in] number_of_parameters_ Number of parameters to keep in the model,
   * for "parameters" pruning.
   * \param[in] pruning_method_ Pruning method, either "cutoff" or "parameters"
   * (which depend on the `number_of_parameters_`)
   */
  PeresShieldsTree(std::string id,
                   seqan3::bitcompressed_vector<alphabet_t> &sequence,
                   size_t max_depth, size_t freq, size_t number_of_parameters,
                   std::string pruning_method, bool multi_core, int paralell_depth)
      : ProbabilisticSuffixTree<alphabet_t>(id, sequence, max_depth, freq,
                                            number_of_parameters,
                                            pruning_method, multi_core,
                                            paralell_depth) {}

protected:
  /**! \brief Removes all nodes from the tree with a delta value below
   * threshold.
   * \details
   * The full tree is iterated to find the leaves in the
   * PST.  These leaves are then iterated, for each leaf the delta value is
   * calculated, and the node is removed if the delta value is below a
   * threshold.  For each removed node, the parent is added to be considered if
   * it now a leaf.
   */
  void cutoff_prune() {
    this->pst_breadth_first_iteration(
        0, 0, [&](int64_t node_index, int64_t level) -> bool {
          int64_t parent_index = this->suffix_links[node_index / 2];

          if (node_index == 0) {
            return true;
          }

          bool prune_tree = this->should_prune_subtree(node_index, level);

          if (prune_tree) {
            this->exclude_subtree(node_index);
          }
          return !prune_tree;
        });
  }

  bool should_prune_subtree(int64_t node_index, int64_t level) {
    float delta = this->calculate_delta(node_index);


    auto deltas = this->get_subtree_level_deltas(node_index, level);
    if (deltas.empty()) {
      return true;
    }


    auto transitions = this->get_transitions(deltas);
    if (transitions.empty()) {
      return true;
    }

    float max_transition =
        *std::max_element(transitions.begin(), transitions.end());

    return max_transition == transitions[0];
  }

  std::vector<float> get_subtree_level_deltas(const int64_t node_index,
                                              const int64_t level) {
    std::vector<float> deltas{};
    float level_max_delta{this->calculate_delta(node_index)};
    int64_t current_level{level};

    this->pst_breadth_first_iteration(
        node_index, level,
        [&](const int64_t child_index, const int64_t child_level) -> bool {
          if (child_level > current_level) {
            deltas.push_back(level_max_delta);


            current_level += 1;
            level_max_delta = 0.0;
          }
          float delta = this->get_subtree_max_delta(child_index, child_level);
          level_max_delta = std::max(delta, level_max_delta);

          return true;
        });
    return deltas;
  }

  float get_subtree_max_delta(const int64_t node_index, const int64_t level) {
    float max_delta{0.0};

    this->pst_breadth_first_iteration(
        node_index, level,
        [&](const int64_t child_index, const int64_t child_level) -> bool {
          if (node_index == child_index) {
            return true;
          }

          float child_delta =
              this->calculate_delta_from(child_index, node_index);
          max_delta = std::max(max_delta, child_delta);

          return true;
        });

    return max_delta;
  }

  std::vector<float> get_transitions(const std::vector<float> &deltas) {
    std::vector<float> transitions{};
    float prev_delta{deltas[0]};
    std::transform(deltas.begin() + 1, deltas.end(),
                   std::back_inserter(transitions),
                   [&](float delta) -> float {
                      float transition{0.0};
                      if (delta != 0) {
                        transition = prev_delta / delta;
                      }
                      prev_delta = delta;
                      return transition;
                    });
    return transitions;
  }

  /**! \brief Mark every node that is a child in the PST of the node as
   * excluded.
   *
   * \param node_index Index of the root of the subtree.
   */
  void exclude_subtree(const int64_t node_index) {
    this->pst_breadth_first_iteration(node_index, 0,
        [&](const int64_t child_index, const int64_t child_level) {
          if (node_index != child_index)
            this->status[child_index / 2] = Status::EXCLUDED;
          return true;
        });
  }

  /**! \brief Peres-Shields delta value for pruning.
   * \details
   * Based on the relative counts between children of the node and its parent.
   * See 10.2202/1544-6115.1214 for details.
   *
   * \param node_index Index to get delta for.
   * \return Delta value.
   */
  float calculate_delta(const int64_t node_index) {
    int64_t parent_index = this->suffix_links[node_index / 2];
    if (parent_index == -1) {
      throw std::invalid_argument(
          "[ps_delta] Given node does not have a parent.");
    }

    return this->calculate_delta_from(node_index, parent_index);
  }

  /**! \brief Peres-Shields delta value for pruning relative to parent node.
   * \details
   * Based on the relative counts between children of the node and its parent.
   * See 10.2202/1544-6115.1214 for details.
   *
   * \param node_index Index to get delta for.
   * \param parent_index Index for the parent to compare node to.
   * \return Delta value.
   */
  float calculate_delta_from(const int64_t node_index, const int64_t parent_index) {
    float delta{0.0};
    float node_count{float(this->get_counts(node_index))};
    float parent_count{float(this->get_counts(parent_index))};

    std::array<int64_t, seqan3::alphabet_size<alphabet_t>> node_child_counts =
        this->get_child_counts(node_index, true);

    std::array<int64_t, seqan3::alphabet_size<alphabet_t>> parent_child_counts =
        this->get_child_counts(parent_index, true);

    for (auto char_rank : this->valid_characters) {
      float new_delta = std::abs(
          float(node_child_counts[char_rank]) -
          (float(parent_child_counts[char_rank]) / parent_count) * node_count);

      delta = std::max(delta, new_delta);
    }

    return delta;
  }
};
} // namespace pst
