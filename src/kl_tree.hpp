#pragma once

#include <limits>
#include <queue>
#include <string>

#include <seqan3/alphabet/concept.hpp>

#include "probabilistic_suffix_tree.hpp"

namespace pst {

template <seqan3::alphabet alphabet_t>
class KullbackLieblerTree : public ProbabilisticSuffixTree<alphabet_t> {
  /**! \brief PST pruned using the Kullback-Liebler estimator.
   *
   */
public:
  KullbackLieblerTree() = default;
  KullbackLieblerTree(KullbackLieblerTree const &) = default;
  ~KullbackLieblerTree() = default;

  /*!\brief Constructor which assumes default values for all parameters.
   * \param[in] id The id of the model.
   * \param[in] sequence The text to construct from.
   */
  KullbackLieblerTree(std::string id,
                      lst::details::sequence_t<alphabet_t> &sequence)
      : ProbabilisticSuffixTree<alphabet_t>(id, sequence), cutoff_value(1.2) {}

  /*!\brief Constructor for parameters pruning.
   * \param[in] id_ The id of the model.
   * \param[in] sequence_ The text to construct from.
   * \param[in] max_depth_ Max length of a branch/context in the tree.
   * \param[in] freq_ Min frequency of each context/node in the tree.
   * \param[in] number_of_parameters_ Number of parameters to keep in the model,
   * for "parameters" pruning.
   * \param[in] multi_core True for parallel execution.
   * \param[in] parallel_depth The maximum depth to spawn new processes, will
   * control task size as `alphabet_size ** depth`.
   */
  KullbackLieblerTree(std::string id,
                      lst::details::sequence_t<alphabet_t> &sequence,
                      size_t max_depth, size_t freq,
                      size_t number_of_parameters, bool multi_core = true,
                      int parallel_depth = 2)
      : ProbabilisticSuffixTree<alphabet_t>(id, sequence, max_depth, freq,
                                            number_of_parameters, "parameters",
                                            multi_core, parallel_depth) {}

  /*!\brief Constructor for threshold pruning.
   * \param[in] id_ The id of the model.
   * \param[in] sequence_ The text to construct from.
   * \param[in] max_depth_ Max length of a branch/context in the tree.
   * \param[in] freq_ Min frequency of each context/node in the tree.
   * \param[in] cutoff_value Cutoff value for the similarity-pruning.
   * \param[in] multi_core True for parallel execution.
   * \param[in] parallel_depth The maximum depth to spawn new processes, will
   * control task size as `alphabet_size ** depth`.
   */
  KullbackLieblerTree(std::string id,
                      lst::details::sequence_t<alphabet_t> &sequence,
                      size_t max_depth, size_t freq, double cutoff_value_,
                      bool multi_core = true, int parallel_depth = 2)
      : ProbabilisticSuffixTree<alphabet_t>(id, sequence, max_depth, freq, 192,
                                            "cutoff", multi_core,
                                            parallel_depth),
        cutoff_value(cutoff_value_) {}

  /*!\brief Constructor for threshold pruning.
   * \param[in] id_ The id of the model.
   * \param[in] sequence_ The text to construct from.
   * \param[in] max_depth_ Max length of a branch/context in the tree.
   * \param[in] freq_ Min frequency of each context/node in the tree.
   * \param[in] number_of_parameters_ Number of parameters to keep in the model,
   * for "parameters" pruning.
   * \param[in] cutoff_value Cutoff value for the similarity-pruning.
   * \param[in] pruning_method_ Pruning method, either "cutoff" or "parameters"
   * (which depend on the `number_of_parameters_`)
   * \param[in] multi_core True for parallel execution.
   * \param[in] parallel_depth The maximum depth to spawn new processes, will
   * control task size as 4 ** depth.
   */
  KullbackLieblerTree(std::string id,
                      lst::details::sequence_t<alphabet_t> &sequence,
                      size_t max_depth, size_t freq, double cutoff_value_,
                      size_t number_of_parameters, std::string pruning_method,
                      bool multi_core, int parallel_depth)
      : ProbabilisticSuffixTree<alphabet_t>(
            id, sequence, max_depth, freq, number_of_parameters, pruning_method,
            multi_core, parallel_depth),
        cutoff_value(cutoff_value_) {}

protected:
  double cutoff_value = 1.2;

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
    auto pst_leaves = this->get_pst_leaves();

    std::queue<size_t> bottom_up{};
    for (auto v : pst_leaves) {
      bottom_up.push(v);
    }

    while (!bottom_up.empty()) {
      auto node_index = bottom_up.front();
      bottom_up.pop();

      if (node_index == 0) {
        continue;
      }

      double delta = calculate_delta(node_index);

      if (delta < this->cutoff_value) {
        this->entries[node_index].included = false;

        auto parent_index = this->suffix_links[node_index];
        if (this->is_pst_leaf(parent_index)) {
          bottom_up.push(parent_index);
        }
      }
    }
  }

  /**! \brief Kullback-Liebler delta value for pruning.
   * \details
   * Based on the relative probabilities between the node and its parent
   * combined with the absolute count of the node.
   *
   * \param node_index Index to get delta for.
   * \return Delta value.
   */
  double calculate_delta(size_t node_index) {
    if (node_index == 0) {
      return std::numeric_limits<double>::max();
    }

    auto parent_index = this->suffix_links[node_index];
    if (parent_index == -1) {
      throw std::invalid_argument(
          "[kl_delta] Given node does not have a parent.");
    }

    double delta = 0;
    for (auto char_rank : this->valid_characters) {
      double prob = this->entries[node_index].probabilities[char_rank];
      double parent_prob = this->entries[parent_index].probabilities[char_rank];

      if (parent_prob == 0.0 || prob == 0.0) {
        continue;
      }

      delta += prob * std::log(prob / parent_prob);
    }
    delta *= this->get_counts(node_index);

    return delta;
  }
};
} // namespace pst
