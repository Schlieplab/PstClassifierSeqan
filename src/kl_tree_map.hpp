#pragma once

#include <limits>
#include <queue>
#include <string>

#include <seqan3/alphabet/concept.hpp>
#include <seqan3/range/container/bitcompressed_vector.hpp>

#include "probabilistic_suffix_tree_map.hpp"

namespace pst {

template <seqan3::alphabet alphabet_t>
class KullbackLieblerTreeMap : public ProbabilisticSuffixTreeMap<alphabet_t> {
  /**! \brief PST pruned using the Kullback-Liebler estimator.
   *
   */
public:
  KullbackLieblerTreeMap() = default;
  KullbackLieblerTreeMap(KullbackLieblerTreeMap const &) = default;
  ~KullbackLieblerTreeMap() = default;

  /*!\brief Constructor which assumes default values for all parameters.
   * \param[in] id The id of the model.
   * \param[in] sequence The text to construct from.
   */
  KullbackLieblerTreeMap(std::string id,
                         lst::details::sequence_t<alphabet_t> &sequence)
      : ProbabilisticSuffixTreeMap<alphabet_t>(id, sequence),
        cutoff_value(3.9075) {}

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
  KullbackLieblerTreeMap(std::string id,
                         lst::details::sequence_t<alphabet_t> &sequence,
                         size_t max_depth, size_t freq,
                         size_t number_of_parameters, bool multi_core = true,
                         int parallel_depth = 2)
      : ProbabilisticSuffixTreeMap<alphabet_t>(
            id, sequence, max_depth, freq, number_of_parameters, "parameters",
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
  KullbackLieblerTreeMap(std::string id,
                         lst::details::sequence_t<alphabet_t> &sequence,
                         size_t max_depth, size_t freq, float cutoff_value_,
                         bool multi_core = true, int parallel_depth = 2)
      : ProbabilisticSuffixTreeMap<alphabet_t>(id, sequence, max_depth, freq,
                                               192, "cutoff", multi_core,
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
  KullbackLieblerTreeMap(std::string id,
                         lst::details::sequence_t<alphabet_t> &sequence,
                         size_t max_depth, size_t freq, float cutoff_value_,
                         size_t number_of_parameters,
                         std::string pruning_method, bool multi_core = true,
                         int parallel_depth = 2)
      : ProbabilisticSuffixTreeMap<alphabet_t>(
            id, sequence, max_depth, freq, number_of_parameters, pruning_method,
            multi_core, parallel_depth),
        cutoff_value(cutoff_value_) {}

  /*!\brief Reads a tree form a file.
   * \param[in] id The id of the model.
   * \param[in] sequence The text to construct from.
   */
  KullbackLieblerTreeMap(std::filesystem::path &filename)
      : ProbabilisticSuffixTreeMap<alphabet_t>(filename), cutoff_value(3.9075) {
  }

protected:
  float cutoff_value = 1.2;

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

    std::queue<std::string> bottom_up{};
    for (auto v : pst_leaves) {
      bottom_up.push(v);
    }

    while (!bottom_up.empty()) {
      auto &node_label = bottom_up.front();

      if (node_label.empty()) {
        bottom_up.pop();
        continue;
      }

      float delta = calculate_delta(node_label);

      if (delta < this->cutoff_value) {
        this->status.erase(node_label);

        auto parent_label = this->get_pst_parent(node_label);
        if (this->is_pst_leaf(parent_label)) {
          bottom_up.push(parent_label);
        }
      }
      bottom_up.pop();
    }
  }

  /**! \brief Kullback-Liebler delta value for pruning.
   * \details
   * Based on the relative probabilities between the node and its parent
   * combined with the absolute count of the node.
   *
   * \param node_label Node to get delta for.
   * \return Delta value.
   */
  float calculate_delta(const std::string &node_label) {
    if (node_label.empty()) {
      return std::numeric_limits<float>::max();
    }

    const std::string parent_label = this->get_pst_parent(node_label);

    float delta = 0;
    for (auto char_rank : this->valid_characters) {
      float prob = this->probabilities[node_label][char_rank];
      float parent_prob = this->probabilities[parent_label][char_rank];

      if (parent_prob == 0 || prob == 0) {
        continue;
      }

      delta += prob * std::log(prob / parent_prob);
    }
    delta *= this->counts[node_label];

    return delta;
  }
};
} // namespace pst
