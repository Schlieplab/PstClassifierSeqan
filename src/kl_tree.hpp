#pragma once

#include <queue>
#include <string>

#include <seqan3/alphabet/concept.hpp>
#include <seqan3/range/container/bitcompressed_vector.hpp>

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
                      seqan3::bitcompressed_vector<alphabet_t> &sequence)
      : KullbackLieblerTree(id, sequence, 15, 100, 1.2, 192, "cutoff") {}

  /*!\brief Constructor for parameters pruning.
   * \param[in] id_ The id of the model.
   * \param[in] sequence_ The text to construct from.
   * \param[in] max_depth_ Max length of a branch/context in the tree.
   * \param[in] freq_ Min frequency of each context/node in the tree.
   * \param[in] number_of_parameters_ Number of parameters to keep in the model,
   * for "parameters" pruning.
   */
  KullbackLieblerTree(std::string id,
                      seqan3::bitcompressed_vector<alphabet_t> &sequence,
                      size_t max_depth, size_t freq,
                      size_t number_of_parameters)
      : ProbabilisticSuffixTree<alphabet_t>(id, sequence, max_depth, freq, 1.2,
                                            number_of_parameters,
                                            "parameters") {}

  /*!\brief Constructor for threshold pruning.
   * \param[in] id_ The id of the model.
   * \param[in] sequence_ The text to construct from.
   * \param[in] max_depth_ Max length of a branch/context in the tree.
   * \param[in] freq_ Min frequency of each context/node in the tree.
   * \param[in] cutoff_value Cutoff value for the similarity-pruning.
   */
  KullbackLieblerTree(std::string id,
                      seqan3::bitcompressed_vector<alphabet_t> &sequence,
                      size_t max_depth, size_t freq, float cutoff_value_)
      : ProbabilisticSuffixTree<alphabet_t>(id, sequence, max_depth, freq, 192,
                                            "cutoff"),
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
   */
  KullbackLieblerTree(std::string id,
                      seqan3::bitcompressed_vector<alphabet_t> &sequence,
                      size_t max_depth, size_t freq, float cutoff_value_,
                      size_t number_of_parameters, std::string pruning_method)
      : ProbabilisticSuffixTree<alphabet_t>(id, sequence, max_depth, freq,
                                            number_of_parameters,
                                            pruning_method),
        cutoff_value(cutoff_value_) {
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

    std::queue<int> bottom_up{};
    for (auto v : pst_leaves) {
      bottom_up.push(v);
    }

    while (!bottom_up.empty()) {
      auto node_index = bottom_up.front();
      bottom_up.pop();

      if (node_index == 0) {
        continue;
      }

      float delta = calculate_delta(node_index);

      if (delta < this->cutoff_value) {
        this->status[node_index / 2] = Status::EXCLUDED;

        int parent_index = this->suffix_links[node_index / 2];
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
  float calculate_delta(int node_index) {
    int parent_index = this->suffix_links[node_index / 2];
    if (parent_index == -1) {
      throw std::invalid_argument(
          "[kl_delta] Given node does not have a parent.");
    }

    float delta = 0;
    for (auto char_rank : this->valid_characters) {
      float prob = this->probabilities[node_index / 2][char_rank];
      float parent_prob = this->probabilities[parent_index / 2][char_rank];

      if (parent_prob == 0 || prob == 0) {
        continue;
      }

      delta += prob * std::log(prob / parent_prob);
    }
    delta *= this->get_counts(node_index);

    return delta;
  }
};
} // namespace pst
