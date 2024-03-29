#pragma once

#include <limits>
#include <queue>
#include <shared_mutex>
#include <string>

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
                         size_t max_depth, size_t freq, double cutoff_value_,
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
                         size_t max_depth, size_t freq, double cutoff_value_,
                         size_t number_of_parameters,
                         const std::string &pruning_method,
                         bool multi_core = true, int parallel_depth = 2,
                         const double pseudo_count_amount = 1.0)
      : ProbabilisticSuffixTreeMap<alphabet_t>(
            id, sequence, max_depth, freq, number_of_parameters, pruning_method,
            multi_core, parallel_depth, pseudo_count_amount),
        cutoff_value(cutoff_value_) {}

  /*!\brief Reads a tree form a file.
   * \param[in] id The id of the model.
   * \param[in] sequence The text to construct from.
   */
  explicit KullbackLieblerTreeMap(const std::filesystem::path &filename)
      : ProbabilisticSuffixTreeMap<alphabet_t>(filename), cutoff_value(3.9075) {
  }

  /**! \brief Prune nodes with new Kullback-Leibler threshold.  New setting has
   * to be larger than previous one for any effect.
   *
   * \param threshold_ new threshold setting.
   */
  void reprune_similarity(double threshold_) {
    this->cutoff_value = threshold_;
    cutoff_prune();
  }

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
    std::shared_mutex remove_mutex{};

    if (this->multi_core) {
      robin_hood::unordered_map<std::string, std::queue<std::string>>
          map_queue{};

      for (std::string &v : pst_leaves) {
        auto key = v.substr(0, this->parallel_depth);
        map_queue[key].push(std::move(v));
      }

      std::vector<std::thread> threads{};
      for (auto &[key, queue] : map_queue) {
        threads.emplace_back(&KullbackLieblerTreeMap::prune_queue, this,
                             std::ref(queue), std::ref(remove_mutex));
      }

      for (auto &thread : threads) {
        if (thread.joinable()) {
          thread.join();
        }
      }
    } else {
      std::queue<std::string> queue{};
      for (std::string &v : pst_leaves) {
        queue.push(std::move(v));
      }
      this->prune_queue(queue, remove_mutex);
    }
  }

  void prune_queue(std::queue<std::string> &queue,
                   std::shared_mutex &remove_mutex) {
    while (!queue.empty()) {
      auto &node_label = queue.front();

      if (node_label.empty()) {
        queue.pop();
        continue;
      }

      double delta = calculate_delta(node_label, remove_mutex);

      if (delta < this->cutoff_value) {
        if (node_label.size() < this->parallel_depth) {
          // These are the only cases with conflicts
          remove_mutex.lock();
        }
        this->counts[node_label].is_included = false;

        auto parent_label = this->get_pst_parent(node_label);
        if (this->is_pst_leaf(parent_label)) {
          queue.push(std::move(parent_label));
        }
        if (node_label.size() < this->parallel_depth) {
          remove_mutex.unlock();
        }
      }
      queue.pop();
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
  double calculate_delta(const std::string &node_label,
                         std::shared_mutex &remove_mutex) {
    if (node_label.empty()) {
      return std::numeric_limits<double>::max();
    }

    const std::string parent_label = this->get_pst_parent(node_label);
    remove_mutex.lock_shared();
    auto node_counts = this->counts[node_label];
    auto parent_counts = this->counts[parent_label];
    remove_mutex.unlock_shared();

    return calculate_delta(node_counts, parent_counts);
  }

  /**! \brief Kullback-Liebler delta value for pruning.
   * \details
   * Based on the relative probabilities between the node and its parent
   * combined with the absolute count of the node.
   *
   * \param node_label Node to get delta for.
   * \return Delta value.
   */
  double calculate_delta(const std::string &node_label) {
    if (node_label.empty()) {
      return std::numeric_limits<double>::max();
    }

    const std::string parent_label = this->get_pst_parent(node_label);
    auto node_counts = this->counts[node_label];
    auto parent_counts = this->counts[parent_label];

    return calculate_delta(node_counts, parent_counts);
  }

  double calculate_delta(hashmap_value<alphabet_t> &node_counts,
                         hashmap_value<alphabet_t> &parent_counts) {
    double delta = 0.0;
    for (auto char_rank : this->valid_characters) {
      double prob, parent_prob;

      prob = node_counts.next_symbol_probabilities[char_rank];
      parent_prob = parent_counts.next_symbol_probabilities[char_rank];

      if (parent_prob == 0.0 || prob == 0.0) {
        continue;
      }

      delta += prob * std::log(prob / parent_prob);
    }
    delta *= node_counts.count;

    return delta;
  }
};

class KullbackLieblerTreeMap5 : public KullbackLieblerTreeMap<seqan3::dna5> {
public:
  KullbackLieblerTreeMap5() = default;
  KullbackLieblerTreeMap5(KullbackLieblerTreeMap5 const &) = default;
  ~KullbackLieblerTreeMap5() = default;

  KullbackLieblerTreeMap5(const std::string &id,
                          lst::details::sequence_t<seqan3::dna5> &sequence,
                          size_t max_depth, size_t freq, double cutoff_value_,
                          size_t number_of_parameters,
                          const std::string &pruning_method,
                          bool multi_core = true, int parallel_depth = 2)
      : KullbackLieblerTreeMap<seqan3::dna5>(
            id, sequence, max_depth, freq, cutoff_value_, number_of_parameters,
            pruning_method, multi_core, parallel_depth, 1.0) {}
};

lst::details::sequence_t<seqan3::dna5>
convert_to_dna(const std::string &sequence_) {
  return sequence_ | seqan3::views::char_to<seqan3::dna5> |
         seqan3::views::to<lst::details::sequence_t<seqan3::dna5>>;
}

KullbackLieblerTreeMap5
create_kl_map(const std::string &id, const std::string &sequence,
              size_t max_depth, size_t freq, double cutoff_value_,
              size_t number_of_parameters, const std::string &pruning_method,
              bool multi_core = true, int parallel_depth = 2) {
  auto seq = convert_to_dna(sequence);
  KullbackLieblerTreeMap5 t{id,
                            seq,
                            max_depth,
                            freq,
                            cutoff_value_,
                            number_of_parameters,
                            pruning_method,
                            multi_core,
                            parallel_depth};
  t.construct_tree();
  return t;
}

} // namespace pst
