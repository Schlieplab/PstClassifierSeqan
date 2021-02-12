#pragma once

#include <algorithm>
#include <array>
#include <atomic>
#include <chrono>
#include <cmath>
#include <ctime>
#include <functional>
#include <mutex>
#include <shared_mutex>
#include <sstream>
#include <stack>
#include <string>
#include <thread>
#include <tuple>
#include <typeinfo>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include <seqan3/alphabet/concept.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/range/views/convert.hpp>
#include <seqan3/range/views/to.hpp>
#include <seqan3/range/views/to_char.hpp>

#include "search/lazy_suffix_tree.hpp"
#include "search/lazy_suffix_tree/iteration.hpp"

namespace pst {

enum Status : unsigned char {
  NONE = 1 << 0,
  INCLUDED = 1 << 1,
  EXCLUDED = 1 << 2,
};

template <seqan3::alphabet alphabet_t> struct PSTEntry {
  Status status;
  size_t count;
  std::array<float, seqan3::alphabet_size<alphabet_t>> probabilities;
};

/*!\brief The probabilistic suffix tree implementation.
 * \tparam alphabet_t   Alphabet type from Seqan3.
 *
 * \details
 *
 * The pst::ProbabilisticSuffixTree is an extension of a lazy suffix tree which
 * is pruned to only contain statistically significant information.
 *
 * ### General information
 *
 * The probabilistic suffix tree is equivalent to a variable length Markov
 * chain.  The implementation follows the AV-2 algorithm from [Schulz et
 * al.](https://doi.org/10.1007/978-3-540-87361-7_26), which is based on a lazy
 * suffix tree.
 *
 * ### Example
 *
 * \include pst-classifier.cpp
 *
 **
 */
template <seqan3::alphabet alphabet_t>
class ProbabilisticSuffixTree : public lst::LazySuffixTree<alphabet_t> {
public:
  ProbabilisticSuffixTree() = default;
  ProbabilisticSuffixTree(ProbabilisticSuffixTree const &) = default;
  ~ProbabilisticSuffixTree() = default;

  /*!\brief Constructor which assumes default values for all parameters.
   * \param[in] id The id of the model.
   * \param[in] sequence The text to construct from.
   */
  ProbabilisticSuffixTree(std::string id,
                          lst::details::sequence_t<alphabet_t> &sequence)

      : ProbabilisticSuffixTree(id, sequence, 15, 100, 1.2, 0, "cutoff", true,
                                2) {} // TODO get number of cores

  /*!\brief Constructor.
   * \param[in] id_ The id of the model.
   * \param[in] sequence_ The text to construct from.
   * \param[in] max_depth_ Max length of a branch/context in the tree.
   * \param[in] freq_ Min frequency of each context/node in the tree.
   * \param[in] number_of_parameters_ Number of parameters to keep in the model,
   * for "parameters" pruning.
   * \param[in] pruning_method_ Pruning method, either "cutoff" or "parameters"
   * (which depend on the `number_of_parameters_`)
   * \param[in] multi_core_ True for parallel execution.
   * \param[in] parallel_depth The maximum depth to spawn new processes, will
   * control task size as `alphabet_size ** depth`.
   */
  ProbabilisticSuffixTree(std::string id_,
                          lst::details::sequence_t<alphabet_t> &sequence_,
                          size_t max_depth_, size_t freq_,
                          size_t number_of_parameters_,
                          std::string pruning_method_, bool multi_core_ = true,
                          int split_depth_ = 2)
      : lst::LazySuffixTree<alphabet_t>(sequence_, multi_core_, split_depth_),
        id(id_), freq(freq_), max_depth(max_depth_),
        number_of_parameters(number_of_parameters_),
        pruning_method(pruning_method_) {
    using seqan3::operator""_dna4;
    seqan3::dna4_vector dna4{"ACGT"_dna4};

    std::vector<seqan3::gapped<alphabet_t>> characters{
        dna4 | seqan3::views::convert<seqan3::gapped<alphabet_t>> |
        seqan3::views::to<std::vector<seqan3::gapped<alphabet_t>>>};

    for (auto c : characters) {
      auto char_rank = seqan3::to_rank(c);
      valid_characters.insert(char_rank);
    }
  }

  /*! \brief Construct tree
   * Constructs the PST by applying the support and similarity pruning
   * steps from the VOMC construction algorithm.
   */
  void construct_tree() {
    this->support_pruning();
    this->add_suffix_links();
    this->add_reverse_suffix_links();
    this->similarity_pruning();
  }

  /**! \brief Support pruning phase of the algorithm.
   * \details
   * Extends a lazy suffix tree as long as the counts of each node is above
   * `freq` and the length is at most `max_depth`.
   *
   * After the full tree is built, it adds implicit nodes and suffix links.
   *
   */
  void support_pruning() {
    this->build_tree();
    entries[0] = {Status::INCLUDED, this->sequence.size(), {}};
  }

  /**! \brief Similarity pruning phase of the algorithm
   * \details
   * Prunes the tree bottom-up (by maintaining a status for each node).  A delta
   * value is assigned to each node, and the nodes with the smallest delta value
   * are removed.  The pruning either stops when a fixed number of parameters
   * are reached, or until a specified threshold value.
   */
  void similarity_pruning() {
    this->compute_probabilities();

    if (this->pruning_method == "cutoff") {
      this->cutoff_prune();
    } else if (this->pruning_method == "parameters") {
      this->parameters_prune();
    }
  }

  /*!\brief Debug printing
   * Prints the tree with the corresponding label, suffix links, delta values
   * etc.
   */
  void print() {
    this->debug_print_node(0, 0, 0);
    std::cout << std::endl;

    this->breadth_first_iteration_sequential(
        0, 0, false,
        [&](int node_index, int lcp, int edge_lcp, int node_count) -> bool {
          if (this->is_excluded(node_index)) {
            return true;
          }
          this->debug_print_node(node_index, lcp, edge_lcp);

          std::cout << std::endl;
          return true;
        });
  }

  void debug_print_node(size_t node_index, size_t lcp, size_t edge_lcp) {
    lst::LazySuffixTree<alphabet_t>::debug_print_node(node_index, lcp,
                                                      edge_lcp);

    if (this->suffix_links.size() > node_index / 2 &&
        this->suffix_links[node_index / 2] != max_size) {
      std::cout << "\tDelta: " << this->calculate_delta(node_index);

      std::cout << "\tPST Leaf: " << this->is_pst_leaf(node_index);

      std::cout << "\tTerminal: " << this->is_terminal(node_index);
    }

    if (this->entries.size() > node_index / 2) {
      std::cout << "\tStatus: " << this->entries[node_index / 2].status;
    }
  }

  /*! \brief Outputs the tree in a .tree format.
   *
   * \details
   * The output string will be in the following format:
   * Name: <id of the tree>
   * Date: <time of generation>
   * Tree: PST
   * Alphabet: <type of alphabet>
   * Number(nodes): <int, number of nodes>
   * Number(parameters). <int, number of parameters (terminal nodes * |alphabet|
   * - 1) Node: <index> <label> [ <counts of children> ] <count of node> [
   * <counts of forward children> ][ <index of children> ]
   * ...
   *
   * \returns std::string with the tree format.
   */
  std::string to_tree() {
    std::ostringstream tree_string;

    auto now = std::chrono::system_clock::now();
    std::time_t time = std::chrono::system_clock::to_time_t(now);

    tree_string << "Name: " << this->id << std::endl;
    tree_string << "Date: " << std::ctime(&time);
    tree_string << "Tree: PST" << std::endl;
    tree_string << "Alphabet: " << lst::get_alphabet_name<alphabet_t>()
                << std::endl;
    tree_string << "Number(nodes): " << nodes_in_tree() << std::endl;

    auto n_parameters =
        this->count_terminal_nodes() * (valid_characters.size() - 1);
    tree_string << "Number(parameters): " << n_parameters << std::endl;

    static std::mutex print_mutex{};

    std::unordered_map<size_t, size_t> lcps{};
    std::unordered_map<size_t, size_t> edge_lcps{};
    lcps[0] = 0;
    edge_lcps[0] = 0;
    this->breadth_first_iteration_sequential([&](size_t node_index, size_t lcp,
                                                 size_t edge_lcp,
                                                 size_t node_count) -> bool {
      std::lock_guard lock{print_mutex};
      lcps[node_index] = lcp;
      edge_lcps[node_index] = edge_lcp;
      return true;
    });

    std::unordered_map<size_t, size_t> iteration_order_indices{};
    size_t i = 0;
    this->pst_breadth_first_iteration(
        [&](size_t node_index, size_t level) -> bool {
          std::lock_guard lock{print_mutex};
          iteration_order_indices[node_index] = i;
          i++;
          return true;
        });

    this->pst_breadth_first_iteration(
        [&](size_t node_index, size_t level) -> bool {
          std::lock_guard lock{print_mutex};
          this->append_node_string(node_index, lcps[node_index],
                                   edge_lcps[node_index],
                                   iteration_order_indices, tree_string);
          return true;
        });

    return tree_string.str();
  }

  /**! \brief Counts the number of terminal nodes in the tree.
   *
   * \return vector of indices to all terminal nodes.
   */
  size_t count_terminal_nodes() {
    size_t n_terminal_nodes = 0;

    this->pst_breadth_first_iteration(
        0, 0, [&](size_t node_index, size_t level) -> bool {
          if (this->is_included(node_index) && this->is_terminal(node_index)) {
            n_terminal_nodes += 1;
          }

          return true;
        });

    return n_terminal_nodes;
  }

  /**! \brief Iterates over the nodes in the PST.
   *
   * \param start_index Index of the node to start iteration at.
   *
   */
  void
  pst_breadth_first_iteration(const std::function<bool(size_t, size_t)> &f) {
    pst_breadth_first_iteration(0, 0, f);
  }

  /**! \brief Iterates over the nodes in the PST.
   *
   * \param start_index Index of the node to start iteration at.
   * \param start_level Starting level for the iteration (0 for root, 1 for
   * ACGT)
   * \param f Callback for each child of start_index in a breadth first
   * fashion.  Is given node_index and depth of the node and should return
   * whether to iterate over the node's children.
   *
   */
  void
  pst_breadth_first_iteration(const size_t start_index,
                              const size_t start_level,
                              const std::function<bool(size_t, size_t)> &f) {
    pst_breadth_first_iteration_(start_index, start_level, f);
  }

  void
  pst_breadth_first_iteration_(const size_t start_index,
                               const size_t start_level,
                               const std::function<bool(size_t, size_t)> &f) {
    std::vector<std::thread> threads{};
    std::queue<std::tuple<size_t, size_t>> queue{};

    queue.emplace(start_index, start_level);

    while (!queue.empty()) {
      auto [node_index, level] = queue.front();
      queue.pop();

      if (f(node_index, level)) {
        if (this->multi_core && this->parallel_depth > level) {
          for (auto child : this->reverse_suffix_links[node_index / 2]) {
            if (child != max_size && !this->skip_node(child)) {
              threads.emplace_back(
                  &ProbabilisticSuffixTree::pst_breadth_first_iteration_, this,
                  child, level + 1, f);
            }
          }
        } else {
          for (auto child : this->reverse_suffix_links[node_index / 2]) {
            if (child != max_size && !this->skip_node(child)) {
              queue.emplace(child, level + 1);
            }
          }
        }
      }
    }

    for (auto &thread : threads) {
      if (thread.joinable()) {
        thread.join();
      }
    }
  }

  /**! \brief Returns the index of the children of in the PST.
   *
   * \param node_index Node to retrieve children of.
   * \return A vector of the children.
   */
  std::vector<size_t> get_pst_children(const size_t node_index) {
    std::vector<size_t> children{};
    for (auto child : this->reverse_suffix_links[node_index / 2]) {
      if (child != max_size && !this->skip_node(child)) {
        children.push_back(child);
      }
    }
    return children;
  }

  /**! \brief Return the parent of node in the PST.
   *
   * \param node_index Node to get parent of.
   * \return Parent index, or -1 if no parent.
   */
  size_t get_pst_parent(const size_t node_index) {
    return this->suffix_links[node_index / 2];
  }

  /**! \brief Returns count of the node
   * \details
   * Finds the count of the node by iterating the tree.  Saves the result in
   * a vector with the size of the tree / 2.
   *
   * \param[in] node_index
   * \return count of the node
   */
  size_t get_counts(size_t node_index) {
    if (node_index == max_size) {
      return 0;
    }
    if (node_index == 0) {
      return this->sequence.size();
    }

    auto c = this->entries[node_index / 2].count;
    if (c == max_size) {
      c = lst::details::node_occurrences(node_index, this->table);
      this->entries[node_index / 2].count = c;
    }
    return c;
  }

  /**! \brief Returns the next-symbol probabilities of the node
   *
   * \param[in] node_index
   * \return next-symbol probabilities of the node.
   */
  std::array<float, seqan3::alphabet_size<alphabet_t>>
  get_probabilities(size_t node_index) {
    if (node_index == max_size) {
      return {};
    }

    auto probs = this->entries[node_index / 2].probabilities;
    if (probs.size() == 0) {
      this->assign_node_probabilities(node_index);
    }
    return this->entries[node_index / 2].probabilities;
  }

  bool is_included(size_t node_index) {
    return (this->entries[node_index / 2].status & Status::INCLUDED) ==
           Status::INCLUDED;
  }

  bool is_excluded(size_t node_index) {
    return (this->entries[node_index / 2].status & Status::EXCLUDED) ==
           Status::EXCLUDED;
  }

  std::string id;

  std::unordered_set<int> valid_characters{};
  size_t freq;
  size_t max_depth;
  size_t number_of_parameters;
  std::string pruning_method;

  std::vector<PSTEntry<alphabet_t>> entries{};

protected:
  friend class ProbabilisticSuffixTreeTest;

  /**! \brief Builds the tree top-down.
   * \details
   * The lazy suffix tree is iterated in a breadth-first fashion and the nodes
   * are saved if the counts of each node is above `freq` and the length is at
   * most `max_depth`.
   *
   * Also saves the count of the node as well as if it is included or
   * excluded.
   */
  void build_tree() {
    std::shared_mutex entries_reallocate_mutex{};

    this->breadth_first_iteration(
        0, 0, true,
        [&](size_t node_index, size_t lcp, size_t &edge_lcp,
            size_t count) -> bool {
          if (this->is_leaf(node_index)) {
            std::lock_guard lock{this->table.mutex};
            expand_implicit_nodes(node_index, lcp, edge_lcp);
            // After adding implicit nodes, we may have to expand the entries
          }

          // If we reach the capacity of the entries vector, we will
          // reallocate it. This can lead to a race condition when a
          // different thread tries to set a value.  Therefore, we lock
          // changes to entries here.
          auto new_table_size = this->table.capacity() / 2;
          if (this->entries.capacity() < new_table_size) {
            std::unique_lock entries_reallocate_lock{entries_reallocate_mutex};
            this->resize_entries();
          }

          std::shared_lock lock{entries_reallocate_mutex};

          return this->check_node(node_index, lcp, edge_lcp, count);
        },
        []() {},
        [&](size_t node_index, size_t lcp, size_t &edge_lcp) {
          expand_implicit_nodes(node_index, lcp, edge_lcp);
        });

    this->entries.resize(this->table.size() / 2);
  }

  void expand_implicit_nodes(size_t node_index, size_t lcp, size_t &edge_lcp) {
    if (edge_lcp > 1) {
      auto max_extension = edge_lcp;
      if (lcp + edge_lcp > this->max_depth) {

        max_extension = this->max_depth - lcp;
      }

      this->add_implicit_nodes(node_index, max_extension);
      edge_lcp = 1;
    }
  }

  void resize_entries() {
    auto new_size = this->table.capacity() / 2;
    this->entries.reserve(new_size);
    this->entries.resize(new_size, {Status::NONE, max_size, {}});
  }

  /**! Check if the node with node_index should be included.
   *
   * \param node_index The index of the node.
   * \param lcp The LCP of the node.
   * \param edge_lcp The edge length (lcp) of the node.
   * \return if the node was included.
   */
  bool check_node(size_t node_index, size_t lcp, size_t edge_lcp,
                  size_t count) {
    auto label_start = this->table[node_index].value - lcp;
    auto label_end = this->table[node_index].value + edge_lcp;

    this->entries[node_index / 2].count = count;

    if (!this->is_leaf(node_index) &&
        this->include_node(label_start, label_end, edge_lcp, count)) {
      this->entries[node_index / 2].status = Status::INCLUDED;
      return true;
    } else {
      this->entries[node_index / 2].status = Status::EXCLUDED;

      // If this node is part of an expanded implicit node, we want to exclude
      // the rest of the implicit nodes as well.
      this->breadth_first_iteration_sequential(
          node_index, lcp, false,
          [&](size_t index, size_t lcp, size_t edge_lcp,
              size_t node_count) -> bool {
            this->entries[index / 2].status = Status::EXCLUDED;
            return true;
          });
      return false;
    }
  }

  /**! \brief Computes and saves the forward probabilities of each node.
   * \details
   * Probabilities for each node are calculated over the sum of the counts
   * of each (forward) child node.  Pseudo counts are added whenever any of
   * the children have a count of 0.
   *
   * Will allocate memory O(number of nodes after support pruning).
   */
  void compute_probabilities() {
    this->assign_node_probabilities(0);

    this->breadth_first_iteration(
        [&](size_t node_index, size_t lcp, size_t edge_lcp, size_t node_count) {
          if (this->is_leaf(node_index) || this->skip_node(node_index)) {
            return false;
          }

          this->assign_node_probabilities(node_index);
          return true;
        });
  }

  /**! \brief Assigns probabilities for the node.
   * \details
   * Sums the counts of all children, with pseudo counts, and assigns the
   * corresponding probabilities.  Iterates over all children to find counts,
   * to sum those counts, and to assign those counts.
   *
   * \param[in] node_index
   */
  void assign_node_probabilities(size_t node_index) {
    auto child_counts = this->get_child_counts(node_index, true);

    float child_sum = 0;
    for (auto c : child_counts) {
      child_sum += c;
    }

    for (size_t i = 0; i < seqan3::alphabet_size<alphabet_t>; i++) {
      this->entries[node_index / 2].probabilities[i] =
          float(child_counts[i]) / child_sum;
    }
  }

  /**! \brief Removes all nodes from the tree with a delta value below
   * threshold.
   * \details
   * The full tree is iterated to find the leaves in the
   * PST.  These leaves are then iterated, for each leaf the delta value is
   * calculated, and the node is removed if the delta value is below a
   * threshold.  For each removed node, the parent is added to be considered
   * if it now a leaf.
   */
  virtual void cutoff_prune() { parameters_prune(); }
  virtual float calculate_delta(size_t node_index) { return 0.0; }

  /**! \brief Removes all nodes until a specified number of parameters are
   * left.
   *
   * \details The full tree is iterated to find the leaves in the
   * PST.  These leaves are then iterated, for each leaf the delta value is
   * calculated.  These leaves are then removed in order of smallest
   * delta value until a given number of parameters are left in the tree.
   * For each removed node, the parent is added to be considered if
   * it now a leaf.
   */
  void parameters_prune() {
    auto pst_leaves = this->get_pst_leaves();

    using q_t = std::tuple<size_t, float>;
    auto cmp = [](q_t left, q_t right) -> bool {
      return std::get<1>(left) < std::get<1>(right);
    };
    std::priority_queue<q_t, std::vector<q_t>, decltype(cmp)> queue{cmp};

    for (auto v : pst_leaves) {
      queue.emplace(v, -this->calculate_delta(v));
    }

    auto n_terminal_nodes = this->count_terminal_nodes();

    auto alphabet_size = this->valid_characters.size() - 1;
    while (!queue.empty() &&
           n_terminal_nodes * alphabet_size > this->number_of_parameters) {
      auto [node_index, delta] = queue.top();
      queue.pop();

      if (node_index == 0) {
        continue;
      }

      auto parent_index = this->suffix_links[node_index / 2];

      this->entries[node_index / 2].status = Status::EXCLUDED;

      if (this->is_pst_leaf(parent_index) && parent_index != 0) {
        queue.emplace(parent_index, -this->calculate_delta(parent_index));
      }

      if (!this->became_terminal(parent_index, node_index)) {
        n_terminal_nodes -= 1;
      }
    }
  }

  /**! \brief Finds the counts of all (forward) children for the node.
   *
   * \param node_index Index to get child counts for.
   * \param with_pseudo_counts Flag for if pseudo counts should be used.
   * \return Array of counts for each children.
   */
  std::array<size_t, seqan3::alphabet_size<alphabet_t>>
  get_child_counts(size_t node_index, bool with_pseudo_counts) {
    std::array<size_t, seqan3::alphabet_size<alphabet_t>> child_counts{};

    int valid_children = 0;

    this->iterate_children(node_index, [&](size_t child_index) {
      auto sequence_index = this->get_sequence_index(child_index);
      auto character = this->get_character(sequence_index);

      if (character == seqan3::gap{}) {
        return;
      }

      auto character_rank = seqan3::to_rank(character);
      child_counts[character_rank] = get_counts(child_index);

      auto char_rank = seqan3::to_rank(character);
      if (this->valid_characters.find(char_rank) ==
          this->valid_characters.end()) {
        return;
      }

      valid_children += 1;
    });

    bool add_pseudo_counts = valid_children != this->valid_characters.size();
    if (with_pseudo_counts && add_pseudo_counts) {
      for (auto char_rank : this->valid_characters) {
        child_counts[char_rank] += 1;
      }
    }

    return child_counts;
  }

  /**! \brief Finds all leaves in the pst.
   *
   * \return vector of indices to all pst leaves.
   */
  std::vector<size_t> get_pst_leaves() {
    std::vector<size_t> bottom_nodes{};

    static std::mutex leaves_mutex{};

    this->pst_breadth_first_iteration(
        0, 0, [&](size_t node_index, size_t level) -> bool {
          if (this->is_included(node_index) && this->is_pst_leaf(node_index)) {
            std::lock_guard lock{leaves_mutex};
            bottom_nodes.emplace_back(node_index);
          }

          return true;
        });

    return bottom_nodes;
  }

  /**! \brief Checks if the node is a leaf in the PST.
   * \details
   * A node is a leaf in the PST if all of the reverse children corresponding
   * to the valid characters are excluded.
   *
   * \param[in] node_index Node to check.
   * \return If all children are missing.
   */
  bool is_pst_leaf(size_t node_index) {
    for (auto char_rank : this->valid_characters) {
      auto child_index = this->reverse_suffix_links[node_index / 2][char_rank];

      if (child_index == 0 || child_index == max_size) {
        continue;
      }

      if (this->is_included(child_index)) {
        return false;
      }
    }

    return true;
  }

  /**! \brief Checks if the node is terminal.
   * \details
   * A node is terminal if any of the reverse children corresponding to
   * the valid characters are excluded.
   *
   * \param[in] node_index Node index to check.
   * \return If the node is terminal (has any missing children).
   */
  bool is_terminal(size_t node_index) {
    for (auto char_rank : this->valid_characters) {
      auto child_index = this->reverse_suffix_links[node_index / 2][char_rank];

      if (child_index == 0) {
        continue;
      }

      if (child_index == max_size || this->is_excluded(child_index)) {
        return true;
      }
    }
    return false;
  }

  /**! \brief Checks if removed_index is the only child that is missing.
   * \details
   * The node_index has to correspond to the parent of removed_index.
   * Iterates through the valid characters of the reverse children of
   * node_index, and checks if any child other than removed_index is excluded.
   *
   * \param[in] node_index Index of node to check.
   * \param[in] removed_index Index of child to node, recently excluded.
   * \return True if removed index is the only missing child of node_index
   */
  bool became_terminal(size_t node_index, size_t removed_index) {
    for (auto char_rank : this->valid_characters) {
      auto child_index = this->reverse_suffix_links[node_index / 2][char_rank];

      if (child_index == 0) {
        continue;
      }

      if (child_index == max_size ||
          (this->is_excluded(child_index) && child_index != removed_index)) {
        return false;
      }
    }

    return is_terminal(node_index);
  }

  /**! \brief Append string for a node to the output stream.
   * \details
   * The format is Node: (index / 2) [ <reverse_children_count> ] count [
   * <children_count>] [ <reverse child indices> ]
   *
   * \param[in] node_index The node index to operate on
   * \param[in] lcp Longest common prefix of the node.
   * \param[in] edge_lcp Length of edge.
   * \param tree_string The output stream to write to.
   */
  void append_node_string(
      size_t node_index, size_t lcp, size_t edge_lcp,
      std::unordered_map<size_t, size_t> &iteration_order_indices,
      std::ostringstream &tree_string) {
    auto label = this->node_label(node_index, lcp, edge_lcp);

    if (node_index == 0) {
      label = "#";
    }

    tree_string << "Node: " << iteration_order_indices[node_index] << " "
                << label << " ";

    append_reverse_child_counts(node_index, tree_string);

    tree_string << " " << get_counts(node_index) << " ";

    append_child_counts(node_index, tree_string);

    //    tree_string << " ";

    append_reverse_children(node_index, iteration_order_indices, tree_string);

    tree_string << std::endl;
  }

  /**! \brief Append the count of forward children to the output stream.
   * \details
   * The format is [ x y z ... ]
   *
   * \param[in] node_index The node index to operate on.
   * \param tree_string The output stream to write to.
   */
  void append_child_counts(size_t node_index, std::ostringstream &tree_string) {

    auto child_counts = this->get_child_counts(node_index, true);

    std::vector<size_t> output(seqan3::alphabet_size<alphabet_t>, max_size);
    for (auto char_rank : this->valid_characters) {
      output[char_rank] = child_counts[char_rank];
    };

    tree_string << "[ ";

    for (auto count : output) {
      if (count == max_size) {
        continue;
      }
      tree_string << count << " ";
    }

    tree_string << "]";
  }

  /**! \brief Append the count of reverse children to the output stream.
   * \details
   * The format is [ x y z ... ]
   *
   * \param[in] node_index The node index to operate on.
   * \param tree_string The output stream to write to.
   */
  void append_reverse_child_counts(size_t node_index,
                                   std::ostringstream &tree_string) {

    std::vector<size_t> output(seqan3::alphabet_size<alphabet_t>, max_size);

    auto reverse_children = this->reverse_suffix_links[node_index / 2];
    for (auto char_rank : this->valid_characters) {

      auto reverse_child = reverse_children[char_rank];
      auto count = get_counts(reverse_child);

      output[char_rank] = count;
    };

    tree_string << "[ ";
    for (auto count : output) {
      if (count == max_size) {
        continue;
      }
      tree_string << count << " ";
    }
    tree_string << "]";
  }

  /**! \brief Append the index of reverse children to the output stream.
   * \details
   * For excluded nodes, a -1 is appended, for included nodes, the
   * internal node index / 2 is used (which to be fair is arbitrary).
   * The format is [ x y z ... ]
   *
   * \param[in] node_index The node index to operate on.
   * \param tree_string The output stream to write to.
   */
  void append_reverse_children(
      size_t node_index,
      std::unordered_map<size_t, size_t> &iteration_order_indices,
      std::ostringstream &tree_string) {

    std::vector<size_t> output(seqan3::alphabet_size<alphabet_t>, -2);

    auto reverse_children = this->reverse_suffix_links[node_index / 2];
    for (auto char_rank : this->valid_characters) {
      auto reverse_child = reverse_children[char_rank];

      if (reverse_child == 0 || reverse_child == max_size ||
          this->is_excluded(reverse_child)) {
        output[char_rank] = max_size;
      } else {
        output[char_rank] = iteration_order_indices[reverse_child];
      }
    };

    tree_string << "[ ";
    for (auto count : output) {
      if (count == -2) {
        continue;
      } else if (count == max_size) {
        tree_string << "-1 ";
      } else {
        tree_string << count << " ";
      }
    }
    tree_string << "]";
  }

  /**! \brief Counts the nodes in the tree.
   *
   * \return The number of nodes in the tree.
   */
  size_t nodes_in_tree() {
    std::atomic_size_t n_nodes = 0;
    this->pst_breadth_first_iteration(
        [&](size_t node_index, size_t level) -> bool {
          if (is_included(node_index)) {
            n_nodes += 1;
          }
          return true;
        });

    return n_nodes;
  }

  /**! \brief Determine if a node should be included in the tree.
   * \details
   * A node is included if the following three conditions are true:
   * - The length of the label is shorter than the max_depth of the tree.
   * - The label occurs at least as often as the specified min frequency.
   * - Every character in the label is included in the valid_characters.
   *
   * \param label_start start index of label in sequence.
   * \param label_end end index of label in sequence.
   * \param count number of times the label occurs in the sequence.
   * \return true if the node should be included.
   */
  bool include_node(size_t label_start, size_t label_end, size_t edge_lcp,
                    size_t count) {
    size_t label_length = label_end - label_start;

    // All characters before label_start + edge_lcp have already been checked.
    return label_length <= this->max_depth && count >= this->freq &&
           label_valid(label_start + edge_lcp, label_end);
  }

  /**! \brief Checks if a label is valid.
   *
   * \param[in] label_start starting index of the label in the sequence.
   * \param[in] label_end ending index of the label in the sequence.
   * \return true if the label is valid, false otherwise.
   */
  bool label_valid(size_t label_start, size_t label_end) {
    for (size_t i = label_start; i < label_end; i++) {
      auto character = this->get_character(i);
      auto char_rank = seqan3::to_rank(character);
      if (this->valid_characters.find(char_rank) ==
          this->valid_characters.end()) {
        return false;
      }
    }

    return true;
  }

  /** \copydoc lst::LazySuffixTree.skip_node()
   */
  bool skip_node(size_t node_index) {
    return this->is_unevaluated(node_index) || this->is_excluded(node_index);
  }
};
} // namespace pst
