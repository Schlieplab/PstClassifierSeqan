#pragma once

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <ctime>
#include <functional>
#include <sstream>
#include <stack>
#include <string>
#include <tuple>
#include <unordered_set>
#include <vector>

#include <seqan3/alphabet/concept.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/range/all.hpp>
#include <seqan3/range/container/bitcompressed_vector.hpp>
#include <seqan3/range/views/to.hpp>

#include "search/lazy_suffix_tree.hpp"
#include "search/lazy_suffix_tree/iteration.hpp"

namespace pst {

enum Status : unsigned char {
  NONE = 1 << 0,
  INCLUDED = 1 << 1,
  EXCLUDED = 1 << 2,
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
  friend class ProbabilisticSuffixTreeTest;

  ProbabilisticSuffixTree() = default;
  ProbabilisticSuffixTree(ProbabilisticSuffixTree const &) = default;
  ~ProbabilisticSuffixTree() = default;

  /*!\brief Constructor which assumes default values for all parameters.
   * \param[in] id The id of the model.
   * \param[in] sequence The text to construct from.
   */
  ProbabilisticSuffixTree(std::string id,
                          seqan3::bitcompressed_vector<alphabet_t> &sequence)
      : ProbabilisticSuffixTree(id, sequence, 15, 100, 1.2, 0, "cutoff") {}

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
  ProbabilisticSuffixTree(std::string id_,
                          seqan3::bitcompressed_vector<alphabet_t> &sequence_,
                          size_t max_depth_, size_t freq_,
                          size_t number_of_parameters_,
                          std::string pruning_method_)
      : lst::LazySuffixTree<alphabet_t>(sequence_), id(id_), freq(freq_),
        max_depth(max_depth_), number_of_parameters(number_of_parameters_),
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
    this->similarity_pruning();
  }

  /*!\brief Debug printing
   * Prints the tree with the corresponding label, suffix links, delta values
   * etc.
   */
  void print() {
    this->debug_print_node(0, 0, 0);
    seqan3::debug_stream << std::endl;

    this->breadth_first_iteration(
        0, 0, false, [&](int node_index, int lcp, int edge_lcp) -> bool {
          if (this->is_excluded(node_index)) {
            return true;
          }
          this->debug_print_node(node_index, lcp, edge_lcp);

          seqan3::debug_stream << std::endl;
          return true;
        });
  }

  void debug_print_node(int node_index, int lcp, int edge_lcp) {
    lst::LazySuffixTree<alphabet_t>::debug_print_node(node_index, lcp,
                                                      edge_lcp);

    if (this->reverse_suffix_links.size() > node_index / 2) {
      seqan3::debug_stream << "\tLeaf: " << this->is_leaf(node_index);
    }

    if (this->suffix_links.size() > node_index / 2 &&
        this->suffix_links[node_index / 2] != -1) {
      seqan3::debug_stream << "\tDelta: " << this->calculate_delta(node_index);

      seqan3::debug_stream << "\tPST Leaf: " << this->is_pst_leaf(node_index);

      seqan3::debug_stream << "\tTerminal: " << this->is_terminal(node_index);
    }

    if (this->status.size() > node_index / 2) {
      seqan3::debug_stream << "\tStatus: " << this->status[node_index / 2];
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
    tree_string << "Alphabet: " << typeid(alphabet_t).name() << std::endl;
    tree_string << "Number(nodes): " << nodes_in_tree() << std::endl;

    auto n_parameters =
        this->count_terminal_nodes() * (valid_characters.size() - 1);
    tree_string << "Number(parameters): " << n_parameters << std::endl;

    this->append_node_string(0, 0, 0, tree_string);

    this->breadth_first_iteration(
        0, 0, false, [&](int node_index, int lcp, int edge_lcp) -> bool {
          if (this->is_included(node_index)) {
            this->append_node_string(node_index, lcp, edge_lcp, tree_string);
          }

          return true;
        });

    return tree_string.str();
  }

  /**! \brief Counts the number of terminal nodes in the tree.
   *
   * \return vector of indices to all terminal nodes.
   */
  int count_terminal_nodes() {
    int n_terminal_nodes = 0;

    this->breadth_first_iteration(
        0, 0, false, [&](int node_index, int lcp, int edge_lcp) -> bool {
          if (this->is_included(node_index) && this->is_terminal(node_index)) {
            n_terminal_nodes += 1;
          }
          return true;
        });

    return n_terminal_nodes;
  }

  void pst_breadth_first_iteration(const int start_index, const int start_level,
                                   const std::function<bool(int, int)> &f) {
    std::queue<std::tuple<int, int>> queue{};

    queue.emplace(start_index, start_level);

    while (!queue.empty()) {
      auto [node_index, level] = queue.front();
      queue.pop();

      if (f(node_index, level)) {
        for (auto child : this->reverse_suffix_links[node_index / 2]) {
          if (child != -1 && !this->skip_node(child)) {
            queue.emplace(child, level + 1);
          }
        }
      }
    }
  }

  std::string id;

  std::unordered_set<int> valid_characters{};
  size_t freq;
  size_t max_depth;
  size_t number_of_parameters;
  std::string pruning_method;

  std::vector<Status> status{};

  std::vector<int> counts{};
  std::vector<std::array<float, seqan3::alphabet_size<alphabet_t>>>
      probabilities{};

protected:
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
    this->expand_implicit_nodes();
    this->add_implicit_node_status();
    this->add_suffix_links();
    this->counts.resize(this->table.size() / 2, -1);
    status[0] = Status::INCLUDED;
  }

  /**! \brief Similarity pruning phase of the algorithm
   * \details
   * Prunes the tree bottom-up (by maintaining a status for each node).  A delta
   * value is assigned to each node, and the nodes with the smallest delta value
   * are removed.  The pruning either stops when a fixed number of parameters
   * are reached, or until a specified threshold value.
   */
  void similarity_pruning() {
    this->add_reverse_suffix_links();
    this->compute_probabilities();
    if (this->pruning_method == "cutoff") {
      this->cutoff_prune();
    } else if (this->pruning_method == "parameters") {
      this->parameters_prune();
    }
  }

  /**! \brief Builds the tree top-down.
   * \details
   * The lazy suffix tree is iterated in a breadth-first fashion and the nodes
   * are saved if the counts of each node is above `freq` and the length is at
   * most `max_depth`.
   *
   * Also saves the count of the node as well as if it is included or excluded.
   */
  void build_tree() {
    this->breadth_first_iteration(
        0, 0, true, [&](int node_index, int lcp, int edge_lcp) -> bool {
          int count = lst::details::node_occurrences(node_index, this->table,
                                                     this->flags);
          if (edge_lcp > 1) {
            int max_extension = edge_lcp;
            if (lcp + edge_lcp > this->max_depth) {
              max_extension = this->max_depth - lcp;
            }
            this->add_implicit_nodes(node_index, max_extension);
            edge_lcp = 1;
          }

          this->counts.resize(this->table.size() / 2, -1);
          this->status.resize(this->table.size() / 2, Status::NONE);

          return this->check_node(node_index, lcp, edge_lcp, count);
        });
  }

  /**! Check if the node with node_index should be included.
   *
   * \param node_index The index of the node.
   * \param lcp The LCP of the node.
   * \param edge_lcp The edge length (lcp) of the node.
   * \return if the node was included.
   */
  bool check_node(int node_index, int lcp, int edge_lcp, int count) {
    this->counts[node_index / 2] = count;

    int label_start = this->table[node_index] - lcp;
    int label_end = this->table[node_index] + edge_lcp;

    if (!this->is_leaf(node_index) &&
        this->include_node(label_start, label_end, count)) {
      this->status[node_index / 2] = Status::INCLUDED;
      return true;
    } else {
      this->status[node_index / 2] = Status::EXCLUDED;
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
    this->probabilities.resize(this->table.size() / 2);
    this->assign_node_probabilities(0);

    this->breadth_first_iteration(
        0, 0, false, [&](int node_index, int lcp, int edge_lcp) {
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
  void assign_node_probabilities(int node_index) {
    auto child_counts = this->get_child_counts(node_index, true);

    float child_sum = 0;
    for (int c : child_counts) {
      child_sum += c;
    }

    for (int i = 0; i < seqan3::alphabet_size<alphabet_t>; i++) {
      this->probabilities[node_index / 2][i] =
          float(child_counts[i]) / child_sum;
    }
  }

  /**! \brief Determines if pseudo counts should be added.
   * \details
   * Checks if any of the children that correspond to valid_characters are
   * missing.  If this is the case, 1 should be added to every child count,
   * to account for the fact that even if we don't observe an event
   * we would expect it to happen with some (small) probability.
   *
   * \param[in] node_index node to check for pseudo counts for.
   * \return true if pseudo counts should be added.
   */
  bool add_pseudo_counts(int node_index) {
    int n_children = 0;
    this->iterate_children(node_index, [&](int child_index) {
      int sequence_index = this->get_sequence_index(child_index);
      auto character = this->get_character(sequence_index);

      auto char_rank = seqan3::to_rank(character);
      if (this->valid_characters.find(char_rank) ==
          this->valid_characters.end()) {
        return;
      }

      n_children += 1;
    });
    return n_children != this->valid_characters.size();
  }

  /**! \brief Assigns node status to implicit nodes.
   * \details
   * Iterates through all nodes and assigns the status for implicit nodes
   * to the status of their parent.
   *
   */
  void add_implicit_node_status() {
    this->status.resize(this->table.size() / 2, Status::NONE);
    std::stack<std::tuple<int, Status>> stack{};

    this->iterate_children(0, [&](int child_index) {
      stack.emplace(child_index, Status::INCLUDED);
    });

    while (!stack.empty()) {
      auto [node_index, parent_status] = stack.top();
      stack.pop();

      Status node_status = this->status[node_index / 2];
      if ((node_status & Status::NONE) == Status::NONE) {
        this->status[node_index / 2] = Status(parent_status);
      }

      iterate_children(
          node_index, this->table, this->flags, [&](int child_index) {
            stack.emplace(child_index, this->status[node_index / 2]);
          });
    }
  }

  /**! \brief Removes all nodes from the tree with a delta value below
   * threshold.
   * \details
   * The full tree is iterated to find the leaves in the
   * PST.  These leaves are then iterated, for each leaf the delta value is
   * calculated, and the node is removed if the delta value is below a
   * threshold.  For each removed node, the parent is added to be considered if
   * it now a leaf.
   */
  virtual void cutoff_prune() { parameters_prune(); }
  virtual float calculate_delta(int node_index) { return 0.0; }

  /**! \brief Removes all nodes until a specified number of parameters are left.
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

    auto cmp = [](std::tuple<int, float> left, std::tuple<int, float> right) {
      return std::get<1>(left) < std::get<1>(right);
    };
    std::priority_queue<std::tuple<int, float>,
                        std::vector<std::tuple<int, float>>, decltype(cmp)>
        queue{cmp};

    for (auto v : pst_leaves) {
      queue.emplace(v, -this->calculate_delta(v));
    }

    auto n_terminal_nodes = this->count_terminal_nodes();

    while (!queue.empty() &&
           n_terminal_nodes * (this->valid_characters.size() - 1) >
               this->number_of_parameters) {
      auto [node_index, delta] = queue.top();
      queue.pop();

      if (node_index == 0) {
        continue;
      }

      int parent_index = this->suffix_links[node_index / 2];

      this->status[node_index / 2] = Status::EXCLUDED;

      if (this->is_pst_leaf(parent_index)) {
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
  std::array<int, seqan3::alphabet_size<alphabet_t>>
  get_child_counts(int node_index, bool with_pseudo_counts) {
    std::array<int, seqan3::alphabet_size<alphabet_t>> child_counts{};

    this->iterate_children(node_index, [&](int child_index) {
      int sequence_index = this->get_sequence_index(child_index);
      auto character = this->get_character(sequence_index);

      if (character == seqan3::gap{}) {
        return;
      }

      int character_rank = seqan3::to_rank(character);

      child_counts[character_rank] = get_counts(child_index);
    });

    bool add_pseudo_counts = this->add_pseudo_counts(node_index);
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
  std::vector<int> get_pst_leaves() {
    std::vector<int> bottom_nodes{};

    this->breadth_first_iteration(
        0, 0, false, [&](int node_index, int lcp, int edge_lcp) -> bool {
          if (this->is_included(node_index) && this->is_pst_leaf(node_index)) {
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
  bool is_pst_leaf(int node_index) {
    for (auto char_rank : this->valid_characters) {
      int child_index = this->reverse_suffix_links[node_index / 2][char_rank];

      if (child_index == 0 || child_index == -1) {
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
  bool is_terminal(int node_index) {
    for (auto char_rank : this->valid_characters) {
      int child_index = this->reverse_suffix_links[node_index / 2][char_rank];

      if (child_index == 0) {
        continue;
      }

      if (child_index == -1 || this->is_excluded(child_index)) {
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
  bool became_terminal(int node_index, int removed_index) {
    for (auto char_rank : this->valid_characters) {
      int child_index = this->reverse_suffix_links[node_index / 2][char_rank];

      if (child_index == 0) {
        continue;
      }

      if (child_index == -1 ||
          (this->is_excluded(child_index) && child_index != removed_index)) {
        return false;
      }
    }

    return is_terminal(node_index);
  }

  /**! \brief Returns count of the node
   * \details
   * Finds the count of the node by iterating the tree.  Saves the result in
   * a vector with the size of the tree / 2.
   *
   * \param[in] node_index
   * \return count of the node
   */
  int get_counts(int node_index) {
    if (node_index == -1) {
      return 0;
    }
    if (node_index == 0) {
      return this->sequence.size();
    }

    int c = this->counts[node_index / 2];
    if (c == -1) {
      c = lst::details::node_occurrences(node_index, this->table, this->flags);
      this->counts[node_index / 2] = c;
    }
    return c;
  }

  bool is_included(int node_index) {
    return (status[node_index / 2] & Status::INCLUDED) == Status::INCLUDED;
  }

  bool is_excluded(int node_index) {
    return (status[node_index / 2] & Status::EXCLUDED) == Status::EXCLUDED;
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
  void append_node_string(int node_index, int lcp, int edge_lcp,
                          std::ostringstream &tree_string) {
    auto label_ = this->node_label(node_index, lcp, edge_lcp);

    if (this->is_leaf(node_index)) {
      label_ = this->leaf_label(node_index, lcp);
    }
    std::string label =
        label_ | seqan3::views::to_char | seqan3::views::to<std::string>;
    if (node_index == 0) {
      label = "#";
    }

    tree_string << "Node: " << node_index / 2 << " " << label << " ";

    append_reverse_child_counts(node_index, tree_string);

    tree_string << " " << get_counts(node_index) << " ";

    append_child_counts(node_index, tree_string);

    tree_string << " ";

    append_reverse_children(node_index, tree_string);

    tree_string << std::endl;
  }

  /**! \brief Append the count of forward children to the output stream.
   * \details
   * The format is [ x y z ... ]
   *
   * \param[in] node_index The node index to operate on.
   * \param tree_string The output stream to write to.
   */
  void append_child_counts(int node_index, std::ostringstream &tree_string) {

    auto child_counts = this->get_child_counts(node_index, true);

    std::vector<int> output(seqan3::alphabet_size<alphabet_t>, -1);
    for (auto char_rank : this->valid_characters) {
      output[char_rank] = child_counts[char_rank];
    };

    tree_string << "[ ";

    for (auto count : output) {
      if (count == -1) {
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
  void append_reverse_child_counts(int node_index,
                                   std::ostringstream &tree_string) {

    std::vector<int> output(seqan3::alphabet_size<alphabet_t>, -1);

    auto reverse_children = this->reverse_suffix_links[node_index / 2];
    for (auto char_rank : this->valid_characters) {

      int reverse_child = reverse_children[char_rank];
      int count = get_counts(reverse_child);

      output[char_rank] = count;
    };

    tree_string << "[ ";
    for (auto count : output) {
      if (count == -1) {
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
  void append_reverse_children(int node_index,
                               std::ostringstream &tree_string) {

    std::vector<int> output(seqan3::alphabet_size<alphabet_t>, -2);

    auto reverse_children = this->reverse_suffix_links[node_index / 2];
    for (auto char_rank : this->valid_characters) {
      int reverse_child = reverse_children[char_rank];

      if (reverse_child == 0 || reverse_child == -1 ||
          this->is_excluded(reverse_child)) {
        output[char_rank] = -1;
      } else {
        output[char_rank] = reverse_child / 2;
      }
    };

    tree_string << "[ ";
    for (auto count : output) {
      if (count == -2) {
        continue;
      }
      tree_string << count << " ";
    }
    tree_string << "]";
  }

  /**! \brief Counts the nodes in the tree.
   *
   * \return The number of nodes in the tree.
   */
  int nodes_in_tree() {
    int n_nodes = 1;
    this->breadth_first_iteration(
        0, 0, false, [&](int node_index, int lcp, int edge_lcp) -> bool {
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
  bool include_node(int label_start, int label_end, int count) {
    int label_length = label_end - label_start;

    return label_length < this->max_depth && count >= this->freq &&
           label_valid(label_start, label_end);
  }

  /**! \brief Checks if a label is valid.
   *
   * \param[in] label_start starting index of the label in the sequence.
   * \param[in] label_end ending index of the label in the sequence.
   * \return true if the label is valid, false otherwise.
   */
  bool label_valid(int label_start, int label_end) {
    for (int i = label_start; i < label_end; i++) {
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
  bool skip_node(int node_index) {
    return this->is_unevaluated(node_index) || this->is_excluded(node_index);
  }
};
} // namespace pst
