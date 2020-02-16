#pragma once

#include <algorithm>
#include <array>
#include <functional>
#include <numeric>
#include <queue>
#include <stack>
#include <tuple>
#include <vector>

#include <seqan3/alphabet/composite/alphabet_variant.hpp>
#include <seqan3/alphabet/concept.hpp>
#include <seqan3/core/concept/tuple.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/range/container/bitcompressed_vector.hpp>
#include <seqan3/range/views/convert.hpp>
#include <seqan3/range/views/to.hpp>

#include "lazy_suffix_tree/construction.hpp"
#include "lazy_suffix_tree/iteration.hpp"
#include "lazy_suffix_tree/suffix_links.hpp"

namespace lst {

/**! \brief Lazy Suffix Tree implementation using the WOTD algorithm.
 *
 * \tparam alphabet_t Alphabet type from Seqan3.
 *
 * \details
 * The WOTD algorithm is described by Giegerich et al. in
 * https://doi.org/10.1002/spe.535. It differs from a suffix tree by only
 * building the parts of the tree that are used.  This massively saves time and
 * memory in cases where we're searching for only a subset of the suffixes.
 *
 */
template <seqan3::alphabet alphabet_t> class LazySuffixTree {

public:
  friend class LazySuffixTreeTest;

  LazySuffixTree() = default;
  LazySuffixTree(LazySuffixTree const &) = default;
  ~LazySuffixTree() = default;

  /**! \brief Constructor.
   * @param sequence_ The sequence to build the tree for.
   */
  LazySuffixTree(seqan3::bitcompressed_vector<alphabet_t> &sequence_)
      : sequence(sequence_) {
    suffixes = std::vector<int>(this->sequence.size() + 1);
    std::iota(this->suffixes.begin(), this->suffixes.end(), 0);

    lst::details::expand_root(this->sequence, this->suffixes, this->table,
                              this->flags);
  }

  /**! \brief Fully builds the lazy suffix tree.
   *
   */
  void expand_all() {
    for (int i = 0; i < this->table.size(); i++) {
      if (this->is_unevaluated(i)) {
        lst::details::expand_node(i, this->sequence, this->suffixes,
                                  this->table, this->flags);
      }
    }
  }

  /**! \brief Finds all labels encoded by the nodes in the tree.
   *
   * \return Vector of tuple of the label and count for each node in the tree.
   */
  std::vector<std::tuple<std::string, int>> get_all_labels() {
    std::vector<std::tuple<std::string, int>> labels{};

    this->breadth_first_iteration(
        0, 0, true, [&](int node_index, int lcp, int edge_lcp) -> bool {
          auto label = node_label(node_index, lcp, edge_lcp);
          int occurrences =
              lst::details::node_occurrences(node_index, table, flags);

          if (is_leaf(node_index)) {
            label = leaf_label(node_index, lcp);
          }
          labels.emplace_back(label, occurrences);
          return true;
        });

    return labels;
  }

  /**! \brief Find the start index of the pattern in the sequence.
   *
   * \param[in] pattern the sequence to search for.
   * \return Vector of indices into the sequence where the pattern is.
   */
  std::vector<int> search(std::vector<alphabet_t> pattern) {
    if (pattern.size() == 0) {
      std::vector<int> all_suffixes{};

      this->iterate_children(0, [&](int index) {
        auto suffixes = this->suffix_indices(index, 0);

        all_suffixes.insert(all_suffixes.end(), suffixes.begin(),
                            suffixes.end());
      });
      return all_suffixes;
    }

    auto [node_index, lcp] = this->find(pattern);

    if (node_index == -1) {
      return std::vector<int>{};
    } else {
      return this->suffix_indices(node_index, lcp);
    }
  }

  /**! \brief Breadth first traversal of the tree.
   * \details
   * Accepts a callback function which for each node gives node index, lcp and
   * edge lcp.
   *
   * \param f callback function which gives the start, end sequence index,
   * edge length and occurrences of each node.  Should return if further
   * iteration of the node is needed.
   */
  void breadth_first_iteration(const std::function<bool(int, int, int)> &f) {
    this->breadth_first_iteration(
        0, 0, false, [&](int node_index, int lcp, int edge_lcp) -> bool {
          if (this->skip_node(node_index)) {
            return true;
          } else {
            return f(node_index, lcp, edge_lcp);
          }
        });
  }

  /**! \brief Expands implicit nodes in the tree.
   * \details
   * I.e. for every edge with an edge length > 1, there is an extra node added.
   * This is a convenience for further implementations and not a requirement
   * for the lazy suffix tree.
   */
  void expand_implicit_nodes() {
    std::queue<int> queue{};
    queue.push(0);

    while (!queue.empty()) {
      int node_index = queue.front();
      queue.pop();

      if (this->skip_node(node_index)) {
        continue;
      }

      if (!is_leaf(node_index)) {
        this->iterate_children(node_index,
                               [&](int index) { queue.push(index); });
      }

      int edge_lcp = this->get_edge_lcp(node_index);
      if (edge_lcp > 1) {
        this->add_implicit_nodes(node_index, edge_lcp);
      }
    }
  }

  /**! \brief Expands the implicit nodes in node_index
   *
   * \param node_index The index of the node to expand implicit nodes of.
   * \param edge_lcp The length of the edge (== number of implicit nodes).
   */
  void add_implicit_nodes(int node_index, int edge_lcp) {
    if (edge_lcp > 1) {
      lst::details::add_implicit_nodes(node_index, edge_lcp, this->table,
                                       this->flags);
    }
  }

  /**! \brief Adds suffix links to the tree.
   *
   * \details
   * Note that this takes a significant amount of extra time and memory.
   *
   */
  void add_suffix_links() {
    suffix_links.resize(table.size() / 2, -1);
    std::fill(suffix_links.begin(), suffix_links.end(), -1);

    lst::details::add_explicit_suffix_links(sequence, suffixes, table, flags,
                                            suffix_links);
    lst::details::add_leaf_suffix_links(sequence, suffixes, table, flags,
                                        suffix_links);
    lst::details::add_implicit_suffix_links(sequence, suffixes, table, flags,
                                            suffix_links);
    suffix_links[0] = -1;
  }

  /**! \brief Adds reverse suffix links to the tree.
   * \details
   * This is done by iterating through the tree in a breadth first fashion
   * and requires (or will compute) the suffix links.
   */
  void add_reverse_suffix_links() {
    if (suffix_links.size() != table.size() / 2) {
      add_suffix_links();
    }

    reverse_suffix_links.resize(table.size() / 2);

    for (auto &reverses : reverse_suffix_links) {
      reverses.fill(-1);
    }

    this->breadth_first_iteration(
        0, 0, false, [&](int node_index, int lcp, int edge_lcp) -> bool {
          if (this->skip_node(node_index)) {
            return true;
          }

          int sequence_index = this->get_sequence_index(node_index);
          int node_start = sequence_index - lcp;

          int suffix_parent = suffix_links[node_index / 2];

          if (suffix_parent == -1) {
            return true;
          }

          auto character = this->get_character(node_start);
          int character_rank = seqan3::to_rank(character);

          if (reverse_suffix_links[suffix_parent / 2][character_rank] == -1) {
            reverse_suffix_links[suffix_parent / 2][character_rank] =
                node_index;
          }
          return true;
        });
  }

  /**! \brief Debug print of the tree.
   *
   */
  virtual void print() {
    this->debug_print_node(0, 0, 0);
    seqan3::debug_stream << std::endl;

    this->breadth_first_iteration(
        0, 0, false, [&](int node_index, int lcp, int edge_lcp) -> bool {
          this->debug_print_node(node_index, lcp, edge_lcp);
          seqan3::debug_stream << std::endl;
          return true;
        });
  }

  virtual void debug_print_node(int node_index, int lcp, int edge_lcp) {
    auto label = node_label(node_index, lcp, edge_lcp);

    if (this->is_leaf(node_index)) {
      label = this->leaf_label(node_index, lcp);
    }

    seqan3::debug_stream << label << "\t" << node_index << "\t"
                         << table[node_index] << "\t" << table[node_index + 1];

    seqan3::debug_stream << "\tLeaf: " << this->is_leaf(node_index);
    seqan3::debug_stream << "\tUnevaluated: "
                         << this->is_unevaluated(node_index);
    seqan3::debug_stream << "\tRightmost child: "
                         << this->is_rightmostchild(node_index);

    if (this->suffix_links.size() > node_index / 2) {
      seqan3::debug_stream << "\tSuffix link: "
                           << this->suffix_links[node_index / 2];
    }

    if (this->reverse_suffix_links.size() > node_index / 2) {
      seqan3::debug_stream << "\tReverse suffix links: "
                           << this->reverse_suffix_links[node_index / 2];
    }
  }

  std::string node_label(int node_index, int lcp, int edge_lcp) {
    int sequence_index = this->get_sequence_index(node_index);

    int node_start = sequence_index - lcp;
    int node_end =
        std::min(sequence_index + edge_lcp, int(this->sequence.size()));

    lst::details::sequence_t<alphabet_t> label(
        this->sequence.begin() + node_start, this->sequence.begin() + node_end);

    std::string label_str =
        label | seqan3::views::to_char | seqan3::views::to<std::string>;

    return label_str;
  }

  lst::details::sequence_t<alphabet_t> sequence;
  std::vector<int> suffixes{};
  std::vector<int> table{0, 2};
  std::vector<lst::details::Flag> flags{lst::details::Flag::RIGHT_MOST_CHILD,
                                        lst::details::Flag::NONE};
  std::vector<int> suffix_links{};
  std::vector<lst::details::alphabet_array<alphabet_t>> reverse_suffix_links{};

protected:
  std::vector<int> suffix_indices(int node_index, int og_lcp) {
    if (node_index >= table.size()) {
      throw std::invalid_argument(
          "[SUFFIX INDICES] Given node index is too large.");
    }

    std::vector<int> start_indices{};
    std::queue<std::tuple<int, int>> queue{};

    queue.emplace(node_index, og_lcp);

    while (!queue.empty()) {
      auto [index, lcp] = queue.front();
      queue.pop();

      if (is_leaf(index)) {
        start_indices.push_back(table[index] - lcp);
      } else if (is_unevaluated(index)) {
        for (int i = table[index]; i < table[index + 1]; i++) {
          start_indices.push_back(suffixes[i] - lcp);
        }
      } else {
        int edge_lcp = this->get_edge_lcp(index);
        int new_lcp = lcp + edge_lcp;
        this->iterate_children(index,
                               [&](int i) { queue.emplace(i, new_lcp); });
      }
    }

    return start_indices;
  }

  std::tuple<int, int> find(std::vector<alphabet_t> pattern) {
    std::queue<std::tuple<int, int>> queue{};

    queue.emplace(0, 0);

    int pattern_lcp = 0;

    while (!queue.empty()) {
      auto [node_index, lcp] = queue.front();
      queue.pop();

      int edge_lcp;
      if (is_unevaluated(node_index)) {
        edge_lcp =
            lst::details::expand_node(node_index, this->sequence,
                                      this->suffixes, this->table, this->flags);
      } else {
        edge_lcp = this->get_edge_lcp(node_index);
      }

      auto edge = edge_label(node_index, edge_lcp);

      bool edge_match = edge_matches(node_index, pattern_lcp, pattern, edge);

      if (!edge_match) {
        continue;
      }

      pattern_lcp += edge.size();
      if (pattern_lcp >= pattern.size()) {
        return std::make_tuple(node_index, lcp);
      }

      empty_queue(queue);

      if (!is_leaf(node_index)) {
        int new_lcp = lcp + edge_lcp;
        this->iterate_children(
            node_index, [&](int index) { queue.emplace(index, new_lcp); });
      }
    }

    return std::make_tuple(-1, -1);
  }

  template <typename T> void empty_queue(std::queue<T> &queue) {
    while (!queue.empty()) {
      queue.pop();
    }
  }

  bool edge_matches(int node_index, int pattern_lcp,
                    std::vector<alphabet_t> &pattern,
                    lst::details::sequence_t<alphabet_t> &edge) {
    for (int i = 0; pattern_lcp + i < pattern.size() && i < edge.size(); i++) {
      int sequence_index = table[node_index] + i;
      auto character = this->get_character(sequence_index);

      if (character != pattern[pattern_lcp + i]) {
        return false;
      }
    }

    return true;
  }

  lst::details::sequence_t<alphabet_t> edge_label(int node_index,
                                                  int edge_lcp) {
    int edge_start = this->get_sequence_index(node_index);
    int edge_end = std::min(edge_start + edge_lcp, int(this->sequence.size()));

    lst::details::sequence_t<alphabet_t> edge(
        this->sequence.begin() + edge_start, this->sequence.begin() + edge_end);

    return edge;
  }

  std::string leaf_label(int node_index, int lcp) {
    int node_start = table[node_index] - lcp;
    lst::details::sequence_t<alphabet_t> label(
        this->sequence.begin() + node_start, this->sequence.end());

    std::string label_str =
        label | seqan3::views::to_char | seqan3::views::to<std::string>;
    return label_str;
  }

  int get_sequence_index(int node_index) {
    return lst::details::get_sequence_index(node_index, this->suffixes,
                                            this->table, this->flags);
  }

  seqan3::gapped<alphabet_t> get_character(size_t index) {
    return lst::details::get_character(this->sequence, index);
  }

  bool is_leaf(int node_index) {
    return lst::details::is_leaf(node_index, flags);
  }

  bool is_unevaluated(int node_index) {
    return lst::details::is_unevaluated(node_index, flags);
  }

  bool is_rightmostchild(int node_index) {
    return lst::details::is_rightmostchild(node_index, flags);
  }

  /** \brief Determines if the node should be skipped during various iterations
   * Useful for subclassing, as further restrictions on which nodes are
   * relevant can be included.
   *
   * \param node_index Node to evaluate if it should be skipped.
   * \return Boolean, true if this node should be skipped.
   */
  virtual bool skip_node(int node_index) {
    return this->is_unevaluated(node_index);
  }

  /**! \brief Iterates through all suffix tree children of the node with
   * node_index.
   *
   * \param[in] node_index index of node to iterate children of.
   * \param f function to call on every child.
   */
  void iterate_children(int node_index, const std::function<void(int)> &f) {
    lst::details::iterate_children(node_index, this->table, this->flags, f);
  }

  /**! \brief Get the edge LCP of the node.
   *
   * \param node_index Index of the node.
   * \return The edge longest common prefix (length of the edge).
   */
  int get_edge_lcp(int node_index) {
    return lst::details::get_edge_lcp(node_index, this->sequence,
                                      this->suffixes, this->table, this->flags);
  }

  void breadth_first_iteration(int node_index, int start_lcp, bool expand_nodes,
                               const std::function<bool(int, int, int)> &f) {
    lst::details::breadth_first_iteration(node_index, start_lcp, this->sequence,
                                          this->suffixes, this->table,
                                          this->flags, expand_nodes, f);
  }
};
} // namespace lst
