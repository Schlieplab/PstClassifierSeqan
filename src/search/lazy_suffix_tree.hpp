#pragma once

#include <algorithm>
#include <array>
#include <functional>
#include <numeric>
#include <queue>
#include <stack>
#include <tuple>
#include <typeinfo>
#include <vector>

#include <seqan3/alphabet/composite/alphabet_variant.hpp>
#include <seqan3/alphabet/concept.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/core/concept/tuple.hpp>
#include <seqan3/range/views/convert.hpp>
#include <seqan3/range/views/to.hpp>
#include <seqan3/range/views/to_char.hpp>

#include "lazy_suffix_tree/construction.hpp"
#include "lazy_suffix_tree/iteration.hpp"
#include "lazy_suffix_tree/suffix_links.hpp"

size_t max_size = (size_t)-1;

namespace lst {

template <typename alph> std::string get_alphabet_name() {
  return typeid(alph).name();
}
template <> std::string get_alphabet_name<seqan3::dna5>() { return "DNA5"; }
template <> std::string get_alphabet_name<seqan3::dna4>() { return "DNA4"; }

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
  friend class ProbabilisticSuffixTreeTest;

  LazySuffixTree() = default;
  LazySuffixTree(LazySuffixTree const &) = default;
  ~LazySuffixTree() = default;

  /**! \brief Constructor.
   * @param sequence_ The sequence to build the tree for.
   * \param[in] multi_core True for parallel execution.
   * \param[in] parallel_depth The maximum depth to spawn new processes, will
   * control task size as `alphabet_size ** depth`.
   */
  LazySuffixTree(lst::details::sequence_t<alphabet_t> &sequence_,
                 bool multi_core_ = true, int parallel_depth = 2)
      : sequence(sequence_), multi_core(multi_core_),
        parallel_depth(parallel_depth) {

    suffixes = std::vector<size_t>(this->sequence.size() + 1);

    lst::details::expand_root(this->sequence, this->suffixes, this->table);
  }

  /**! \brief Fully builds the lazy suffix tree.
   *
   */
  void expand_all() {
    for (size_t i = 0; i < this->table.size(); i++) {
      if (this->is_unevaluated(i)) {
        lst::details::expand_node(i, this->sequence, this->suffixes,
                                  this->table, [](size_t n, size_t &e) {});
      }
    }
  }

  /**! \brief Finds all labels encoded by the nodes in the tree.
   *
   * \return Vector of tuple of the label and count for each node in the tree.
   */
  std::vector<std::tuple<std::string, size_t>> get_all_labels() {
    std::vector<std::tuple<std::string, size_t>> labels{};
    static std::mutex labels_mutex{};

    this->breadth_first_iteration(
        0, 0, true,
        [&](size_t node_index, size_t lcp, size_t edge_lcp, size_t occurrences,
            lst::details::alphabet_array<size_t, alphabet_t> &child_counts)
            -> bool {
          occurrences = node_occurrences(node_index);
          auto label = node_label(node_index, lcp, edge_lcp);

          std::lock_guard labels_lock{labels_mutex};
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
  std::vector<size_t> search(std::vector<alphabet_t> pattern) {
    if (pattern.size() == 0) {
      std::vector<size_t> all_suffixes{};

      this->iterate_children(0, [&](size_t index) {
        auto suffixes_ = this->suffix_indices(index, 0);

        all_suffixes.insert(all_suffixes.end(), suffixes_.begin(),
                            suffixes_.end());
      });
      return all_suffixes;
    }

    auto [node_index, lcp] = this->find(pattern);

    if (node_index == -1) {
      return std::vector<size_t>{};
    } else {
      return this->suffix_indices(node_index, lcp);
    }
  }

  /**! \brief Breadth first traversal of the tree.
   * \details
   * Accepts a callback function which for each node gives node index, lcp,
   * edge lcp, and node count.
   *
   * Does not expand the tree.
   *
   * \param f callback function which gives the start, end sequence index,
   * edge length and occurrences of each node.  Should return if further
   * iteration of the node is needed.
   */
  void breadth_first_iteration(
      const std::function<
          bool(size_t, size_t, size_t, size_t,
               lst::details::alphabet_array<size_t, alphabet_t> &)> &f) {
    this->breadth_first_iteration(
        0, 0, false,
        [&](size_t node_index, size_t lcp, size_t edge_lcp, size_t node_count,
            lst::details::alphabet_array<size_t, alphabet_t> &child_counts)
            -> bool {
          if (this->skip_node(node_index)) {
            return true;
          } else {
            return f(node_index, lcp, edge_lcp, node_count, child_counts);
          }
        });
  }

  /**! \brief Breadth first traversal of the tree.
   * \details
   * Accepts a callback function which for each node gives node index, lcp,
   * edge lcp, and node count.
   *
   * Does not expand the tree.
   *
   * \param f callback function which gives the start, end sequence index,
   * edge length and occurrences of each node.  Should return if further
   * iteration of the node is needed.
   * \param done callback function to signal that no more nodes can be iterated.
   * Useful for synchronisation of accumulated results in e.g. threads.
   * \param[in] locked_callback Callback for when the tree array has grown and a
   * user wants to do something in a synchronisation lock.
   */
  void breadth_first_iteration(
      const std::function<
          bool(size_t, size_t, size_t, size_t,
               lst::details::alphabet_array<size_t, alphabet_t> &)> &f,
      const std::function<void()> &done,
      const std::function<void()> &locked_callback) {
    this->breadth_first_iteration(
        0, 0, false,
        [&](size_t node_index, size_t lcp, size_t edge_lcp, size_t node_count,
            lst::details::alphabet_array<size_t, alphabet_t> &child_counts)
            -> bool {
          if (this->skip_node(node_index)) {
            return true;
          } else {
            return f(node_index, lcp, edge_lcp, node_count);
          }
        },
        done, locked_callback);
  }

  /**! \brief Sequential breadth first traversal of the tree.
   * \details
   * Accepts a callback function which for each node gives node index, lcp,
   * edge lcp, and node count.
   *
   * Does not expand the tree.
   *
   * \param f callback function which gives the start, end sequence index,
   * edge length and occurrences of each node.  Should return if further
   * iteration of the node is needed.
   */
  void breadth_first_iteration_sequential(
      const std::function<
          bool(size_t, size_t, size_t, size_t,
               lst::details::alphabet_array<size_t, alphabet_t> &)> &f) {
    this->breadth_first_iteration_sequential(
        0, 0, false,
        [&](size_t node_index, size_t lcp, size_t edge_lcp, size_t node_count,
            lst::details::alphabet_array<size_t, alphabet_t> &child_counts)
            -> bool {
          if (this->skip_node(node_index)) {
            return true;
          } else {
            return f(node_index, lcp, edge_lcp, node_count, child_counts);
          }
        });
  }

  /**! \brief Breadth first traversal of the tree, without constructing the
   * table (which contains the tree structure).
   *
   * \details Accepts a callback function which for each node gives node
   * index, lcp and edge lcp.
   *
   * \param f callback function which gives the start, end sequence index,
   * edge length and occurrences of each node.  Should return if further
   * iteration of the node is needed.
   * \param done callback function to signal that no more nodes can be iterated.
   */
  void breadth_first_iteration_table_less(
      const std::function<
          bool(size_t, size_t, size_t, size_t,
               lst::details::alphabet_array<size_t, alphabet_t> &, bool)> &f,
      const std::function<void()> &done) {
    lst::details::breadth_first_iteration_table_less(
        this->parallel_depth, this->table, this->sequence, this->suffixes, f, done);
  }

  /**! \brief Expands implicit nodes in the tree.
   * \details
   * I.e. for every edge with an edge length > 1, there is an extra node added.
   * This is a convenience for further implementations and not a requirement
   * for the lazy suffix tree.
   */
  void expand_implicit_nodes() {
    std::queue<size_t> queue{};
    queue.push(0);

    while (!queue.empty()) {
      auto node_index = queue.front();
      queue.pop();

      if (this->skip_node(node_index)) {
        continue;
      }

      if (!is_leaf(node_index)) {
        this->iterate_children(node_index,
                               [&](size_t index) { queue.push(index); });
      }

      auto edge_lcp = this->get_edge_lcp(node_index);
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
  void add_implicit_nodes(size_t node_index, size_t edge_lcp) {
    if (edge_lcp > 1) {
      lst::details::add_implicit_nodes(node_index, edge_lcp, this->table);
    }
  }

  /**! \brief Adds suffix links to the tree.
   *
   * \details
   * Note that this takes a significant amount of extra time and memory.
   *
   */
  void add_suffix_links() {
    suffix_links.resize(table.size(), (size_t)-1);
    std::fill(suffix_links.begin(), suffix_links.end(), -1);

    lst::details::add_explicit_suffix_links(sequence, suffixes, table,
                                            suffix_links, this->multi_core,
                                            this->parallel_depth);

    lst::details::add_leaf_suffix_links(sequence, suffixes, table, suffix_links,
                                        this->multi_core, this->parallel_depth);

    lst::details::add_implicit_suffix_links(sequence, suffixes, table,
                                            suffix_links, this->multi_core,
                                            this->parallel_depth);

    suffix_links[0] = -1;
  }

  /**! \brief Adds reverse suffix links to the tree.
   * \details
   * This is done by iterating through the tree in a breadth first fashion
   * and requires (or will compute) the suffix links.
   */
  void add_reverse_suffix_links() {
    if (suffix_links.size() != table.size()) {
      add_suffix_links();
    }

    reverse_suffix_links.resize(table.size());

    for (auto &reverses : reverse_suffix_links) {
      reverses.fill(-1);
    }

    this->breadth_first_iteration(
        0, 0, false,
        [&](size_t node_index, size_t lcp, size_t edge_lcp, size_t node_count,
            lst::details::alphabet_array<size_t, alphabet_t> &child_counts)
            -> bool {
          if (this->skip_node(node_index)) {
            return true;
          }

          auto sequence_index = this->get_sequence_index(node_index);
          auto node_start = sequence_index - lcp;

          auto suffix_parent = suffix_links[node_index];

          if (suffix_parent == -1) {
            return true;
          }

          auto character_rank = this->get_character_rank(node_start);

          if (reverse_suffix_links[suffix_parent][character_rank] == -1) {
            reverse_suffix_links[suffix_parent][character_rank] =
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
    std::cout << std::endl;

    this->breadth_first_iteration_sequential(
        0, 0, false,
        [&](size_t node_index, size_t lcp, size_t edge_lcp, size_t node_count,
            lst::details::alphabet_array<size_t, alphabet_t> &child_counts)
            -> bool {
          this->debug_print_node(node_index, lcp, edge_lcp);
          std::cout << std::endl;
          return true;
        });
  }

  virtual void debug_print_node(size_t node_index, size_t lcp,
                                size_t edge_lcp) {
    auto label = node_label(node_index, lcp, edge_lcp);

    if (this->is_leaf(node_index)) {
      label = this->leaf_label(node_index, lcp);
    }

    std::cout << label << "\t" << node_index << "\t" << table[node_index].first
              << "\t" << table[node_index].second;

    std::cout << "\tLeaf: " << this->is_leaf(node_index);
    std::cout << "\tUnevaluated: " << this->is_unevaluated(node_index);

    std::cout << "\tRightmost child: " << this->is_rightmostchild(node_index);

    //    if (this->suffix_links.size() > node_index) {
    //      std::cout << "\tSuffix link: " << this->suffix_links[node_index /
    //      2];
    //    }

    //    if (this->reverse_suffix_links.size() > node_index) {
    //      std::cout << "\tReverse suffix links: " <<
    //      this->reverse_suffix_links[node_index];
    //    }
  }

  std::string node_label(size_t node_index, size_t lcp, size_t edge_lcp) {
    if (this->is_leaf(node_index)) {
      return this->leaf_label(node_index, lcp);
    }
    auto sequence_index = this->get_sequence_index(node_index);

    auto node_start = sequence_index - lcp;
    auto node_end = std::min(sequence_index + edge_lcp, this->sequence.size());

    lst::details::sequence_t<alphabet_t> label(
        this->sequence.begin() + node_start, this->sequence.begin() + node_end);

    std::string label_str =
        label | seqan3::views::to_char | seqan3::views::to<std::string>;

    return label_str;
  }

  size_t node_occurrences(size_t node_index) {
    return lst::details::node_occurrences(node_index, this->table,
                                          this->sequence, this->suffixes);
  }

  lst::details::sequence_t<alphabet_t> sequence;
  std::vector<size_t> suffixes{};
  lst::details::Table<> table{{0, lst::details::Flag::RIGHT_MOST_CHILD},
                              {2, lst::details::Flag::NONE}};
  std::vector<size_t> suffix_links{};
  std::vector<lst::details::alphabet_array<size_t, alphabet_t>>
      reverse_suffix_links{};

  bool multi_core;
  int parallel_depth;

protected:
  std::vector<size_t> suffix_indices(size_t node_index, size_t og_lcp) {
    if (node_index >= table.size()) {
      throw std::invalid_argument(
          "[SUFFIX INDICES] Given node index is too large.");
    }

    std::vector<size_t> start_indices{};
    std::queue<std::tuple<size_t, size_t>> queue{};

    queue.emplace(node_index, og_lcp);

    while (!queue.empty()) {
      auto [index, lcp] = queue.front();
      queue.pop();

      if (is_leaf(index)) {
        start_indices.push_back(table[index].first - lcp);
      } else if (is_unevaluated(index)) {
        for (auto i = table[index].first; i < table[index].second; i++) {
          start_indices.push_back(suffixes[i] - lcp);
        }
      } else {
        auto edge_lcp = this->get_edge_lcp(index);
        auto new_lcp = lcp + edge_lcp;
        this->iterate_children(index,
                               [&](size_t i) { queue.emplace(i, new_lcp); });
      }
    }

    return start_indices;
  }

  std::tuple<size_t, size_t> find(std::vector<alphabet_t> pattern) {
    std::queue<std::tuple<size_t, size_t>> queue{};

    queue.emplace(0, 0);

    size_t pattern_lcp = 0;

    while (!queue.empty()) {
      auto [node_index, lcp] = queue.front();
      queue.pop();

      size_t edge_lcp;
      if (is_unevaluated(node_index)) {
        auto [edge_lcp_, _1, _2] = lst::details::expand_node(
            node_index, this->sequence, this->suffixes, this->table,
            [](size_t n, size_t &e) {});
        edge_lcp = edge_lcp_;
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
        size_t new_lcp = lcp + edge_lcp;
        this->iterate_children(
            node_index, [&](size_t index) { queue.emplace(index, new_lcp); });
      }
    }

    return std::make_tuple(-1, -1);
  }

  template <typename T> void empty_queue(std::queue<T> &queue) {
    while (!queue.empty()) {
      queue.pop();
    }
  }

  bool edge_matches(size_t node_index, size_t pattern_lcp,
                    std::vector<alphabet_t> &pattern,
                    lst::details::sequence_t<alphabet_t> &edge) {
    for (size_t i = 0; pattern_lcp + i < pattern.size() && i < edge.size();
         i++) {
      auto sequence_index = table[node_index].first + i;
      auto character = this->get_character(sequence_index);

      if (character != pattern[pattern_lcp + i]) {
        return false;
      }
    }

    return true;
  }

  lst::details::sequence_t<alphabet_t> edge_label(size_t node_index,
                                                  size_t edge_lcp) {
    auto edge_start = this->get_sequence_index(node_index);
    auto edge_end = std::min(edge_start + edge_lcp, this->sequence.size());

    lst::details::sequence_t<alphabet_t> edge(
        this->sequence.begin() + edge_start, this->sequence.begin() + edge_end);

    return edge;
  }

  std::string leaf_label(size_t node_index, size_t lcp) {
    auto node_start = table[node_index].first - lcp;
    lst::details::sequence_t<alphabet_t> label(
        this->sequence.begin() + node_start, this->sequence.end());

    std::string label_str =
        label | seqan3::views::to_char | seqan3::views::to<std::string>;
    return label_str;
  }

  size_t get_sequence_index(size_t node_index) {
    return lst::details::get_sequence_index(node_index, this->suffixes,
                                            this->table);
  }

  seqan3::gapped<alphabet_t> get_character(size_t index) {
    return lst::details::get_character(this->sequence, index);
  }

  size_t get_character_rank(size_t index) {
    return lst::details::get_character_rank(this->sequence, index);
  }

  bool is_leaf(size_t node_index) {
    return lst::details::is_leaf(node_index, this->table);
  }

  bool is_unevaluated(size_t node_index) {
    return lst::details::is_unevaluated(node_index, this->table);
  }

  bool is_rightmostchild(size_t node_index) {
    return lst::details::is_rightmostchild(node_index, this->table);
  }

  /** \brief Determines if the node should be skipped during various iterations
   * Useful for subclassing, as further restrictions on which nodes are
   * relevant can be included.
   *
   * \param node_index Node to evaluate if it should be skipped.
   * \return Boolean, true if this node should be skipped.
   */
  virtual bool skip_node(size_t node_index) {
    return this->is_unevaluated(node_index);
  }

  /**! \brief Iterates through all suffix tree children of the node with
   * node_index.
   *
   * \param[in] node_index index of node to iterate children of.
   * \param f function to call on every child.
   */
  void iterate_children(size_t node_index,
                        const std::function<void(size_t)> &f) {
    lst::details::iterate_children(node_index, this->table, f);
  }

  /**! \brief Get the edge LCP of the node.
   *
   * \param node_index Index of the node.
   * \return The edge longest common prefix (length of the edge).
   */
  size_t get_edge_lcp(size_t node_index) {
    return lst::details::get_edge_lcp(node_index, this->sequence,
                                      this->suffixes, this->table);
  }

  void breadth_first_iteration(
      size_t node_index, size_t start_lcp, bool expand_nodes,
      const std::function<
          bool(size_t, size_t, size_t &, size_t,
               lst::details::alphabet_array<size_t, alphabet_t> &)> &f) {
    breadth_first_iteration(
        node_index, start_lcp, expand_nodes, f, []() {},
        [](size_t n, size_t l, size_t &e) {});
  }

  void breadth_first_iteration(
      size_t node_index, size_t start_lcp, bool expand_nodes,
      const std::function<
          bool(size_t, size_t, size_t &, size_t,
               lst::details::alphabet_array<size_t, alphabet_t> &)> &f,
      const std::function<void()> &done,
      const std::function<void(size_t, size_t, size_t &)> &locked_callback) {
    if (this->multi_core) {
      lst::details::breadth_first_iteration_parallel(
          node_index, start_lcp, this->sequence, this->suffixes, this->table,
          expand_nodes, f, done, this->parallel_depth, locked_callback);
    } else {
      lst::details::breadth_first_iteration(
          node_index, start_lcp, this->sequence, this->suffixes, this->table,
          expand_nodes, f, locked_callback);
      done();
    }
  }

  void breadth_first_iteration_sequential(
      size_t node_index, size_t start_lcp, bool expand_nodes,
      const std::function<
          bool(size_t, size_t, size_t &, size_t,
               lst::details::alphabet_array<size_t, alphabet_t> &)> &f) {
    lst::details::breadth_first_iteration(
        node_index, start_lcp, this->sequence, this->suffixes, this->table,
        expand_nodes, f, [](size_t n, size_t l, size_t &e) {});
  }
};
} // namespace lst
