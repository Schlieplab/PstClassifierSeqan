#pragma once

#include <queue>
#include <stack>
#include <vector>

#include <seqan3/alphabet/concept.hpp>

#include "construction.hpp"
#include "iteration.hpp"

namespace lst::details {

/**! \brief Get number of children for node.
 *
 * \param[in] node_index  Index to get number of children for.
 * \param[in] table Table of the tree.
 * \param[in] flags Flags of the tree.
 * \return
 */
int get_number_of_children(int node_index, std::vector<int> &table,
                           std::vector<Flag> &flags) {
  int n_children = 0;
  iterate_children(node_index, table, flags, [&](int index) { n_children++; });

  return n_children;
}

/**! \brief Calculates the height/depth of the tree and of every node.
 *
 * @tparam alphabet_t Type of alphabet used (from seqan3)
 * \param[in] sequence Sequence of the tree.
 * \param[in] suffixes Suffixes of the tree.
 * \param[in] table Table of the tree.
 * \param[in] flags Flags of the tree.
 * \return the height of the tree (longest node contained).
 */
template <seqan3::alphabet alphabet_t>
int tree_height(sequence_t<alphabet_t> &sequence, std::vector<int> &suffixes,
                std::vector<int> &table, std::vector<Flag> &flags) {
  std::queue<std::tuple<int, int>> queue{};
  queue.emplace(0, 0);

  int tree_height = 0;

  while (!queue.empty()) {
    auto [node_index, parent_depth] = queue.front();
    queue.pop();

    if (is_leaf(node_index, flags)) {
      continue;
    }

    int edge_lcp = get_edge_lcp(node_index, sequence, suffixes, table, flags);
    int node_depth = parent_depth + edge_lcp;

    iterate_children(node_index, table, flags,
                     [&](int index) { queue.emplace(index, node_depth); });

    tree_height = std::max(tree_height, node_depth);
  }

  return tree_height;
}

/**! \brief Gets the "leaf" index of the node.
 * This is the index the label of the leaf starts at in the sequence.
 *
 * \param node_index Index of the node.
 * \param lcp Longest common prefix of the node.
 * \param suffixes Suffixes of the tree.
 * \param table Table of the tree.
 * \param flags Flags of the tree.
 * \return The index at which the leaf starts in the sequence or -1 if not leaf.
 */
int get_leaf_index(int node_index, int lcp, std::vector<int> &suffixes,
                   std::vector<int> &table, std::vector<Flag> &flags) {
  if (is_leaf(node_index, flags)) {
    return table[node_index] - lcp;
  } else if (is_unevaluated(node_index, flags)) {
    return suffixes[table[node_index]] - lcp;
  } else {
    return -1;
  }
}

/**! \brief Adds suffix links to all explicit nodes in the tree.
 *
 * This follows the algorithm described by Maaß in
 *  doi:10.1016/j.ipl.2005.12.012
 *  With modifications to allow for unevaluated nodes.
 *
 *  Takes O(sequence) memory and time.
 *
 * \tparam alphabet_t seqan3 alphabet for the tree.
 * \param sequence Sequence of the tree
 * \param suffixes Suffixes of the tree.
 * \param table Table of the tree.
 * \param flags Flags of the tree.
 * \param[out] suffix_links Suffix links of each explicit node in the tree.
 */
template <seqan3::alphabet alphabet_t>
void add_explicit_suffix_links(sequence_t<alphabet_t> &sequence,
                               std::vector<int> &suffixes,
                               std::vector<int> &table,
                               std::vector<Flag> &flags,
                               std::vector<int> &suffix_links) {
  std::vector<std::tuple<int, int>> cause(suffixes.size());

  prepare_suffix_links(0, 0, cause, sequence, suffixes, table, flags);

  int height = tree_height(sequence, suffixes, table, flags);

  std::vector<int> branch(height + 1, -1);

  compute_suffix_links(cause, branch, sequence, suffixes, table, flags,
                       suffix_links);
}

/**! \brief For each node in the tree, computes which node caused that node.
 *
 * This is the prepare step of the algorithm by Maaß.
 *
 * \tparam alphabet_t seqan3 alphabet type.
 * \param node_index Index of the node.
 * \param lcp Longest common prefix of the node.
 * \param[out] cause At each index the node that was caused by that index is
 * saved.
 * \param sequence Sequence of the tree
 * \param suffixes Suffixes of the tree.
 * \param table Table of the tree.
 * \param flags Flags of the tree.
 * \return The smallest index of the node.
 */
template <seqan3::alphabet alphabet_t>
int prepare_suffix_links(int node_index, int lcp,
                         std::vector<std::tuple<int, int>> &cause,
                         sequence_t<alphabet_t> &sequence,
                         std::vector<int> &suffixes, std::vector<int> &table,
                         std::vector<Flag> &flags) {
  if (is_leaf(node_index, flags) || is_unevaluated(node_index, flags)) {
    return get_leaf_index(node_index, lcp, suffixes, table, flags);
  } else {
    int edge_lcp = get_edge_lcp(node_index, sequence, suffixes, table, flags);

    int smallest_child = suffixes.size();
    int second_smallest_child = suffixes.size();

    iterate_children(node_index, table, flags, [&](int index) {
      int child = prepare_suffix_links(index, lcp + edge_lcp, cause, sequence,
                                       suffixes, table, flags);

      if (child < smallest_child) {
        second_smallest_child = smallest_child;
        smallest_child = child;
      } else if (child < second_smallest_child) {
        second_smallest_child = child;
      }
    });

    int n_children = get_number_of_children(node_index, table, flags);
    // If the node only has one child it is an implicit node
    // and we don't want to add it as a causing node.
    if (n_children == 1) {
      return smallest_child;
    }

    if (second_smallest_child < suffixes.size()) {
      cause[second_smallest_child + 1] =
          std::make_tuple(node_index, lcp + edge_lcp);
    }

    return smallest_child;
  }
}

/**! \brief Assign suffix link corresponding to the leaf index.
 *
 * \param leaf_index Leaf index to set suffix link for.
 * \param[in] cause At each index the node that was caused by that index is
 * saved.
 * \param[in] branch Contains the branching node at each depth.
 * \param[in] table Table of the tree.
 * \param[in] flags Flags of the tree.
 * \param[out] suffix_links The suffix link for each node.
 */
void assign_link(int leaf_index, std::vector<std::tuple<int, int>> &cause,
                 std::vector<int> &branch, std::vector<int> &table,
                 std::vector<Flag> &flags, std::vector<int> &suffix_links) {
  auto &[caused, depth] = cause[leaf_index];

  if (caused != -1 && depth != -1) {
    suffix_links[caused / 2] = branch[depth - 1];
  }
}

/**! \brief Adds suffix links to all explicit nodes.
 *
 * \tparam alphabet_t seqan3 alphabet type.
 * \param[in] cause At each index the node that was caused by that index is
 * saved.
 * \param[in] branch Contains the branching node at each depth.
 * \param sequence Sequence of the tree
 * \param suffixes Suffixes of the tree.
 * \param table Table of the tree.
 * \param flags Flags of the tree.
 * \param[out] suffix_links The suffix link for each node.
 */
template <seqan3::alphabet alphabet_t>
void compute_suffix_links(std::vector<std::tuple<int, int>> &cause,
                          std::vector<int> &branch,
                          sequence_t<alphabet_t> &sequence,
                          std::vector<int> &suffixes, std::vector<int> &table,
                          std::vector<Flag> &flags,
                          std::vector<int> &suffix_links) {
  std::stack<std::tuple<int, int>> stack{};
  stack.emplace(0, 0);

  while (!stack.empty()) {
    auto [node_index, lcp] = stack.top();
    stack.pop();

    if (is_leaf(node_index, flags)) {
      int leaf_index = get_leaf_index(node_index, lcp, suffixes, table, flags);

      assign_link(leaf_index, cause, branch, table, flags, suffix_links);
    } else if (is_unevaluated(node_index, flags)) {
      int edge_lcp = get_edge_lcp(node_index, sequence, suffixes, table, flags);

      int height = lcp + edge_lcp;
      branch[height] = node_index;
      std::fill(branch.begin() + height + 1, branch.end(), -1);

      for (auto i = table[node_index]; i < table[node_index + 1]; i++) {
        auto leaf_index = suffixes[i] - lcp;
        assign_link(leaf_index, cause, branch, table, flags, suffix_links);
      }
    } else {
      int edge_lcp = get_edge_lcp(node_index, sequence, suffixes, table, flags);
      int height = lcp + edge_lcp;

      int n_children = get_number_of_children(node_index, table, flags);
      // If the node only has one child it is an implicit node
      // and we don't want to add it as a branching point.
      if (n_children != 1) {
        branch[height] = node_index;
      }

      iterate_children(node_index, table, flags,
                       [&](int index) { stack.emplace(index, height); });
    }
  }
}

/**! \brief Adds suffix links to all leaves nodes in the tree.
 *
 * These are easy to compute as they follow directly from the next index
 * of the leaf in the sequence.
 *  Takes O(sequence) memory and time.
 *
 * \tparam alphabet_t seqan3 alphabet for the tree.
 * \param sequence Sequence of the tree
 * \param suffixes Suffixes of the tree.
 * \param table Table of the tree.
 * \param flags Flags of the tree.
 * \param[out] suffix_links Suffix links of each explicit node in the tree.
 */
template <seqan3::alphabet alphabet_t>
void add_leaf_suffix_links(sequence_t<alphabet_t> &sequence,
                           std::vector<int> &suffixes, std::vector<int> &table,
                           std::vector<Flag> &flags,
                           std::vector<int> &suffix_links) {

  std::vector<int> leaf_indices(suffixes.size() + 1, -1);
  leaf_indices[suffixes.size()] = 0;

  prepare_leaf_suffix_links(leaf_indices, sequence, suffixes, table, flags);

  compute_leaf_suffix_links(leaf_indices, sequence, suffixes, table, flags,
                            suffix_links);
}

/**! \brief For each node in the tree, computes which node caused that node.
 *
 * This is the prepare step of the algorithm by Maaß.
 *
 * \tparam alphabet_t seqan3 alphabet type.
 * \param node_index Index of the node.
 * \param lcp Longest common prefix of the node.
 * \param[out] leaf_indices At each index the previous leaf is stored.
 * \param sequence Sequence of the tree
 * \param suffixes Suffixes of the tree.
 * \param table Table of the tree.
 * \param flags Flags of the tree.
 * \return The smallest index of the node.
 */
template <seqan3::alphabet alphabet_t>
void prepare_leaf_suffix_links(std::vector<int> &leaf_indices,
                               sequence_t<alphabet_t> &sequence,
                               std::vector<int> &suffixes,
                               std::vector<int> &table,
                               std::vector<Flag> &flags) {
  std::stack<std::tuple<int, int>> stack{};
  stack.emplace(0, 0);

  while (!stack.empty()) {
    auto [node_index, lcp] = stack.top();
    stack.pop();

    if (is_leaf(node_index, flags)) {
      int leaf_index = get_leaf_index(node_index, lcp, suffixes, table, flags);

      leaf_indices[leaf_index] = node_index;
    }
    int edge_lcp = get_edge_lcp(node_index, sequence, suffixes, table, flags);
    iterate_children(node_index, table, flags,
                     [&](int index) { stack.emplace(index, lcp + edge_lcp); });
  }
}

void assign_leaf_link(int node_index, int leaf_index,
                      std::vector<int> &suffix_links,
                      std::vector<int> &leaf_indices) {
  if (leaf_indices[leaf_index + 1] != -1) {
    suffix_links[node_index / 2] = leaf_indices[leaf_index + 1];
  }
}

/**! \brief Adds suffix links to all explicit nodes.
 *
 * \tparam alphabet_t seqan3 alphabet type.
 * \param[in] leaf_indices At each index the previous leaf is stored.
 * \param sequence Sequence of the tree
 * \param suffixes Suffixes of the tree.
 * \param table Table of the tree.
 * \param flags Flags of the tree.
 * \param[out] suffix_links The suffix link for each node.
 */
template <seqan3::alphabet alphabet_t>
void compute_leaf_suffix_links(std::vector<int> &leaf_indices,
                               sequence_t<alphabet_t> &sequence,
                               std::vector<int> &suffixes,
                               std::vector<int> &table,
                               std::vector<Flag> &flags,
                               std::vector<int> &suffix_links) {
  std::stack<std::tuple<int, int>> stack{};
  stack.emplace(0, 0);

  while (!stack.empty()) {
    auto [node_index, lcp] = stack.top();
    stack.pop();

    if (is_leaf(node_index, flags)) {
      int leaf_index = get_leaf_index(node_index, lcp, suffixes, table, flags);

      assign_leaf_link(node_index, leaf_index, suffix_links, leaf_indices);
    }

    int edge_lcp = get_edge_lcp(node_index, sequence, suffixes, table, flags);

    iterate_children(node_index, table, flags,
                     [&](int index) { stack.emplace(index, lcp + edge_lcp); });
  }
}

/** \brief Checks if the sequences corresponding to node_index and
 * suffix_link_child_index matches.
 *
 * We only have to check the last edge_lcp characters of the nodes because of
 * the suffix link between the parent of node_index and a parent of the
 * suffix_link_child_index.
 *
 * \tparam alphabet_t seqan3 alphabet for the tree.
 * \param node_index The index of the node to find the suffix link for.
 * \param edge_lcp The edge length of node_index.
 * \param suffix_link_child_index The potential suffix link of node_index.
 * \param suffix_link_edge_lcp edge lcp of suffix_link_child_index
 * \param sequence Sequence of the tree
 * \param suffixes Suffixes of the tree.
 * \param table Table of the tree.
 * \param flags Flags of the tree.
 * \return True if the two sequences matches.
 */
template <seqan3::alphabet alphabet_t>
bool sequences_match(int node_index, int edge_lcp, int suffix_link_child_index,
                     int suffix_link_edge_lcp, sequence_t<alphabet_t> &sequence,
                     std::vector<int> &suffixes, std::vector<int> &table,
                     std::vector<Flag> &flags) {
  int node_start = get_sequence_index(node_index, suffixes, table, flags);
  int suffix_link_child_end = std::min(
      int(sequence.size()),
      get_sequence_index(suffix_link_child_index, suffixes, table, flags) +
          suffix_link_edge_lcp);

  int suffix_link_child_start = suffix_link_child_end - edge_lcp;

  for (int i = 0; i < edge_lcp; i++) {
    if (sequence[node_start + i] != sequence[suffix_link_child_start + i]) {
      return false;
    }
  }

  return true;
}

/** \brief Returns the correct suffix link for node_index if possible.
 *
 * Iterates through all of the children of parent_suffix_link up to a maximum
 * depth of edge_lcp.  For each child with we check if it is long enough and
 * if it matches the sequence of node_index.  If so, return as the suffix link
 * destination of node_index.  Otherwise, return -1.
 *
 * \tparam alphabet_t seqan3 alphabet for the tree.
 * \param node_index The index of the node to find the suffix link for.
 * \param edge_lcp The edge length of node_index.
 * \param parent_suffix_link The suffix link of the parent of node_index.
 * \param sequence Sequence of the tree
 * \param suffixes Suffixes of the tree.
 * \param table Table of the tree.
 * \param flags Flags of the tree.
 * \return suffix link destination of node_index, or -1 if none found.
 */
template <seqan3::alphabet alphabet_t>
int find_suffix_match(int node_index, int edge_lcp, int parent_suffix_link,
                      sequence_t<alphabet_t> &sequence,
                      std::vector<int> &suffixes, std::vector<int> &table,
                      std::vector<Flag> &flags) {
  std::queue<std::tuple<int, int>> suffix_link_queue{};
  suffix_link_queue.emplace(parent_suffix_link, 0);

  while (!suffix_link_queue.empty()) {
    auto &[suffix_link_child, suffix_link_lcp] = suffix_link_queue.front();
    suffix_link_queue.pop();
    auto suffix_link_edge_lcp =
        get_edge_lcp(suffix_link_child, sequence, suffixes, table, flags);

    if (suffix_link_lcp == edge_lcp) {
      bool match = sequences_match(node_index, edge_lcp, suffix_link_child,
                                   suffix_link_edge_lcp, sequence, suffixes,
                                   table, flags);
      if (match) {
        return suffix_link_child;
      }
    } else if (suffix_link_lcp < edge_lcp) {

      iterate_children(suffix_link_child, table, flags, [&](int index) {
        suffix_link_queue.emplace(index,
                                  suffix_link_lcp + suffix_link_edge_lcp);
      });
    }
  }
  return -1;
}

/**! \brief Adds suffix links for all (extended) implicit nodes in the tree.
 *
 * Finds nodes which have not been assigned suffix links.  If the parent
 * has a suffix link, we can iterate from the parent's suffix link to find
 * the suffix link of the child, if it exists.
 *
 * \tparam alphabet_t seqan3 alphabet for the tree.
 * \param sequence Sequence of the tree
 * \param suffixes Suffixes of the tree.
 * \param table Table of the tree.
 * \param flags Flags of the tree.
 * \param[out] suffix_links Suffix links of each explicit node in the tree.
 */
template <seqan3::alphabet alphabet_t>
void add_implicit_suffix_links(sequence_t<alphabet_t> &sequence,
                               std::vector<int> &suffixes,
                               std::vector<int> &table,
                               std::vector<Flag> &flags,
                               std::vector<int> &suffix_links) {

  std::queue<std::tuple<int, int>> queue{};
  queue.emplace(0, 0);

  while (!queue.empty()) {
    auto [node_index, parent_index] = queue.front();
    queue.pop();

    auto edge_lcp = get_edge_lcp(node_index, sequence, suffixes, table, flags);

    if (suffix_links[node_index / 2] == -1 &&
        suffix_links[parent_index / 2] != -1) {
      auto parent_suffix_link = suffix_links[parent_index / 2];
      int suffix_link_destination =
          find_suffix_match(node_index, edge_lcp, parent_suffix_link, sequence,
                            suffixes, table, flags);

      suffix_links[node_index / 2] = suffix_link_destination;
    }

    iterate_children(node_index, table, flags,
                     [&](int index) { queue.emplace(index, node_index); });
  }
}

} // namespace lst::details
