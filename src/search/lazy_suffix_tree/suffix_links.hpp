#pragma once

#include <queue>
#include <stack>
#include <vector>

#include <seqan3/alphabet/concept.hpp>

#include "construction.hpp"
#include "iteration.hpp"

namespace lst::details {

/**! \brief Calculates the height/depth of the tree and of every node.
 *
 * @tparam alphabet_t Type of alphabet used (from seqan3)
 * \param[out] depths Depths of each node in the tree.
 * \param[in] sequence Sequence of the tree.
 * \param[in] suffixes Suffixes of the tree.
 * \param[in] table Table of the tree.
 * \param[in] flags Flags of the tree.
 * \return the heigh of the tree (longest node contained).
 */
template <seqan3::alphabet alphabet_t>
int tree_height(std::vector<int> &depths, sequence_t<alphabet_t> &sequence,
                std::vector<int> &suffixes, std::vector<int> &table,
                std::vector<Flag> &flags) {
  std::queue<std::tuple<int, int>> queue{};
  queue.emplace(0, 0);

  int tree_height = 0;

  while (!queue.empty()) {
    auto [node_index, parent_depth] = queue.front();
    queue.pop();

    int node_depth = parent_depth +
                     get_edge_lcp(node_index, sequence, suffixes, table, flags);

    depths[node_index / 2] = node_depth;

    iterate_children(node_index, table, flags,
                     [&](int index) { queue.emplace(index, node_depth); });

    if (node_depth > tree_height) {
      tree_height = node_depth;
    }
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
  std::vector<int> cause(suffixes.size() + 2, -1);

  std::vector<int> leaf_indices(suffixes.size() + 1, -1);
  leaf_indices[suffixes.size()] = 0;

  prepare_suffix_links(0, 0, cause, leaf_indices, sequence, suffixes, table,
                       flags);

  std::vector<int> depths(table.size() / 2, -1);
  int height = tree_height(depths, sequence, suffixes, table, flags);

  std::vector<int> branch(height + 1, -1);

  compute_suffix_links(cause, branch, depths, leaf_indices, sequence, suffixes,
                       table, flags, suffix_links);
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
 * \param[out] leaf_indices At each index the previous leaf is stored.
 * \param sequence Sequence of the tree
 * \param suffixes Suffixes of the tree.
 * \param table Table of the tree.
 * \param flags Flags of the tree.
 * \return The smallest index of the node.
 */
template <seqan3::alphabet alphabet_t>
int prepare_suffix_links(int node_index, int lcp, std::vector<int> &cause,
                         std::vector<int> &leaf_indices,
                         sequence_t<alphabet_t> &sequence,
                         std::vector<int> &suffixes, std::vector<int> &table,
                         std::vector<Flag> &flags) {
  if (is_leaf(node_index, flags)) {
    int leaf_index = get_leaf_index(node_index, lcp, suffixes, table, flags);

    leaf_indices[leaf_index] = node_index;

    return leaf_index;
  } else if (is_unevaluated(node_index, flags)) {
    // This isn't technically correct.  The following gives the second smallest
    // suffix _overall_ and the algorithm calls for the second smallest suffix
    // from the immediate children.  Not entirely sure that it is important
    // since the values will still be unique.

    int second_smallest_child = suffixes[table[node_index] + 1] - lcp;
    cause[second_smallest_child + 1] = node_index;

    return get_leaf_index(node_index, lcp, suffixes, table, flags);
  } else {
    int edge_lcp = get_edge_lcp(node_index, sequence, suffixes, table, flags);

    int smallest_child = suffixes.size();
    int second_smallest_child = suffixes.size();

    iterate_children(node_index, table, flags, [&](int index) {
      int child =
          prepare_suffix_links(index, lcp + edge_lcp, cause, leaf_indices,
                               sequence, suffixes, table, flags);

      if (child < smallest_child) {
        second_smallest_child = smallest_child;
        smallest_child = child;
      } else if (child < second_smallest_child) {
        second_smallest_child = child;
      }
    });

    cause[second_smallest_child + 1] = node_index;

    return smallest_child;
  }
}

void assign_link(int leaf_index, std::vector<int> &cause,
                 std::vector<int> &branch, std::vector<int> &depths,
                 std::vector<int> &suffix_links) {
  int caused = cause[leaf_index];
  if (caused != -1 && depths[caused / 2] != 0) {
    suffix_links[caused / 2] = branch[depths[caused / 2] - 1];
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
 * \param[in] cause At each index the node that was caused by that index is
 * saved.
 * \param[in] branch Contains the branching node at each depth.
 * \param[in] depths Depths of each node in the tree
 * \param[in] leaf_indices At each index the previous leaf is stored.
 * \param sequence Sequence of the tree
 * \param suffixes Suffixes of the tree.
 * \param table Table of the tree.
 * \param flags Flags of the tree.
 * \param[out] suffix_links The suffix link for each node.
 */
template <seqan3::alphabet alphabet_t>
void compute_suffix_links(std::vector<int> &cause, std::vector<int> &branch,
                          std::vector<int> &depths,
                          std::vector<int> &leaf_indices,
                          sequence_t<alphabet_t> &sequence,
                          std::vector<int> &suffixes, std::vector<int> &table,
                          std::vector<Flag> &flags,
                          std::vector<int> &suffix_links) {
  std::stack<std::tuple<int, int>> stack{};
  stack.emplace(0, 0);

  while (!stack.empty()) {
    auto [node_index, lcp] = stack.top();
    stack.pop();

    if (!is_leaf(node_index, flags)) {
      branch[depths[node_index / 2]] = node_index;
    }

    if (is_leaf(node_index, flags)) {
      int leaf_index = get_leaf_index(node_index, lcp, suffixes, table, flags);

      assign_link(leaf_index, cause, branch, depths, suffix_links);
      assign_leaf_link(node_index, leaf_index, suffix_links, leaf_indices);
    }

    if (is_unevaluated(node_index, flags)) {
      // Maybe possible to stop this iteration early.
      for (int i = table[node_index]; i < table[node_index + 1]; i++) {
        int leaf_index = suffixes[i] - lcp;
        assign_link(leaf_index, cause, branch, depths, suffix_links);
      }
    }

    int edge_lcp = get_edge_lcp(node_index, sequence, suffixes, table, flags);

    iterate_children(node_index, table, flags,
                     [&](int index) { stack.emplace(index, lcp + edge_lcp); });
  }
}

/**! \brief Computes the parents and the closest parent with a suffix link.
 *
 * \tparam alphabet_t seqan3 alphabet type.
 * \param[out] closest_suffix_link_destination Saves the closest parent with a
 * suffix link (and the distance in characters to that parent).
 * \param[out] parent_links Contains the parent for each node for a bottom-up
 * traversal.
 * \param[in] sequence Sequence of the tree
 * \param[in] suffixes Suffixes of the tree.
 * \param[in] table Table of the tree.
 * \param[in] flags Flags of the tree.
 * \param[in] suffix_links The suffix link for each node.
 */
template <seqan3::alphabet alphabet_t>
void prepare_implicit_suffix_links(
    std::vector<std::tuple<int, int>> &closest_suffix_link_destination,
    std::vector<int> &parent_links, sequence_t<alphabet_t> &sequence,
    std::vector<int> &suffixes, std::vector<int> &table,
    std::vector<Flag> &flags, std::vector<int> &suffix_links) {

  std::stack<std::tuple<int, int, int>> stack{};
  stack.emplace(0, 0, 0);

  std::vector<std::tuple<int, int>> missing_suffix_links{};
  while (!stack.empty()) {
    auto [node_index, parent, parent_lcp] = stack.top();
    stack.pop();

    int edge_lcp = get_edge_lcp(node_index, sequence, suffixes, table, flags);
    int lcp = parent_lcp + edge_lcp;

    parent_links[node_index / 2] = parent;

    if (!is_leaf(node_index, flags) && suffix_links[node_index / 2] == -1) {
      missing_suffix_links.emplace_back(node_index, lcp);

    } else if (missing_suffix_links.size() > 0) {
      int suffix_link_destination = suffix_links[node_index / 2];

      for (int i = 0; i < missing_suffix_links.size(); i++) {
        auto &[missing, missing_lcp] = missing_suffix_links[i];

        closest_suffix_link_destination[missing] =
            std::make_tuple(suffix_link_destination, lcp - missing_lcp);
      }

      missing_suffix_links.resize(0);
    }

    iterate_children(node_index, table, flags,
                     [&](int child) { stack.emplace(child, node_index, lcp); });
  }
}

/**! \brief Sets suffix links for implcit (expanded) nodes with the help
 * of the distance to the closest parent with a suffix link and the parent
 * links.
 *
 * \tparam alphabet_t seqan3 alphabet type.
 * \param[in] closest_suffix_link_destination Saves the closest parent with a
 * suffix link (and the distance in characters to that parent).
 * \param[in] parent_links Contains the parent for each node for a bottom-up
 * traversal.
 * \param[in] sequence Sequence of the tree
 * \param[in] suffixes Suffixes of the tree.
 * \param[in] table Table of the tree.
 * \param[in] flags Flags of the tree.
 * \param[out] suffix_links The suffix link for each node.
 */
template <seqan3::alphabet alphabet_t>
void compute_implicit_suffix_links(
    std::vector<std::tuple<int, int>> &closest_suffix_link_destination,
    std::vector<int> &parent_links, sequence_t<alphabet_t> &sequence,
    std::vector<int> &suffixes, std::vector<int> &table,
    std::vector<Flag> &flags, std::vector<int> &suffix_links) {
  std::queue<int> queue{};
  queue.push(0);

  while (!queue.empty()) {
    int node_index = queue.front();
    queue.pop();

    auto [suffix_link_destination, distance] =
        closest_suffix_link_destination[node_index];

    if (suffix_link_destination == -1) {
      suffix_links[node_index / 2] = -1;
    } else if (suffix_link_destination != 0) {
      int destination_parent = suffix_link_destination;
      int edge_lcp =
          get_edge_lcp(destination_parent, sequence, suffixes, table, flags);
      distance -= edge_lcp;
      while (distance > 0) {
        destination_parent = parent_links[destination_parent / 2];
        edge_lcp =
            get_edge_lcp(destination_parent, sequence, suffixes, table, flags);
        distance -= edge_lcp;
      }
      destination_parent = parent_links[destination_parent / 2];
      suffix_links[node_index / 2] = destination_parent;
    }

    iterate_children(node_index, table, flags,
                     [&](int index) { queue.push(index); });
  }
}

/**! \brief Adds suffix links for all (extended) implicit nodes in the tree.
 *
 * Prepare the two vectors by storing the distance to the closest
 * child which has a suffix link, as well as all parents links.
 *
 * When these vectors are filled, it allows us to assign the correct
 * suffix links of the implicit nodes by determining
 * which other node is at an equal distance to a child with a suffix link.
 * That node must therefore have a suffix link to the current node.
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
  std::vector<std::tuple<int, int>> closest_suffix_link_destination(
      table.size());

  std::vector<int> parent_links(table.size() / 2, -1);

  prepare_implicit_suffix_links(closest_suffix_link_destination, parent_links,
                                sequence, suffixes, table, flags, suffix_links);
  compute_implicit_suffix_links(closest_suffix_link_destination, parent_links,
                                sequence, suffixes, table, flags, suffix_links);
}
} // namespace lst::details
