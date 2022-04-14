#pragma once

#include <algorithm>
#include <chrono>
// #include <execution>
#include <future>
#include <mutex>
#include <queue>
#include <stack>
#include <thread>
#include <vector>

#include "seqan3/alphabet/concept.hpp"

#include "construction.hpp"
#include "iteration.hpp"

namespace lst::details {

thread_local size_t max_val = 0;

/**! \brief Get number of children for node.
 *
 * \param[in] node_index  Index to get number of children for.
 * \param[in] table Table of the tree.
 * \return
 */
size_t get_number_of_children(size_t node_index, const Table<> &table) {
  size_t n_children = 0;
  iterate_children(node_index, table, [&](size_t index) { n_children++; });

  return n_children;
}

/**! \brief Calculates the internal height/depth of the tree and of every node.
 * @tparam alphabet_t Type of alphabet used (from seqan3)
 * \param[in] sequence Sequence of the tree.
 * \param[in] suffixes Suffixes of the tree.
 * \param[in] table Table of the tree.
 * \return the height of the tree (longest node contained, not leaf).
 */
template <seqan3::alphabet alphabet_t>
size_t tree_height(const sequence_t<alphabet_t> &sequence,
                   const std::vector<size_t> &suffixes, const Table<> &table) {

  std::queue<std::tuple<size_t, size_t>> queue{};
  queue.emplace(0, 0);

  size_t tree_height = 0;

  while (!queue.empty()) {
    auto [node_index, parent_depth] = queue.front();
    queue.pop();

    if (is_leaf(node_index, table)) {
      continue;
    }

    auto edge_lcp = get_edge_lcp(node_index, sequence, suffixes, table);
    auto node_depth = parent_depth + edge_lcp;

    tree_height = std::max(tree_height, node_depth);

    iterate_children(node_index, table,
                     [&](size_t index) { queue.emplace(index, node_depth); });
  }
  return tree_height;
}

/**! \brief Calculates the height/depth of the tree and of every node.
 *
 * @tparam alphabet_t Type of alphabet used (from seqan3)
 * \param[in] sequence Sequence of the tree.
 * \param[in] suffixes Suffixes of the tree.
 * \param[in] table Table of the tree.
 * \return the height of the tree (longest node contained, not leaf).
 */
template <seqan3::alphabet alphabet_t>
size_t tree_height_parallel(const sequence_t<alphabet_t> &sequence,
                            const std::vector<size_t> &suffixes,
                            const Table<> &table, int parallel_depth) {
  size_t tree_height_found = 0;
  max_val = 0;

  std::mutex height_lock;

  breadth_first_iteration_parallel<alphabet_t>(
      sequence, const_cast<std::vector<size_t> &>(suffixes),
      const_cast<Table<> &>(table), false,
      [&](size_t node_index, size_t lcp, size_t edge_lcp, size_t node_count,
          lst::details::alphabet_array<size_t, alphabet_t> &child_counts)
          -> bool {
        if (is_leaf(node_index, table)) {
          return false;
        }

        max_val = std::max(max_val, lcp + edge_lcp);
        return true;
      },
      [&]() {
        std::lock_guard lock{height_lock};
        tree_height_found = std::max(max_val, tree_height_found);
      },
      parallel_depth, [](size_t n, size_t l, size_t &e) {});

  return tree_height_found;
}

/**! \brief Gets the "leaf" index of the node.
 * This is the index the label of the leaf starts at in the sequence.
 *
 * \param node_index Index of the node.
 * \param lcp Longest common prefix of the node.
 * \param suffixes Suffixes of the tree.
 * \param table Table of the tree.
 * \return The index at which the leaf starts in the sequence or -1 if not leaf.
 */
size_t get_leaf_index(size_t node_index, size_t lcp,
                      const std::vector<size_t> &suffixes,
                      const Table<> &table) {
  if (is_leaf(node_index, table)) {
    return table[node_index].first - lcp;
  } else if (is_unevaluated(node_index, table)) {
    return suffixes[table[node_index].first] - lcp;
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
 * \param[out] suffix_links Suffix links of each explicit node in the tree.
 */
template <seqan3::alphabet alphabet_t>
void add_explicit_suffix_links(const sequence_t<alphabet_t> &sequence,
                               const std::vector<size_t> &suffixes,
                               const Table<> &table,
                               std::vector<size_t> &suffix_links,
                               bool multi_core, int parallel_depth) {
  std::vector<std::tuple<size_t, unsigned short>> cause(suffixes.size(),
                                                        {-1, -1});

  if (multi_core) {
    prepare_suffix_links(0, 0, cause, sequence, suffixes, table,
                         parallel_depth);

    auto height =
        tree_height_parallel(sequence, suffixes, table, parallel_depth);

    std::vector<size_t> branch(height + 1, -1);

    compute_suffix_links_parallel(0, 0, 0, cause, branch, sequence, suffixes,
                                  table, suffix_links, parallel_depth);
  } else {
    prepare_suffix_links(0, 0, cause, sequence, suffixes, table, 0);

    auto height = tree_height(sequence, suffixes, table);

    std::vector<size_t> branch(height + 1, -1);

    compute_suffix_links(cause, branch, sequence, suffixes, table,
                         suffix_links);
  }
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
 * \return The smallest index of the node.
 */
template <seqan3::alphabet alphabet_t>
size_t
prepare_suffix_links(size_t node_index, size_t lcp,
                     std::vector<std::tuple<size_t, unsigned short>> &cause,
                     const sequence_t<alphabet_t> &sequence,
                     const std::vector<size_t> &suffixes, const Table<> &table,
                     int parallel_depth) {

  if (is_leaf(node_index, table) || is_unevaluated(node_index, table)) {
    return get_leaf_index(node_index, lcp, suffixes, table);
  } else {
    auto edge_lcp = get_edge_lcp(node_index, sequence, suffixes, table);

    std::vector<size_t> child_values{};

    if (lcp + edge_lcp < parallel_depth) {
      std::vector<std::future<size_t>> child_futures{};

      iterate_children(node_index, table, [&](size_t index) {
        std::future<size_t> child = std::async(
            std::launch::async, prepare_suffix_links<alphabet_t>, index,
            lcp + edge_lcp, std::ref(cause), std::ref(sequence),
            std::ref(suffixes), std::ref(table), parallel_depth);
        child_futures.push_back(std::move(child));
      });

      for (auto &f : child_futures) {
        child_values.push_back(f.get());
      }
    } else {
      iterate_children(node_index, table, [&](size_t index) {
        auto child =
            prepare_suffix_links(index, lcp + edge_lcp, cause, sequence,
                                 suffixes, table, parallel_depth);
        child_values.push_back(child);
      });
    }

    auto smallest_child = suffixes.size();
    auto second_smallest_child = suffixes.size();

    for (auto &child : child_values) {
      if (child < smallest_child) {
        second_smallest_child = smallest_child;
        smallest_child = child;
      } else if (child < second_smallest_child) {
        second_smallest_child = child;
      }
    }

    auto n_children = child_values.size();
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
 * \param[out] suffix_links The suffix link for each node.
 */
inline void
assign_link(size_t leaf_index,
            const std::vector<std::tuple<size_t, unsigned short>> &cause,
            const std::vector<size_t> &branch,
            std::vector<size_t> &suffix_links) {
  auto &[caused, depth] = cause[leaf_index];

  if (caused != std::numeric_limits<size_t>::max() &&
      depth != std::numeric_limits<unsigned short>::max()) {
    suffix_links[caused] = branch[depth - 1];
  }
}

void assign_leaf_link(size_t node_index, size_t leaf_index,
                      std::vector<size_t> &suffix_links,
                      const std::vector<size_t> &leaf_indices) {
  if (leaf_indices[leaf_index + 1] != -1) {
    suffix_links[node_index] = leaf_indices[leaf_index + 1];
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
 * \param[out] suffix_links The suffix link for each node.
 */
template <seqan3::alphabet alphabet_t>
void compute_suffix_links(
    std::vector<std::tuple<size_t, unsigned short>> &cause,
    std::vector<size_t> &branch, const sequence_t<alphabet_t> &sequence,
    const std::vector<size_t> &suffixes, const Table<> &table,
    std::vector<size_t> &suffix_links) {

  std::stack<std::tuple<size_t, size_t>> stack{};
  stack.emplace(0, 0);

  while (!stack.empty()) {
    auto [node_index, lcp] = stack.top();
    stack.pop();

    auto height =
        compute_single_suffix_link(node_index, lcp, cause, branch, sequence,
                                   suffixes, table, suffix_links);

    iterate_children(node_index, table,
                     [&](size_t index) { stack.emplace(index, height); });
  }
}

/**! \brief Adds suffix links to all explicit nodes in parallel.
 *
 * \tparam alphabet_t seqan3 alphabet type.
 * \param[in] node_index The node index to start execution at.
 * \param[in] lcp The lcp of the node of node_index.
 * \param[in] depth Depth of the node in the tree (to control parallelisation).
 * \param[in] cause At each index the node that was caused by that index is
 * saved.
 * \param[in] branch Contains the branching node at each depth.
 * \param sequence Sequence of the tree
 * \param suffixes Suffixes of the tree.
 * \param table Table of the tree.
 * \param[out] suffix_links The suffix link for each node.
 * \param[in] parallel_depth The depth to stop spawning threads at.
 */
template <seqan3::alphabet alphabet_t>
void compute_suffix_links_parallel(
    size_t node_index, size_t lcp, size_t depth,
    const std::vector<std::tuple<size_t, unsigned short>> &cause,
    std::vector<size_t> &branch, const sequence_t<alphabet_t> &sequence,
    const std::vector<size_t> &suffixes, const Table<> &table,
    std::vector<size_t> &suffix_links, int parallel_depth) {

  std::vector<std::thread> threads{};
  std::vector<std::vector<size_t>> branch_copies(
      seqan3::alphabet_size<alphabet_t> + 1, std::vector<size_t>());
  size_t branch_copy_idx = 0;

  std::stack<std::tuple<size_t, size_t, size_t>> stack{};
  stack.emplace(node_index, lcp, depth);

  while (!stack.empty()) {
    auto [node_index, lcp, depth] = stack.top();
    stack.pop();

    auto height =
        compute_single_suffix_link(node_index, lcp, cause, branch, sequence,
                                   suffixes, table, suffix_links);

    if (depth < parallel_depth) {
      iterate_children(node_index, table, [&, depth = depth](size_t index) {
        branch_copies[branch_copy_idx] = std::vector<size_t>(branch);

        threads.push_back(std::thread{
            compute_suffix_links_parallel<alphabet_t>, index, height, depth + 1,
            std::ref(cause), std::ref(branch_copies[branch_copy_idx]),
            std::ref(sequence), std::ref(suffixes), std::ref(table),
            std::ref(suffix_links), parallel_depth});
        branch_copy_idx++;
      });

    } else {
      iterate_children(node_index, table, [&, depth = depth](size_t index) {
        stack.emplace(index, height, depth + 1);
      });
    }
  }

  for (auto &thread : threads) {
    if (thread.joinable()) {
      thread.join();
    }
  }
}

template <seqan3::alphabet alphabet_t>
size_t compute_single_suffix_link(
    size_t node_index, size_t lcp,
    const std::vector<std::tuple<size_t, unsigned short>> &cause,
    std::vector<size_t> &branch, const sequence_t<alphabet_t> &sequence,
    const std::vector<size_t> &suffixes, const Table<> &table,
    std::vector<size_t> &suffix_links) {
  auto edge_lcp = get_edge_lcp(node_index, sequence, suffixes, table);
  auto height = lcp + edge_lcp;

  if (is_leaf(node_index, table)) {
    auto leaf_index = get_leaf_index(node_index, lcp, suffixes, table);

    assign_link(leaf_index, cause, branch, suffix_links);

    return 0;
  } else if (is_unevaluated(node_index, table)) {
    branch[height] = node_index;

    // For a correct, general, lazy suffix tree.
    //    deal_with_implicit_nodes(node_index, lcp, cause, branch, sequence,
    //    suffixes, table, suffix_links);
  } else {
    auto n_children = get_number_of_children(node_index, table);
    // If the node only has one child it is an implicit node
    // and we don't want to add it as a branching point.
    if (n_children != 1) {
      branch[height] = node_index;
    }
  }

  return height;
}

template <seqan3::alphabet alphabet_t>
void deal_with_implicit_nodes(
    size_t node_index, size_t lcp,
    const std::vector<std::tuple<size_t, unsigned short>> &cause,
    std::vector<size_t> &branch, const sequence_t<alphabet_t> &sequence,
    const std::vector<size_t> &suffixes, const Table<> &table,
    std::vector<size_t> &suffix_links) {
  for (auto i = table[node_index].first; i < table[node_index].second; i++) {
    auto leaf_index = suffixes[i] - lcp;
    assign_link(leaf_index, cause, branch, suffix_links);
  }
}

/**! \brief For each node in the tree, computes which node caused that node.
 *
 * \param[out] leaf_indices At each index the previous leaf is stored.
 * \param sequence Sequence of the tree
 * \param suffixes Suffixes of the tree.
 * \param table Table of the tree.
 * \param leaves Leaf and lcp of leaf in the tree.
 */
void prepare_leaf_suffix_links(
    std::vector<size_t> &leaf_indices, const std::vector<size_t> &suffixes,
    const Table<> &table,
    const std::vector<std::tuple<size_t, size_t>> &leaves) {

  std::for_each(
      leaves.begin(), leaves.end(), [&](std::tuple<size_t, size_t> leaf) {
        auto &[node_index, lcp] = leaf;
        auto leaf_index = get_leaf_index(node_index, lcp, suffixes, table);

        leaf_indices[leaf_index] = node_index;
      });
}

/**! \brief For each node in the tree, computes which node caused that node.
 *
 *
 * \param[out] leaf_indices At each index the previous leaf is stored.
 * \param sequence Sequence of the tree
 * \param suffixes Suffixes of the tree.
 * \param table Table of the tree.
 * \param leaves Leaf and lcp of leaf in the tree.
 */
void prepare_leaf_suffix_links_parallel(
    std::vector<size_t> &leaf_indices, const std::vector<size_t> &suffixes,
    const Table<> &table,
    const std::vector<std::tuple<size_t, size_t>> &leaves) {

  std::for_each(
      leaves.begin(), leaves.end(), [&](std::tuple<size_t, size_t> leaf) {
        auto &[node_index, lcp] = leaf;
        auto leaf_index = get_leaf_index(node_index, lcp, suffixes, table);

        leaf_indices[leaf_index] = node_index;
      });
}

/**! \brief Adds suffix links to all explicit nodes.
 *
 * \param[in] leaf_indices At each index the previous leaf is stored.
 * \param sequence Sequence of the tree
 * \param suffixes Suffixes of the tree.
 * \param table Table of the tree.
 * \param leaves Leaf and lcp of leaf in the tree.
 * \param[out] suffix_links The suffix link for each node.
 */
void compute_leaf_suffix_links(
    const std::vector<size_t> &leaf_indices,
    const std::vector<size_t> &suffixes, const Table<> &table,
    const std::vector<std::tuple<size_t, size_t>> &leaves,
    std::vector<size_t> &suffix_links) {

  std::for_each(
      leaves.begin(), leaves.end(), [&](std::tuple<size_t, size_t> leaf) {
        auto &[node_index, lcp] = leaf;
        auto leaf_index = get_leaf_index(node_index, lcp, suffixes, table);

        assign_leaf_link(node_index, leaf_index, suffix_links, leaf_indices);
      });
}

/**! \brief Adds suffix links to all explicit nodes in parallel.
 *
 * \param[in] leaf_indices At each index the previous leaf is stored.
 * \param sequence Sequence of the tree
 * \param suffixes Suffixes of the tree.
 * \param table Table of the tree.
 * \param leaves Leaf and lcp of leaf in the tree
 * \param[out] suffix_links The suffix link for each node.
 */
void compute_leaf_suffix_links_parallel(
    std::vector<size_t> &leaf_indices, const std::vector<size_t> &suffixes,
    const Table<> &table, const std::vector<std::tuple<size_t, size_t>> &leaves,
    std::vector<size_t> &suffix_links) {

  std::for_each(
      leaves.begin(), leaves.end(), [&](std::tuple<size_t, size_t> leaf) {
        auto &[node_index, lcp] = leaf;
        auto leaf_index = get_leaf_index(node_index, lcp, suffixes, table);

        assign_leaf_link(node_index, leaf_index, suffix_links, leaf_indices);
      });
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
 * \param[out] suffix_links Suffix links of each explicit node in the tree.
 * \param[in] parallel_depth The depth to stop spawning threads at.
 */
template <seqan3::alphabet alphabet_t>
void add_leaf_suffix_links(const sequence_t<alphabet_t> &sequence,
                           const std::vector<size_t> &suffixes,
                           const Table<> &table,
                           std::vector<size_t> &suffix_links, bool multi_core,
                           int parallel_depth) {

  std::vector<size_t> leaf_indices(suffixes.size() + 1, -1);
  leaf_indices[suffixes.size()] = 0;

  auto leaves =
      get_leaves(sequence, suffixes, table, multi_core, parallel_depth);

  if (multi_core) {
    prepare_leaf_suffix_links_parallel(leaf_indices, suffixes, table, leaves);

    compute_leaf_suffix_links_parallel(leaf_indices, suffixes, table, leaves,
                                       suffix_links);

  } else {
    prepare_leaf_suffix_links(leaf_indices, suffixes, table, leaves);

    compute_leaf_suffix_links(leaf_indices, suffixes, table, leaves,
                              suffix_links);
  }
}

template <seqan3::alphabet alphabet_t>
std::vector<std::tuple<size_t, size_t>>
get_leaves(const sequence_t<alphabet_t> &sequence,
           const std::vector<size_t> &suffixes, const Table<> &table,
           bool multi_core, int parallel_depth) {
  std::vector<std::tuple<size_t, size_t>> leaves{};

  if (multi_core) {
    std::mutex leaves_mutex{};
    thread_local std::vector<std::tuple<size_t, size_t>> local_leaves{};

    breadth_first_iteration_parallel<alphabet_t>(
        sequence, const_cast<std::vector<size_t> &>(suffixes),
        const_cast<Table<> &>(table), false,
        [&](size_t node_index, size_t lcp, size_t edge_lcp, size_t node_count,
            lst::details::alphabet_array<size_t, alphabet_t> &child_counts)
            -> bool {
          if (is_leaf(node_index, table)) {
            local_leaves.emplace_back(node_index, lcp);
          }

          return true;
        },
        [&]() {
          std::lock_guard<std::mutex> lock(leaves_mutex);
          std::move(local_leaves.begin(), local_leaves.end(),
                    std::back_inserter(leaves));
        },
        parallel_depth, [](size_t n, size_t l, size_t &e) {});

  } else {
    std::queue<std::tuple<size_t, size_t>> queue{};
    queue.emplace(0, 0);
    while (!queue.empty()) {
      auto [node_index, lcp] = queue.front();
      queue.pop();

      if (is_leaf(node_index, table)) {
        leaves.emplace_back(node_index, lcp);
        continue;
      }

      auto edge_lcp = get_edge_lcp(node_index, sequence, suffixes, table);

      auto new_lcp = lcp + edge_lcp;

      iterate_children(node_index, table,
                       [&](size_t index) { queue.emplace(index, new_lcp); });
    }
  }
  return leaves;
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
 * \return True if the two sequences matches.
 */
template <seqan3::alphabet alphabet_t>
bool sequences_match(size_t node_index, size_t edge_lcp,
                     size_t suffix_link_child_index,
                     size_t suffix_link_edge_lcp,
                     const sequence_t<alphabet_t> &sequence,
                     const std::vector<size_t> &suffixes,
                     const Table<> &table) {
  auto node_start = get_sequence_index(node_index, suffixes, table);
  size_t suffix_link_child_end =
      std::min(sequence.size(),
               get_sequence_index(suffix_link_child_index, suffixes, table) +
                   suffix_link_edge_lcp);

  size_t suffix_link_child_start = suffix_link_child_end - edge_lcp;

  for (size_t i = 0; i < edge_lcp; i++) {
    if (sequence[node_start + i] != sequence[suffix_link_child_start + i]) {
      return false;
    }
  }

  return true;
}

/** \brief Returns the correct suffix link for node_index if possible.
 *
 * Iterates through all of the children of parent_suffix_link up to a maximum
 * depth of edge_lcp.  For each child we check if it is long enough and
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
 * \return suffix link destination of node_index, or -1 if none found.
 */
template <seqan3::alphabet alphabet_t>
size_t
find_suffix_match(size_t node_index, size_t edge_lcp, size_t parent_suffix_link,
                  const sequence_t<alphabet_t> &sequence,
                  const std::vector<size_t> &suffixes, const Table<> &table) {
  std::queue<std::tuple<size_t, size_t>> suffix_link_queue{};
  iterate_children(parent_suffix_link, table,
                   [&](size_t index) { suffix_link_queue.emplace(index, 0); });

  while (!suffix_link_queue.empty()) {
    auto [suffix_link_child, suffix_link_parent_lcp] =
        suffix_link_queue.front();
    suffix_link_queue.pop();

    // If we're dealing with implicit nodes, all edges are 1 long.
    //    size_t suffix_link_edge_lcp = get_edge_lcp(suffix_link_child,
    //    sequence, suffixes, table);
    size_t suffix_link_edge_lcp = 1;

    size_t suffix_link_lcp = suffix_link_parent_lcp + suffix_link_edge_lcp;

    if (suffix_link_lcp == edge_lcp) {
      bool match =
          sequences_match(node_index, edge_lcp, suffix_link_child,
                          suffix_link_edge_lcp, sequence, suffixes, table);
      if (match) {
        return suffix_link_child;
      }
    } else if (suffix_link_lcp < edge_lcp) {
      iterate_children(suffix_link_child, table, [&](size_t index) {
        suffix_link_queue.emplace(index, suffix_link_lcp);
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
 * \param[out] suffix_links Suffix links of each explicit node in the tree.
 */
template <seqan3::alphabet alphabet_t>
void add_implicit_suffix_links(const sequence_t<alphabet_t> &sequence,
                               const std::vector<size_t> &suffixes,
                               const Table<> &table,
                               std::vector<size_t> &suffix_links,
                               bool multi_core, int parallel_depth) {
  if (multi_core) {
    add_implicit_suffix_links_parallel(0, 0, 0, sequence, suffixes, table,
                                       suffix_links, parallel_depth);
  } else {
    add_implicit_suffix_links_sequential(sequence, suffixes, table,
                                         suffix_links);
  }
}

/**! \brief Adds suffix links for all (extended) implicit nodes in the tree
 * sequentially.
 *
 * Finds nodes which have not been assigned suffix links.  If the parent
 * has a suffix link, we can iterate from the parent's suffix link to find
 * the suffix link of the child, if it exists.
 *
 * \tparam alphabet_t seqan3 alphabet for the tree.
 * \param sequence Sequence of the tree
 * \param suffixes Suffixes of the tree.
 * \param table Table of the tree.
 * \param[out] suffix_links Suffix links of each explicit node in the tree.
 */
template <seqan3::alphabet alphabet_t>
void add_implicit_suffix_links_sequential(
    const sequence_t<alphabet_t> &sequence, const std::vector<size_t> &suffixes,
    const Table<> &table, std::vector<size_t> &suffix_links) {

  std::queue<std::tuple<size_t, size_t>> queue{};
  queue.emplace(0, 0);

  while (!queue.empty()) {
    auto [node_index, parent_index] = queue.front();
    queue.pop();

    visit_implicit_suffix_links_node(node_index, parent_index, sequence,
                                     suffixes, table, suffix_links);

    iterate_children(node_index, table,
                     [&, node_index = node_index](size_t index) {
                       queue.emplace(index, node_index);
                     });
  }
}

template <seqan3::alphabet alphabet_t>
void visit_implicit_suffix_links_node(size_t node_index, size_t parent_index,
                                      const sequence_t<alphabet_t> &sequence,
                                      const std::vector<size_t> &suffixes,
                                      const Table<> &table,
                                      std::vector<size_t> &suffix_links) {
  size_t max_size = (size_t)-1;
  // If we're dealing with implicit nodes, all edges are 1 long.
  //  size_t edge_lcp = get_edge_lcp(node_index, sequence, suffixes, table);
  size_t edge_lcp = 1;

  if (suffix_links[node_index] == max_size && parent_index == 0) {
    suffix_links[node_index] = 0;
  } else if (suffix_links[node_index] == max_size &&
             suffix_links[parent_index] != max_size) {
    auto parent_suffix_link = suffix_links[parent_index];

    auto suffix_link_destination = find_suffix_match(
        node_index, edge_lcp, parent_suffix_link, sequence, suffixes, table);

    suffix_links[node_index] = suffix_link_destination;
  }
}

/**! \brief Adds suffix links for all (extended) implicit nodes in the tree
 * sequentially.
 *
 * Finds nodes which have not been assigned suffix links.  If the parent
 * has a suffix link, we can iterate from the parent's suffix link to find
 * the suffix link of the child, if it exists.
 *
 * \tparam alphabet_t seqan3 alphabet for the tree.
 * \param sequence Sequence of the tree
 * \param suffixes Suffixes of the tree.
 * \param table Table of the tree.
 * \param[out] suffix_links Suffix links of each explicit node in the tree.
 */
template <seqan3::alphabet alphabet_t>
void add_implicit_suffix_links_parallel(size_t start_index, size_t start_parent,
                                        size_t start_depth,
                                        const sequence_t<alphabet_t> &sequence,
                                        const std::vector<size_t> &suffixes,
                                        const Table<> &table,
                                        std::vector<size_t> &suffix_links,
                                        int parallel_depth) {

  std::queue<std::tuple<size_t, size_t, size_t>> queue{};
  std::vector<std::thread> threads{};
  queue.emplace(start_index, start_parent, start_depth);

  while (!queue.empty()) {
    auto [node_index, parent_index, depth] = queue.front();
    queue.pop();

    visit_implicit_suffix_links_node(node_index, parent_index, sequence,
                                     suffixes, table, suffix_links);

    if (depth < parallel_depth) {
      iterate_children(
          node_index, table,
          [&, node_index = node_index, depth = depth](size_t index) {
            threads.push_back(std::thread{
                add_implicit_suffix_links_parallel<alphabet_t>, index,
                node_index, depth + 1, std::ref(sequence), std::ref(suffixes),
                std::ref(table), std::ref(suffix_links), parallel_depth});
          });
    } else {
      iterate_children(
          node_index, table,
          [&, node_index = node_index, depth = depth](size_t index) {
            queue.emplace(index, node_index, depth + 1);
          });
    }
  }
  for (auto &thread : threads) {
    if (thread.joinable()) {
      thread.join();
    }
  }
}

} // namespace lst::details
