#pragma once

#include <queue>
#include <stack>
#include <vector>

#include <seqan3/alphabet/concept.hpp>

#include "construction.hpp"
#include "iteration.hpp"
#include <chrono>
#include <mutex>
#include <thread>

std::mutex height_lock;
using namespace std::chrono;
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

/**! \brief Calculates the internal height/depth of the tree and of every node.
 * @tparam alphabet_t Type of alphabet used (from seqan3)
 * \param[in] sequence Sequence of the tree.
 * \param[in] suffixes Suffixes of the tree.
 * \param[in] table Table of the tree.
 * \param[in] flags Flags of the tree.
 * \return the height of the tree (longest node contained, not leaf).
 */
template <seqan3::alphabet alphabet_t>
void tree_height_parallel(sequence_t<alphabet_t> &sequence,
                          std::vector<int> &suffixes, std::vector<int> &table,
                          std::vector<Flag> &flags, int &tree_height,
                          int node_index, int parent_depth) {

  std::queue<std::tuple<int, int>> queue{};
  queue.emplace(node_index, parent_depth);

  while (!queue.empty()) {
    auto [node_index, parent_depth] = queue.front();
    queue.pop();

    if (is_leaf(node_index, flags)) {
      continue;
    }

    int node_depth = parent_depth +
                     get_edge_lcp(node_index, sequence, suffixes, table, flags);

    iterate_children(node_index, table, flags,
                     [&](int index) { queue.emplace(index, node_depth); });

    // Has to be before the check to prevent read before write.
    std::lock_guard<std::mutex> lock(height_lock);
    if (node_depth > tree_height) {
      tree_height = node_depth;
    }
  }
}

/**! \brief Calculates the height/depth of the tree and of every node.
 *
 * @tparam alphabet_t Type of alphabet used (from seqan3)
 * \param[in] sequence Sequence of the tree.
 * \param[in] suffixes Suffixes of the tree.
 * \param[in] table Table of the tree.
 * \param[in] flags Flags of the tree.
 * \return the height of the tree (longest node contained, not leaf).
 */
template <seqan3::alphabet alphabet_t>
int tree_height(sequence_t<alphabet_t> &sequence, std::vector<int> &suffixes,
                std::vector<int> &table, std::vector<Flag> &flags,
                bool &multi_core, int parallel_depth) {
  std::queue<std::tuple<int, int, int>> queue{};
  queue.emplace(0, 0, 0);

  int tree_height = 0;

  // TODO refactor this
  while (!queue.empty()) {
    auto [node_index, parent_depth, level] = queue.front();
    if (parallel_depth == level && multi_core) {
      std::thread threads[queue.size()];
      seqan3::debug_stream << "Starting upto " << queue.size() << "threads."
                           << std::endl;
      int thread_index = 0;
      while (!queue.empty()) {
        auto [node_index, parent_depth, level] = queue.front();
        queue.pop();

        if (sequence[table[node_index]].to_rank() != 3 &&
            sequence[table[node_index]].to_rank() != 5) {
          threads[thread_index] = std::thread(
              tree_height_parallel<alphabet_t>, std::ref(sequence),
              std::ref(suffixes), std::ref(table), std::ref(flags),
              std::ref(tree_height), node_index, parent_depth);
          thread_index++;
        }
      }
      for (int i = 0; i < thread_index; ++i) {
        threads[i].join();
      }
      seqan3::debug_stream << "All threads returned." << std::endl;
      return tree_height;
    }

    queue.pop();

    if (is_leaf(node_index, flags)) {
      continue;
    }

    int edge_lcp = get_edge_lcp(node_index, sequence, suffixes, table, flags);
    int node_depth = parent_depth + edge_lcp;

    level++;
    iterate_children(node_index, table, flags, [&](int index) {
      queue.emplace(index, node_depth, level);
    });

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
                               std::vector<int> &suffix_links, bool &multi_core,
                               int parallel_depth) {

  std::vector<std::tuple<int, int>> cause(suffixes.size());

  seqan3::debug_stream << "    Preparing Suffix Links..." << std::endl;

  auto start = std::chrono::high_resolution_clock::now();

  prepare_suffix_links(0, 0, cause, sequence, suffixes, table, flags);

  seqan3::debug_stream << "    Calculating tree height..." << std::endl;

  auto t1 = std::chrono::high_resolution_clock::now();

  int height =
      tree_height(sequence, suffixes, table, flags, multi_core, parallel_depth);

  std::vector<int> branch(height + 1, -1);

  auto t2 = std::chrono::high_resolution_clock::now();

  seqan3::debug_stream << "    Computing Suffix Links..." << std::endl;

  compute_suffix_links(cause, branch, sequence, suffixes, table, flags,
                       suffix_links, multi_core, parallel_depth);
  auto stop = std::chrono::high_resolution_clock::now();

  auto duration = std::chrono::duration_cast<seconds>(stop - start);
  auto compute = std::chrono::duration_cast<seconds>(stop - t2);
  auto tree_height_time = std::chrono::duration_cast<seconds>(t2 - t1);
  auto prepare_time = std::chrono::duration_cast<seconds>(t1 - start);
  seqan3::debug_stream << "    Duration: " << duration.count() << std::endl;
  seqan3::debug_stream << "        Prepare: " << prepare_time.count()
                       << std::endl;
  seqan3::debug_stream << "        Tree Height: " << tree_height_time.count()
                       << std::endl;
  seqan3::debug_stream << "        Compute Suffix Links: " << compute.count()
                       << std::endl;
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
 * \param[out] suffix_links The suffix link for each node.
 */
void assign_link(int leaf_index, std::vector<std::tuple<int, int>> &cause,
                 std::vector<int> &branch,
                  std::vector<int> &suffix_links) {
  auto &[caused, depth] = cause[leaf_index];

  if (caused != -1 && depth != -1) {
    suffix_links[caused / 2] = branch[depth - 1];
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
 * \param sequence Sequence of the tree
 * \param suffixes Suffixes of the tree.
 * \param table Table of the tree.
 * \param flags Flags of the tree.
 * \param[out] suffix_links The suffix link for each node.
 */
template <seqan3::alphabet alphabet_t>
void compute_suffix_links_parallel(
    std::vector<std::tuple<int, int>> &cause, std::vector<int> &branch,
    sequence_t<alphabet_t> &sequence, std::vector<int> &suffixes,
    std::vector<int> &table, std::vector<Flag> &flags,
    std::vector<int> &suffix_links, int node_index, int lcp) {

  std::stack<std::tuple<int, int>> stack{};
  stack.emplace(node_index, lcp);

  while (!stack.empty()) {
    auto [node_index, lcp] = stack.top();

    stack.pop();
    if (is_leaf(node_index, flags)) {
      int leaf_index = get_leaf_index(node_index, lcp, suffixes, table, flags);

      assign_link(leaf_index, cause, branch, suffix_links);
    } else {
      int height =
          lcp + get_edge_lcp(node_index, sequence, suffixes, table, flags);
      branch[height] = node_index;

      if (is_unevaluated(node_index, flags)) {
        for (int i = table[node_index]; i < table[node_index + 1]; i++) {
          int leaf_index = suffixes[i] - lcp;
          assign_link(leaf_index, cause, branch, suffix_links);
        }
      } else {
        iterate_children(node_index, table, flags,
                         [&](int index) { stack.emplace(index, height); });
      }
    }
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
                          std::vector<int> &suffix_links, bool &multi_core,
                          int parallel_depth) {

  std::stack<std::tuple<int, int, int>> stack{};
  std::vector<std::thread> threads{};
  stack.emplace(0, 0, 0);

  // TODO Refactor this, combine with _parallel
  while (!stack.empty()) {
    auto [node_index, lcp, level] = stack.top();
    stack.pop();
    if (level == parallel_depth && multi_core) {
      if (sequence[table[node_index]].to_rank() != 3 &&
          sequence[table[node_index]].to_rank() != 5) {
        seqan3::debug_stream
            << "Node_ID: " << std::setw(4) << node_index
            << " | Table[Node_id]: " << table[node_index]
            << " | Sequence[Table[node_id]]: " << sequence[table[node_index]]
            << " " << sequence[table[node_index]].to_rank() << std::endl;

        threads.push_back(
            std::thread(compute_suffix_links_parallel<alphabet_t>,
                        std::ref(cause), std::ref(branch), std::ref(sequence),
                        std::ref(suffixes), std::ref(table), std::ref(flags),
                        std::ref(suffix_links), node_index, lcp));
      }
      continue;
    } else {

      if (is_leaf(node_index, flags)) {
        int leaf_index =
            get_leaf_index(node_index, lcp, suffixes, table, flags);

        assign_link(leaf_index, cause, branch, suffix_links);
      } else {

        int height =
            lcp + get_edge_lcp(node_index, sequence, suffixes, table, flags);
        branch[height] = node_index;

        if (is_unevaluated(node_index, flags)) {
          for (int i = table[node_index]; i < table[node_index + 1]; i++) {
            int leaf_index = suffixes[i] - lcp;
            assign_link(leaf_index, cause, branch, suffix_links);
          }
        } else {
          level++;
          iterate_children(node_index, table, flags, [&](int index) {
            stack.emplace(index, height, level);
          });
        }
      }
    }
  }
  for (int i = 0; i < threads.size(); ++i) {
    threads[i].join();
  }
  seqan3::debug_stream << "All threads returned." << std::endl;
  return;
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
                           std::vector<int> &suffix_links, bool &multi_core,
                           int parallel_depth) {

  auto start = std::chrono::high_resolution_clock::now();

  std::vector<int> leaf_indices(suffixes.size() + 1, -1);
  leaf_indices[suffixes.size()] = 0;

  prepare_leaf_suffix_links(leaf_indices, sequence, suffixes, table, flags,
                            multi_core, parallel_depth);
  auto t1 = std::chrono::high_resolution_clock::now();

  compute_leaf_suffix_links(leaf_indices, sequence, suffixes, table, flags,
                            suffix_links, multi_core, parallel_depth);

  auto stop = std::chrono::high_resolution_clock::now();

  auto duration = std::chrono::duration_cast<seconds>(stop - start);
  auto compute = std::chrono::duration_cast<seconds>(stop - t1);
  auto prepare_time = std::chrono::duration_cast<seconds>(t1 - start);
  std::cout << "    Duration: " << duration.count() << std::endl;
  std::cout << "        Prepare: " << prepare_time.count() << std::endl;
  std::cout << "        Compute Suffix Links: " << compute.count() << std::endl;
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
void prepare_leaf_suffix_links_p(std::vector<int> &leaf_indices,
                                 sequence_t<alphabet_t> &sequence,
                                 std::vector<int> &suffixes,
                                 std::vector<int> &table,
                                 std::vector<Flag> &flags, int node_index,
                                 int lcp) {
  std::stack<std::tuple<int, int>> stack{};
  stack.emplace(node_index, lcp);

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
template <seqan3::alphabet alphabet_t>
void prepare_leaf_suffix_links(std::vector<int> &leaf_indices,
                               sequence_t<alphabet_t> &sequence,
                               std::vector<int> &suffixes,
                               std::vector<int> &table,
                               std::vector<Flag> &flags, bool multi_core,
                               int parallel_depth) {
  std::stack<std::tuple<int, int, int>> stack{};
  stack.emplace(0, 0, 0);
  std::vector<std::thread> threads{};

  while (!stack.empty()) {
    auto [node_index, lcp, level] = stack.top();
    stack.pop();
    if (level == parallel_depth && multi_core) {
      if (sequence[table[node_index]].to_rank() != 3 &&
          sequence[table[node_index]].to_rank() != 5) {
        threads.push_back(std::thread(
            prepare_leaf_suffix_links_p<alphabet_t>, std::ref(leaf_indices),
            std::ref(sequence), std::ref(suffixes), std::ref(table),
            std::ref(flags), node_index, lcp));
      }
      continue;
    } else {
      if (is_leaf(node_index, flags)) {
        int leaf_index =
            get_leaf_index(node_index, lcp, suffixes, table, flags);

        leaf_indices[leaf_index] = node_index;
      }
      int edge_lcp = get_edge_lcp(node_index, sequence, suffixes, table, flags);
      level++;
      iterate_children(node_index, table, flags, [&](int index) {
        stack.emplace(index, lcp + edge_lcp, level);
      });
    }
  }
  for (int i = 0; i < threads.size(); ++i) {
    threads[i].join();
  }
  seqan3::debug_stream << "All threads returned." << std::endl;
  return;
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
void compute_leaf_suffix_links_p(std::vector<int> &leaf_indices,
                                 sequence_t<alphabet_t> &sequence,
                                 std::vector<int> &suffixes,
                                 std::vector<int> &table,
                                 std::vector<Flag> &flags,
                                 std::vector<int> &suffix_links, int node_index,
                                 int lcp) {

  std::stack<std::tuple<int, int>> stack{};
  stack.emplace(node_index, lcp);

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
                               std::vector<int> &suffix_links, bool multi_core,
                               int parallel_depth) {

  std::stack<std::tuple<int, int, int>> stack{};
  stack.emplace(0, 0, 0);
  std::vector<std::thread> threads{};

  while (!stack.empty()) {
    auto [node_index, lcp, level] = stack.top();
    stack.pop();
    if (level == parallel_depth && multi_core) {
      if (sequence[table[node_index]].to_rank() != 3 &&
          sequence[table[node_index]].to_rank() != 5) {
        threads.push_back(std::thread(
            compute_leaf_suffix_links_p<alphabet_t>, std::ref(leaf_indices),
            std::ref(sequence), std::ref(suffixes), std::ref(table),
            std::ref(flags), std::ref(suffix_links), node_index, lcp));
      }
      continue;
    } else {
      if (is_leaf(node_index, flags)) {
        int leaf_index =
            get_leaf_index(node_index, lcp, suffixes, table, flags);

        assign_leaf_link(node_index, leaf_index, suffix_links, leaf_indices);
      }

      int edge_lcp = get_edge_lcp(node_index, sequence, suffixes, table, flags);
      level++;
      iterate_children(node_index, table, flags, [&](int index) {
        stack.emplace(index, lcp + edge_lcp, level);
      });
    }
  }
  for (int i = 0; i < threads.size(); ++i) {
    threads[i].join();
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
  iterate_children(parent_suffix_link, table, flags,
                   [&](int index) { suffix_link_queue.emplace(index, 0); });

  while (!suffix_link_queue.empty()) {
    auto &[suffix_link_child, suffix_link_parent_lcp] =
        suffix_link_queue.front();
    suffix_link_queue.pop();
    auto suffix_link_edge_lcp =
        get_edge_lcp(suffix_link_child, sequence, suffixes, table, flags);
    auto suffix_link_lcp = suffix_link_parent_lcp + suffix_link_edge_lcp;

    if (suffix_link_lcp == edge_lcp) {
      bool match = sequences_match(node_index, edge_lcp, suffix_link_child,
                                   suffix_link_edge_lcp, sequence, suffixes,
                                   table, flags);
      if (match) {
        return suffix_link_child;
      }
    } else if (suffix_link_lcp < edge_lcp) {

      iterate_children(suffix_link_child, table, flags, [&](int index) {
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
 * \param flags Flags of the tree.
 * \param[out] suffix_links Suffix links of each explicit node in the tree.
 */
template <seqan3::alphabet alphabet_t>
void add_implicit_suffix_links_p(sequence_t<alphabet_t> &sequence,
                                 std::vector<int> &suffixes,
                                 std::vector<int> &table,
                                 std::vector<Flag> &flags,
                                 std::vector<int> &suffix_links, int node_index,
                                 int parent_index) {

  std::queue<std::tuple<int, int>> queue{};
  queue.emplace(node_index, parent_index);

  while (!queue.empty()) {
    auto [node_index, parent_index] = queue.front();
    queue.pop();

    auto edge_lcp = get_edge_lcp(node_index, sequence, suffixes, table, flags);

    if (suffix_links[node_index / 2] == -1 && parent_index == 0 &&
        edge_lcp == 1) {
      suffix_links[node_index / 2] = 0;
    } else if (suffix_links[node_index / 2] == -1 &&
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

template <seqan3::alphabet alphabet_t>
void add_implicit_suffix_links(sequence_t<alphabet_t> &sequence,
                               std::vector<int> &suffixes,
                               std::vector<int> &table,
                               std::vector<Flag> &flags,
                               std::vector<int> &suffix_links, bool multi_core,
                               int parallel_depth) {

  std::queue<std::tuple<int, int, int>> queue{};
  std::vector<std::thread> threads{};
  queue.emplace(0, 0, 0);

  // TODO refactor this with above
  while (!queue.empty()) {
    auto [node_index, parent_index, level] = queue.front();
    queue.pop();
    if (parallel_depth == level && multi_core) {
      if (sequence[table[node_index]].to_rank() != 3 &&
          sequence[table[node_index]].to_rank() != 5) {

        threads.push_back(std::thread(
            add_implicit_suffix_links_p<alphabet_t>, std::ref(sequence),
            std::ref(suffixes), std::ref(table), std::ref(flags),
            std::ref(suffix_links), std::ref(node_index), std::ref(parent_index)));
      }

      continue;
    } else {
      auto edge_lcp =
          get_edge_lcp(node_index, sequence, suffixes, table, flags);

      if (suffix_links[node_index / 2] == -1 && parent_index == 0 &&
          edge_lcp == 1) {
        suffix_links[node_index / 2] = 0;
      } else if (suffix_links[node_index / 2] == -1 &&
                 suffix_links[parent_index / 2] != -1) {
        auto parent_suffix_link = suffix_links[parent_index / 2];
        int suffix_link_destination =
            find_suffix_match(node_index, edge_lcp, parent_suffix_link,
                              sequence, suffixes, table, flags);

        suffix_links[node_index / 2] = suffix_link_destination;
      }

      level++;
      iterate_children(node_index, table, flags, [&](int index) {
        queue.emplace(index, node_index, level);
      });
    }
  }
  for (int i = 0; i < threads.size(); ++i) {
    threads[i].join();
  }
}

} // namespace lst::details
