#pragma once

#include <queue>
#include <stack>
#include <vector>

#include <seqan3/alphabet/concept.hpp>

#include "construction.hpp"
#include "iteration.hpp"
#include <thread>
#include <mutex>
#include <chrono>

std::mutex hight_lock;
using namespace std::chrono;
namespace lst::details {

template <seqan3::alphabet alphabet_t>
void tree_height_paralell( sequence_t<alphabet_t> &sequence,
                           std::vector<int64_t> &suffixes,
                           std::vector<int64_t> &table,
                           std::vector<Flag> &flags, int64_t &tree_height,
                           int64_t node_index, int64_t parent_depth) {

  std::queue<std::tuple<int64_t, int64_t>> queue{};
  queue.emplace(node_index, parent_depth);

  while (!queue.empty()) {
    auto [node_index, parent_depth] = queue.front();
    queue.pop();

    if (is_leaf(node_index, flags)) {
      continue;
    }

    int64_t node_depth = parent_depth +
              get_edge_lcp(node_index, sequence, suffixes, table, flags);

    iterate_children(node_index, table, flags,
                     [&](int64_t index) {
                         queue.emplace(index, node_depth);
                     });

    // Has to be before the check to prevent read before write.
    std::lock_guard<std::mutex> lock(hight_lock);
    if (node_depth > tree_height) {
      tree_height = node_depth;
    }
  }

}
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
int64_t tree_height(sequence_t<alphabet_t> &sequence,
                    std::vector<int64_t> &suffixes,
                    std::vector<int64_t> &table,
                    std::vector<Flag> &flags,
                    bool &multi_core, int paralell_depth) {

  std::queue<std::tuple<int64_t, int64_t, int>> queue{};
  queue.emplace(0, 0, 0);
  int64_t tree_height = 0;
  while (!queue.empty()) {

    auto [node_index, parent_depth, level] = queue.front();


    if (paralell_depth == level && multi_core){
      std::thread threads[queue.size()];
      seqan3::debug_stream << "Starting upto " << queue.size() << "threads." << std::endl;
      int thread_index = 0;
      while (!queue.empty()) {
        auto [node_index, parent_depth, level] = queue.front();
        queue.pop();
        if (sequence[table[node_index]].to_rank() != 3 && sequence[table[node_index]].to_rank() != 5) {
          threads[thread_index] =
                  std::thread(tree_height_paralell<seqan3::dna5>,
                              std::ref(sequence),
                              std::ref(suffixes),
                              std::ref(table),
                              std::ref(flags),
                              std::ref(tree_height),
                              node_index, parent_depth);
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
    int64_t node_depth = parent_depth +
                            get_edge_lcp(node_index, sequence, suffixes, table, flags);

    level++;
    iterate_children(node_index, table, flags,
                     [&](int64_t index) {
                       queue.emplace(index, node_depth, level);
                     });

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
int64_t get_leaf_index(int64_t node_index, int64_t lcp,
                            std::vector<int64_t> &suffixes,
                            std::vector<int64_t> &table,
                            std::vector<Flag> &flags) {

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
                               std::vector<int64_t> &suffixes,
                               std::vector<int64_t> &table,
                               std::vector<Flag> &flags,
                               std::vector<int64_t> &suffix_links,
                               bool &multi_core, int paralell_depth) {

  std::vector<std::tuple<int64_t, int64_t>> cause(suffixes.size());

  seqan3::debug_stream << "    Preparing Suffix Links..." <<  std::endl;

  auto start = high_resolution_clock::now();

  prepare_suffix_links(0, 0, cause, sequence, suffixes, table, flags);

  seqan3::debug_stream << "    Calculating tree height..." <<  std::endl;

  auto t1 = high_resolution_clock::now();

  int64_t height = tree_height(sequence, suffixes, table, flags, multi_core, paralell_depth);

  std::vector<int64_t> branch(height + 1, -1);

  auto t2 = high_resolution_clock::now();

  seqan3::debug_stream << "    Computing Suffix Links..." <<  std::endl;

  compute_suffix_links(cause, branch, sequence, suffixes,
                       table, flags, suffix_links, multi_core,
                       paralell_depth);
  multi_core = false;
  auto stop = high_resolution_clock::now();

  auto duration = duration_cast<seconds>(stop - start);
  auto compute = duration_cast<seconds>(stop - t2);
  auto tree_height_time = duration_cast<seconds>(t2 - t1);
  auto prepare_time = duration_cast<seconds>(t1 - start);
      std::cout << "    Duration: " << duration.count() << std::endl;
      std::cout << "        Prepare: " << prepare_time.count() << std::endl;
      std::cout << "        Tree Height: " << tree_height_time.count() << std::endl;
      std::cout << "        Compute Suffix Links: " << compute.count() << std::endl;


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
int64_t prepare_suffix_links(int64_t node_index, int64_t lcp,
                            std::vector<std::tuple<int64_t, int64_t>> &cause,
                            sequence_t<alphabet_t> &sequence,
                            std::vector<int64_t> &suffixes,
                            std::vector<int64_t> &table,
                            std::vector<Flag> &flags) {

  if (is_leaf(node_index, flags)) {
    return get_leaf_index(node_index, lcp, suffixes, table, flags);
  } else if (is_unevaluated(node_index, flags)) {
    // This isn't technically correct.  The following gives the second smallest
    // suffix _overall_ and the algorithm calls for the second smallest suffix
    // from the immediate children.  Not entirely sure that it is important
    // since the values will still be unique.

//    auto second_smallest_child = suffixes[table[node_index] + 1] - lcp;

//    if (second_smallest_child < suffixes.size()) {
//      cause[second_smallest_child + 1] =
//          std::make_tuple(node_index, lcp + edge_lcp);
//    }

    return get_leaf_index(node_index, lcp, suffixes, table, flags);
  } else {

    int64_t edge_lcp =
            get_edge_lcp(node_index, sequence, suffixes, table, flags);

    int64_t smallest_child        = suffixes.size();
    int64_t second_smallest_child = suffixes.size();

    iterate_children(node_index, table,  flags, [&](int64_t index) {
      int64_t child = prepare_suffix_links(index, lcp + edge_lcp, cause, sequence,
                                       suffixes, table, flags);

      if (child < smallest_child) {
        second_smallest_child = smallest_child;
        smallest_child = child;
      } else if (child < second_smallest_child) {
        second_smallest_child = child;
      }
    });

    cause[second_smallest_child + 1] =
          std::make_tuple(node_index, lcp + edge_lcp);

    return smallest_child;
  }
}

void assign_link(int64_t leaf_index,
                 std::vector<std::tuple<int64_t, int64_t>> &cause,
                 std::vector<int64_t> &branch,
                 std::vector<int64_t> &suffix_links) {
  auto &[caused, depth] = cause[leaf_index];
  if (caused != -1 && depth != -1) {
    suffix_links[caused / 2] = branch[depth - 1];
  }
}

void assign_leaf_link(int64_t node_index, int64_t leaf_index,
                      std::vector<int64_t> &suffix_links,
                      std::vector<int64_t> &leaf_indices) {
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
void compute_suffix_links_paralell(std::vector<std::tuple<int64_t, int64_t>> &cause,
                                   std::vector<int64_t> &branch,
                                   sequence_t<alphabet_t> &sequence,
                                   std::vector<int64_t> &suffixes,
                                   std::vector<int64_t> &table,
                                   std::vector<Flag> &flags,
                                   std::vector<int64_t> &suffix_links,
                                   int64_t node_index, int64_t lcp) {

 std::stack<std::tuple<int64_t, int64_t>> stack{};
 stack.emplace(node_index, lcp);

 while (!stack.empty()) {
   auto [node_index, lcp] = stack.top();

   stack.pop();
    if (is_leaf(node_index, flags)) {
      int64_t leaf_index =
               get_leaf_index(node_index, lcp, suffixes, table, flags);

      assign_link(leaf_index, cause, branch, suffix_links);
    } else {
      int64_t height = lcp +
             get_edge_lcp(node_index, sequence, suffixes, table, flags);
      branch[height] = node_index;

      if (is_unevaluated(node_index, flags)) {
        for (int64_t i = table[node_index]; i < table[node_index + 1]; i++) {
          int64_t leaf_index = suffixes[i] - lcp;
          assign_link(leaf_index, cause, branch, suffix_links);
        }
      } else {
        iterate_children(node_index, table, flags,
                 [&](int64_t index) { stack.emplace(index, height); });
     }
   }
 }
}

template <seqan3::alphabet alphabet_t>
void compute_suffix_links(std::vector<std::tuple<int64_t, int64_t>> &cause,
                          std::vector<int64_t> &branch,
                          sequence_t<alphabet_t> &sequence,
                          std::vector<int64_t> &suffixes,
                          std::vector<int64_t> &table,
                          std::vector<Flag> &flags,
                          std::vector<int64_t> &suffix_links,
                          bool &multi_core, int paralell_depth) {

  std::stack<std::tuple<int64_t, int64_t>> stack{};
  stack.emplace(0, 0);
  int level = 0;
  while (!stack.empty()) {
    auto [node_index, lcp] = stack.top();

    if (level == 1 && multi_core){
      std::thread threads[stack.size()];
      int thread_index = 0;
      while (!stack.empty()) {
        auto [node_index, lcp] = stack.top();
        stack.pop();
        if (sequence[table[node_index]].to_rank() != 3 && sequence[table[node_index]].to_rank() != 5) {
          seqan3::debug_stream << "Node_ID: " << std::setw(4) << node_index << " | Table[Node_id]: " << table[node_index]<< " | Sequence[Table[node_id]]: " << sequence[table[node_index]] << " " << sequence[table[node_index]].to_rank() << std::endl;
          threads[thread_index] =
                  std::thread(compute_suffix_links_paralell<seqan3::dna5>,
                                                            std::ref(cause),
                                                            std::ref(branch),
                                                            std::ref(sequence),
                                                            std::ref(suffixes),
                                                            std::ref(table),
                                                            std::ref(flags),
                                                            std::ref(suffix_links),
                                                            node_index, lcp);
          thread_index++;
        }
      }
      for (int i = 0; i < thread_index; ++i) {
        threads[i].join();
      }
      seqan3::debug_stream << "All threads returned." << std::endl;
      return;
    }

    stack.pop();
    if (is_leaf(node_index, flags)) {
      int64_t leaf_index =
                get_leaf_index(node_index, lcp, suffixes, table, flags);

      assign_link(leaf_index, cause, branch, suffix_links);
    } else {

      int64_t height = lcp +
              get_edge_lcp(node_index, sequence, suffixes, table, flags);
      branch[height] = node_index;

      if (is_unevaluated(node_index, flags)) {
        for (int64_t i = table[node_index]; i < table[node_index + 1]; i++) {
          int64_t leaf_index = suffixes[i] - lcp;
          assign_link(leaf_index, cause, branch, suffix_links);
        }
      } else {
        level++;
        iterate_children(node_index, table, flags,
                    [&](int64_t index) { stack.emplace(index, height); });
      }
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
                           std::vector<int64_t> &suffixes,
                           std::vector<int64_t> &table,
                           std::vector<Flag> &flags,
                           std::vector<int64_t> &suffix_links) {

  std::vector<int64_t> leaf_indices(suffixes.size() + 1, -1);
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
void prepare_leaf_suffix_links(std::vector<int64_t> &leaf_indices,
                               sequence_t<alphabet_t> &sequence,
                               std::vector<int64_t> &suffixes,
                               std::vector<int64_t> &table,
                               std::vector<Flag> &flags) {
  std::stack<std::tuple<int64_t, int64_t>> stack{};
  stack.emplace(0, 0);

  while (!stack.empty()) {
    auto [node_index, lcp] = stack.top();
    stack.pop();

    if (is_leaf(node_index, flags)) {
      int64_t leaf_index = get_leaf_index(node_index, lcp, suffixes, table, flags);

      leaf_indices[leaf_index] = node_index;
    }
    int64_t edge_lcp = get_edge_lcp(node_index, sequence, suffixes, table, flags);
    iterate_children(node_index, table, flags,
                     [&](int64_t index) { stack.emplace(index, lcp + edge_lcp); });
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
void compute_leaf_suffix_links(std::vector<int64_t> &leaf_indices,
                               sequence_t<alphabet_t> &sequence,
                               std::vector<int64_t> &suffixes,
                               std::vector<int64_t> &table,
                               std::vector<Flag> &flags,
                               std::vector<int64_t> &suffix_links) {
  std::stack<std::tuple<int64_t, int64_t>> stack{};
  stack.emplace(0, 0);

  while (!stack.empty()) {
    auto [node_index, lcp] = stack.top();
    stack.pop();

    if (is_leaf(node_index, flags)) {
      int64_t leaf_index = get_leaf_index(node_index, lcp, suffixes, table, flags);

      assign_leaf_link(node_index, leaf_index, suffix_links, leaf_indices);
    }

    int64_t edge_lcp = get_edge_lcp(node_index, sequence, suffixes, table, flags);

    iterate_children(node_index, table, flags,
                     [&](int64_t index) { stack.emplace(index, lcp + edge_lcp); });
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
    std::vector<std::tuple<int64_t, int64_t>> &closest_suffix_link_destination,
    std::vector<int64_t> &parent_links, sequence_t<alphabet_t> &sequence,
    std::vector<int64_t> &suffixes, std::vector<int64_t> &table,
    std::vector<Flag> &flags, std::vector<int64_t> &suffix_links) {

  std::stack<std::tuple<int64_t, int64_t, int64_t>> stack{};
  stack.emplace(0, 0, 0);

  std::vector<std::tuple<int64_t, int64_t>> missing_suffix_links{};
  while (!stack.empty()) {
    auto [node_index, parent, parent_lcp] = stack.top();
    stack.pop();

    int64_t edge_lcp = get_edge_lcp(node_index, sequence, suffixes, table, flags);
    int64_t lcp = parent_lcp + edge_lcp;

    parent_links[node_index / 2] = parent;

    if (!is_leaf(node_index, flags) && suffix_links[node_index / 2] == -1) {
      missing_suffix_links.emplace_back(node_index, lcp);

    } else if (missing_suffix_links.size() > 0) {
      int64_t suffix_link_destination = suffix_links[node_index / 2];

      for (int64_t i = 0; i < missing_suffix_links.size(); i++) {
        auto &[missing, missing_lcp] = missing_suffix_links[i];

        closest_suffix_link_destination[missing] =
            std::make_tuple(suffix_link_destination, lcp - missing_lcp);
      }

      missing_suffix_links.resize(0);
    }

    iterate_children(node_index, table, flags,
                     [&](int64_t child) { stack.emplace(child, node_index, lcp); });
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
    std::vector<std::tuple<int64_t, int64_t>> &closest_suffix_link_destination,
    std::vector<int64_t> &parent_links, sequence_t<alphabet_t> &sequence,
    std::vector<int64_t> &suffixes, std::vector<int64_t> &table,
    std::vector<Flag> &flags, std::vector<int64_t> &suffix_links) {
  std::queue<int64_t> queue{};
  queue.push(0);

  while (!queue.empty()) {
    int64_t node_index = queue.front();
    queue.pop();

    auto [suffix_link_destination, distance] =
        closest_suffix_link_destination[node_index];

    if (suffix_link_destination == -1) {
      suffix_links[node_index / 2] = -1;
    } else if (suffix_link_destination != 0) {
      int64_t destination_parent = suffix_link_destination;

      int64_t edge_lcp =
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
                     [&](int64_t index) { queue.push(index); });
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
                               std::vector<int64_t> &suffixes,
                               std::vector<int64_t> &table,
                               std::vector<Flag> &flags,
                               std::vector<int64_t> &suffix_links) {
  std::vector<std::tuple<int64_t, int64_t>> closest_suffix_link_destination(
      table.size());

  std::vector<int64_t> parent_links(table.size() / 2, -1);

  prepare_implicit_suffix_links(closest_suffix_link_destination, parent_links,
                                sequence, suffixes, table, flags, suffix_links);
  compute_implicit_suffix_links(closest_suffix_link_destination, parent_links,
                                sequence, suffixes, table, flags, suffix_links);
}
} // namespace lst::details
