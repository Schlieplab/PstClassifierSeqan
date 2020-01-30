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

  std::stack<std::tuple<int64_t, int64_t, int>> stack{};
  std::vector<std::thread> threads{};
  stack.emplace(0, 0, 0);
  while (!stack.empty()) {
    auto [node_index, lcp, level] = stack.top();
    stack.pop();
    if (level == paralell_depth && multi_core){
      if (sequence[table[node_index]].to_rank() != 3 && sequence[table[node_index]].to_rank() != 5) {
        seqan3::debug_stream << "Node_ID: " << std::setw(4) << node_index << " | Table[Node_id]: "
                             << table[node_index] << " | Sequence[Table[node_id]]: " << sequence[table[node_index]]
                             << " " << sequence[table[node_index]].to_rank() << std::endl;
        threads.push_back(std::thread(compute_suffix_links_paralell<seqan3::dna5>,
                                      std::ref(cause),
                                      std::ref(branch),
                                      std::ref(sequence),
                                      std::ref(suffixes),
                                      std::ref(table),
                                      std::ref(flags),
                                      std::ref(suffix_links),
                                      node_index, lcp));

      }
      continue;
    } else {

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
                           [&](int64_t index) { stack.emplace(index, height, level); });
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
                           std::vector<int64_t> &suffixes,
                           std::vector<int64_t> &table,
                           std::vector<Flag> &flags,
                           std::vector<int64_t> &suffix_links,
                           bool &multi_core, int paralell_depth) {

  auto start = high_resolution_clock::now();

  std::vector<int64_t> leaf_indices(suffixes.size() + 1, -1);
  leaf_indices[suffixes.size()] = 0;

  prepare_leaf_suffix_links(leaf_indices, sequence, suffixes, table, flags, multi_core, paralell_depth);
  auto t1 = high_resolution_clock::now();

  compute_leaf_suffix_links(leaf_indices, sequence, suffixes, table, flags,
                            suffix_links, multi_core, paralell_depth);

  auto stop = high_resolution_clock::now();

  auto duration = duration_cast<seconds>(stop - start);
  auto compute = duration_cast<seconds>(stop - t1);
  auto prepare_time = duration_cast<seconds>(t1 - start);
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
void prepare_leaf_suffix_links_p(std::vector<int64_t> &leaf_indices,
                               sequence_t<alphabet_t> &sequence,
                               std::vector<int64_t> &suffixes,
                               std::vector<int64_t> &table,
                               std::vector<Flag> &flags,
                               int64_t node_index, int64_t lcp) {
  std::stack<std::tuple<int64_t, int64_t>> stack{};
  stack.emplace(node_index, lcp);

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
template <seqan3::alphabet alphabet_t>
void prepare_leaf_suffix_links(std::vector<int64_t> &leaf_indices,
                               sequence_t<alphabet_t> &sequence,
                               std::vector<int64_t> &suffixes,
                               std::vector<int64_t> &table,
                               std::vector<Flag> &flags,
                               bool multi_core, int paralell_depth) {
  std::stack<std::tuple<int64_t, int64_t, int>> stack{};
  stack.emplace(0, 0, 0);
  std::vector<std::thread> threads{};

  while (!stack.empty()) {
    auto [node_index, lcp, level] = stack.top();
    stack.pop();
    if (level == paralell_depth && multi_core){
      if (sequence[table[node_index]].to_rank() != 3 && sequence[table[node_index]].to_rank() != 5) {
        threads.push_back(std::thread(prepare_leaf_suffix_links_p<seqan3::dna5>,
                                      std::ref(leaf_indices),
                                      std::ref(sequence),
                                      std::ref(suffixes),
                                      std::ref(table),
                                      std::ref(flags),
                                      node_index, lcp));
      }
      continue;
    } else {
      if (is_leaf(node_index, flags)) {
        int64_t leaf_index = get_leaf_index(node_index, lcp, suffixes, table, flags);

        leaf_indices[leaf_index] = node_index;
      }
      int64_t edge_lcp = get_edge_lcp(node_index, sequence, suffixes, table, flags);
      level++;
      iterate_children(node_index, table, flags,
                       [&](int64_t index) { stack.emplace(index, lcp + edge_lcp, level); });
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
void compute_leaf_suffix_links_p(std::vector<int64_t> &leaf_indices,
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
     int64_t leaf_index = get_leaf_index(node_index, lcp, suffixes, table, flags);

     assign_leaf_link(node_index, leaf_index, suffix_links, leaf_indices);
   }

   int64_t edge_lcp = get_edge_lcp(node_index, sequence, suffixes, table, flags);

   iterate_children(node_index, table, flags,
                    [&](int64_t index) { stack.emplace(index, lcp + edge_lcp); });
 }
}

template <seqan3::alphabet alphabet_t>
void compute_leaf_suffix_links(std::vector<int64_t> &leaf_indices,
                               sequence_t<alphabet_t> &sequence,
                               std::vector<int64_t> &suffixes,
                               std::vector<int64_t> &table,
                               std::vector<Flag> &flags,
                               std::vector<int64_t> &suffix_links,
                               bool multi_core, int paralell_depth) {

  std::stack<std::tuple<int64_t, int64_t, int>> stack{};
  stack.emplace(0, 0, 0);
  std::vector<std::thread> threads{};

  while (!stack.empty()) {
    auto [node_index, lcp, level] = stack.top();
    stack.pop();
    if (level == paralell_depth && multi_core){
      if (sequence[table[node_index]].to_rank() != 3 && sequence[table[node_index]].to_rank() != 5) {
        threads.push_back(std::thread(compute_leaf_suffix_links_p<seqan3::dna5>,
                                      std::ref(leaf_indices),
                                      std::ref(sequence),
                                      std::ref(suffixes),
                                      std::ref(table),
                                      std::ref(flags),
                                      std::ref(suffix_links),
                                      node_index, lcp));
      }
      continue;
    } else {
      if (is_leaf(node_index, flags)) {
        int64_t leaf_index = get_leaf_index(node_index, lcp, suffixes, table, flags);

        assign_leaf_link(node_index, leaf_index, suffix_links, leaf_indices);
      }

      int64_t edge_lcp = get_edge_lcp(node_index, sequence, suffixes, table, flags);
      level++;
      iterate_children(node_index, table, flags,
                       [&](int64_t index) { stack.emplace(index, lcp + edge_lcp, level); });
    }
  }
  for (int i = 0; i < threads.size(); ++i) {
    threads[i].join();
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
void prepare_implicit_suffix_links_p(
   std::vector<std::tuple<int64_t, int64_t>> &closest_suffix_link_destination,
   std::vector<int64_t> &parent_links, sequence_t<alphabet_t> &sequence,
   std::vector<int64_t> &suffixes, std::vector<int64_t> &table,
   std::vector<Flag> &flags, std::vector<int64_t> &suffix_links,
   int64_t node_index, int64_t parent, int64_t parent_lcp) {

 std::stack<std::tuple<int64_t, int64_t, int64_t>> stack{};
 stack.emplace(node_index, parent, parent_lcp);

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
bool sequences_match(int64_t node_index, int64_t edge_lcp, int64_t suffix_link_child_index,
                     int64_t suffix_link_edge_lcp, sequence_t<alphabet_t> &sequence,
                     std::vector<int64_t> &suffixes, std::vector<int64_t> &table,
                     std::vector<Flag> &flags) {
    int64_t node_start = get_sequence_index(node_index, suffixes, table, flags);
    int64_t suffix_link_child_end = std::min(
            int64_t(sequence.size()),
          get_sequence_index(suffix_link_child_index, suffixes, table, flags) +
          suffix_link_edge_lcp);
    int64_t suffix_link_child_start = suffix_link_child_end - edge_lcp;

  for (int64_t i = 0; i < edge_lcp; i++) {
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
int64_t find_suffix_match(int64_t node_index, int64_t edge_lcp, int64_t parent_suffix_link,
                      sequence_t<alphabet_t> &sequence,
                      std::vector<int64_t> &suffixes, std::vector<int64_t> &table,
                      std::vector<Flag> &flags) {
  std::queue<std::tuple<int64_t, int64_t>> suffix_link_queue{};
  iterate_children(parent_suffix_link, table, flags,
                   [&](int64_t index) { suffix_link_queue.emplace(index, 0); });

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

      iterate_children(suffix_link_child, table, flags, [&](int64_t index) {
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
                               std::vector<int64_t> &suffixes,
                               std::vector<int64_t> &table,
                               std::vector<Flag> &flags,
                               std::vector<int64_t> &suffix_links,
                               int64_t node_index, int64_t parent_index) {

      std::queue<std::tuple<int64_t, int64_t>> queue{};
      queue.emplace(node_index, parent_index);

      while (!queue.empty()) {
        auto[node_index, parent_index] = queue.front();
        queue.pop();

        auto edge_lcp = get_edge_lcp(node_index, sequence, suffixes, table, flags);

        if (suffix_links[node_index / 2] == -1 && parent_index == 0 &&
            edge_lcp == 1) {
          suffix_links[node_index / 2] = 0;
        } else if (suffix_links[node_index / 2] == -1 &&
                   suffix_links[parent_index / 2] != -1) {
          auto parent_suffix_link = suffix_links[parent_index / 2];
          int64_t suffix_link_destination =
                  find_suffix_match(node_index, edge_lcp, parent_suffix_link, sequence,
                                    suffixes, table, flags);

          suffix_links[node_index / 2] = suffix_link_destination;
        }

        iterate_children(node_index, table, flags,
                         [&](int64_t index) { queue.emplace(index, node_index); });
      }
}
template <seqan3::alphabet alphabet_t>
void add_implicit_suffix_links(sequence_t<alphabet_t> &sequence,
                               std::vector<int64_t> &suffixes,
                               std::vector<int64_t> &table,
                               std::vector<Flag> &flags,
                               std::vector<int64_t> &suffix_links,
                               bool multi_core, int paralell_depth) {

  std::queue<std::tuple<int64_t, int64_t, int>> queue{};
  std::vector<std::thread> threads{};
  queue.emplace(0, 0, 0);

  while (!queue.empty()) {
    auto[node_index, parent_index, level] = queue.front();
    queue.pop();
    if (paralell_depth == level && multi_core) {
      if (sequence[table[node_index]].to_rank() != 3 && sequence[table[node_index]].to_rank() != 5) {

        threads.push_back(std::thread(add_implicit_suffix_links_p<seqan3::dna5>,
                                      std::ref(sequence),
                                      std::ref(suffixes),
                                      std::ref(table),
                                      std::ref(flags),
                                      std::ref(suffix_links),
                                      node_index, parent_index));

      }

      continue;
    }else{
      auto edge_lcp = get_edge_lcp(node_index, sequence, suffixes, table, flags);

      if (suffix_links[node_index / 2] == -1 && parent_index == 0 &&
          edge_lcp == 1) {
        suffix_links[node_index / 2] = 0;
      } else if (suffix_links[node_index / 2] == -1 &&
                 suffix_links[parent_index / 2] != -1) {
        auto parent_suffix_link = suffix_links[parent_index / 2];
        int64_t suffix_link_destination =
                find_suffix_match(node_index, edge_lcp, parent_suffix_link, sequence,
                                  suffixes, table, flags);

        suffix_links[node_index / 2] = suffix_link_destination;
      }

      level++;
      iterate_children(node_index, table, flags,
                       [&](int64_t index) { queue.emplace(index, node_index, level); });
    }
  }
  for (int i = 0; i < threads.size(); ++i) {
    threads[i].join();
  }
}

} // namespace lst::details
