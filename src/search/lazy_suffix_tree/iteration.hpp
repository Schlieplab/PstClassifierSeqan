#pragma once

#include <functional>
#include <map>
#include <mutex>
#include <queue>
#include <thread>
#include <vector>

#include <seqan3/alphabet/concept.hpp>

#include "construction.hpp"

namespace lst::details {

int next_child_index(int node_index, const Table<> &table) {
  if (is_leaf(node_index, table)) {
    // Should be 1, but I've added a value to leaves to allow for
    // explicit nodes.
    return node_index + 2;
  } else {
    return node_index + 2;
  }
}

void iterate_children(int node_index, const Table<> &table,
                      const std::function<void(int)> &f) {
  if (is_leaf(node_index, table) || is_unevaluated(node_index, table)) {
    return;
  }

  int first_child = table[node_index + 1].value;

  for (int i = first_child; i <= table.size();) {
    f(i);

    if (is_rightmostchild(i, table)) {
      break;
    }
    i = next_child_index(i, table);
  }
}

int node_occurrences(int node_index, const Table<> &table) {
  assert(node_index <= table.size());

  int occurrences = 0;
  std::queue<int> queue{};

  queue.push(node_index);

  while (!queue.empty()) {
    int index = queue.front();
    queue.pop();

    if (is_leaf(index, table)) {
      occurrences += 1;
    } else if (is_unevaluated(index, table)) {
      occurrences += table[index + 1].value - table[index].value;
    } else {
      iterate_children(index, table, [&](int i) { queue.push(i); });
    }
  }

  return occurrences;
}

/**
 * Iterates over the tree in a breadth-first fashion without constructing the
 * tree.
 *
 * \tparam alphabet_t seqan3::alphabet type
 * \param start_lower_bound lower bound of unevaluated node
 * \param start_upper_bound upper bound of unevaluated node
 * \param start_lcp Longest common prefix of the node.
 * \param start_level Depth in tree.
 * \param parallel_depth Depth to spawn threads to.
 * \param sequence Sequence the tree is built over.
 * \param suffixes Suffix vector of the tree.
 * \param table Table vector of the tree.
 * \param f Callback for the visited nodes.
 * \param done Callback to signal that there are no more nodes to iterate.
 */
template <seqan3::alphabet alphabet_t>
void breadth_first_iteration_table_less_(
    int start_lower_bound, int start_upper_bound, int start_lcp,
    int start_level, int parallel_depth, const sequence_t<alphabet_t> &sequence,
    std::vector<int> &suffixes,
    const std::function<bool(int, int, int, int, bool)> &f,
    const std::function<void()> &done) {
  std::vector<std::thread> threads{};

  std::queue<std::tuple<int, int, int, int>> queue{};
  queue.emplace(start_lower_bound, start_upper_bound, start_lcp, start_level);

  while (!queue.empty()) {
    auto [lower, upper, lcp, level] = queue.front();
    queue.pop();
    if (lower == 0 && upper == suffixes.size()) {
    }

    if (upper == -1) {
      f(lower, lcp, sequence.size() - lower + 1, 1, true);
      continue;
    }

    auto [sequence_index, edge_lcp, node_count, children] =
        evaluate_node(lower, upper, sequence, suffixes);

    bool consider_children =
        f(sequence_index, lcp, edge_lcp, node_count, false);

    if (!consider_children) {
      continue;
    }
    int new_lcp = lcp + edge_lcp;

    if (level < parallel_depth) {
      for (auto &[child_lower, child_upper] : children) {
        threads.emplace_back(breadth_first_iteration_table_less_<alphabet_t>,
                             child_lower, child_upper, new_lcp, level + 1,
                             parallel_depth, std::ref(sequence),
                             std::ref(suffixes), f, done);
      }
    } else {
      for (auto &[child_lower, child_upper] : children) {
        queue.emplace(child_lower, child_upper, new_lcp, level + 1);
      }
    }
  }

  done();

  for (auto &thread : threads) {
    if (thread.joinable()) {
      thread.join();
    }
  }
}

/**
 * Iterates over the tree in a breadth-first fashion without constructing the
 * tree.
 *
 * \tparam alphabet_t seqan3::alphabet type
 * \param start_lower_bound lower bound of unevaluated node
 * \param start_upper_bound upper bound of unevaluated node
 * \param start_lcp Longest common prefix of the node.
 * \param start_level Depth in tree.
 * \param parallel_depth Depth to spawn threads to.
 * \param sequence Sequence the tree is built over.
 * \param suffixes Suffix vector of the tree.
 * \param table Table vector of the tree.
 * \param f Callback for the visited nodes.
 * \param done Callback to signal that there are no more nodes to iterate.
 */
template <seqan3::alphabet alphabet_t>
void breadth_first_iteration_table_less(
    int parallel_depth, const sequence_t<alphabet_t> &sequence,
    std::vector<int> &suffixes,
    const std::function<bool(int, int, int, int, bool)> &f,
    const std::function<void()> &done) {
  breadth_first_iteration_table_less_<alphabet_t>(
      0, suffixes.size(), 0, 0, parallel_depth, sequence, suffixes, f, done);
}

/**
 * Iterates over the structure in a breadth-first fashion.
 *
 * \tparam alphabet_t seqan3::alphabet type
 * \param sequence Sequence the tree is built over.
 * \param suffixes Suffix vector of the tree.
 * \param table Table vector of the tree.
 * \param expand_nodes Determines if nodes should be explored (construction) or
 * if we only want to visit nodes.
 * \param f Callback for the visited nodes.
 * \param[in] locked_callback Callback for when the tree array has grown and a
 * user wants to do something in a synchronisation lock.
 */
template <seqan3::alphabet alphabet_t>
void breadth_first_iteration(
    const sequence_t<alphabet_t> &sequence, std::vector<int> &suffixes,
    Table<> &table, bool expand_nodes,
    const std::function<bool(int, int, int &, int)> &f,
    const std::function<void(int, int, int &)> &locked_callback) {
  breadth_first_iteration(0, 0, sequence, suffixes, table, expand_nodes, f,
                          locked_callback);
}

/**
 * Iterates over the structure in a breadth-first fashion.
 *
 * \tparam alphabet_t seqan3::alphabet type
 * \param start_index Node index to start iteration from.
 * \param start_lcp Longest common prefix of the start_index node.
 * \param sequence Sequence the tree is built over.
 * \param suffixes Suffix vector of the tree.
 * \param table Table vector of the tree.
 * \param expand_nodes Determines if nodes should be explored (construction) or
 * if we only want to visit nodes.
 * \param f Callback for the visited nodes.
 * \param[in] locked_callback Callback for when the tree array has grown and a
 * user wants to do something in a synchronisation lock.
 */
template <seqan3::alphabet alphabet_t>
void breadth_first_iteration(
    int start_index, int start_lcp, const sequence_t<alphabet_t> &sequence,
    std::vector<int> &suffixes, Table<> &table, bool expand_nodes,
    const std::function<bool(int, int, int &, int)> &f,
    const std::function<void(int, int, int &)> &locked_callback) {
  std::queue<std::tuple<int, int>> queue{};
  queue.emplace(start_index, start_lcp);

  while (!queue.empty()) {
    auto [node_index, lcp] = queue.front();
    queue.pop();

    auto [new_lcp, consider_children] =
        visit_top_node(node_index, lcp, sequence, suffixes, table, expand_nodes,
                       f, locked_callback);

    if (!consider_children) {
      continue;
    }

    iterate_children(node_index, table, [&, new_lcp = new_lcp](int index) {
      queue.emplace(index, new_lcp);
    });
  }
}

/**
 * Iterates over the structure a breadth-first fashion in parallel.
 *
 * \tparam alphabet_t seqan3::alphabet type
 * \param sequence Sequence the tree is built over.
 * \param suffixes Suffix vector of the tree.
 * \param table Table vector of the tree.
 * \param expand_nodes Determines if nodes should be explored (construction) or
 * if we only want to visit nodes. \param f Callback for the visited nodes.
 * \param parallel_depth Number of levels to spawn new processes for.
 * \param[in] locked_callback Callback for when the tree array has grown and a
 * user wants to do something in a synchronisation lock.
 */
template <seqan3::alphabet alphabet_t>
void breadth_first_iteration_parallel(
    const sequence_t<alphabet_t> &sequence, std::vector<int> &suffixes,
    Table<> &table, bool expand_nodes,
    const std::function<bool(int, int, int &, int)> &f,
    const std::function<void()> &done, int parallel_depth,
    const std::function<void(int, int, int &)> &locked_callback) {
  breadth_first_iteration_parallel_(0, 0, 0, sequence, suffixes, table,
                                    expand_nodes, f, done, parallel_depth,
                                    locked_callback);
}

/**
 * Iterates over the structure a breadth-first fashion in parallel.
 *
 * \tparam alphabet_t seqan3::alphabet type
 * \param sequence Sequence the tree is built over.
 * \param suffixes Suffix vector of the tree.
 * \param table Table vector of the tree.
 * \param expand_nodes Determines if nodes should be explored (construction) or
 * if we only want to visit nodes.
 * \param f Callback for the visited nodes.
 * \param parallel_depth Number of levels to spawn new processes for.
 * \param[in] locked_callback Callback for when the tree array has grown and a
 * user wants to do something in a synchronisation lock.
 */
template <seqan3::alphabet alphabet_t>
void breadth_first_iteration_parallel(
    int start_index, int start_lcp, const sequence_t<alphabet_t> &sequence,
    std::vector<int> &suffixes, Table<> &table, bool expand_nodes,
    const std::function<bool(int, int, int &, int)> &f,
    const std::function<void()> &done, int parallel_depth,
    const std::function<void(int, int, int &)> &locked_callback) {
  breadth_first_iteration_parallel_(start_index, start_lcp, 0, sequence,
                                    suffixes, table, expand_nodes, f, done,
                                    parallel_depth, locked_callback);
}

template <seqan3::alphabet alphabet_t>
void breadth_first_iteration_parallel_(
    int start_index, int start_lcp, int start_depth,
    const sequence_t<alphabet_t> &sequence, std::vector<int> &suffixes,
    Table<> &table, bool expand_nodes,
    const std::function<bool(int, int, int &, int)> &f,
    const std::function<void()> &done, int parallel_depth,
    const std::function<void(int, int, int &)> &locked_callback) {

  std::vector<std::thread> threads{};

  std::queue<std::tuple<int, int, int>> queue{};
  queue.emplace(start_index, start_lcp, start_depth);

  while (!queue.empty()) {
    auto [node_index, lcp, depth] = queue.front();
    queue.pop();

    auto [new_lcp, consider_children] =
        visit_top_node(node_index, lcp, sequence, suffixes, table, expand_nodes,
                       f, locked_callback);

    if (!consider_children) {
      continue;
    }

    if (depth < parallel_depth) {
      // Spawn threads per children
      iterate_children(
          node_index, table, [&, new_lcp = new_lcp, depth = depth](int index) {
            threads.emplace_back(breadth_first_iteration_parallel_<alphabet_t>,
                                 index, new_lcp, depth + 1, std::ref(sequence),
                                 std::ref(suffixes), std::ref(table),
                                 expand_nodes, f, done, parallel_depth,
                                 locked_callback);
          });

    } else {
      // Continue working in this thread.
      iterate_children(node_index, table,
                       [&, new_lcp = new_lcp, depth = depth](int index) {
                         queue.emplace(index, new_lcp, depth + 1);
                       });
    }
  }

  done();

  for (auto &thread : threads) {
    if (thread.joinable()) {
      thread.join();
    }
  }
}

using namespace std::placeholders;

template <seqan3::alphabet alphabet_t>
std::tuple<int, bool>
visit_top_node(int node_index, int lcp, const sequence_t<alphabet_t> &sequence,
               std::vector<int> &suffixes, Table<> &table, bool expand_nodes,
               const std::function<bool(int, int, int &, int)> &f,
               const std::function<void(int, int, int &)> &locked_callback) {

  if (node_index == 0) {
    return {get_edge_lcp(node_index, sequence, suffixes, table), true};
  }

  auto modified_locked_callback = std::bind(locked_callback, _1, lcp, _2);

  int edge_lcp, node_count;
  if (is_unevaluated(node_index, table) && expand_nodes) {
    auto [edge_lcp_, node_count_] = expand_node(
        node_index, sequence, suffixes, table, modified_locked_callback);
    edge_lcp = edge_lcp_;
    node_count = node_count_;
  } else {
    edge_lcp = get_edge_lcp(node_index, sequence, suffixes, table);
    node_count = node_occurrences(node_index, table);
  }

  bool consider_children = f(node_index, lcp, edge_lcp, node_count);
  if (!consider_children) {
    return {-1, false};
  }

  // It is possible that the call to f expands implicit nodes, may need to
  // recalculate the edge_lcp.
  // int old_edge_lcp = edge_lcp;
  //  edge_lcp = get_edge_lcp(node_index, sequence, suffixes, table);

  int new_lcp = lcp + edge_lcp;
  return {new_lcp, true};
}

template <seqan3::alphabet alphabet_t>
int get_edge_lcp(int node_index, const sequence_t<alphabet_t> &sequence,
                 const std::vector<int> &suffixes, const Table<> &table) {
  if (node_index == 0) {
    return 0;
  }

  if (is_leaf(node_index, table)) {
    return suffixes.size() - table[node_index].value;
  }

  if (is_unevaluated(node_index, table)) {
    return lst::details::longest_common_prefix(table[node_index].value,
                                               table[node_index + 1].value,
                                               sequence, suffixes);
  }

  int smallest_child_index = suffixes.size();

  iterate_children(node_index, table, [&](int index) {
    int sequence_index = get_sequence_index(index, suffixes, table);
    smallest_child_index = std::min(smallest_child_index, sequence_index);
  });

  assert(smallest_child_index > table[node_index].value);

  int edge_lcp = smallest_child_index - table[node_index].value;

  return edge_lcp;
}

} // namespace lst::details
