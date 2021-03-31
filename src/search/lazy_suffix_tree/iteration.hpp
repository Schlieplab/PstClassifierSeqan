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

size_t next_child_index(size_t node_index, const Table<> &table) {
  return node_index + 1;
}

void iterate_children(size_t node_index, const Table<> &table,
                      const std::function<void(size_t)> &f) {
  if (is_leaf(node_index, table) || is_unevaluated(node_index, table)) {
    return;
  }

  auto first_child = table[node_index].second;

  for (size_t i = first_child; i <= table.size();) {
    f(i);

    if (is_rightmostchild(i, table)) {
      break;
    }
    i = next_child_index(i, table);
  }
}

template <seqan3::alphabet alphabet_t>
size_t node_occurrences(size_t node_index, const Table<> &table,
                        const sequence_t<alphabet_t> &sequence,
                        std::vector<size_t> &suffixes) {
  assert(node_index <= table.size());

  size_t occurrences = 0;
  std::queue<size_t> queue{};

  queue.push(node_index);

  while (!queue.empty()) {
    size_t index = queue.front();
    queue.pop();

    if (is_leaf(index, table)) {
      occurrences += 1;
    } else if (is_unevaluated(index, table)) {
      auto lower_bound = table[index].first;
      auto upper_bound = table[index].second;

      auto [counts, count] =
          count_suffixes(lower_bound, upper_bound, sequence, suffixes);

      occurrences += count;
    } else {
      iterate_children(index, table, [&](size_t i) { queue.push(i); });
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
    size_t start_lower_bound, size_t start_upper_bound, size_t start_lcp,
    int level, int parallel_depth, const sequence_t<alphabet_t> &sequence,
    std::vector<size_t> &suffixes,
    const std::function<bool(size_t, size_t, size_t, size_t,
                             alphabet_array<size_t, alphabet_t> &, bool)> &f,
    const std::function<void()> &done) {
  std::vector<std::thread> threads{};

  std::queue<std::tuple<size_t, size_t, size_t>> queue{};
  queue.emplace(start_lower_bound, start_upper_bound, start_lcp);

  while (!queue.empty()) {
    auto [lower, upper, lcp] = queue.front();
    queue.pop();

    if (upper == (size_t)-1) {
      alphabet_array<size_t, alphabet_t> no_counts{};
      f(lower, lcp, sequence.size() - lower + 1, 1, no_counts, true);
      continue;
    }

    auto [sequence_index, edge_lcp, node_count, child_counts, children] =
        evaluate_node(lower, upper, sequence, suffixes);

    bool consider_children =
        f(sequence_index, lcp, edge_lcp, node_count, child_counts, false);

    if (!consider_children) {
      continue;
    }
    size_t new_lcp = lcp + edge_lcp;

    if (level < parallel_depth) {
      for (auto &[child_lower, child_upper] : children) {
        threads.emplace_back(breadth_first_iteration_table_less_<alphabet_t>,
                             child_lower, child_upper, new_lcp, level + 1,
                             parallel_depth, std::ref(sequence),
                             std::ref(suffixes), f, done);
      }
    } else {
      for (auto &[child_lower, child_upper] : children) {
        queue.emplace(child_lower, child_upper, new_lcp);
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
    std::vector<size_t> &suffixes,
    const std::function<bool(size_t, size_t, size_t, size_t,
                             alphabet_array<size_t, alphabet_t> &, bool)> &f,
    const std::function<void()> &done) {
  breadth_first_iteration_table_less_<alphabet_t>(
      0, suffixes.size(), 0, 0, parallel_depth, sequence, suffixes, f, done);
}

/**
 * Iterates over the tree in a breadth-first fashion without constructing the
 * tree, but starts with an root-expanded tree.
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
    int parallel_depth, const Table<> &table,
    const sequence_t<alphabet_t> &sequence, std::vector<size_t> &suffixes,
    const std::function<bool(size_t, size_t, size_t, size_t,
                             alphabet_array<size_t, alphabet_t> &, bool)> &f,
    const std::function<void()> &done) {
  if (parallel_depth > 0) {
    std::vector<std::thread> threads{};
    iterate_children(0, table, [&](size_t index) {
      size_t lower_bound, upper_bound;

      if (is_unevaluated(index, table)) {
        lower_bound = table[index].first;
        upper_bound = table[index].second;
      } else if (is_leaf(index, table)) {
        lower_bound = table[index].first;
        upper_bound = (size_t)-1;
      }

      threads.emplace_back(breadth_first_iteration_table_less_<alphabet_t>,
                           lower_bound, upper_bound, 0, 1, parallel_depth,
                           std::ref(sequence), std::ref(suffixes), f, done);
    });

    for (auto &thread : threads) {
      if (thread.joinable()) {
        thread.join();
      }
    }

  } else {
    iterate_children(0, table, [&](size_t index) {
      size_t lower_bound, upper_bound;

      if (is_unevaluated(index, table)) {
        lower_bound = table[index].first;
        upper_bound = table[index].second;
      } else if (is_leaf(index, table)) {
        lower_bound = table[index].first;
        upper_bound = -1;
      }

      breadth_first_iteration_table_less_<alphabet_t>(
          lower_bound, upper_bound, 0, 1, parallel_depth, sequence, suffixes, f,
          []() {});
    });
  }

  done();
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
    const sequence_t<alphabet_t> &sequence, std::vector<size_t> &suffixes,
    Table<> &table, bool expand_nodes,
    const std::function<bool(size_t, size_t, size_t &, size_t,
                             alphabet_array<size_t, alphabet_t> &)> &f,
    const std::function<void(size_t, size_t, size_t &)> &locked_callback) {
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
    size_t start_index, size_t start_lcp,
    const sequence_t<alphabet_t> &sequence, std::vector<size_t> &suffixes,
    Table<> &table, bool expand_nodes,
    const std::function<bool(size_t, size_t, size_t &, size_t,
                             alphabet_array<size_t, alphabet_t> &)> &f,
    const std::function<void(size_t, size_t, size_t &)> &locked_callback) {
  std::queue<std::tuple<size_t, size_t>> queue{};
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

    iterate_children(node_index, table, [&, new_lcp = new_lcp](size_t index) {
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
    const sequence_t<alphabet_t> &sequence, std::vector<size_t> &suffixes,
    Table<> &table, bool expand_nodes,
    const std::function<bool(size_t, size_t, size_t &, size_t,
                             alphabet_array<size_t, alphabet_t> &)> &f,
    const std::function<void()> &done, size_t parallel_depth,
    const std::function<void(size_t, size_t, size_t &)> &locked_callback) {
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
    size_t start_index, size_t start_lcp,
    const sequence_t<alphabet_t> &sequence, std::vector<size_t> &suffixes,
    Table<> &table, bool expand_nodes,
    const std::function<bool(size_t, size_t, size_t &, size_t,
                             alphabet_array<size_t, alphabet_t> &)> &f,
    const std::function<void()> &done, size_t parallel_depth,
    const std::function<void(size_t, size_t, size_t &)> &locked_callback) {
  breadth_first_iteration_parallel_(start_index, start_lcp, 0, sequence,
                                    suffixes, table, expand_nodes, f, done,
                                    parallel_depth, locked_callback);
}

template <seqan3::alphabet alphabet_t>
void breadth_first_iteration_parallel_(
    size_t start_index, size_t start_lcp, int depth,
    const sequence_t<alphabet_t> &sequence, std::vector<size_t> &suffixes,
    Table<> &table, bool expand_nodes,
    const std::function<bool(size_t, size_t, size_t &, size_t,
                             alphabet_array<size_t, alphabet_t> &)> &f,
    const std::function<void()> &done, int parallel_depth,
    const std::function<void(size_t, size_t, size_t &)> &locked_callback) {

  std::vector<std::thread> threads{};

  std::queue<std::tuple<size_t, size_t>> queue{};
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

    if (depth < parallel_depth) {
      // Spawn threads per children
      iterate_children(node_index, table, [&, new_lcp = new_lcp](size_t index) {
        threads.emplace_back(breadth_first_iteration_parallel_<alphabet_t>,
                             index, new_lcp, depth + 1, std::ref(sequence),
                             std::ref(suffixes), std::ref(table), expand_nodes,
                             f, done, parallel_depth, locked_callback);
      });

    } else {
      // Continue working in this thread.
      iterate_children(node_index, table, [&, new_lcp = new_lcp](size_t index) {
        queue.emplace(index, new_lcp);
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

template <seqan3::alphabet alphabet_t>
std::tuple<size_t, bool> visit_top_node(
    size_t node_index, size_t lcp, const sequence_t<alphabet_t> &sequence,
    std::vector<size_t> &suffixes, Table<> &table, bool expand_nodes,
    const std::function<bool(size_t, size_t, size_t &, size_t,
                             alphabet_array<size_t, alphabet_t> &)> &f,
    const std::function<void(size_t, size_t, size_t &)> &locked_callback) {

  if (node_index == 0) {
    return {get_edge_lcp(node_index, sequence, suffixes, table), true};
  }

  auto modified_locked_callback = [&](size_t node_index, size_t &edge_lcp) {
    return locked_callback(node_index, lcp, edge_lcp);
  };

  size_t edge_lcp, node_count;
  alphabet_array<size_t, alphabet_t> child_counts;
  if (is_unevaluated(node_index, table) && expand_nodes) {
    auto [edge_lcp_, node_count_, child_counts_] = expand_node(
        node_index, sequence, suffixes, table, modified_locked_callback);
    edge_lcp = edge_lcp_;
    node_count = node_count_;
    child_counts = child_counts_;
  } else {
    edge_lcp = get_edge_lcp(node_index, sequence, suffixes, table);
    node_count = (size_t)-1; // node_occurrences(node_index, table);
    child_counts = {};
  }

  bool consider_children =
      f(node_index, lcp, edge_lcp, node_count, child_counts);
  if (!consider_children) {
    return {-1, false};
  }

  size_t new_lcp = lcp + edge_lcp;
  return {new_lcp, true};
}

template <seqan3::alphabet alphabet_t>
size_t get_edge_lcp(size_t node_index, const sequence_t<alphabet_t> &sequence,
                    const std::vector<size_t> &suffixes, const Table<> &table) {
  if (node_index == 0) {
    return 0;
  }

  if (is_leaf(node_index, table)) {
    return suffixes.size() - table[node_index].first;
  }

  if (is_unevaluated(node_index, table)) {
    return lst::details::longest_common_prefix(
        table[node_index].first, table[node_index].second, sequence, suffixes);
  }

  auto smallest_child_index = suffixes.size();

  iterate_children(node_index, table, [&](size_t index) {
    auto sequence_index = get_sequence_index(index, suffixes, table);
    smallest_child_index = std::min(smallest_child_index, sequence_index);
  });

  assert(smallest_child_index > table[node_index].first);

  size_t edge_lcp = smallest_child_index - table[node_index].first;

  return edge_lcp;
}

} // namespace lst::details
