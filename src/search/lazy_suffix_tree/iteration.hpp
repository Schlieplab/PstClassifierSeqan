#pragma once

#include <functional>
#include <mutex>
#include <queue>
#include <seqan3/alphabet/concept.hpp>
#include <thread>
#include <vector>

#include "construction.hpp"

std::mutex expand_lock;
std::mutex child_lock;
namespace lst::details {

int next_child_index(int node_index, std::vector<Flag> &flags) {
  if (is_leaf(node_index, flags)) {
    // Should be 1, but I've added a value to leaves to allow for
    // explicit nodes.
    return node_index + 2;
  } else {
    return node_index + 2;
  }
}

void iterate_children(int node_index, std::vector<int> &table,
                      std::vector<Flag> &flags,
                      const std::function<void(int)> &f) {
  if (is_leaf(node_index, flags) || is_unevaluated(node_index, flags)) {
    return;
  }

  int first_child = table[node_index + 1];

  for (int i = first_child; i <= table.size();) {
    f(i);

    if (is_rightmostchild(i, flags)) {
      break;
    }
    i = next_child_index(i, flags);
  }
}

int node_occurrences(int node_index, std::vector<int> &table,
                         std::vector<Flag> &flags) {
  if (node_index > table.size()) {
    throw std::invalid_argument(
        "[NODE OCCURRENCES] Given node index is too large.");
  }

  int occurrences = 0;
  std::queue<int> queue{};

  queue.push(node_index);

  while (!queue.empty()) {
    int index = queue.front();
    queue.pop();

    if (is_leaf(index, flags)) {
      occurrences += 1;
    } else if (is_unevaluated(index, flags)) {
      occurrences += table[index + 1] - table[index];
    } else {
      iterate_children(index, table, flags, [&](int i) { queue.push(i); });
    }
  }

  return occurrences;
}

template <seqan3::alphabet alphabet_t>
void BFIp(sequence_t<alphabet_t> &sequence, std::vector<int> &suffixes,
          std::vector<int> &table, std::vector<Flag> &flags,
          bool expand_nodes,
          const std::function<bool(int, int, int)> &f,
          int start_node, int initial_lcp) {

  std::queue<std::tuple<int, int>> queue{};
  iterate_children(start_node, table, flags,
                   [&](int index) { queue.emplace(index, initial_lcp); });

  while (!queue.empty()) {
    auto [node_index, lcp] = queue.front();
    queue.pop();

    int edge_lcp;
    // expand_lock.lock();
    if (is_unevaluated(node_index, flags) && expand_nodes) {
      edge_lcp = lst::details::expand_node(node_index, sequence, suffixes,
                                           table, flags);
    } else {
      edge_lcp = get_edge_lcp(node_index, sequence, suffixes, table, flags);
    }
    // expand_lock.unlock();

    child_lock.lock();
    bool consider_children = f(node_index, lcp, edge_lcp);
    child_lock.unlock();
    if (!consider_children) {
      continue;
    }

    int new_lcp = lcp + edge_lcp;

    iterate_children(node_index, table, flags,
                     [&](int index) { queue.emplace(index, new_lcp); });
  }
}

template <seqan3::alphabet alphabet_t>
void breadth_first_iteration(
    sequence_t<alphabet_t> &sequence, std::vector<int> &suffixes,
    std::vector<int> &table, std::vector<Flag> &flags, bool expand_nodes,
    const std::function<bool(int, int, int)> &f) {
  breadth_first_iteration(0, 0, sequence, suffixes, table, flags, expand_nodes,
                          f, false, 0);
}

template <seqan3::alphabet alphabet_t>
void breadth_first_iteration(
    int start_index, int start_lcp, sequence_t<alphabet_t> &sequence,
    std::vector<int> &suffixes, std::vector<int> &table,
    std::vector<Flag> &flags, bool expand_nodes,
    const std::function<bool(int, int, int)> &f, bool multi_core,
    int parallel_depth) {
  // TODO Refactor this function.

  std::queue<std::tuple<int, int>> queue{};
  iterate_children(start_index, table, flags,
                   [&](int index) { queue.emplace(index, start_lcp); });

  // Calculate lower and upper bound for nodes index counting
  // TODO after merge: can probably skip the if statement with a simple std::max(parallel_depth - 1, 1)
  int node_lower_bound = 2;
  int node_upper_bound = pow(5, parallel_depth) * 2;

  if (parallel_depth > 1) {
    node_lower_bound = pow(5, (parallel_depth - 1)) * 2 + 2;
  }

  while (!queue.empty()) {
    auto [node_index, lcp] = queue.front();
    // Enter if all serial nodes are expanded.
    if (node_index > node_upper_bound && multi_core) {
      int number_of_threads = 4;
      int thread_index = 0;
      // Spawn threads
      if (parallel_depth > 1) {
        number_of_threads = 0;
        for (int j = node_lower_bound; j < node_upper_bound + 1; j += 2) {
          // Check if the character is a gap ('_') or an 'N'.
          // TODO replace with alphabet_t('N').to_rank() or some such.
          if (sequence[table[j]].to_rank() != 3 &&
              sequence[table[j]].to_rank() != 5) {
            number_of_threads++;
          }
        }
      }
      seqan3::debug_stream << "Started " << number_of_threads << "threads."
                           << std::endl;
      std::thread threads[number_of_threads];
      for (int j = node_lower_bound; j < node_upper_bound + 1; j += 2) {
        // Filter out all nodes with N as label
        // Check if the character is a gap ('_') or an 'N'.
        // TODO replace with alphabet_t('N').to_rank() or some such.
        if (sequence[table[j]].to_rank() != 3 &&
            sequence[table[j]].to_rank() != 5) {
          int edge_lcp = get_edge_lcp(j, sequence, suffixes, table, flags);
          threads[thread_index] = std::thread(
              BFIp<alphabet_t>, std::ref(sequence), std::ref(suffixes),
              std::ref(table), std::ref(flags), expand_nodes, std::ref(f), j, edge_lcp);
          thread_index++;
        }
      }

      for (int i = 0; i < number_of_threads; ++i) {
        seqan3::debug_stream << i << std::endl;
        threads[i].join();
      }
      seqan3::debug_stream << "All threads returned." << std::endl;
      return;
    }
    queue.pop();

    int edge_lcp;

    if (is_unevaluated(node_index, flags) && expand_nodes) {
      edge_lcp = lst::details::expand_node(node_index, sequence, suffixes,
                                           table, flags);
    } else {
      edge_lcp = get_edge_lcp(node_index, sequence, suffixes, table, flags);
    }

    bool consider_children = f(node_index, lcp, edge_lcp);

    // It is possible that the call to f expands implicit nodes, need to
    // recalculate the edge_lcp.
    edge_lcp = get_edge_lcp(node_index, sequence, suffixes, table, flags);

    if (!consider_children) {
      continue;
    }

    int new_lcp = lcp + edge_lcp;
    iterate_children(node_index, table, flags,
                     [&](int index) { queue.emplace(index, new_lcp); });
  }
}

template <seqan3::alphabet alphabet_t>
int get_edge_lcp(int node_index, sequence_t<alphabet_t> &sequence,
                 std::vector<int> &suffixes, std::vector<int> &table,
                 std::vector<Flag> &flags) {
  if (is_leaf(node_index, flags)) {
    return suffixes.size() - table[node_index];
  }

  if (is_unevaluated(node_index, flags)) {
    return lst::details::longest_common_prefix(
        table[node_index], table[node_index + 1], sequence, suffixes);
  }

  int smallest_child_index = suffixes.size();

  iterate_children(node_index, table, flags, [&](int index) {
    int table_index = table[index];

    if (is_unevaluated(index, flags)) {
      table_index = suffixes[table[index]];
    }

    if (table_index < smallest_child_index) {
      smallest_child_index = table_index;
    }
  });

  return smallest_child_index - table[node_index];
}

} // namespace lst::details
