#pragma once

#include <functional>
#include <queue>
#include <vector>

#include <seqan3/alphabet/concept.hpp>

#include "construction.hpp"

namespace lst::details {

int next_child_index(int node_index, std::vector<Flag> &flags) {
  if (is_leaf(node_index, flags)) {
    // Should be 1, but I've addded a value to leaves to allow for
    // explicit nodes.
    node_index += 2;

  } else {
    node_index += 2;
  }

  return node_index;
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
void breadth_first_iteration(sequence_t<alphabet_t> &sequence,
                             std::vector<int> &suffixes,
                             std::vector<int> &table, std::vector<Flag> &flags,
                             const std::function<bool(int, int, int)> &f) {
  breadth_first_iteration(sequence, suffixes, table, flags, true, f);
}

template <seqan3::alphabet alphabet_t>
void breadth_first_iteration(sequence_t<alphabet_t> &sequence,
                             std::vector<int> &suffixes,
                             std::vector<int> &table, std::vector<Flag> &flags,
                             bool expand_nodes,
                             const std::function<bool(int, int, int)> &f) {
  std::queue<std::tuple<int, int>> queue{};
  iterate_children(0, table, flags,
                   [&](int index) { queue.emplace(index, 0); });

  while (!queue.empty()) {
    auto [node_index, lcp] = queue.front();
    queue.pop();

    int edge_lcp;

    if (is_unevaluated(node_index, flags) && expand_nodes) {
      edge_lcp = lst::details::expand_node(node_index, sequence, suffixes,
                                           table, flags);
    } else {
      edge_lcp = get_edge_lcp(node_index, sequence, suffixes, table, flags);
    }

    bool consider_children = f(node_index, lcp, edge_lcp);

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
