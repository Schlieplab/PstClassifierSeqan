#pragma once

#include <functional>
#include <queue>
#include <vector>

#include <seqan3/alphabet/all.hpp>

#include "construction.hpp"

int biggest_index = 0;
namespace lst::details {

int next_child_index(int node_index, std::vector<Flag> &flags) {
  if (node_index < biggest_index){
    node_index = biggest_index;
  }
  if (is_leaf(node_index, flags)) {
    // Should be 1, but I've addded a value to leaves to allow for
    // explicit nodes.
    node_index += 2;

  } else {
    node_index += 2;
  }

  biggest_index = node_index;
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

template <seqan3::Alphabet alphabet_t>
void breadth_first_iteration(sequence_t<alphabet_t> &sequence,
                             std::vector<int> &suffixes,
                             std::vector<int> &table, std::vector<Flag> &flags,
                             const std::function<bool(int, int, int)> &f) {
  breadth_first_iteration(sequence, suffixes, table, flags, true, f);
}

template <seqan3::Alphabet alphabet_t>
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
    if (node_index > 8 ){
      std::vector<int> aSuffix;
      std::vector<int> cSuffix;
      std::vector<int> gSuffix;
      std::vector<int> tSuffix;

      for (int j : suffixes) {
        int letterCode = sequence[j].to_rank();
        if (letterCode == 0){
          aSuffix.push_back(j+1);
        }
        else if (letterCode == 1){
          cSuffix.push_back(j+1);
        }
        else if (letterCode == 2){
          gSuffix.push_back(j+1);
        }
        else if (letterCode == 3){
          tSuffix.push_back(j+1);
        }
      }
      for (int j : cSuffix) {
        seqan3::debug_stream << suffixes[j] << std::endl;
      }
      breadth_first_iteration_paralell(sequence, aSuffix, table, flags, true, 2, f);
      breadth_first_iteration_paralell(sequence, cSuffix, table, flags, true, 4, f);
      breadth_first_iteration_paralell(sequence, gSuffix, table, flags, true, 6, f);
      breadth_first_iteration_paralell(sequence, tSuffix, table, flags, true, 8, f);

      break;
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

    if (!consider_children) {
      continue;
    }

    int new_lcp = lcp + edge_lcp;
    iterate_children(node_index, table, flags,
                     [&](int index) { queue.emplace(index, new_lcp); });
  }
}

template <seqan3::Alphabet alphabet_t>
void breadth_first_iteration_paralell(sequence_t<alphabet_t> &sequence,
                             std::vector<int> &suffixes,
                             std::vector<int> &table, std::vector<Flag> &flags,
                             bool expand_nodes, int start_node,
                             const std::function<bool(int, int, int)> &f){

    std::queue<std::tuple<int, int>> queue{};
    iterate_children(start_node, table, flags,
                     [&](int index) { queue.emplace(index, start_node); });

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


template <seqan3::Alphabet alphabet_t>
int get_edge_lcp(int node_index, sequence_t<alphabet_t> &sequence,
                 std::vector<int> &suffixes, std::vector<int> &table,
                 std::vector<Flag> &flags) {
  if (is_leaf(node_index, flags)) {
    return sequence.size() - table[node_index];
  }

  if (is_unevaluated(node_index, flags)) {
    return lst::details::longest_common_prefix(
        table[node_index], table[node_index + 1], sequence, suffixes);
  }

  int smallest_child_index = sequence.size();

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
