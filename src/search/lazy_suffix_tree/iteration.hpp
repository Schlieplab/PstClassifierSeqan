#pragma once

#include <functional>
#include <queue>
#include <vector>
#include <mutex>
#include <thread>
#include <seqan3/alphabet/concept.hpp>

#include "construction.hpp"

std::mutex expand_lock;
std::mutex child_lock;
namespace lst::details {

int64_t next_child_index(int64_t node_index,
                              std::vector<Flag> &flags) {

  if (is_leaf(node_index, flags)) {
    // Should be 1, but I've addded a value to leaves to allow for
    // explicit nodes.
    node_index += 2;

  } else {
    node_index += 2;
  }
  return node_index;
}

void iterate_children(int64_t node_index, std::vector<int64_t> &table,
                      std::vector<Flag> &flags,
                      const std::function<void(int64_t)> &f) {
  if (is_leaf(node_index, flags) || is_unevaluated(node_index, flags)) {
    return;
  }

  int64_t first_child = table[node_index + 1];

  for (int64_t i = first_child; i <= table.size();) {
    f(i);

    if (is_rightmostchild(i, flags)) {
      break;
    }
    i = next_child_index(i, flags);
  }
}

int64_t node_occurrences(int64_t node_index,
                              std::vector<int64_t> &table,
                     std::vector<Flag> &flags) {
  if (node_index > table.size()) {
    throw std::invalid_argument(
        "[NODE OCCURRENCES] Given node index is too large.");
  }

  int64_t occurrences = 0;
  std::queue<int64_t> queue{};

  queue.push(node_index);

  while (!queue.empty()) {
    int64_t index = queue.front();
    queue.pop();

    if (is_leaf(index, flags)) {
      occurrences += 1;
    } else if (is_unevaluated(index, flags)) {
      occurrences += table[index + 1] - table[index];
    } else {
      iterate_children(index, table, flags,
        [&](int64_t i) {
          queue.push(i);
        });
    }
  }
  return occurrences;
}
bool label_valid(int64_t label_start, int64_t label_end) {
  return true;
}
bool include_node(int64_t label_start, int64_t label_end,
                  int64_t count) {
  int64_t label_length = label_end - label_start;

  return label_length < 15 && count >= 100 &&
         label_valid(label_start, label_end);
}


template <seqan3::alphabet alphabet_t>
void BFIp(sequence_t<alphabet_t> &sequence,
          std::vector<int64_t> &suffixes,
          std::vector<int64_t> &table,
          std::vector<Flag> &flags,
          bool expand_nodes,
          const std::function<bool(int64_t, int64_t, int64_t)> &f,
          int64_t start_node, int64_t initial_lcp){

  std::queue<std::tuple<int64_t, int64_t>> queue{};
  iterate_children(start_node, table, flags,
                [&](int64_t index) { queue.emplace(index, initial_lcp); });

  while (!queue.empty()) {
    auto [node_index, lcp] = queue.front();
    queue.pop();

    int64_t edge_lcp;
    //expand_lock.lock();
    if (is_unevaluated(node_index, flags) && expand_nodes) {
      edge_lcp = lst::details::expand_node(node_index, sequence, suffixes,
                                           table, flags);
    } else {
      edge_lcp = get_edge_lcp(node_index, sequence, suffixes, table, flags);
    }
    //expand_lock.unlock();

    child_lock.lock();
    bool consider_children = f(node_index, lcp, edge_lcp);
    child_lock.unlock();
    if (!consider_children) {
      continue;
    }

    int64_t new_lcp = lcp + edge_lcp;

    iterate_children(node_index, table, flags,
                   [&](int64_t index) { queue.emplace(index, new_lcp); });
  }


}

template <seqan3::alphabet alphabet_t>
void breadth_first_iteration(sequence_t<alphabet_t> &sequence,
     std::vector<int64_t> &suffixes,
     std::vector<int64_t> &table,
     std::vector<Flag> &flags, bool expand_nodes,
     const std::function<bool(int64_t, int64_t, int64_t)> &f) {
        breadth_first_iteration(0, 0, sequence, suffixes, table, flags,
                                expand_nodes, f, false, 0);
}

template <seqan3::alphabet alphabet_t>
void breadth_first_iteration(int64_t start_index, int64_t start_lcp,
       sequence_t<alphabet_t> &sequence,
       std::vector<int64_t> &suffixes,
       std::vector<int64_t> &table,
       std::vector<Flag> &flags, bool expand_nodes,
       const std::function<bool(int64_t, int64_t, int64_t)> &f,
       bool &multi_core, int paralell_depth) {


  std::queue<std::tuple<int64_t, int64_t>> queue{};
  iterate_children(start_index, table, flags,
                   [&](int64_t index) { queue.emplace(index, start_lcp); });

   // Calculate lower and upper bound for nodes index counting
  int node_lower_bound  = 2;
  int node_upper_bound  = pow(5,paralell_depth)*2;


  if (paralell_depth > 1){
    node_lower_bound = pow(5,(paralell_depth-1))*2+2;
  }

  while (!queue.empty()) {
    auto [node_index, lcp] = queue.front();
    // Enter if all serial nodes are expanded.
    if (node_index > node_upper_bound && multi_core){
      int number_of_threads = 4;
      int thread_index      = 0;
      // Spawn threads
      if (paralell_depth > 1){
        number_of_threads = 0;
        for (int j = node_lower_bound; j < node_upper_bound+1; j+=2) {
          if (sequence[table[j]].to_rank() != 3 && sequence[table[j]].to_rank() != 5) {
            number_of_threads++;
          }
        }
      }
      seqan3::debug_stream << "Started " << number_of_threads << "threads." << std::endl;
      std::thread threads[number_of_threads];
      for (int j = node_lower_bound; j < node_upper_bound+1; j+=2) {
        // Filter out all nodes with N as label
        if (sequence[table[j]].to_rank() != 3 && sequence[table[j]].to_rank() != 5) {
          int64_t edge_lcp = get_edge_lcp(j, sequence, suffixes, table, flags);
          threads[thread_index] = std::thread(BFIp<seqan3::dna5>,
                                              std::ref(sequence),
                                              std::ref(suffixes),
                                              std::ref(table),
                                              std::ref(flags),
                                              expand_nodes, f, j,
                                              edge_lcp);
          thread_index++;
        }
      }

      multi_core = false;
      for (int i = 0; i < number_of_threads; ++i) {
        seqan3::debug_stream << i << std::endl;
        threads[i].join();
      }
      seqan3::debug_stream << "All threads returned." << std::endl;
      multi_core = true;
      return;
    }
    queue.pop();

    int64_t edge_lcp;

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

    int64_t new_lcp = lcp + edge_lcp;
    iterate_children(node_index, table, flags,
                     [&](int64_t index) { queue.emplace(index, new_lcp); });
  }
}

template <seqan3::alphabet alphabet_t>
int64_t get_edge_lcp(int64_t node_index,
                          sequence_t<alphabet_t> &sequence,
                           std::vector<int64_t> &suffixes,
                           std::vector<int64_t> &table,
                           std::vector<Flag> &flags) {
  if (is_leaf(node_index, flags)) {
    return suffixes.size() - table[node_index];
  }

  if (is_unevaluated(node_index, flags)) {
    return lst::details::longest_common_prefix(
        table[node_index], table[node_index + 1], sequence, suffixes);
  }

  int64_t smallest_child_index = suffixes.size();

  iterate_children(node_index, table, flags, [&](int64_t index) {
    int64_t table_index = table[index];

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
