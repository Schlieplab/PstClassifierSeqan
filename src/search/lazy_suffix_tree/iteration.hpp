#pragma once

#include <functional>
#include <queue>
#include <vector>
#include <mutex>
#include <thread>
#include <seqan3/alphabet/all.hpp>

#include "construction.hpp"

std::mutex expand_lock;
std::mutex child_lock;
namespace lst::details {

int next_child_index(int node_index, std::vector<Flag> &flags) {
//if (node_index < biggest_index){
 //   node_index = biggest_index;
  //}
  if (is_leaf(node_index, flags)) {
    // Should be 1, but I've addded a value to leaves to allow for
    // explicit nodes.
    node_index += 2;

  } else {
    node_index += 2;
  }

  //biggest_index = node_index;
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
bool label_valid(int label_start, int label_end) {
  return true;
}
bool include_node(int label_start, int label_end, int count) {
  int label_length = label_end - label_start;

  return label_length < 15 && count >= 100 &&
         label_valid(label_start, label_end);
}

void wrapper(auto &sequence,
             std::vector<uint> &suffixes,
             std::vector<int> &table, std::vector<Flag> &flags,
             bool expand_nodes, int start_node) {
  BFIp(&sequence, &suffixes, &table, &flags,
        expand_nodes, start_node);
}

template <seqan3::Alphabet alphabet_t>
void BFIp(sequence_t<alphabet_t> &sequence,
                                      std::vector<uint> &suffixes,
                                      std::vector<int> &table, std::vector<Flag> &flags,
                                      bool expand_nodes,  const std::function<bool(int, int, int)> &f,
                                      int start_node, int initial_lcp){

  std::queue<std::tuple<int, int>> queue{};
  iterate_children(start_node, table, flags,
                   [&](int index) { queue.emplace(index, initial_lcp); });

  while (!queue.empty()) {
    auto [node_index, lcp] = queue.front();
    queue.pop();

    int edge_lcp;
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

    int new_lcp = lcp + edge_lcp;

    iterate_children(node_index, table, flags,
                     [&](int index) { queue.emplace(index, new_lcp); });
  }


}


template <seqan3::Alphabet alphabet_t>
void breadth_first_iteration(sequence_t<alphabet_t> &sequence,
                             std::vector<uint> &suffixes,
                             std::vector<int> &table, std::vector<Flag> &flags,
                             const std::function<bool(int, int, int)> &f) {
  breadth_first_iteration(sequence, suffixes, table, flags, false, f, false, 0);
}



template <seqan3::Alphabet alphabet_t>
void breadth_first_iteration(sequence_t<alphabet_t> &sequence,    
                             std::vector<uint> &suffixes,
                             std::vector<int> &table, std::vector<Flag> &flags,
                             bool expand_nodes,
                             const std::function<bool(int, int, int)> &f, bool build, int paralell_depth) {


  std::queue<std::tuple<int, int>> queue{};
  iterate_children(0, table, flags,
                   [&](int index) { queue.emplace(index, 0); });

   // Calculate lower and upper bound for nodes index counting
  int node_lower_bound = 2;
  if (paralell_depth > 1){
    node_lower_bound = pow(5,(paralell_depth-1))*2+2;
  }
  int node_upper_bound = pow(5,paralell_depth)*2;

  while (!queue.empty()) {
    auto [node_index, lcp] = queue.front();

    // Enter if all serial nodes are expanded.
    if (node_index > node_upper_bound && build){
      // Spawn threads
      int number_of_threads = 4;
      if (paralell_depth > 1){
        number_of_threads = 0;
        for (int j = node_lower_bound; j < node_upper_bound+1; j+=2) {
          if (sequence[table[j]].to_rank() != 3 && sequence[table[j]].to_rank() != 5) {
            number_of_threads++;
          }
        }
      }
      std::thread threads[number_of_threads];
      int thread_index = 0;
      for (int j = node_lower_bound; j < node_upper_bound+1; j+=2) {
        // Filter out all nodes with N as label
        if (sequence[table[j]].to_rank() != 3 && sequence[table[j]].to_rank() != 5) {
          seqan3::debug_stream << "Node_ID: " << j << " |Table[Node_id]: " << table[j]<< " |Sequence[Table[node_id]]: " << sequence[table[j]] << " " << sequence[table[j]].to_rank() << std::endl;
          int edge_lcp = get_edge_lcp(j, sequence, suffixes, table, flags);
          threads[thread_index] = std::thread(BFIp<seqan3::dna5>, std::ref(sequence), std::ref(suffixes),
                                              std::ref(table), std::ref(flags), expand_nodes, f, j, edge_lcp);
          thread_index++;
        }
      }
      for (int i = 0; i < number_of_threads; ++i) {
        seqan3::debug_stream << i << std::endl;
        threads[i].join();
      }
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
                 std::vector<uint> &suffixes, std::vector<int> &table,
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
