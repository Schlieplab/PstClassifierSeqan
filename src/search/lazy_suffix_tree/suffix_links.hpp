#pragma once

#include <queue>
#include <stack>
#include <vector>

#include "construction.hpp"
#include "iteration.hpp"

namespace lst::details {

template <seqan3::Alphabet alphabet_t>
int tree_height(std::vector<int> &depths, sequence_t<alphabet_t> &sequence,
                std::vector<int> &suffixes, std::vector<int> &table,
                std::vector<Flag> &flags) {
  std::queue<std::tuple<int, int>> queue{};
  queue.emplace(0, 0);

  int tree_height = 0;

  while (!queue.empty()) {
    auto [node_index, parent_depth] = queue.front();
    queue.pop();

    int node_depth = parent_depth +
                     get_edge_lcp(node_index, sequence, suffixes, table, flags);

    depths[node_index / 2] = node_depth;

    iterate_children(node_index, table, flags,
                     [&](int index) { queue.emplace(index, node_depth); });

    if (node_depth > tree_height) {
      tree_height = node_depth;
    }
  }

  return tree_height;
}

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

template <seqan3::Alphabet alphabet_t>
void add_explicit_suffix_links(sequence_t<alphabet_t> &sequence,
                               std::vector<int> &suffixes,
                               std::vector<int> &table,
                               std::vector<Flag> &flags,
                               std::vector<int> &suffix_links) {
  // This follows the algorithm described by Maaß in
  // doi:10.1016/j.ipl.2005.12.012
  // With modifications to allow for unevaluated nodes.

  std::vector<int> cause(sequence.size() + 1, -1);

  std::vector<int> leaf_indicies(sequence.size() + 1, -1);
  leaf_indicies[sequence.size()] = 0;

  prepare_suffix_links(cause, leaf_indicies, sequence, suffixes, table, flags);

  std::vector<int> depths(table.size() / 2, -1);
  int height = tree_height(depths, sequence, suffixes, table, flags);

  std::vector<int> branch(height, -1);

  compute_suffix_links(cause, branch, depths, leaf_indicies, sequence, suffixes,
                       table, flags, suffix_links);
}

template <seqan3::Alphabet alphabet_t>
void prepare_suffix_links(std::vector<int> &cause,
                          std::vector<int> &leaf_indicies,
                          sequence_t<alphabet_t> &sequence,
                          std::vector<int> &suffixes, std::vector<int> &table,
                          std::vector<Flag> &flags) {
  prepare_suffix_links(0, 0, cause, leaf_indicies, sequence, suffixes, table,
                       flags);
}

template <seqan3::Alphabet alphabet_t>
int prepare_suffix_links(int node_index, int lcp, std::vector<int> &cause,
                         std::vector<int> &leaf_indicies,
                         sequence_t<alphabet_t> &sequence,
                         std::vector<int> &suffixes, std::vector<int> &table,
                         std::vector<Flag> &flags) {
  if (is_leaf(node_index, flags)) {
    int leaf_index = get_leaf_index(node_index, lcp, suffixes, table, flags);

    leaf_indicies[leaf_index] = node_index;

    return leaf_index;
  } else if (is_unevaluated(node_index, flags)) {
    // This isn't technically correct.  The following gives the second smallest
    // suffix _overall_ and the algorithm calls for the second smallest suffix
    // from the immediate children.  Not entirely sure that it is important,
    // since the values will still be unique.

    int second_smallest_child = suffixes[table[node_index] + 1] - lcp;
    cause[second_smallest_child + 1] = node_index;

    return get_leaf_index(node_index, lcp, suffixes, table, flags);
  } else {
    int edge_lcp = get_edge_lcp(node_index, sequence, suffixes, table, flags);

    int smallest_child = sequence.size();
    int second_smallest_child = sequence.size();

    iterate_children(node_index, table, flags, [&](int index) {
      int child =
          prepare_suffix_links(index, lcp + edge_lcp, cause, leaf_indicies,
                               sequence, suffixes, table, flags);

      if (child < smallest_child) {
        second_smallest_child = smallest_child;
        smallest_child = child;
      } else if (child < second_smallest_child) {
        second_smallest_child = child;
      }
    });

    cause[second_smallest_child + 1] = node_index;

    return smallest_child;
  }
}

template <seqan3::Alphabet alphabet_t>
void compute_suffix_links(std::vector<int> &cause, std::vector<int> &branch,
                          std::vector<int> &depths,
                          std::vector<int> &leaf_indicies,
                          sequence_t<alphabet_t> &sequence,
                          std::vector<int> &suffixes, std::vector<int> &table,
                          std::vector<Flag> &flags,
                          std::vector<int> &suffix_links) {
  compute_suffix_links(0, 0, cause, branch, depths, leaf_indicies, sequence,
                       suffixes, table, flags, suffix_links);
}

void assign_link(int leaf_index, std::vector<int> &cause,
                 std::vector<int> &branch, std::vector<int> &depths,
                 std::vector<int> &suffix_links) {
  if (cause[leaf_index] != -1) {
    int caused = cause[leaf_index];

    suffix_links[caused / 2] = branch[depths[caused / 2] - 1];
  }
}

void assign_leaf_link(int node_index, int leaf_index,
                      std::vector<int> &suffix_links,
                      std::vector<int> &leaf_indicies) {
  if (leaf_indicies[leaf_index + 1] != -1) {
    suffix_links[node_index / 2] = leaf_indicies[leaf_index + 1];
  }
}

template <seqan3::Alphabet alphabet_t>
void compute_suffix_links(int node_index, int lcp, std::vector<int> &cause,
                          std::vector<int> &branch, std::vector<int> &depths,
                          std::vector<int> &leaf_indicies,
                          sequence_t<alphabet_t> &sequence,
                          std::vector<int> &suffixes, std::vector<int> &table,
                          std::vector<Flag> &flags,
                          std::vector<int> &suffix_links) {
  if (!is_leaf(node_index, flags)) {
    branch[depths[node_index / 2]] = node_index;
  }

  if (is_leaf(node_index, flags)) {
    int leaf_index = get_leaf_index(node_index, lcp, suffixes, table, flags);

    assign_link(leaf_index, cause, branch, depths, suffix_links);
    assign_leaf_link(node_index, leaf_index, suffix_links, leaf_indicies);
  }

  if (is_unevaluated(node_index, flags)) {
    // Maybe possible to stop this iteration early.
    for (int i = table[node_index]; i < table[node_index + 1]; i++) {
      int leaf_index = suffixes[i] - lcp;
      assign_link(leaf_index, cause, branch, depths, suffix_links);
    }
  }

  int edge_lcp = get_edge_lcp(node_index, sequence, suffixes, table, flags);

  iterate_children(node_index, table, flags, [&](int index) {
    compute_suffix_links(index, lcp + edge_lcp, cause, branch, depths,
                         leaf_indicies, sequence, suffixes, table, flags,
                         suffix_links);
  });
}

void prepare_implicit_suffix_links(
    std::vector<std::tuple<int, int>> &closest_suffix_link_destination,
    std::vector<int> &parent_links, std::vector<int> &table,
    std::vector<Flag> &flags, std::vector<int> &suffix_links) {

  std::stack<std::tuple<int, int>> stack{};
  stack.emplace(0, 0);

  std::vector<int> missing_suffix_links{};
  while (!stack.empty()) {
    auto [node_index, parent] = stack.top();
    stack.pop();

    parent_links[node_index / 2] = parent;

    if (!is_leaf(node_index, flags) && suffix_links[node_index / 2] == -1) {
      missing_suffix_links.push_back(node_index);
    } else if (missing_suffix_links.size() > 0) {
      int suffix_link_destination = suffix_links[node_index / 2];

      for (int i = 0; i < missing_suffix_links.size(); i++) {
        int missing = missing_suffix_links[i];

        closest_suffix_link_destination[missing] = std::make_tuple(
            suffix_link_destination, missing_suffix_links.size() - 1 - i);
      }

      missing_suffix_links.resize(0);
    }

    iterate_children(node_index, table, flags,
                     [&](int child) { stack.emplace(child, node_index); });
  }
}

void compute_implicit_suffix_links(
    std::vector<std::tuple<int, int>> &closest_suffix_link_destination,
    std::vector<int> &parent_links, std::vector<int> &table,
    std::vector<Flag> &flags, std::vector<int> &suffix_links) {
  std::queue<int> queue{};
  queue.push(0);

  while (!queue.empty()) {
    int node_index = queue.front();
    queue.pop();

    auto [suffix_link_destination, distance] =
        closest_suffix_link_destination[node_index];

    if (suffix_link_destination == -1) {
      suffix_links[node_index / 2] = -1;
    } else if (suffix_link_destination != 0) {
      int destination_parent = parent_links[suffix_link_destination / 2];
      for (int i = 0; i < distance; i++) {
        destination_parent = parent_links[destination_parent / 2];
      }

      suffix_links[node_index / 2] = destination_parent;
    }

    iterate_children(node_index, table, flags,
                     [&](int index) { queue.push(index); });
  }
}

void add_implicit_suffix_links(std::vector<int> &table,
                               std::vector<Flag> &flags,
                               std::vector<int> &suffix_links) {
  // Prepare the two vectors by storing the distance to the closest
  // child which has a suffix link, as well as all parents links.

  // When these vectors are filled, it allows us to assign the correct
  // suffix links of the implicit nodes by determining
  // which other node is at an equal distance to a child with a suffix link.
  // That node must therefore have a suffix link to the current node.

  std::vector<std::tuple<int, int>> closest_suffix_link_destination(
      table.size());

  std::vector<int> parent_links(table.size() / 2, -1);

  prepare_implicit_suffix_links(closest_suffix_link_destination, parent_links,
                                table, flags, suffix_links);
  compute_implicit_suffix_links(closest_suffix_link_destination, parent_links,
                                table, flags, suffix_links);
}
} // namespace lst::details
