#pragma once

#include <array>
#include <functional>
#include <numeric>
#include <queue>
#include <tuple>
#include <vector>

#include <seqan3/alphabet/all.hpp>
#include <seqan3/alphabet/composite/alphabet_variant.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/range/view/convert.hpp>

#include "lazy_suffix_tree_construction.hpp"

using namespace lst::details;

namespace lst {

template <seqan3::Alphabet alphabet_t = seqan3::dna5> class LazySuffixTree {

public:
  friend class LazySuffixTreeTest;

  LazySuffixTree(std::vector<alphabet_t> &sequence_) {
    sequence = sequence_ | seqan3::view::convert<seqan3::gapped<alphabet_t>>;
    sequence.push_back(seqan3::gap{});

    suffixes = std::vector<int>(sequence.size());
    std::iota(suffixes.begin(), suffixes.end(), 0);

    lst::details::expand_root(sequence, suffixes, table, flags);
  }

  void expand_all() {
    for (int i = 0; i < table.size(); i++) {
      if ((flags[i] & Flag::Unevaluated) == Flag::Unevaluated) {
        lst::details::expand_node(i, sequence, suffixes, table, flags);
      }
    }
  }

  std::vector<std::tuple<sequence_t<alphabet_t>, int>> get_all_labels() {
    std::vector<std::tuple<sequence_t<alphabet_t>, int>> labels{};

    breadth_first_iteration(
        [&](int node_index, int lcp, int edge_lcp, int occurrences) -> bool {
          auto label = node_label(node_index, lcp, edge_lcp);
          if ((flags[node_index] & Flag::Leaf) == Flag::Leaf) {
            label = leaf_label(node_index, lcp);
          }

          labels.emplace_back(label, occurrences);
          return true;
        });

    return labels;
  }

  std::vector<int> search(std::vector<alphabet_t> pattern) {
    if (pattern.size() == 0) {
      std::vector<int> all_suffixes{};
      iterate_children(0, [&](int index) {
        auto suffixes = suffix_indicies(index, 0);
        all_suffixes.insert(all_suffixes.end(), suffixes.begin(),
                            suffixes.end());
      });
      return all_suffixes;
    }

    auto [node_index, lcp] = find(pattern);

    if (node_index == -1) {
      return std::vector<int>{};
    } else {
      return suffix_indicies(node_index, lcp);
    }
  }

  void
  breadth_first_traversal(const std::function<bool(int, int, int, int)> &f) {
    breadth_first_iteration(
        [&](int node_index, int lcp, int edge_lcp, int occurrences) -> bool {
          int node_start = table[node_index] - lcp;
          int node_end = table[node_index] + edge_lcp;

          if ((flags[node_index] & Flag::Leaf) == Flag::Leaf) {
            node_end = suffixes.size() - 1;
          }

          return f(node_start, node_end, edge_lcp, occurrences);
        });
  }

  void expand_implicit_nodes() {
    std::queue<int> queue{};
    queue.push(0);

    while (!queue.empty()) {
      int node_index = queue.front();
      queue.pop();

      if ((flags[node_index] & Flag::Unevaluated) == Flag::Unevaluated) {
        continue;
      }

      if ((flags[node_index] & Flag::Leaf) != Flag::Leaf) {
        iterate_children(node_index, [&](int index) { queue.push(index); });
      }

      int edge_lcp = get_edge_lcp(node_index);
      if (edge_lcp > 1) {
        lst::details::add_implicit_nodes(node_index, edge_lcp, table, flags);
      }
    }
  }

private:
  // TODO: Can we make the vectors bitcompressed, as we rarely need the full
  // bytes per input??
  std::vector<seqan3::gapped<alphabet_t>> sequence{};
  std::vector<int> suffixes{};
  std::vector<int> table{0, 2};
  std::vector<Flag> flags{Flag::RightMostChild, Flag::None};

  void
  breadth_first_iteration(const std::function<bool(int, int, int, int)> &f) {
    std::queue<std::tuple<int, int>> queue{};
    queue.emplace(0, 0);

    while (!queue.empty()) {
      auto [node_index, lcp] = queue.front();
      queue.pop();

      int edge_lcp, occurrences;

      if ((flags[node_index] & Flag::Unevaluated) == Flag::Unevaluated) {
        std::tie(occurrences, edge_lcp) = lst::details::expand_node(
            node_index, sequence, suffixes, table, flags);
      } else {
        edge_lcp = get_edge_lcp(node_index);
        occurrences = node_occurrences(node_index);
      }

      bool consider_children = f(node_index, lcp, edge_lcp, occurrences);

      if (!consider_children) {
        continue;
      }

      if ((flags[node_index] & Flag::Leaf) != Flag::Leaf) {
        int new_lcp = lcp + edge_lcp;
        iterate_children(node_index,
                         [&](int index) { queue.emplace(index, new_lcp); });
      }
    }
  }

  std::vector<int> suffix_indicies(int node_index, int og_lcp) {
    if (node_index >= table.size()) {
      throw std::invalid_argument(
          "[SUFFIX INDICIES] Given node index is too large.");
    }

    std::vector<int> start_indicies{};
    std::queue<std::tuple<int, int>> queue{};

    queue.emplace(node_index, og_lcp);

    while (!queue.empty()) {
      auto [index, lcp] = queue.front();
      queue.pop();

      if ((flags[index] & Flag::Leaf) == Flag::Leaf) {
        start_indicies.push_back(table[index] - lcp);
      } else if ((flags[index] & Flag::Unevaluated) == Flag::Unevaluated) {
        for (int i = table[index]; i < table[index + 1]; i++) {
          start_indicies.push_back(suffixes[i] - lcp);
        }
      } else {
        int edge_lcp = get_edge_lcp(index);
        int new_lcp = lcp + edge_lcp;
        iterate_children(index, [&](int i) { queue.emplace(i, new_lcp); });
      }
    }

    return start_indicies;
  }

  std::tuple<int, int> find(std::vector<alphabet_t> pattern) {
    std::queue<std::tuple<int, int>> queue{};

    queue.emplace(0, 0);

    int pattern_lcp = 0;

    while (!queue.empty()) {
      auto [node_index, lcp] = queue.front();
      queue.pop();

      int edge_lcp, occurrences;
      if ((flags[node_index] & Flag::Unevaluated) == Flag::Unevaluated) {
        std::tie(occurrences, edge_lcp) = lst::details::expand_node(
            node_index, sequence, suffixes, table, flags);
      } else {
        edge_lcp = get_edge_lcp(node_index);
      }

      sequence_t<alphabet_t> edge = edge_label(node_index, edge_lcp);

      bool edge_match = edge_matches(node_index, pattern_lcp, pattern, edge);

      if (!edge_match) {
        continue;
      }

      pattern_lcp += edge.size();
      if (pattern_lcp >= pattern.size()) {
        return std::make_tuple(node_index, lcp);
      }

      empty_queue(queue);

      if ((flags[node_index] & Flag::Leaf) != Flag::Leaf) {
        int new_lcp = lcp + edge_lcp;
        iterate_children(node_index,
                         [&](int index) { queue.emplace(index, new_lcp); });
      }
    }

    return std::make_tuple(-1, -1);
  }

  int node_occurrences(int node_index) {
    if (node_index > table.size()) {
      throw std::invalid_argument(
          "[NODE OCCURRENCES] Given node index is too large.");
    }

    int occurrences = 0;
    std::queue<int> queue{};

    queue.push(node_index);

    while (!queue.empty()) {
      auto index = queue.front();
      queue.pop();

      if ((flags[index] & Flag::Leaf) == Flag::Leaf) {
        occurrences += 1;
      } else if ((flags[index] & Flag::Unevaluated) == Flag::Unevaluated) {
        occurrences += table[index + 1] - table[index];
      } else {
        iterate_children(index, [&](int i) { queue.push(i); });
      }
    }

    return occurrences;
  }

  int get_edge_lcp(int node_index) {
    if ((flags[node_index] & Flag::Leaf) == Flag::Leaf) {
      return sequence.size() - table[node_index];
    }

    if ((flags[node_index] & Flag::Unevaluated) == Flag::Unevaluated) {
      return lst::details::longest_common_prefix(
          table[node_index], table[node_index + 1], sequence, suffixes);
    }

    iterate_children(node_index, [&](int index) {
      if ((flags[index] & Flag::Unevaluated) == Flag::Unevaluated) {
        lst::details::expand_node(index, sequence, suffixes, table, flags);
      }
    });

    int first_child = table[node_index + 1];
    int smallest_child_index = table[first_child];
    iterate_children(node_index, [&](int index) {
      if (table[index] < smallest_child_index) {
        smallest_child_index = table[index];
      }
    });

    return smallest_child_index - table[node_index];
  }

  void iterate_children(int node_index, const std::function<void(int)> &f) {
    if ((flags[node_index] & Flag::Leaf) != Flag::Leaf) {
      int first_child = table[node_index + 1];

      for (int i = first_child; i <= table.size();) {
        f(i);

        if ((flags[i] & Flag::RightMostChild) == Flag::RightMostChild) {
          break;
        }
        i = next_child_index(i);
      }
    }
  }

  sequence_t<alphabet_t> edge_label(int node_index, int edge_lcp) {
    int edge_start = table[node_index];
    sequence_t<alphabet_t> edge(sequence.begin() + edge_start,
                                sequence.begin() + edge_start + edge_lcp);

    return edge;
  }

  sequence_t<alphabet_t> node_label(int node_index, int lcp, int edge_lcp) {
    int node_start = table[node_index] - lcp;
    int node_end = table[node_index] + edge_lcp;
    sequence_t<alphabet_t> label(sequence.begin() + node_start,
                                 sequence.begin() + node_end);

    return label;
  }

  sequence_t<alphabet_t> leaf_label(int node_index, int lcp) {
    int node_start = table[node_index] - lcp;
    sequence_t<alphabet_t> label(sequence.begin() + node_start, sequence.end());

    return label;
  }

  int next_child_index(int node_index) {
    if ((flags[node_index] & Flag::Leaf) == Flag::Leaf) {
      // Should be 1, but I've addded a value to leaves to allow for
      // explicit nodes.
      node_index += 2;

    } else {
      node_index += 2;
    }

    return node_index;
  }

  bool edge_matches(int node_index, int pattern_lcp,
                    std::vector<alphabet_t> &pattern,
                    sequence_t<alphabet_t> &edge) {
    for (int i = 0; pattern_lcp + i < pattern.size() && i < edge.size(); i++) {
      int sequence_index = table[node_index] + i;

      if (sequence[sequence_index] != pattern[pattern_lcp + i]) {
        return false;
      }
    }

    return true;
  }

  template <typename T> void empty_queue(std::queue<T> &queue) {
    while (!queue.empty()) {
      queue.pop();
    }
  }
};
} // namespace lst
