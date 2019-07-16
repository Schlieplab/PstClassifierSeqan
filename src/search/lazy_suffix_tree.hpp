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

namespace lst {

enum Flag : unsigned char {
  None = 0,
  RightMostChild = 1 << 0,
  Leaf = 1 << 1,
  Unevaluated = 1 << 2,
};

template <seqan3::Alphabet alphabet_t = seqan3::dna5>
using alphabet_count =
    std::array<int, seqan3::alphabet_size<seqan3::gapped<alphabet_t>>>;

template <seqan3::Alphabet alphabet_t = seqan3::dna5>
using sequence_t = std::vector<seqan3::gapped<alphabet_t>>;

template <seqan3::Alphabet alphabet_t = seqan3::dna5> class LazySuffixTree {

public:
  friend class LazySuffixTreeTest;

  LazySuffixTree(std::vector<alphabet_t> sequence_) {
    sequence = sequence_ | seqan3::view::convert<seqan3::gapped<alphabet_t>>;
    sequence.push_back(seqan3::gap{});

    suffixes = std::vector<int>(sequence.size());
    std::iota(suffixes.begin(), suffixes.end(), 0);
  }

  void expand_all() {
    if (table.size() == 0) {
      expand_root();
    }

    int table_size = table.size();
    for (int i = 0; i <= table_size; i++) {

      if ((flags[i] & Flag::Unevaluated) == Flag::Unevaluated) {
        expand_node(i);
      }

      table_size = table.size();
    }
  }

  std::vector<std::tuple<sequence_t<alphabet_t>, int>> get_all_labels() {
    std::vector<std::tuple<sequence_t<alphabet_t>, int>> labels{};

    breadth_first_iteration([&](int node_index, int lcp, int edge_lcp) -> bool {
      auto label = node_label(node_index, lcp, edge_lcp);
      int counts = suffix_indicies(node_index, lcp).size();
      labels.emplace_back(label, counts);
      return false;
    });

    return labels;
  }

  std::vector<int> search(std::vector<alphabet_t> pattern) {
    if (pattern.size() == 0) {
      std::vector<int> all_suffixes{};
      iterate_root_children([&](int index) {
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

  void breadth_first_traversal(
      const std::function<bool(std::vector<alphabet_t>, int)> &f) {
    breadth_first_iteration([&](int node_index, int lcp, int edge_lcp) -> bool {
      auto label = node_label(node_index, lcp, edge_lcp);

      int counts = count(node_index);

      std::vector<alphabet_t> out_label{};
      for (auto c : label) {
        if (c != seqan3::gap{}) {
          out_label.push_back(c.template convert_to<alphabet_t>());
        }
      }

      return f(out_label, counts);
    });
  }

private:
  // TODO: Can we make the vectors bitcompressed, as we rarely need the full
  // bytes per input??
  std::vector<seqan3::gapped<alphabet_t>> sequence{};
  std::vector<int> suffixes{};
  std::vector<int> table{};
  std::vector<Flag> flags{};

  void expand_root() {
    int lower_bound = 0;
    int upper_bound = sequence.size();

    auto counts = count_suffixes(lower_bound, upper_bound);

    sort_suffixes(counts, lower_bound, upper_bound);

    add_children(counts, lower_bound);
  }

  void expand_node(int node_index) {
    if ((flags[node_index] & Flag::Unevaluated) != Flag::Unevaluated) {
      throw std::invalid_argument("Given node is already expanded");
    }

    int lower_bound = table[node_index];
    int upper_bound = table[node_index + 1];

    table[node_index] = suffixes[lower_bound];
    table[node_index + 1] = table.size();

    add_lcp_to_suffixes(lower_bound, upper_bound);

    auto counts = count_suffixes(lower_bound, upper_bound);

    sort_suffixes(counts, lower_bound, upper_bound);

    add_children(counts, lower_bound);

    flags[node_index] = Flag(flags[node_index] & ~Flag::Unevaluated);
  }

  void add_children(alphabet_count<alphabet_t> counts, int lower_bound) {
    bool last_added_leaf = false;
    int index = lower_bound;

    for (auto count : counts) {
      if (count == 0) {
        continue;
      } else if (count == 1) {
        add_leaf(index);
        last_added_leaf = true;
      } else {
        add_branching_node(index, count);
        last_added_leaf = false;
      }

      index += count;
    }

    int right_most_child_index = flags.size() - 2;

    if (last_added_leaf) {
      right_most_child_index = flags.size() - 1;
    }

    flags[right_most_child_index] =
        Flag(flags[right_most_child_index] | Flag::RightMostChild);
  }

  void sort_suffixes(alphabet_count<alphabet_t> counts, int lower_bound,
                     int upper_bound) {
    std::vector<int> temp_suffixes(suffixes.begin() + lower_bound,
                                   suffixes.begin() + upper_bound);

    auto pointers = suffix_pointers(counts);

    for (int i = lower_bound; i < upper_bound; i++) {
      int character_rank =
          seqan3::to_rank(sequence[temp_suffixes[i - lower_bound]]);

      int suffix_index = pointers[character_rank] + lower_bound;

      suffixes[suffix_index] = temp_suffixes[i - lower_bound];

      pointers[character_rank] += 1;
    }
  }

  alphabet_count<alphabet_t>
  suffix_pointers(alphabet_count<alphabet_t> counts) {
    alphabet_count<alphabet_t> pointers{};

    int counter = 0;
    for (int i = 0; i < counts.size(); i++) {
      pointers[i] = counter;
      counter += counts[i];
    }

    return pointers;
  }

  void add_lcp_to_suffixes(int lower_bound, int upper_bound) {
    int lcp = longest_common_prefix(lower_bound, upper_bound);

    for (int i = lower_bound; i < upper_bound; i++) {
      suffixes[i] += lcp;
    }
  }

  alphabet_count<alphabet_t> count_suffixes(int lower_bound, int upper_bound) {
    if (upper_bound > suffixes.size())
      throw std::invalid_argument(
          "Upper bound is larger than the size of the sequence.");

    if (lower_bound < 0)
      throw std::invalid_argument("Lower bound is less than 0.");

    alphabet_count<alphabet_t> count{};

    std::vector<int> group(suffixes.begin() + lower_bound,
                           suffixes.begin() + upper_bound);

    for (auto character_i : group) {
      int character_rank = seqan3::to_rank(sequence[character_i]);
      count[character_rank] += 1;
    }

    return count;
  }

  int longest_common_prefix(int lower_bound, int upper_bound) {
    std::vector<int> group(suffixes.begin() + lower_bound,
                           suffixes.begin() + upper_bound);
    for (int prefix_length = 0;; prefix_length++) {
      if (prefix_length + group.back() >= sequence.size()) {
        return prefix_length - 1;
      }

      seqan3::gapped<alphabet_t> character = sequence[group[0] + prefix_length];

      for (int i = 1; i < group.size(); i++) {
        if (sequence[group[i] + prefix_length] != character) {
          return prefix_length;
        }
      }
    }

    return -1;
  }

  void add_branching_node(int index, int count) {
    table.push_back(index);
    table.push_back(index + count);

    flags.push_back(Flag::Unevaluated);
    flags.push_back(Flag::None);
  }

  void add_leaf(int index) {
    table.push_back(suffixes[index]);
    flags.push_back(Flag::Leaf);
  }

  void breadth_first_iteration(const std::function<bool(int, int, int)> &f) {
    std::queue<std::tuple<int, int>> queue{};
    iterate_root_children([&](int index) { queue.emplace(index, 0); });

    while (queue.size() != 0) {
      auto [node_index, lcp] = queue.front();
      queue.pop();

      if ((flags[node_index] & Flag::Unevaluated) == Flag::Unevaluated) {
        expand_node(node_index);
      }
      int edge_lcp = get_edge_lcp(node_index);

      bool consider_children = f(node_index, lcp, edge_lcp);

      if (!consider_children) {
        continue;
      }

      if ((flags[node_index] & Flag::Leaf) != Flag::Leaf) {
        int new_lcp = lcp + edge_lcp - table[node_index];
        iterate_children(node_index,
                         [&](int index) { queue.emplace(index, new_lcp); });
      }
    }
  }

  std::vector<int> suffix_indicies(int node_index, int og_lcp) {
    if (node_index >= table.size()) {
      throw std::invalid_argument("Given node index is too large.");
    }

    std::vector<int> start_indicies{};
    std::queue<std::tuple<int, int>> queue{};

    queue.emplace(node_index, og_lcp);

    while (queue.size() != 0) {
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
        int new_lcp = lcp + edge_lcp - table[index];
        iterate_children(index, [&](int i) { queue.emplace(i, new_lcp); });
      }
    }

    return start_indicies;
  }

  int count(int node_index) {
    if (node_index >= table.size()) {
      throw std::invalid_argument("Given node index is too large.");
    }

    int occurrences = 0;
    std::queue<int> queue{};

    queue.push(node_index);

    while (queue.size() != 0) {
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

  std::tuple<int, int> find(std::vector<alphabet_t> pattern) {
    std::queue<std::tuple<int, int>> queue{};

    iterate_root_children([&](int index) { queue.emplace(index, 0); });

    int pattern_lcp = 0;

    while (!queue.empty()) {
      auto [node_index, lcp] = queue.front();
      queue.pop();

      if ((flags[node_index] & Flag::Unevaluated) == Flag::Unevaluated) {
        expand_node(node_index);
      }
      int edge_lcp = get_edge_lcp(node_index);

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
        int new_lcp = lcp + edge_lcp - table[node_index];
        iterate_children(node_index,
                         [&](int index) { queue.emplace(index, new_lcp); });
      }
    }

    return std::make_tuple(-1, -1);
  }

  int get_edge_lcp(int node_index) {
    if ((flags[node_index] & Flag::Leaf) == Flag::Leaf) {
      return sequence.size();
    }

    if ((flags[node_index] & Flag::Unevaluated) == Flag::Unevaluated) {
      return longest_common_prefix(table[node_index], table[node_index + 1]);
    }

    iterate_children(node_index, [&](int index) {
      if ((flags[index] & Flag::Unevaluated) == Flag::Unevaluated) {
        expand_node(index);
      }
    });

    int first_child = table[node_index + 1];
    int edge_lcp = table[first_child];
    iterate_children(node_index, [&](int index) {
      if (table[index] < edge_lcp) {
        edge_lcp = table[index];
      }
    });

    return edge_lcp;
  }

  void iterate_children(int node_index, const std::function<void(int)> &f) {
    int first_child = table[node_index + 1];

    iterate(first_child, f);
  }

  void iterate_root_children(const std::function<void(int)> f) {
    if (table.size() == 0) {
      expand_root();
    }
    iterate(0, f);
  }

  void iterate(int start_index, const std::function<void(int)> f) {
    for (int i = start_index; i <= table.size();) {
      f(i);

      if ((flags[i] & Flag::RightMostChild) == Flag::RightMostChild) {
        break;
      }
      i = next_child_index(i);
    }
  }

  sequence_t<alphabet_t> edge_label(int node_index, int edge_lcp) {

    int edge_start = table[node_index];
    sequence_t<alphabet_t> edge(sequence.begin() + edge_start,
                                sequence.begin() + edge_lcp);

    return edge;
  }

  sequence_t<alphabet_t> node_label(int node_index, int lcp, int end_index) {

    int node_start = table[node_index] - lcp;
    sequence_t<alphabet_t> label(sequence.begin() + node_start,
                                 sequence.begin() + end_index);

    return label;
  }

  int next_child_index(int node_index) {
    if ((flags[node_index] & Flag::Leaf) == Flag::Leaf) {
      node_index += 1;
    } else {
      node_index += 2;
    }

    return node_index;
  }

  bool edge_matches(int node_index, int pattern_lcp,
                    std::vector<alphabet_t> pattern,
                    sequence_t<alphabet_t> edge) {
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
