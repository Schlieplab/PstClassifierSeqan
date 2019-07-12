#include <array>
#include <functional>
#include <numeric>
#include <queue>
#include <tuple>
#include <vector>

#include <seqan3/alphabet/all.hpp>
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

  LazySuffixTree(std::vector<alphabet_t> sequence_) {
    sequence = sequence_ | seqan3::view::convert<seqan3::gapped<alphabet_t>>;
    sequence.push_back(seqan3::gap{});

    suffixes = std::vector<int>(sequence.size());
    std::iota(suffixes.begin(), suffixes.end(), 0);
  }

  std::vector<sequence_t<alphabet_t>> get_all_labels() {
    std::vector<sequence_t<alphabet_t>> labels{};

    breadth_first_iteration(
        [&](int node_index, int lcp, int smallest_child_index) -> bool {
          auto label = node_label(node_index, lcp, smallest_child_index);
          labels.push_back(label);
          return false;
        });

    return labels;
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

    std::vector<int> group(suffixes.begin() + lower_bound,
                           suffixes.begin() + upper_bound);

    int lcp = longest_common_prefix(group);

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

  int longest_common_prefix(std::vector<int> group) {
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
    if (table.size() == 0) {
      expand_root();
    }

    std::queue<std::tuple<int, int>> queue{};

    for (int i = 0; i <= table.size();) {
      queue.push(std::make_tuple(i, 0));

      if ((flags[i] & Flag::RightMostChild) == Flag::RightMostChild) {
        break;
      }

      i = next_child_index(i);
    }

    bool stop_looping = false;

    while (queue.size() != 0 && !stop_looping) {
      auto [node_index, lcp] = queue.front();
      queue.pop();

      int smallest_child_index = evaluate_node(node_index);

      stop_looping = f(node_index, lcp, smallest_child_index);

      if ((flags[node_index] & Flag::Leaf) != Flag::Leaf) {
        int new_lcp = lcp + smallest_child_index - table[node_index];
        child_iteration(node_index, [&](int index) {
          queue.push(std::make_tuple(index, new_lcp));
        });
      }
    }
  }

  int evaluate_node(int node_index) {
    if ((flags[node_index] & Flag::Unevaluated) == Flag::Unevaluated) {
      expand_node(node_index);
    }

    if ((flags[node_index] & Flag::Leaf) == Flag::Leaf) {
      return sequence.size();
    }

    child_iteration(node_index, [&](int index) {
      if ((flags[index] & Flag::Unevaluated) == Flag::Unevaluated) {
        expand_node(index);
      }
    });

    int first_child = table[node_index + 1];
    int smallest_child_index = table[first_child];

    child_iteration(node_index, [&](int index) {
      if (table[index] < smallest_child_index) {
        smallest_child_index = table[index];
      }
    });

    return smallest_child_index;
  }

  void child_iteration(int node_index, const std::function<void(int)> &f) {
    int first_child = table[node_index + 1];

    for (int i = first_child; i <= table.size();) {
      f(i);

      if ((flags[i] & Flag::RightMostChild) == Flag::RightMostChild) {
        break;
      }
      i = next_child_index(i);
    }
  }

  sequence_t<alphabet_t> edge_label(int node_index, int smallest_child_index) {

    int edge_start = table[node_index];
    int edge_end = smallest_child_index;
    sequence_t<alphabet_t> edge(sequence.begin() + edge_start,
                                sequence.begin() + edge_end);

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
};
} // namespace lst
