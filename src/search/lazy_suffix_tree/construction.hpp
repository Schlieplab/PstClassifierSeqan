#pragma once

#include <mutex>
#include <tuple>
#include <vector>

#include <seqan3/alphabet/concept.hpp>
#include <seqan3/alphabet/gap/all.hpp>
#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/range/container/bitcompressed_vector.hpp>

namespace lst::details {

std::mutex expand_table_mutex;

template <seqan3::alphabet alphabet_t>
using sequence_t = std::vector<alphabet_t>;

template <seqan3::alphabet alphabet_t = seqan3::dna5>
using alphabet_array =
    std::array<int, seqan3::alphabet_size<seqan3::gapped<alphabet_t>>>;

enum Flag : unsigned char {
  NONE = 0,
  RIGHT_MOST_CHILD = 1 << 0,
  LEAF = 1 << 1,
  UNEVALUATED = 1 << 2,
};

bool is_leaf(int node_index, const std::vector<Flag> &flags) {
  return (flags[node_index] & Flag::LEAF) == Flag::LEAF;
}

bool is_unevaluated(int node_index, const std::vector<Flag> &flags) {
  return (flags[node_index] & Flag::UNEVALUATED) == Flag::UNEVALUATED;
}

bool is_rightmostchild(int node_index, const std::vector<Flag> &flags) {
  return (flags[node_index] & Flag::RIGHT_MOST_CHILD) == Flag::RIGHT_MOST_CHILD;
}

void add_branching_node(int index, int count, std::vector<int> &table,
                        std::vector<Flag> &flags) {
  table.push_back(index);
  table.push_back(index + count);

  flags.push_back(Flag::UNEVALUATED);
  flags.push_back(Flag::NONE);
}

void add_leaf(int index, std::vector<int> &table, std::vector<Flag> &flags,
              const std::vector<int> &suffixes) {
  table.push_back(suffixes[index]);
  // Add extra index for leaves to allow for explicit nodes.
  table.push_back(0);

  flags.push_back(Flag::LEAF);
  flags.push_back(Flag::NONE);
}

void add_lcp_to_suffixes(int lower_bound, int upper_bound, int lcp,
                         std::vector<int> &suffixes) {
  static std::mutex lcp_mutex{};
  std::lock_guard<std::mutex> lock{lcp_mutex};
  for (int i = lower_bound; i < upper_bound; i++) {
    suffixes[i] += lcp;
  }
}

template <seqan3::alphabet alphabet_t>
seqan3::gapped<alphabet_t> get_character(const sequence_t<alphabet_t> &sequence,
                                         size_t index) {
  if (index >= sequence.size()) {
    return seqan3::gapped<alphabet_t>(seqan3::gap{});
  } else {
    alphabet_t character = sequence[index];
    return seqan3::gapped<alphabet_t>(character);
  }
}

template <seqan3::alphabet alphabet_t>
int longest_common_prefix(int lower_bound, int upper_bound,
                          const sequence_t<alphabet_t> &sequence,
                          const std::vector<int> &suffixes) {
  for (int prefix_length = 0;; prefix_length++) {
    if (prefix_length + suffixes[upper_bound - 1] >= suffixes.size()) {
      return prefix_length - 1;
    }

    auto character =
        get_character(sequence, suffixes[lower_bound] + prefix_length);

    for (int i = lower_bound + 1; i < upper_bound; i++) {
      if (get_character(sequence, suffixes[i] + prefix_length) != character) {
        return prefix_length;
      }
    }
  }

  return -1;
}
template <seqan3::alphabet alphabet_t>
void add_children(const alphabet_array<alphabet_t> &counts, int lower_bound,
                  const std::vector<int> &suffixes, std::vector<int> &table,
                  std::vector<Flag> &flags) {
  bool last_added_leaf = false;
  int index = lower_bound;

  for (auto count : counts) {
    if (count == 0) {
      continue;
    } else if (count == 1) {
      add_leaf(index, table, flags, suffixes);
      last_added_leaf = true;
    } else {
      add_branching_node(index, count, table, flags);
      last_added_leaf = false;
    }

    index += count;
  }

  int right_most_child_index = flags.size() - 2;

  if (last_added_leaf) {
    // Should be 1, but I've added a value to leaves to allow for
    // explicit nodes.
    right_most_child_index = flags.size() - 2;
  }

  flags[right_most_child_index] =
      Flag(flags[right_most_child_index] | Flag::RIGHT_MOST_CHILD);
}

template <seqan3::alphabet alphabet_t>
alphabet_array<alphabet_t>
suffix_pointers(const alphabet_array<alphabet_t> &counts) {
  alphabet_array<alphabet_t> pointers{};

  int counter = 0;
  for (int i = 0; i < counts.size(); i++) {
    pointers[i] = counter;
    counter += counts[i];
  }

  return pointers;
}

template <seqan3::alphabet alphabet_t>
void sort_suffixes(alphabet_array<alphabet_t> counts, int lower_bound,
                   int upper_bound, const sequence_t<alphabet_t> &sequence,
                   std::vector<int> &suffixes) {
  std::vector<int> temp_suffixes(suffixes.begin() + lower_bound,
                                 suffixes.begin() + upper_bound);

  auto pointers = suffix_pointers<alphabet_t>(counts);

  for (int i = lower_bound; i < upper_bound; i++) {
    auto character = get_character(sequence, temp_suffixes[i - lower_bound]);
    int character_rank = seqan3::to_rank(character);

    int suffix_index = pointers[character_rank] + lower_bound;
    assert(suffix_index < upper_bound);

    suffixes[suffix_index] = temp_suffixes[i - lower_bound];

    pointers[character_rank] += 1;
  }
}

template <seqan3::alphabet alphabet_t>
alphabet_array<alphabet_t>
count_suffixes(int lower_bound, int upper_bound,
               const sequence_t<alphabet_t> &sequence,
               const std::vector<int> &suffixes) {
  assert(upper_bound <= suffixes.size());
  assert(lower_bound >= 0);

  alphabet_array<alphabet_t> count{};
  count.fill(0);

  for (int i = lower_bound; i < upper_bound; i++) {
    auto character = get_character(sequence, suffixes[i]);
    int character_rank = seqan3::to_rank(character);
    count[character_rank] += 1;
  }

  return count;
}

template <seqan3::alphabet alphabet_t>
void expand_root(const sequence_t<alphabet_t> &sequence,
                 std::vector<int> &suffixes, std::vector<int> &table,
                 std::vector<Flag> &flags) {
  int lower_bound = 0;
  int upper_bound = suffixes.size();

  auto counts = count_suffixes(lower_bound, upper_bound, sequence, suffixes);

  sort_suffixes(counts, lower_bound, upper_bound, sequence, suffixes);

  add_children<alphabet_t>(counts, lower_bound, suffixes, table, flags);
}

template <seqan3::alphabet alphabet_t>
int expand_node(int node_index, const sequence_t<alphabet_t> &sequence,
                std::vector<int> &suffixes, std::vector<int> &table,
                std::vector<Flag> &flags) {
  assert(is_unevaluated(node_index, flags));

  int lower_bound = table[node_index];
  int upper_bound = table[node_index + 1];
  int suffix_lower_bound = suffixes[lower_bound];

  int lcp = longest_common_prefix(lower_bound, upper_bound, sequence, suffixes);
  assert(lcp > 0);
  add_lcp_to_suffixes(lower_bound, upper_bound, lcp, suffixes);

  alphabet_array<alphabet_t> counts =
      count_suffixes(lower_bound, upper_bound, sequence, suffixes);

  sort_suffixes(counts, lower_bound, upper_bound, sequence, suffixes);

  std::lock_guard<std::mutex> unevaluated_lock{expand_table_mutex};
  table[node_index] = suffix_lower_bound;
  table[node_index + 1] = table.size();

  add_children<alphabet_t>(counts, lower_bound, suffixes, table, flags);
  flags[node_index] = Flag(flags[node_index] & ~Flag::UNEVALUATED);
  return lcp;
}

void add_implicit_nodes(int node_index, int edge_lcp, std::vector<int> &table,
                        std::vector<Flag> &flags) {
  assert(!is_unevaluated(node_index, flags));
  std::lock_guard<std::mutex> expand_table_lock{expand_table_mutex};

  int previous_child = table[node_index + 1];
  table[node_index + 1] = table.size();

  int start = table[node_index];
  for (auto i = start + 1; i < start + edge_lcp; i++) {
    table.push_back(i);
    table.push_back(table.size() + 1);

    flags.push_back(Flag::RIGHT_MOST_CHILD);
    flags.push_back(Flag::NONE);
  }

  if (is_leaf(node_index, flags)) {
    flags[node_index] = Flag(flags[node_index] & ~Flag::LEAF);
    flags[flags.size() - 2] = Flag(flags[flags.size() - 2] | Flag::LEAF);

    table[table.size() - 1] = 0;
  } else {
    table[table.size() - 1] = previous_child;
  }
}

int get_sequence_index(int node_index, const std::vector<int> &suffixes,
                       const std::vector<int> &table,
                       const std::vector<Flag> &flags) {
  if (is_unevaluated(node_index, flags)) {
    return suffixes[table[node_index]];
  } else {
    return table[node_index];
  }
}

} // namespace lst::details
