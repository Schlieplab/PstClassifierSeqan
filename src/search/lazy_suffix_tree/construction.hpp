#pragma once

#include <tuple>
#include <vector>

#include <mutex>
#include <seqan3/alphabet/concept.hpp>
#include <seqan3/alphabet/gap/all.hpp>
#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/range/container/bitcompressed_vector.hpp>

std::mutex l1;
std::mutex l2;
std::mutex l3;
std::mutex l4;


namespace lst::details {
    std::mutex sort_suffix_lock;
    std::mutex add_child_lock;
    std::mutex lcp_lock;
template <seqan3::alphabet alphabet_t>
using sequence_t = seqan3::bitcompressed_vector<alphabet_t>;

template <seqan3::alphabet alphabet_t = seqan3::dna5>
using alphabet_array =
    std::array<int64_t, seqan3::alphabet_size<seqan3::gapped<alphabet_t>>>;

enum Flag : unsigned char {
  NONE = 0,
  RIGHT_MOST_CHILD = 1 << 0,
  LEAF = 1 << 1,
  UNEVALUATED = 1 << 2,
};

bool is_leaf(int64_t node_index, std::vector<Flag> &flags) {
  return (flags[node_index] & Flag::LEAF) == Flag::LEAF;
}

bool is_unevaluated(int64_t node_index, std::vector<Flag> &flags) {
  return (flags[node_index] & Flag::UNEVALUATED) == Flag::UNEVALUATED;
}

bool is_rightmostchild(int64_t node_index, std::vector<Flag> &flags) {
  return (flags[node_index] & Flag::RIGHT_MOST_CHILD) == Flag::RIGHT_MOST_CHILD;
}

void add_branching_node(int64_t index, int64_t count,
                        std::vector<int64_t> &table,
                        std::vector<Flag> &flags) {
  table.push_back(index);
  table.push_back(index + count);

  flags.push_back(Flag::UNEVALUATED);
  flags.push_back(Flag::NONE);
}

void add_leaf(int64_t index,
              std::vector<int64_t> &table,
              std::vector<Flag> &flags,
              std::vector<int64_t> &suffixes) {
  table.push_back(suffixes[index]);
  // Add extra index for leaves to allow for explicit nodes.
  // It's not technically needed.
  table.push_back(0);

  flags.push_back(Flag::LEAF);
  flags.push_back(Flag::NONE);
}


template <seqan3::alphabet alphabet_t>
seqan3::gapped<alphabet_t> get_character(sequence_t<alphabet_t> &sequence,
                                         int64_t index) {
  if (index >= sequence.size()) {
    return seqan3::gapped<alphabet_t>(seqan3::gap{});
  } else {
    return seqan3::gapped<alphabet_t>(sequence[index]);
  }
}

template <seqan3::alphabet alphabet_t>
int64_t longest_common_prefix(int64_t lower_bound,
                                   int64_t upper_bound,
                                   sequence_t<alphabet_t> &sequence,
                                   std::vector<int64_t> &suffixes) {
  for (int64_t prefix_length = 0;; prefix_length++) {
    if (prefix_length + suffixes[upper_bound - 1] >= suffixes.size()) {
      return prefix_length - 1;
    }

    auto character =
        get_character(sequence, suffixes[lower_bound] + prefix_length);

    for (int64_t i = lower_bound + 1; i < upper_bound; i++) {
      if (get_character(sequence, suffixes[i] + prefix_length) != character) {
        return prefix_length;
      }
    }
  }

  return -1;
}

template <seqan3::alphabet alphabet_t>
int64_t add_lcp_to_suffixes(int64_t lower_bound,
                                 int64_t upper_bound,
                                 sequence_t<alphabet_t> &sequence,
                                 std::vector<int64_t> &suffixes) {

  int64_t lcp = longest_common_prefix(lower_bound,
                                           upper_bound,
                                           sequence,
                                           suffixes);

  std::lock_guard<std::mutex> lock(lcp_lock);
  for (int64_t i = lower_bound; i < upper_bound; i++) {
    suffixes[i] += lcp;
  }
  return lcp;
}



template <seqan3::alphabet alphabet_t>
void add_children(alphabet_array<alphabet_t> &counts,
                  int64_t lower_bound,
                  std::vector<int64_t> &suffixes,
                  std::vector<int64_t> &table,
                  std::vector<Flag> &flags) {
  bool last_added_leaf = false;
  int64_t index = lower_bound;

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

  int64_t right_most_child_index = flags.size() - 2;

  if (last_added_leaf) {
    // Should be 1, but I've addded a value to leaves to allow for
    // explicit nodes.
    right_most_child_index = flags.size() - 2;
  }

  flags[right_most_child_index] =
      Flag(flags[right_most_child_index] | Flag::RIGHT_MOST_CHILD);
}

template <seqan3::alphabet alphabet_t>
alphabet_array<alphabet_t> suffix_pointers(alphabet_array<alphabet_t> &counts) {
  alphabet_array<alphabet_t> pointers{};

  int64_t counter = 0;
  for (int64_t i = 0; i < counts.size(); i++) {
    pointers[i] = counter;
    counter += counts[i];
  }

  return pointers;
}

template <seqan3::alphabet alphabet_t>
void sort_suffixes_root(alphabet_array<alphabet_t> counts,
                        int64_t lower_bound, int64_t upper_bound,
                        sequence_t<alphabet_t> &sequence,
                        std::vector<int64_t > &suffixes) {
  std::vector<int64_t> temp_suffixes(suffixes.begin() + lower_bound,
                                     suffixes.begin() + upper_bound);

  auto pointers = suffix_pointers<alphabet_t>(counts);

  for (int64_t i = lower_bound; i < upper_bound; i++) {
    auto character = get_character(sequence, temp_suffixes[i - lower_bound]);
    int character_rank = seqan3::to_rank(character);

    int64_t suffix_index = pointers[character_rank] + lower_bound;

    suffixes[suffix_index] = temp_suffixes[i - lower_bound];

    pointers[character_rank] += 1;
  }
}

template <seqan3::alphabet alphabet_t>
void sort_suffixes(alphabet_array<alphabet_t> counts,
                   int64_t lower_bound,
                   int64_t upper_bound,
                   sequence_t<alphabet_t> &sequence,
                   std::vector<int64_t> &suffixes) {

  std::vector<int64_t> temp_suffixes(suffixes.begin() + lower_bound,
                                          suffixes.begin() + upper_bound);
  std::vector<int64_t> temp_storage;
  auto pointers = suffix_pointers<alphabet_t>(counts);

  for (int64_t i = lower_bound; i < upper_bound; i++) {
    auto character = get_character(sequence, temp_suffixes[i - lower_bound]);
    int character_rank = seqan3::to_rank(character);

    int64_t suffix_index = pointers[character_rank] + lower_bound;

    temp_storage.push_back(suffix_index);
    temp_storage.push_back(temp_suffixes[i - lower_bound]);
    pointers[character_rank] += 1;
  }

  std::lock_guard<std::mutex> lock(sort_suffix_lock);
  for (auto j = 0; j < temp_storage.size(); j += 2) {
    suffixes[temp_storage[j]] = temp_storage[j + 1];
  }

}
template <seqan3::alphabet alphabet_t>
alphabet_array<alphabet_t> count_suffixes(int64_t lower_bound,
                                          int64_t upper_bound,
                                          sequence_t<alphabet_t> &sequence,
                                          std::vector<int64_t> &suffixes) {
  if (upper_bound > suffixes.size()) {
    seqan3::debug_stream << "Upper Bound " << upper_bound << " Suffix Size " << suffixes.size() << std::endl;
    throw std::invalid_argument("[COUNT SUFFIXES] Upper bound is larger than "
                                "the size of the sequence.");
  }
  if (lower_bound < 0) {
    throw std::invalid_argument("[COUNT SUFFIXES] Lower bound is less than 0.");
  }
  alphabet_array<alphabet_t> count{};

  for (int64_t i = lower_bound; i < upper_bound; i++) {

    auto character     = get_character(sequence, suffixes[i]);
    int character_rank = seqan3::to_rank(character);
    count[character_rank] += 1;
  }

  return count;
}

template <seqan3::alphabet alphabet_t>
void expand_root(sequence_t<alphabet_t> &sequence,
                 std::vector<int64_t> &suffixes,
                 std::vector<int64_t> &table,
                 std::vector<Flag> &flags) {
  int64_t lower_bound = 0;
  int64_t upper_bound = suffixes.size();

  auto counts = count_suffixes(lower_bound, upper_bound, sequence, suffixes);

  sort_suffixes_root(counts, lower_bound, upper_bound, sequence, suffixes);

  add_children<alphabet_t>(counts, lower_bound, suffixes, table, flags);
}

template <seqan3::alphabet alphabet_t>
int64_t expand_node(int64_t node_index,
                sequence_t<alphabet_t> &sequence,
                std::vector<int64_t> &suffixes,
                std::vector<int64_t> &table,
                std::vector<Flag> &flags) {

  if (!is_unevaluated(node_index, flags)) {
    throw std::invalid_argument("[EXPAND NODE] Given node is already expanded");
  }

  int64_t lower_bound  = table[node_index];
  int64_t upper_bound  = table[node_index + 1];
  table[node_index]= suffixes[lower_bound];
  int64_t lcp = add_lcp_to_suffixes(lower_bound,
                                         upper_bound,
                                         sequence,
                                         suffixes);

  alphabet_array<alphabet_t> counts = count_suffixes(lower_bound,
                                                     upper_bound,
                                                     sequence,
                                                     suffixes);

  sort_suffixes(counts, lower_bound, upper_bound, sequence, suffixes);

  std::lock_guard<std::mutex> lock(add_child_lock);
  table[node_index + 1] = table.size();
  add_children<alphabet_t>(counts, lower_bound, suffixes, table, flags);
  flags[node_index]     = Flag(flags[node_index] & ~Flag::UNEVALUATED);

  return lcp;
}

void add_implicit_nodes(int64_t node_index, int64_t edge_lcp,
                        std::vector<int64_t> &table,
                        std::vector<Flag> &flags) {

  std::lock_guard<std::mutex> lock(add_child_lock);
  int64_t previous_child = table[node_index + 1];
  table[node_index + 1] = table.size();

  int64_t start = table[node_index];
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

int64_t get_sequence_index(int64_t node_index, std::vector<int64_t> &suffixes,
                           std::vector<int64_t> &table, std::vector<Flag> &flags) {
  if (is_unevaluated(node_index, flags)) {
    return suffixes[table[node_index]];
  } else {
    return table[node_index];
  }
}

} // namespace lst::details
