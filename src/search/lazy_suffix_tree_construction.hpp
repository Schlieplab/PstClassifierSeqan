#pragma once

#include <tuple>
#include <vector>

#include <seqan3/alphabet/all.hpp>

namespace lst::details {

template <seqan3::Alphabet alphabet_t = seqan3::dna5>
using sequence_t = std::vector<seqan3::gapped<alphabet_t>>;

template <seqan3::Alphabet alphabet_t = seqan3::dna5>
using alphabet_count =
    std::array<int, seqan3::alphabet_size<seqan3::gapped<alphabet_t>>>;

enum Flag : unsigned char {
  None = 0,
  RightMostChild = 1 << 0,
  Leaf = 1 << 1,
  Unevaluated = 1 << 2,
};

void add_branching_node(int index, int count, std::vector<int> &table,
                        std::vector<Flag> &flags) {
  table.push_back(index);
  table.push_back(index + count);

  flags.push_back(Flag::Unevaluated);
  flags.push_back(Flag::None);
}

void add_leaf(int index, std::vector<int> &table, std::vector<Flag> &flags,
              std::vector<int> &suffixes) {
  table.push_back(suffixes[index]);
  // Add extra index for leaves to allow for explicit nodes.
  // It's not technically needed.
  table.push_back(0);

  flags.push_back(Flag::Leaf);
  flags.push_back(Flag::None);
}

void add_lcp_to_suffixes(int lower_bound, int upper_bound, int lcp,
                         std::vector<int> &suffixes) {
  for (int i = lower_bound; i < upper_bound; i++) {
    suffixes[i] += lcp;
  }
}

template <seqan3::Alphabet alphabet_t>
int longest_common_prefix(int lower_bound, int upper_bound,
                          sequence_t<alphabet_t> &sequence,
                          std::vector<int> &suffixes) {
  for (int prefix_length = 0;; prefix_length++) {
    if (prefix_length + suffixes[upper_bound - 1] >= sequence.size()) {
      return prefix_length - 1;
    }

    seqan3::gapped<alphabet_t> &character =
        sequence[suffixes[lower_bound] + prefix_length];

    for (int i = lower_bound + 1; i < upper_bound; i++) {
      if (sequence[suffixes[i] + prefix_length] != character) {
        return prefix_length;
      }
    }
  }

  return -1;
}

template <seqan3::Alphabet alphabet_t>
void add_children(alphabet_count<alphabet_t> &counts, int lower_bound,
                  std::vector<int> &suffixes, std::vector<int> &table,
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
    // Should be 1, but I've addded a value to leaves to allow for
    // explicit nodes.
    right_most_child_index = flags.size() - 2;
  }

  flags[right_most_child_index] =
      Flag(flags[right_most_child_index] | Flag::RightMostChild);
}

template <seqan3::Alphabet alphabet_t>
alphabet_count<alphabet_t> suffix_pointers(alphabet_count<alphabet_t> counts) {
  alphabet_count<alphabet_t> pointers{};

  int counter = 0;
  for (int i = 0; i < counts.size(); i++) {
    pointers[i] = counter;
    counter += counts[i];
  }

  return pointers;
}

template <seqan3::Alphabet alphabet_t>
void sort_suffixes(alphabet_count<alphabet_t> counts, int lower_bound,
                   int upper_bound, sequence_t<alphabet_t> &sequence,
                   std::vector<int> &suffixes) {
  std::vector<int> temp_suffixes(suffixes.begin() + lower_bound,
                                 suffixes.begin() + upper_bound);

  auto pointers = suffix_pointers<alphabet_t>(counts);

  for (int i = lower_bound; i < upper_bound; i++) {
    int character_rank =
        seqan3::to_rank(sequence[temp_suffixes[i - lower_bound]]);

    int suffix_index = pointers[character_rank] + lower_bound;

    suffixes[suffix_index] = temp_suffixes[i - lower_bound];

    pointers[character_rank] += 1;
  }
}

template <seqan3::Alphabet alphabet_t>
alphabet_count<alphabet_t> count_suffixes(int lower_bound, int upper_bound,
                                          sequence_t<alphabet_t> &sequence,
                                          std::vector<int> &suffixes) {
  if (upper_bound > suffixes.size())
    throw std::invalid_argument(
        "Upper bound is larger than the size of the sequence.");

  if (lower_bound < 0)
    throw std::invalid_argument("Lower bound is less than 0.");

  alphabet_count<alphabet_t> count{};

  for (int i = lower_bound; i < upper_bound; i++) {
    int character_i = suffixes[i];
    int character_rank = seqan3::to_rank(sequence[character_i]);
    count[character_rank] += 1;
  }

  return count;
}

template <seqan3::Alphabet alphabet_t>
void expand_root(sequence_t<alphabet_t> &sequence, std::vector<int> &suffixes,
                 std::vector<int> &table, std::vector<Flag> &flags) {
  int lower_bound = 0;
  int upper_bound = sequence.size();

  auto counts = count_suffixes(lower_bound, upper_bound, sequence, suffixes);

  sort_suffixes(counts, lower_bound, upper_bound, sequence, suffixes);

  add_children<alphabet_t>(counts, lower_bound, suffixes, table, flags);
}

template <seqan3::Alphabet alphabet_t>
std::tuple<int, int>
expand_node(int node_index, sequence_t<alphabet_t> &sequence,
            std::vector<int> &suffixes, std::vector<int> &table,
            std::vector<Flag> &flags) {
  if ((flags[node_index] & Flag::Unevaluated) != Flag::Unevaluated) {
    throw std::invalid_argument("Given node is already expanded");
  }

  int lower_bound = table[node_index];
  int upper_bound = table[node_index + 1];

  table[node_index] = suffixes[lower_bound];
  table[node_index + 1] = table.size();

  int lcp = longest_common_prefix(lower_bound, upper_bound, sequence, suffixes);

  add_lcp_to_suffixes(lower_bound, upper_bound, lcp, suffixes);

  alphabet_count<alphabet_t> counts =
      count_suffixes(lower_bound, upper_bound, sequence, suffixes);

  sort_suffixes(counts, lower_bound, upper_bound, sequence, suffixes);

  add_children<alphabet_t>(counts, lower_bound, suffixes, table, flags);

  flags[node_index] = Flag(flags[node_index] & ~Flag::Unevaluated);

  int occurrences = upper_bound - lower_bound;

  return std::make_tuple(occurrences, lcp);
}

} // namespace lst::details
