#pragma once

#include <mutex>
#include <tuple>
#include <vector>

#include <seqan3/alphabet/concept.hpp>
#include <seqan3/alphabet/gap/all.hpp>
#include <seqan3/alphabet/nucleotide/dna5.hpp>

namespace lst::details {

template <seqan3::alphabet alphabet_t>
using sequence_t = std::vector<alphabet_t>;

template <class T = int, seqan3::alphabet alphabet_t = seqan3::dna5>
using alphabet_array =
    std::array<T, seqan3::alphabet_size<seqan3::gapped<alphabet_t>>>;

enum Flag : unsigned char {
  NONE = 0,
  RIGHT_MOST_CHILD = 1 << 0,
  LEAF = 1 << 1,
  UNEVALUATED = 1 << 2,
};

template <class size_t = int, class flag_t = Flag> struct Entry {
  size_t value;
  flag_t flag;
};

template <class size_t = int, class flag_t = Flag> struct Table {
  std::vector<Entry<size_t, flag_t>> table;
  mutable std::mutex mutex;

  using size_type = std::vector<Entry<>>::size_type;

  Table(std::initializer_list<Entry<size_t, flag_t>> args...) {
    table = std::vector<Entry<size_t, flag_t>>{args};
  }

  Table(Table &&other) noexcept {
    std::lock_guard<std::mutex> lock(other.mutex);
    table = std::move(other.table);
  }

  Table(const Table &other) {
    std::lock_guard<std::mutex> lock(other.mutex);
    table = other.table;
  }

  // Move assignment
  Table &operator=(Table &&other) noexcept {
    std::lock(mutex, other.mutex);
    std::lock_guard<std::mutex> self_lock(mutex, std::adopt_lock);
    std::lock_guard<std::mutex> other_lock(other.mutex, std::adopt_lock);
    table = std::move(other.table);
    return *this;
  }

  // Copy assignment
  Table &operator=(const Table &other) {
    std::lock(mutex, other.mutex);
    std::lock_guard<std::mutex> self_lock(mutex, std::adopt_lock);
    std::lock_guard<std::mutex> other_lock(other.mutex, std::adopt_lock);
    table = other.table;
    return *this;
  }

  Entry<> &operator[](size_type pos) { return table[pos]; }

  Entry<> operator[](size_type pos) const { return table[pos]; }
  size_type size() const { return table.size(); }
};

bool is_leaf(int node_index, const Table<> &table) {
  return (table[node_index].flag & Flag::LEAF) == Flag::LEAF;
}

bool is_unevaluated(int node_index, const Table<> &table) {
  return (table[node_index].flag & Flag::UNEVALUATED) == Flag::UNEVALUATED;
}

bool is_rightmostchild(int node_index, const Table<> &table) {
  return (table[node_index].flag & Flag::RIGHT_MOST_CHILD) ==
         Flag::RIGHT_MOST_CHILD;
}

void add_branching_node(int index, int count, Table<> &table) {
  table.table.push_back({index, Flag::UNEVALUATED});
  table.table.push_back({index + count, Flag::NONE});
}

void add_leaf(int index, Table<> &table, const std::vector<int> &suffixes) {
  table.table.push_back({suffixes[index], Flag::LEAF});
  // Add extra index for leaves to allow for explicit nodes.
  table.table.push_back({0, Flag::NONE});
}

void add_lcp_to_suffixes(int lower_bound, int upper_bound, int lcp,
                         std::vector<int> &suffixes) {
  for (int i = lower_bound; i < upper_bound; i++) {
    suffixes[i] += lcp;
  }
}

template <seqan3::alphabet alphabet_t>
inline seqan3::gapped<alphabet_t>
get_character(const sequence_t<alphabet_t> &sequence, int index) {
  // In testing, this is faster than manually checking bounds
  try {
    alphabet_t character = sequence.at(index);
    return seqan3::gapped<alphabet_t>{character};
  } catch (std::out_of_range const &error) {
    return seqan3::gapped<alphabet_t>{seqan3::gap{}};
  }
}

template <seqan3::alphabet alphabet_t>
inline int get_character_rank(const sequence_t<alphabet_t> &sequence,
                              int index) {
  // In testing, this is faster than manually checking bounds
  try {
    return sequence.at(index).to_rank();
  } catch (std::out_of_range const &error) {
    return seqan3::alphabet_size<alphabet_t>;
  }
}

template <seqan3::alphabet alphabet_t>
int longest_common_prefix(int lower_bound, int upper_bound,
                          const sequence_t<alphabet_t> &sequence,
                          const std::vector<int> &suffixes) {
  for (int prefix_length = 1;; prefix_length++) {
    if (prefix_length + suffixes[upper_bound - 1] >= suffixes.size()) {
      return prefix_length - 1;
    }

    auto character_rank =
        get_character_rank(sequence, suffixes[lower_bound] + prefix_length);

    for (int i = lower_bound + 1; i < upper_bound; i++) {
      if (get_character_rank(sequence, suffixes[i] + prefix_length) !=
          character_rank) {
        return prefix_length;
      }
    }
  }

  return -1;
}
template <seqan3::alphabet alphabet_t>
void add_children(const alphabet_array<int, alphabet_t> &counts,
                  int lower_bound, const std::vector<int> &suffixes,
                  Table<> &table) {
  bool last_added_leaf = false;
  int index = lower_bound;

  for (auto count : counts) {
    if (count == 0) {
      continue;
    } else if (count == 1) {
      add_leaf(index, table, suffixes);
      last_added_leaf = true;
    } else {
      add_branching_node(index, count, table);
      last_added_leaf = false;
    }

    index += count;
  }

  auto right_most_child_index = table.size() - 2;

  //  if (last_added_leaf) {
  //    // Should be 1, but I've added a value to leaves to allow for
  //    // explicit nodes.
  //    right_most_child_index = flags.size() - 2;
  //  }

  table[right_most_child_index].flag =
      Flag(table[right_most_child_index].flag | Flag::RIGHT_MOST_CHILD);
}

template <seqan3::alphabet alphabet_t>
alphabet_array<int, alphabet_t>
suffix_pointers(const alphabet_array<int, alphabet_t> &counts) {
  alphabet_array<int, alphabet_t> pointers{};

  int counter = 0;
  for (int i = 0; i < counts.size(); i++) {
    pointers[i] = counter;
    counter += counts[i];
  }

  return pointers;
}

template <seqan3::alphabet alphabet_t>
void sort_suffixes(const alphabet_array<int, alphabet_t> &counts,
                   int lower_bound, int upper_bound,
                   const sequence_t<alphabet_t> &sequence,
                   std::vector<int> &suffixes) {
  std::vector<int> temp_suffixes(suffixes.begin() + lower_bound,
                                 suffixes.begin() + upper_bound);

  auto pointers = suffix_pointers<alphabet_t>(counts);
  for (auto &p : pointers) {
    p += lower_bound;
  }

  for (int suffix : temp_suffixes) {
    auto character_rank = get_character_rank(sequence, suffix);

    int suffix_index = pointers[character_rank]++;

    suffixes[suffix_index] = suffix;
  }
}

template <seqan3::alphabet alphabet_t>
void sort_suffixes_root(const alphabet_array<int, alphabet_t> &counts,
                        const sequence_t<alphabet_t> &sequence,
                        std::vector<int> &suffixes) {
  // Haven't got suffixes before root is created.
  auto pointers = suffix_pointers<alphabet_t>(counts);

  for (int i = 0; i < suffixes.size(); i++) {
    auto character_rank = get_character_rank(sequence, i);

    int suffix_index = pointers[character_rank];

    suffixes[suffix_index] = i;

    pointers[character_rank] += 1;
  }
}

template <seqan3::alphabet alphabet_t>
alphabet_array<int, alphabet_t>
count_suffixes(int lower_bound, int upper_bound,
               const sequence_t<alphabet_t> &sequence,
               const std::vector<int> &suffixes) {
  assert(upper_bound <= suffixes.size());
  assert(lower_bound >= 0);

  alphabet_array<int, alphabet_t> count{};

  for (int i = lower_bound; i < upper_bound; i++) {
    auto character_rank = get_character_rank(sequence, suffixes[i]);
    count[character_rank] += 1;
  }

  return count;
}

template <seqan3::alphabet alphabet_t>
alphabet_array<int, alphabet_t>
count_suffixes_root(const sequence_t<alphabet_t> &sequence) {
  // Haven't got suffixes before root is created.
  alphabet_array<int, alphabet_t> count{};

  for (int i = 0; i < sequence.size() + 1; i++) {
    auto character_rank = get_character_rank(sequence, i);
    count[character_rank] += 1;
  }

  return count;
}

template <seqan3::alphabet alphabet_t>
void expand_root(const sequence_t<alphabet_t> &sequence,
                 std::vector<int> &suffixes, Table<> &table) {
  int lower_bound = 0;

  auto counts = count_suffixes_root(sequence);

  sort_suffixes_root(counts, sequence, suffixes);

  add_children<alphabet_t>(counts, lower_bound, suffixes, table);
}

template <seqan3::alphabet alphabet_t>
std::tuple<int, int> expand_node(int node_index,
                                 const sequence_t<alphabet_t> &sequence,
                                 std::vector<int> &suffixes, Table<> &table) {
  assert(is_unevaluated(node_index, table));

  int lower_bound = table[node_index].value;
  int upper_bound = table[node_index + 1].value;
  int suffix_lower_bound = suffixes[lower_bound];

  int lcp = longest_common_prefix(lower_bound, upper_bound, sequence, suffixes);
  assert(lcp > 0);
  add_lcp_to_suffixes(lower_bound, upper_bound, lcp, suffixes);

  auto counts = count_suffixes(lower_bound, upper_bound, sequence, suffixes);

  sort_suffixes(counts, lower_bound, upper_bound, sequence, suffixes);

  std::lock_guard<std::mutex> unevaluated_lock{table.mutex};

  table[node_index].value = suffix_lower_bound;
  table[node_index + 1].value = table.size();

  add_children<alphabet_t>(counts, lower_bound, suffixes, table);
  table[node_index].flag = Flag(table[node_index].flag & ~Flag::UNEVALUATED);

  int count = upper_bound - lower_bound;
  return {lcp, count};
}

void add_implicit_nodes(int node_index, int edge_lcp, Table<> &table) {
  assert(!is_unevaluated(node_index, table));

  std::lock_guard<std::mutex> expand_table_lock{table.mutex};

  int previous_child = table[node_index + 1].value;
  table[node_index + 1].value = table.size();

  int start = table[node_index].value;
  for (auto i = start + 1; i < start + edge_lcp; i++) {
    table.table.push_back({i, Flag::RIGHT_MOST_CHILD});
    table.table.push_back({static_cast<int>(table.size() + 1), Flag::NONE});
  }

  if (is_leaf(node_index, table)) {
    table[node_index].flag = Flag(table[node_index].flag & ~Flag::LEAF);
    table[table.size() - 2].flag =
        Flag(table[table.size() - 2].flag | Flag::LEAF);

    table[table.size() - 1].value = 0;
  } else {
    table[table.size() - 1].value = previous_child;
  }
}

int get_sequence_index(int node_index, const std::vector<int> &suffixes,
                       const Table<> &table) {
  if (is_unevaluated(node_index, table)) {
    return suffixes[table[node_index].value];
  } else {
    return table[node_index].value;
  }
}

} // namespace lst::details
