#pragma once

#include <array>
#include <chrono>
#include <cmath>
#include <ctime>
#include <functional>
#include <sstream>
#include <string>
#include <tuple>
#include <unordered_set>
#include <vector>

#include <seqan3/alphabet/all.hpp>
#include <seqan3/range/view/to_char.hpp>

using seqan3::operator""_dna4;

#include "search/lazy_suffix_tree.hpp"
#include "search/lazy_suffix_tree/iteration.hpp"

namespace pst {

enum Status : unsigned char {
  None = 0,
  Included = 1 << 0,
  Excluded = 1 << 1,
};

template <seqan3::Alphabet alphabet_t = seqan3::dna5>
class ProbabilisticSuffixTree : public lst::LazySuffixTree<alphabet_t> {
  friend class ProbabilisticSuffixTreeTest;

  ProbabilisticSuffixTree() {}

  ProbabilisticSuffixTree(std::string id, std::vector<alphabet_t> &sequence)
      : ProbabilisticSuffixTree(id, sequence, 15, 100, 1.2, 0, "cutoff", "PS") {
  }

  ProbabilisticSuffixTree(std::string id, std::vector<alphabet_t> &sequence,
                          size_t max_depth, size_t freq, float cutoff_value)
      : ProbabilisticSuffixTree(id, sequence, max_depth, freq, cutoff_value, 0,
                                "cutoff", "KL") {}

  ProbabilisticSuffixTree(std::string id, std::vector<alphabet_t> &sequence,
                          size_t max_depth, size_t freq)
      : ProbabilisticSuffixTree(id, sequence, max_depth, freq, cutoff_value, 0,
                                "cutoff", "PS") {}

  ProbabilisticSuffixTree(std::string id, std::vector<alphabet_t> &sequence,
                          size_t max_depth, size_t freq, float cutoff_value,
                          std::string estimator)
      : ProbabilisticSuffixTree(id, sequence, max_depth, freq, cutoff_value, 0,
                                "cutoff", estimator) {}

  ProbabilisticSuffixTree(std::string id_, std::vector<alphabet_t> &sequence_,
                          size_t max_depth_, size_t freq_, float cutoff_value_,
                          size_t number_of_parameters_,
                          std::string pruning_method_, std::string estimator_)
      : lst::LazySuffixTree<alphabet_t>(sequence_), id(id_), freq(freq_),
        max_depth(max_depth_), cutoff_value(cutoff_value_),
        number_of_parameters(number_of_parameters_),
        pruning_method(pruning_method_), estimator(estimator_) {

    seqan3::dna4_vector dna4{"ACGT"_dna4};
    std::vector<seqan3::gapped<alphabet_t>> characters{
        dna4 | seqan3::view::convert<seqan3::gapped<alphabet_t>>};

    for (auto c : characters) {
      valid_characters.insert(c);
    }

    if (estimator_ == "PS") {
      this->cutoff_value = std::pow(float(sequence_.size()), 3.0 / 4.0);
    }

    this->support_pruning();
    this->similarity_pruning();
  }

  void print() {
    lst::details::breadth_first_iteration(
        this->sequence, this->suffixes, this->table, this->flags, false,
        [&](int node_index, int lcp, int edge_lcp) -> bool {
          if (this->is_excluded(node_index)) {
            return false;
          }

          auto label = this->node_label(node_index, lcp, edge_lcp);

          if (this->is_leaf(node_index)) {
            label = this->leaf_label(node_index, lcp);
          }

          seqan3::debug_stream << label << "\t" << node_index << "\t"
                               << this->table[node_index] << "\t"
                               << this->table[node_index + 1];

          if (this->suffix_links.size() > node_index / 2) {
            seqan3::debug_stream << "\t" << this->suffix_links[node_index / 2];
          }

          if (this->reverse_suffix_links.size() > node_index / 2) {
            seqan3::debug_stream << "\t"
                                 << this->reverse_suffix_links[node_index / 2];
          }

          seqan3::debug_stream << std::endl;
          return true;
        });
  }

  std::string to_tree() {
    std::ostringstream tree_string;

    auto now = std::chrono::system_clock::now();
    std::time_t time = std::chrono::system_clock::to_time_t(now);

    tree_string << "Name: " << this->id << std::endl;
    tree_string << "Date: " << std::ctime(&time);
    tree_string << "Tree: PST" << std::endl;
    tree_string << "Alphabet: " << typeid(alphabet_t).name() << std::endl;
    tree_string << "Number(nodes): " << nodes_in_tree() << std::endl;

    int n_parameters =
        get_terminal_nodes().size() * (seqan3::alphabet_size<alphabet_t> - 1);
    tree_string << "Number(parameters): " << n_parameters << std::endl;

    this->append_node_string(0, 0, 0, tree_string);

    lst::details::breadth_first_iteration(
        this->sequence, this->suffixes, this->table, this->flags, false,
        [&](int node_index, int lcp, int edge_lcp) -> bool {
          if (this->is_excluded(node_index)) {
            return false;
          }
          this->append_node_string(node_index, lcp, edge_lcp, tree_string);

          return true;
        });

    return tree_string.str();
  }

protected:
  std::string id;

  std::unordered_set<seqan3::gapped<alphabet_t>> valid_characters{};
  size_t freq;
  size_t max_depth;
  float cutoff_value;
  size_t number_of_parameters;
  std::string pruning_method;
  std::string estimator;

  std::vector<Status> status{};

  std::vector<int> counts{};
  std::vector<std::array<float, seqan3::alphabet_size<alphabet_t>>>
      probabilities{};

  void support_pruning() {
    this->build_tree();
    this->expand_implicit_nodes();
    this->add_suffix_links();
    this->counts.resize(this->table.size() / 2);
    this->status.resize(this->table.size() / 2);
    status[0] = Status::Included;
  }

  void similarity_pruning() {
    this->add_reverse_suffix_links();
    this->compute_probabilities();
    if (this->pruning_method == "cutoff") {
      this->cutoff_prune();
    } else if (this->pruning_method == "parameters") {
      this->parameters_prune();
    }
  }

  void build_tree() {
    lst::details::breadth_first_iteration(
        this->sequence, this->suffixes, this->table, this->flags,
        [&](int node_index, int lcp, int edge_lcp) -> bool {
          int count = lst::details::node_occurrences(node_index, this->table,
                                                     this->flags);

          this->counts.resize(this->table.size() / 2);
          this->counts[node_index / 2] = count;

          int label_start = this->table[node_index] - lcp;
          int label_end = this->table[node_index] + edge_lcp;
          int label_length = label_end - label_start;

          this->status.resize(this->table.size() / 2);
          if (label_length < this->max_depth && count >= this->freq &&
              label_valid(label_start, label_end)) {
            this->status[node_index / 2] = Status::Included;
            return true;
          }
          this->status[node_index / 2] = Status::Excluded;
          return false;
        });
  }

  void compute_probabilities() {
    this->probabilities.resize(this->table.size() / 2);

    std::queue<int> queue{};
    queue.push(0);

    while (!queue.empty()) {
      int node_index = queue.front();
      queue.pop();

      if (this->is_leaf(node_index) || this->is_unevaluated(node_index) ||
          this->is_excluded(node_index)) {
        continue;
      }

      int child_sum = sum_child_counts(node_index);

      lst::details::iterate_children(
          node_index, this->table, this->flags, [&](int index) {
            assign_probabilities(node_index, index, child_sum);
          });

      lst::details::iterate_children(node_index, this->table, this->flags,
                                     [&](int index) { queue.push(index); });
    }
  }

  void assign_probabilities(int node_index, int child_index, int child_sum) {

    int sequence_index = this->get_sequence_index(child_index);
    if (this->sequence[sequence_index] == seqan3::gap{}) {
      return;
    }

    int character_rank = seqan3::to_rank(this->sequence[sequence_index]);

    float child_count = get_counts(child_index);

    this->probabilities[node_index / 2][character_rank] =
        child_count / child_sum;
  }

  int sum_child_counts(int node_index) {
    int child_sum = 0;
    lst::details::iterate_children(
        node_index, this->table, this->flags, [&](int index) {
          int sequence_index = this->get_sequence_index(index);
          if (this->sequence[sequence_index] == seqan3::gap{}) {
            return;
          }

          int child_count = get_counts(index);
          child_sum += child_count;
        });
    return child_sum;
  }

  bool label_valid(int label_start, int label_end) {
    for (int i = label_start; i < label_end; i++) {
      if (valid_characters.find(this->sequence[i]) == valid_characters.end()) {
        return false;
      }
    }

    return true;
  }

  void expand_implicit_nodes() {
    std::queue<int> queue{};
    queue.push(0);

    while (!queue.empty()) {
      int node_index = queue.front();
      queue.pop();

      if (this->is_unevaluated(node_index) || this->is_excluded(node_index)) {
        continue;
      }

      lst::details::iterate_children(node_index, this->table, this->flags,
                                     [&](int index) { queue.push(index); });

      int edge_lcp = lst::details::get_edge_lcp(
          node_index, this->sequence, this->suffixes, this->table, this->flags);
      if (edge_lcp > 1) {
        lst::details::add_implicit_nodes(node_index, edge_lcp, this->table,
                                         this->flags);
      }
    }
  }

  void cutoff_prune() {
    std::vector<int> terminal_nodes = get_terminal_nodes();
    std::queue<int> bottom_up{};
    for (int v : terminal_nodes) {
      bottom_up.push(v);
    }

    while (!bottom_up.empty()) {
      int node_index = bottom_up.front();
      bottom_up.pop();

      if (node_index == 0 || !this->is_terminal(node_index)) {
        continue;
      }

      float delta = calculate_delta(node_index);

      if (delta < this->cutoff_value) {
        this->status[node_index / 2] = Status::Excluded;

        int parent_index = this->suffix_links[node_index / 2];
        bottom_up.push(parent_index);
      }
    }
  }

  void parameters_prune() {
    std::vector<int> terminal_nodes = get_terminal_nodes();
    auto cmp = [](std::tuple<int, float> left, std::tuple<int, float> right) {
      return std::get<1>(left) < std::get<1>(right);
    };
    std::priority_queue<std::tuple<int, float>,
                        std::vector<std::tuple<int, float>>, decltype(cmp)>
        queue{cmp};
    for (int v : terminal_nodes) {
      queue.emplace(v, -this->calculate_delta(v));
    }

    int current_number_of_parameters =
        terminal_nodes.size() * (seqan3::alphabet_size<alphabet_t> - 1);

    seqan3::debug_stream << current_number_of_parameters << std::endl;

    while (!queue.empty() &&
           current_number_of_parameters >= this->number_of_parameters) {
      auto [node_index, delta] = queue.top();
      queue.pop();

      if (node_index == 0) {
        continue;
      }

      current_number_of_parameters -= (seqan3::alphabet_size<alphabet_t> - 1);

      this->status[node_index / 2] = Status::Excluded;

      int parent_index = this->suffix_links[node_index / 2];
      if (this->is_terminal(parent_index)) {
        queue.emplace(parent_index, -this->calculate_delta(parent_index));
      }
    }
  }

  float calculate_delta(int node_index) {
    if (this->estimator == "PS") {
      return ps_delta(node_index);
    } else if (this->estimator == "KL") {
      return kl_delta(node_index);
    } else {
      return 0.0;
    }
  }

  float kl_delta(int node_index) {
    int parent_index = this->suffix_links[node_index / 2];
    if (parent_index == -1) {
      throw std::invalid_argument(
          "[CALCULATE DELTA] Given node does not have a parent.");
    }

    int size = seqan3::alphabet_size<alphabet_t>;

    float delta = 0;
    for (int i = 0; i < size; i++) {
      float prob = this->probabilities[node_index / 2][i];
      float parent_prob = this->probabilities[parent_index / 2][i];

      if (parent_prob == 0 || prob == 0) {
        continue;
      }

      delta += prob * std::log(prob / parent_prob);
    }
    delta *= this->get_counts(node_index);

    return delta;
  }

  float ps_delta(int node_index) {
    int parent_index = this->suffix_links[node_index / 2];
    if (parent_index == -1) {
      throw std::invalid_argument(
          "[CALCULATE DELTA] Given node does not have a parent.");
    }

    float delta = 0;
    float node_count = this->get_counts(node_index);
    float parent_count = this->get_counts(parent_index);

    std::array<int, seqan3::alphabet_size<alphabet_t>> node_child_counts =
        get_child_counts(node_index);

    std::array<int, seqan3::alphabet_size<alphabet_t>> parent_child_counts =
        get_child_counts(parent_index);

    for (int i = 0; i < seqan3::alphabet_size<alphabet_t>; i++) {
      float new_delta =
          std::abs(float(node_child_counts[i]) -
                   (float(parent_child_counts[i]) / parent_count) * node_count);

      if (delta < new_delta) {
        delta = new_delta;
      }
    }

    return delta;
  }

  std::array<int, seqan3::alphabet_size<alphabet_t>>
  get_child_counts(int node_index) {
    std::array<int, seqan3::alphabet_size<alphabet_t>> child_counts{};

    lst::details::iterate_children(
        node_index, this->table, this->flags, [&](int child_index) {
          int sequence_index = this->get_sequence_index(child_index);

          if (this->sequence[sequence_index] == seqan3::gap{}) {
            return;
          }

          int character_rank = seqan3::to_rank(this->sequence[sequence_index]);

          child_counts[character_rank] = get_counts(child_index);
        });

    return child_counts;
  }

  std::vector<int> get_terminal_nodes() {
    std::queue<int> top_down{};
    std::vector<int> terminal_nodes{};
    lst::details::iterate_children(0, this->table, this->flags,
                                   [&](int index) { top_down.push(index); });

    while (!top_down.empty()) {
      int node_index = top_down.front();
      top_down.pop();

      if (this->is_unevaluated(node_index) || this->is_excluded(node_index)) {
        continue;
      }
      if (this->is_terminal(node_index)) {
        terminal_nodes.push_back(node_index);
      }

      lst::details::iterate_children(node_index, this->table, this->flags,
                                     [&](int index) { top_down.push(index); });
    }

    return terminal_nodes;
  }

  bool is_terminal(int node_index) {
    for (int i = 0; i < seqan3::alphabet_size<alphabet_t>; i++) {
      int child_index = this->reverse_suffix_links[node_index / 2][i];
      if (child_index == 0) {
        continue;
      }

      if (this->is_included(child_index)) {
        return false;
      }
    }
    return true;
  }

  int get_counts(int node_index) {
    int c = this->counts[node_index / 2];
    if (c == 0) {
      c = lst::details::node_occurrences(node_index, this->table, this->flags);
      this->counts[node_index / 2] = c;
    }
    return c;
  }

  bool is_included(int node_index) {
    return (status[node_index / 2] & Status::Included) == Status::Included;
  }

  bool is_excluded(int node_index) {
    return (status[node_index / 2] & Status::Excluded) == Status::Excluded;
  }

  void append_node_string(int node_index, int lcp, int edge_lcp,
                          std::ostringstream &tree_string) {
    auto label_ = this->node_label(node_index, lcp, edge_lcp);

    if (this->is_leaf(node_index)) {
      label_ = this->leaf_label(node_index, lcp);
    }
    std::string label = label_ | seqan3::view::to_char;
    if (node_index == 0) {
      label = "#";
    }

    tree_string << "Node: " << node_index / 2 << " " << label << " ";

    append_reverse_child_counts(node_index, tree_string);

    tree_string << " " << get_counts(node_index) << " ";

    append_child_counts(node_index, tree_string);

    tree_string << " ";

    append_reverse_children(node_index, tree_string);

    tree_string << std::endl;
  }

  void append_child_counts(int node_index, std::ostringstream &tree_string) {
    tree_string << "[ ";

    auto child_counts = get_child_counts(node_index);
    for (auto c : child_counts) {
      tree_string << c << " ";
    }

    tree_string << "]";
  }

  void append_reverse_child_counts(int node_index,
                                   std::ostringstream &tree_string) {
    tree_string << "[ ";
    auto reverse_children = this->reverse_suffix_links[node_index / 2];
    for (int i = 0; i < seqan3::alphabet_size<alphabet_t>; i++) {
      int reverse_child = reverse_children[i];
      int count = get_counts(reverse_child);

      if (reverse_child == 0) {
        count = 0;
      }

      tree_string << count << " ";
    }
    tree_string << "]";
  }

  void append_reverse_children(int node_index,
                               std::ostringstream &tree_string) {
    tree_string << "[ ";
    auto reverse_children = this->reverse_suffix_links[node_index / 2];
    for (int i = 0; i < seqan3::alphabet_size<alphabet_t>; i++) {
      int reverse_child = reverse_children[i];

      if (reverse_child == 0 || this->is_excluded(reverse_child)) {
        reverse_child = -2;
      }

      tree_string << reverse_child / 2 << " ";
    }
    tree_string << "]";
  }

  int nodes_in_tree() {
    int n_nodes = 0;
    lst::details::breadth_first_iteration(
        this->sequence, this->suffixes, this->table, this->flags, false,
        [&](int node_index, int lcp, int edge_lcp) -> bool {
          if (is_excluded(node_index)) {
            return false;
          } else {
            n_nodes += 1;
            return true;
          }
        });

    return n_nodes;
  }
};
} // namespace pst
