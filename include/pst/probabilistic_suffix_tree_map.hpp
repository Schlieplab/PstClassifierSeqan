#pragma once

#include <algorithm>
#include <array>
#include <atomic>
#include <chrono>
#include <cmath>
#include <ctime>
#include <fstream>
#include <functional>
#include <iostream>
#include <mutex>
#include <queue>
#include <random>
#include <shared_mutex>
#include <sstream>
#include <stack>
#include <string>
#include <string_view>
#include <thread>
#include <tuple>
#include <utility>
#include <vector>

#include "seqan3/alphabet/concept.hpp"
#include "seqan3/alphabet/nucleotide/dna4.hpp"
#include "seqan3/alphabet/nucleotide/dna5.hpp"
#include "seqan3/alphabet/views/all.hpp"
#include <filesystem>
#include "seqan3/std/ranges"
#include "seqan3/utility/views/convert.hpp"

#include "robin_hood.h"

#include "cereal/archives/binary.hpp"
#include "cereal/cereal.hpp"

#include "vlmc_from_kmers/kmer.hpp"

#include "search/lazy_suffix_tree.hpp"
#include "search/lazy_suffix_tree/iteration.hpp"

namespace pst {

template <seqan3::alphabet alphabet_t>
struct hashmap_value {
  size_t count;
  std::array<double, seqan3::alphabet_size<alphabet_t>> next_symbol_probabilities;
  bool is_included;
};

template <seqan3::alphabet alphabet_t>
using vec_t =
    std::deque<std::tuple<std::string, size_t,
                          lst::details::alphabet_array<size_t, alphabet_t>>>;

/*!\brief The probabilistic suffix tree implementation.
 * \tparam alphabet_t   Alphabet type from Seqan3.
 *
 * \details
 *
 * The pst::ProbabilisticSuffixTree is an extension of a lazy suffix tree which
 * is pruned to only contain statistically significant information.
 *
 * ### General information
 *
 * The probabilistic suffix tree is equivalent to a variable length Markov
 * chain.  The implementation follows the AV-2 algorithm from [Schulz et
 * al.](https://doi.org/10.1007/978-3-540-87361-7_26), which is based on a lazy
 * suffix tree.
 *
 * ### Example
 *
 * \include pst-classifier.cpp
 *
 **
 */
template <seqan3::alphabet alphabet_t = seqan3::dna5>
class ProbabilisticSuffixTreeMap : public lst::LazySuffixTree<alphabet_t> {
public:
  ProbabilisticSuffixTreeMap() = default;
  ProbabilisticSuffixTreeMap(ProbabilisticSuffixTreeMap const &) = default;
  ~ProbabilisticSuffixTreeMap() = default;

  /*!\brief Constructor which assumes default values for all parameters.
   * \param[in] id The id of the model.
   * \param[in] sequence The text to construct from.
   */
  ProbabilisticSuffixTreeMap(std::string id,
                             lst::details::sequence_t<alphabet_t> &sequence)

      : ProbabilisticSuffixTreeMap(id, sequence, 15, 100, 1.2, 0, "cutoff",
                                   true, 2) {}

  /*!\brief Constructor.
   * \param[in] id_ The id of the model.
   * \param[in] sequence_ The text to construct from.
   * \param[in] max_depth_ Max length of a branch/context in the tree.
   * \param[in] freq_ Min frequency of each context/node in the tree.
   * \param[in] number_of_parameters_ Number of parameters to keep in the model,
   * for "parameters" pruning.
   * \param[in] pruning_method_ Pruning method, either "cutoff" or "parameters"
   * (which depend on the `number_of_parameters_`)
   * \param[in] multi_core_ True for parallel execution.
   * \param[in] parallel_depth The maximum depth to spawn new processes, will
   * control task size as `alphabet_size ** depth`.
   */
  ProbabilisticSuffixTreeMap(std::string id_,
                             lst::details::sequence_t<alphabet_t> &sequence_,
                             size_t max_depth_, size_t freq_,
                             size_t number_of_parameters_,
                             std::string pruning_method_,
                             bool multi_core_ = true, int parallel_depth = 2,
                             const double pseudo_count_amount_ = 1.0)
      : lst::LazySuffixTree<alphabet_t>(sequence_, multi_core_, parallel_depth),
        id(std::move(id_)), freq(freq_), max_depth(max_depth_),
        number_of_parameters(number_of_parameters_),
        pruning_method(std::move(pruning_method_)),
        pseudo_count_amount(pseudo_count_amount_) {
    auto characters = this->get_valid_characters();
    for (auto &c : characters) {
      auto char_rank = seqan3::to_rank(c);
      valid_characters.insert(char_rank);
      valid_character_chars.push_back(c.to_char());
    }

    if (multi_core_ == false) {
      this->parallel_depth = 0;
    }
  }

  /*!\brief Reads a tree form a file.
   * \param[in] id The id of the model.
   * \param[in] sequence The text to construct from.
   */
  ProbabilisticSuffixTreeMap(const std::filesystem::path &filename,
                             const double pseudo_count_amount_ = 1.0)
      : pseudo_count_amount(pseudo_count_amount_) {
    auto characters = this->get_valid_characters();
    for (auto &c : characters) {
      auto char_rank = seqan3::to_rank(c);
      valid_characters.insert(char_rank);
      valid_character_chars.push_back(c.to_char());
    }

    if (filename.extension() == ".tree") {
      std::ifstream input(filename);
      if (!input.is_open()) {
        throw std::invalid_argument{"Failed to open file."};
      }

      for (std::string line; std::getline(input, line, '\n');) {
        parse_line(line, characters);
      }
    } else if (filename.extension() == ".bintree") {
      this->id = filename.stem().string();
      std::ifstream file_stream(filename, std::ios::binary);
      cereal::BinaryInputArchive iarchive(file_stream);

      vlmc::VLMCKmer kmer{};

      while (file_stream.peek() != EOF) {
        try {
          iarchive(kmer);
          double child_count = std::accumulate(kmer.next_symbol_counts.begin(),
                                               kmer.next_symbol_counts.end(),
                                               pseudo_count_amount_ * 4);
          this->counts[kmer.to_string()] = {
              kmer.count,
              {(double(kmer.next_symbol_counts[0]) + pseudo_count_amount_) /
                   child_count,
               (double(kmer.next_symbol_counts[1]) + pseudo_count_amount_) /
                   child_count,
               (double(kmer.next_symbol_counts[2]) + pseudo_count_amount_) /
                   child_count,
               0,
               (double(kmer.next_symbol_counts[3]) + pseudo_count_amount_) /
                   child_count},
              true};

        } catch (const cereal::Exception &e) {
          std::cout << (file_stream.peek() == EOF) << std::endl;
          std::cout << e.what() << std::endl;
        }
      }
      file_stream.close();
    }
    this->set_root_state();
  }

  /*!\brief Reads a from the tree format.
   * \param[in] sequence The text to construct from.
   */
  ProbabilisticSuffixTreeMap(const std::string &tree,
                             const double pseudo_count_amount_ = 1.0)
      : pseudo_count_amount(pseudo_count_amount_) {
    auto characters = this->get_valid_characters();

    for (auto &c : characters) {
      auto char_rank = seqan3::to_rank(c);
      valid_characters.insert(char_rank);
      valid_character_chars.push_back(c.to_char());
    }

    std::stringstream tree_stream{tree};

    for (std::string line; std::getline(tree_stream, line, '\n');) {
      this->parse_line(line, characters);
    }
    this->set_root_state();
  }

  void set_root_state() {
    auto root_iter = this->counts.find("");
    if (root_iter != this->counts.end()) {
      this->root_state = root_iter->second;
    }
  }

  /*! \brief Construct tree
   * Constructs the PST by applying the support and similarity pruning
   * steps from the VOMC construction algorithm.
   */
  void construct_tree() {
    if (this->sequence.size() <= 0) {
      throw std::invalid_argument{"Construct tree called without a sequence."};
    }
    this->support_pruning();
    this->similarity_pruning();
    this->set_root_state();
  }

  /*!\brief Debug printing
   * Prints the tree with the corresponding label, suffix links, delta values
   * etc.
   */
  void print() {
    std::string root{""};
    this->debug_print_node(root);
    std::cout << std::endl;

    this->pst_breadth_first_iteration(
        root, 0, [&](const std::string label, size_t level) -> bool {
          if (this->is_excluded(label)) {
            return true;
          }
          this->debug_print_node(label);

          std::cout << std::endl;
          return true;
        });
  }

  virtual void debug_print_node(const std::string &node_label) {
    std::cout << node_label << " ";

    std::cout << "\tDelta: " << this->calculate_delta(node_label);

    std::cout << "\tPST Leaf: " << this->is_pst_leaf(node_label);

    std::cout << "\tTerminal: " << this->is_terminal(node_label);

    std::cout << "\tStatus: " << this->is_included(node_label);
  }

  /*! \brief Outputs the tree in a .tree format.
   *
   * \details
   * The output string will be in the following format:
   * Name: <id of the tree>
   * Date: <time of generation>
   * Tree: PST
   * Alphabet: <type of alphabet>
   * Number(nodes): <int, number of nodes>
   * Number(parameters). <int, number of parameters (terminal nodes * |alphabet|
   * - 1) Node: <index> <label> [ <counts of children> ] <count of node> [
   * <counts of forward children> ][ <index of children> ]
   * ...
   *
   * \returns std::string with the tree format.
   */
  std::string to_tree() {
    std::ostringstream tree_string;

    auto now = std::chrono::system_clock::now();
    std::time_t time = std::chrono::system_clock::to_time_t(now);

    tree_string << "Name: " << this->id << std::endl;
    tree_string << "Date: " << std::ctime(&time);
    tree_string << "Tree: PST" << std::endl;
    tree_string << "Alphabet: " << lst::get_alphabet_name<alphabet_t>()
                << std::endl;
    tree_string << "Number(nodes): " << this->nodes_in_tree() << std::endl;

    auto n_parameters =
        this->count_terminal_nodes() * (valid_characters.size() - 1);
    tree_string << "Number(parameters): " << n_parameters << std::endl;

    robin_hood::unordered_map<std::string, size_t> iteration_order_indices{};
    size_t i = 0;
    // old pst-classifier compatibility:
    //    this->pst_breadth_first_iteration(
    //        [&](const std::string child_label, size_t level) -> bool {
    //          iteration_order_indices[child_label] = i;
    //          i++;
    //          return true;
    //        });

    this->pst_breadth_first_iteration(
        [&](const std::string child_label, size_t level) -> bool {
          this->append_node_string(child_label, iteration_order_indices,
                                   tree_string);
          return true;
        });

    return tree_string.str();
  }

  /**! \brief Counts the number of terminal nodes in the tree.
   *
   * \return vector of indices to all terminal nodes.
   */
  size_t count_terminal_nodes() {
    std::atomic<size_t> n_terminal_nodes = 0;

    for (auto &[node, v] : this->counts) {
      bool included = v.is_included;
      if (included && this->is_terminal(node)) {
        n_terminal_nodes += 1;
      }
    }

    return n_terminal_nodes;
  }

  std::vector<std::string> terminal_nodes{};
  /**! \brief Get the terminal nodes in the tree.
   *
   * \return vector of all terminal nodes, sorted.
   */
  std::vector<std::string> get_terminal_nodes() {
    static std::shared_mutex terminal_mutex{};
    terminal_mutex.lock_shared();
    if (!terminal_nodes.empty()) {
      terminal_mutex.unlock_shared();
      return terminal_nodes;
    }
    terminal_mutex.unlock_shared();

    std::vector<std::string> nodes{};
    for (const auto &[child_label, v] : this->counts) {
      bool included = v.is_included;
      if (included && this->is_terminal(child_label)) {
        nodes.push_back(child_label);
      }
    }

    std::sort(nodes.begin(), nodes.end());

    std::lock_guard terminal_lock{terminal_mutex};
    terminal_nodes = nodes;

    return nodes;
  }

  std::vector<std::string> sorted_contexts{};
  /**! \brief Get the terminal nodes in the tree.
   *
   * \return vector of all terminal nodes, sorted.
   */
  std::vector<std::string> get_sorted_contexts() {
    static std::shared_mutex context_mutex{};
    context_mutex.lock_shared();
    if (!sorted_contexts.empty()) {
      context_mutex.unlock_shared();
      return sorted_contexts;
    }
    context_mutex.unlock_shared();

    std::vector<std::string> nodes{};
    for (const auto &[child_label, v] : this->counts) {
      bool included = v.is_included;
      if (included) {
        nodes.push_back(child_label);
      }
    }

    std::sort(nodes.begin(), nodes.end());

    std::lock_guard terminal_lock{context_mutex};
    sorted_contexts = nodes;

    return nodes;
  }

  std::string generate_sequence(size_t length) {
    std::random_device rd;
    std::mt19937 gen = std::mt19937{rd()};

    using seqan3::operator""_dna5;

    static constexpr char char_codes[4] = {'A', 'C', 'G', 'T'};

    std::uniform_int_distribution<> distrib(0, 3);

    size_t order_max = this->get_max_order();
    std::stringstream ss;

    std::string subsequence{};
    subsequence.resize(order_max);
    for (size_t i = 0; i < length; i++) {
      size_t start_index = i - order_max;
      size_t end_index = i;
      if (order_max > i) {
        start_index = 0;
        end_index = i;
      }

      auto [_context, val] = this->get_closest_state(subsequence);
      auto probs = val.next_symbol_probabilities;
      std::discrete_distribution<> d({probs[0], probs[1], probs[2], probs[4]});
      auto rand = d(gen);

      char c = char_codes[rand];

      ss << c;
      subsequence.push_back(c);
      if (subsequence.length() > order_max) {
        subsequence.erase(0, 1);
      }
    }

    return ss.str();
  }

  /**! \brief Iterates over the nodes in the PST.
   *
   * \param start_index Index of the node to start iteration at.
   *
   */
  void pst_breadth_first_iteration(
      const std::function<bool(const std::string, size_t)> &f) {
    pst_breadth_first_iteration("", 0, f);
  }

  /**! \brief Iterates over the nodes in the PST.
   *
   * \param start_index Index of the node to start iteration at.
   * \param start_level Starting level for the iteration (0 for root, 1 for
   * ACGT) \param f Callback for each child of start_index in a breadth first
   * fashion.  Is given node_index and depth of the node and should return
   * wether to iterate over the node's children.
   *
   */
  void pst_breadth_first_iteration(
      const std::string &start_label, const size_t start_level,
      const std::function<bool(const std::string, size_t)> &f) {
    std::queue<std::tuple<std::string, size_t>> queue{};

    queue.emplace(start_label, start_level);

    while (!queue.empty()) {
      auto &[label, level] = queue.front();

      std::vector<std::string> children{};
      this->iterate_pst_children(label, [&](const std::string child_label) {
        children.emplace_back(std::move(child_label));
      });

      if (f(std::move(label), level)) {
        for (auto &child : children) {
          queue.emplace(std::move(child), level + 1);
        }
      }
      queue.pop();
    }
  }

  void breadth_first_iteration_p(
      const std::function<bool(const std::string, size_t)> &f) {
    std::string root{""};
    this->breadth_first_iteration_p(root, 0, f);
  }

  void breadth_first_iteration_p(
      const std::string &start_label, const size_t start_level,
      const std::function<bool(const std::string, size_t)> &f) {
    breadth_first_iteration_p_(start_label, start_level, f);
  }

  void breadth_first_iteration_p_(
      const std::string start_label, const size_t start_level,
      const std::function<bool(const std::string, size_t)> &f) {
    std::vector<std::thread> threads{};

    std::queue<std::tuple<std::string, size_t>> queue{};

    queue.emplace(start_label, start_level);

    while (!queue.empty()) {
      auto &[label, level] = queue.front();

      std::vector<std::string> children{};
      this->iterate_pst_children(label, [&](std::string child_label) {
        if (this->is_included(child_label)) {
          children.emplace_back(std::move(child_label));
        }
      });

      if (f(std::move(label), level)) {
        if (this->multi_core && this->parallel_depth > level) {
          for (auto &child : children) {
            threads.emplace_back(
                &ProbabilisticSuffixTreeMap::breadth_first_iteration_p_, this,
                std::move(child), level + 1, f);
          }
        } else {
          for (auto &child : children) {
            queue.emplace(std::move(child), level + 1);
          }
        }
      }
      queue.pop();
    }

    for (auto &thread : threads) {
      if (thread.joinable()) {
        thread.join();
      }
    }
  }

  /**! \brief Returns the index of the children of in the PST.
   *
   * \param node_index Node to retrieve children of.
   * \return A vector of the children.
   */
  std::vector<std::string> get_pst_children(const std::string &label) {
    std::vector<std::string> children{};
    this->iterate_pst_children(label, [&](const std::string child_label) {
      children.push_back(std::move(child_label));
    });
    return children;
  }

  /**! \brief Return the parent of node in the PST.
   *
   * \param node_index Node to get parent of.
   * \return Parent index, or -1 if no parent.
   */
  std::string get_pst_parent(const std::string &label) {
    return label.substr(1);
  }

  /**! \brief Remove nodes that are less frequent that min_count and longer than
   * max_depth.
   *
   * For any effect, the new settings have to be more stringent than the
   * previous settings.
   *
   * \param min_count_ new min count setting.
   * \param max_depth_ new max depth setting.
   */
  void reprune_support(size_t min_count_, int max_depth_) {
    this->pst_breadth_first_iteration(
        [&](const std::string context, size_t level) -> bool {
          if (this->get_count(context) < min_count_ || level > max_depth_) {
            this->counts[context].is_included = false;
          }
          return true;
        });
  }


  std::tuple<std::string, hashmap_value<alphabet_t>&> get_closest_state(const std::string_view &label) {
    if (this->root_state.count == 0) {
      this->set_root_state();
    }

    if (label.empty()) {
      return {"", this->root_state};
    }

    std::string sublabel{label};
    for (size_t i = 0; i < label.size(); i++) {
      auto [is_included, val] = this->is_included_get_state(sublabel);
      if (is_included) {
        return {sublabel, val};
      }
      sublabel.erase(0, 1);
    }

    return {"", this->root_state};
  }

  std::tuple<std::string, hashmap_value<alphabet_t>&> get_closest_state(const std::string &label) {
    if (this->root_state.count == 0) {
      this->set_root_state();
    }

    if (label.empty()) {
      return {"", this->root_state};
    }

    auto sublabel = label;
    for (size_t i = 0; i < label.size(); i++) {
      auto [is_included, val] = this->is_included_get_state(sublabel);
      if (is_included) {
        return {sublabel, val};
      }
      sublabel.erase(0, 1);
    }

    return {"", this->root_state};
  }

  double get_transition_probability(const hashmap_value<alphabet_t> &val,
                                    const size_t &char_rank) {
    return val.next_symbol_probabilities[char_rank];
  }

  double get_transition_probability(const std::string &label,
                                    const size_t &char_rank) {
    return this->counts[label].next_symbol_probabilities[char_rank];
  }


  double get_transition_probability(const hashmap_value<alphabet_t> &val,
                                    const char &character) {
    auto c = seqan3::assign_char_to(character, alphabet_t{});
    auto char_rank = c.to_rank();
    return val.next_symbol_probabilities[char_rank];
  }

  double get_transition_probability(const std::string &label,
                                    const char &character) {
    auto c = seqan3::assign_char_to(character, alphabet_t{});
    auto char_rank = c.to_rank();
    return this->counts[label].next_symbol_probabilities[char_rank];
  }


  size_t get_max_order() {
    if (this->max_order != max_size) {
      return this->max_order;
    }

    size_t max_order_ = 0;
    for (const auto &[label, v] : this->counts) {
      bool included = v.is_included;
      if (included) {
        max_order_ = std::max(label.size(), max_order_);
      }
    }

    this->max_order = max_order_;

    return max_order_;
  }

  /**! \brief Support pruning phase of the algorithm.
   * \details
   * Extends a lazy suffix tree as long as the counts of each node is above
   * `freq` and the length is at most `max_depth`.
   *
   * After the full tree is built, it adds implicit nodes and suffix links.
   *
   */
  void support_pruning() {
    std::string root{""};
    this->counts[root] = {this->sequence.size(), {}, true};
    this->build_tree();
  }

  /**! \brief Similarity pruning phase of the algorithm
   * \details
   * Prunes the tree bottom-up (by maintaining a status for each node).  A delta
   * value is assigned to each node, and the nodes with the smallest delta value
   * are removed.  The pruning either stops when a fixed number of parameters
   * are reached, or until a specified threshold value.
   */
  void similarity_pruning() {
    // Probabilities are computed during construction.
    //    this->compute_probabilities();

    if (this->pruning_method == "cutoff") {
      this->cutoff_prune();
    } else if (this->pruning_method == "parameters") {
      this->parameters_prune();
    }
  }

  std::string id;

  robin_hood::unordered_map<std::string, hashmap_value<alphabet_t>> counts{};
  hashmap_value<alphabet_t> root_state{};

  robin_hood::unordered_set<size_t> valid_characters{};
  std::vector<char> valid_character_chars{};
  size_t max_order = max_size;

  friend class ProbabilisticSuffixTreeTest;
  size_t freq;
  size_t max_depth;
  size_t number_of_parameters;

  double pseudo_count_amount = 1.0;

  std::string pruning_method;

  thread_local static inline vec_t<alphabet_t> thread_vec{};

  /**! \brief Set the characters that are considered valid.
   * This needs to be modified to accept a custom set of valid characters.
   * \return vector of valid characters.
   */
  std::vector<seqan3::gapped<alphabet_t>> get_valid_characters() const {
    using seqan3::operator""_dna4;
    seqan3::dna4_vector dna4{"ACGT"_dna4};

    std::vector<seqan3::gapped<alphabet_t>> characters{
        dna4 | seqan3::views::convert<seqan3::gapped<alphabet_t>> |
        ranges::_to_::to<std::vector<seqan3::gapped<alphabet_t>>>};

    return characters;
  }

  bool build_tree_callback(
      size_t sequence_index, size_t lcp, size_t edge_lcp, size_t node_count,
      lst::details::alphabet_array<size_t, alphabet_t> &child_counts,
      bool is_leaf) {
    if (is_leaf && edge_lcp == 1) {
      return false;
    }

    if (edge_lcp == 0) {
      return true;
    }

    bool include_node = true;

    for (size_t i = 1; i <= edge_lcp && include_node; i++) {
      bool include_sub_node =
          this->include_node(sequence_index, lcp, i, node_count);

      auto label = get_label(sequence_index, lcp, i);

      // If we're expanding the node, the child counts will be the count
      // of this node, but at the next character
      lst::details::alphabet_array<size_t, alphabet_t> cc{};
      auto next_node_end = std::min(sequence_index + i, this->sequence.size());
      auto character_rank = this->get_character_rank(next_node_end);
      cc[character_rank] = node_count;

      if (include_sub_node) {
        thread_vec.emplace_back(std::move(label), node_count, std::move(cc));
      }

      include_node = include_sub_node;
    }

    if (include_node) {
      // If we added nodes, add the correct counts to the last one.
      std::get<2>(thread_vec.back()) = std::move(child_counts);
    }

    return include_node;
  }

  std::string get_label(size_t sequence_index, size_t lcp, size_t edge_lcp) {
    auto node_start = sequence_index - lcp;
    auto node_end = std::min(sequence_index + edge_lcp, this->sequence.size());

    lst::details::sequence_t<alphabet_t> label_dna(
        this->sequence.begin() + node_start, this->sequence.begin() + node_end);

    std::string label =
        label_dna | seqan3::views::to_char | seqan3::views::to<std::string>;

    return label;
  }

  void merge_tree_callback(std::shared_mutex &counts_mutex) {
    {
      std::unique_lock lock{counts_mutex};
      for (auto &[label, node_count, _cc] : thread_vec) {
        this->counts[label] = {node_count, {}, true};
      }
    }

    for (auto &[label, _nc, child_counts] : thread_vec) {
      size_t child_sum =
          std::accumulate(child_counts.begin(), child_counts.end() - 1, 0);
      if (child_sum == 0) {
        std::shared_lock count_lock{counts_mutex};
        this->assign_node_probabilities(label);
      } else {
        this->assign_node_probabilities(label, child_counts, counts_mutex);
      }
    }
  }

  /**! \brief Builds the tree top-down.
   * \details
   * The lazy suffix tree is iterated in a breadth-first fashion and the nodes
   * are saved if the counts of each node is above `freq` and the length is at
   * most `max_depth`.
   *
   * Also saves the count of the node as well as if it is included or excluded.
   */
  void build_tree() {
    // Reset main thread.
    thread_vec = vec_t<alphabet_t>{};
    std::shared_mutex counts_mutex{};

    using namespace std::placeholders;
    auto build_callback = [&](auto &&...args) -> bool {
      return this->build_tree_callback(args...);
    };

    auto merge_callback = [&]() { merge_tree_callback(counts_mutex); };

    this->breadth_first_iteration_table_less(build_callback, merge_callback);
    this->assign_node_probabilities("");
  }

  /**! \brief Computes and saves the forward probabilities of each node.
   * \details
   * Probabilities for each node are calculated over the sum of the counts
   * of each (forward) child node.  Pseudo counts are added to all counts.
   *
   */
  void compute_probabilities() {
    this->assign_node_probabilities("");

    this->breadth_first_iteration_p([&](const std::string label, size_t level) {
      this->assign_node_probabilities(label);
      return true;
    });
  }

  /**! \brief Assigns probabilities for the node.
   * \details
   * Sums the counts of all children, with pseudo counts, and assigns the
   * corresponding probabilities.  Iterates over all children to find counts,
   * to sum those counts, and to assign those counts.
   *
   * \param[in] label
   */
  void assign_node_probabilities(const std::string &label) {
    auto child_counts = this->get_child_counts(label);

    double child_sum =
        std::accumulate(child_counts.begin(), child_counts.end(), 0.0);

    for (size_t i = 0; i < seqan3::alphabet_size<alphabet_t>; i++) {
      this->counts[label].next_symbol_probabilities[i] = double(child_counts[i]) / child_sum;
    }
  }

  /**! \brief Assigns probabilities for the node.
   * \details
   * Sums the counts of all children, with pseudo counts, and assigns the
   * corresponding probabilities.  Iterates over all children to find counts,
   * to sum those counts, and to assign those counts.
   *
   * \param[in] label
   */
  void assign_node_probabilities(
      const std::string &label,
      lst::details::alphabet_array<size_t, alphabet_t> &child_counts,
      std::shared_mutex &counts_mutex) {

    double child_sum = 0.0;
    for (auto char_rank : this->valid_characters) {
      child_sum += child_counts[char_rank];
    }
    child_sum += 4;

    std::shared_lock lock{counts_mutex};
    for (auto char_rank : this->valid_characters) {
      this->counts[label].next_symbol_probabilities[char_rank] =
          double(child_counts[char_rank] + 1) / child_sum;
    }
  }

  /**! \brief Removes all nodes from the tree with a delta value below
   * threshold.
   * \details
   * The full tree is iterated to find the leaves in the
   * PST.  These leaves are then iterated, for each leaf the delta value is
   * calculated, and the node is removed if the delta value is below a
   * threshold.  For each removed node, the parent is added to be considered if
   * it now a leaf.
   */
  virtual void cutoff_prune() { parameters_prune(); }
  virtual double calculate_delta(const std::string &node_label) { return 0.0; }

  /**! \brief Removes all nodes until a specified number of parameters are left.
   *
   * \details The full tree is iterated to find the leaves in the
   * PST.  These leaves are then iterated, for each leaf the delta value is
   * calculated.  These leaves are then removed in order of smallest
   * delta value until a given number of parameters are left in the tree.
   * For each removed node, the parent is added to be considered if
   * it now a leaf.
   */
  void parameters_prune() {
    auto pst_leaves = this->get_pst_leaves();

    using q_t = std::tuple<std::string, double>;
    auto cmp = [](q_t left, q_t right) -> bool {
      return std::get<1>(left) < std::get<1>(right);
    };
    std::priority_queue<q_t, std::vector<q_t>, decltype(cmp)> queue{cmp};

    for (auto &v : pst_leaves) {
      queue.emplace(v, -this->calculate_delta(v));
    }

    auto n_terminal_nodes = this->count_terminal_nodes();

    auto alphabet_size = this->valid_characters.size() - 1;
    while (!queue.empty() &&
           n_terminal_nodes * alphabet_size > this->number_of_parameters) {
      auto &[node_label, delta] = queue.top();

      if (node_label.empty()) {
        queue.pop();
        continue;
      }
      this->counts[node_label].is_included = false;

      const std::string parent_label = this->get_pst_parent(node_label);

      if (this->is_pst_leaf(parent_label) && !parent_label.empty()) {
        queue.emplace(parent_label, -this->calculate_delta(parent_label));
      }

      if (!this->became_terminal(parent_label, node_label)) {
        n_terminal_nodes -= 1;
      }
      queue.pop();
    }
  }

  /**! \brief Finds the counts of all (forward) children for the node.
   *
   * \param label std::string to get child counts for.
   * \return Array of counts for each children.
   */
  std::array<size_t, seqan3::alphabet_size<alphabet_t>>
  get_child_counts(const std::string &label) {
    std::array<size_t, seqan3::alphabet_size<alphabet_t>> child_counts{};

    std::string child_label = label + ' ';
    alphabet_t c{};
    for (auto char_rank : this->valid_characters) {
      seqan3::assign_rank_to(char_rank, c);
      child_label[child_label.size() - 1] = c.to_char();

      auto count = this->get_count(child_label);
      child_counts[char_rank] = count + this->pseudo_count_amount;
    }

    return child_counts;
  }

  /**! \brief Finds all leaves in the pst.
   *
   * \return vector of indices to all pst leaves.
   */
  std::vector<std::string> get_pst_leaves() {
    std::vector<std::string> bottom_nodes{};

    static std::mutex leaves_mutex{};
    this->breadth_first_iteration_p([&](const std::string label, size_t level) {
      if (this->is_pst_leaf(label)) {
        std::lock_guard lock{leaves_mutex};
        bottom_nodes.emplace_back(std::move(label));
      }
      return true;
    });

    return bottom_nodes;
  }

  /**! \brief Finds the counts of all (forward) children for the node.
   *
   * \param label std::string to get child counts for.
   * \return Array of counts for each children.
   */
  std::array<size_t, seqan3::alphabet_size<alphabet_t>>
  get_prob_estimated_child_counts(const std::string &label) {
    std::array<size_t, seqan3::alphabet_size<alphabet_t>> child_counts{};

    auto &[count, probs, included] = this->counts[label];
    size_t count_with_pseudo = count + 4 * this->pseudo_count_amount;
    for (auto char_rank : this->valid_characters) {
      child_counts[char_rank] = std::round(
          probs[char_rank] * count_with_pseudo - this->pseudo_count_amount);
    }

    return child_counts;
  }

  /**! \brief Checks if the node is a leaf in the PST.
   * \details
   * A node is a leaf in the PST if all of the reverse children corresponding
   * to the valid characters are excluded.
   *
   * \param[in] node_index Node to check.
   * \return If all children are missing.
   */
  bool is_pst_leaf(const std::string &node_label) {
    if (this->is_excluded(node_label)) {
      return false;
    }

    std::string child_label{' ' + node_label};

    for (auto char_rank : this->valid_characters) {
      alphabet_t c = seqan3::assign_rank_to(char_rank, alphabet_t{});
      child_label[0] = c.to_char();

      if (this->is_included(child_label)) {
        return false;
      }
    }

    return true;
  }

  /**! \brief Checks if the node is terminal.
   * \details
   * A node is terminal if any of the reverse children corresponding to
   * the valid characters are excluded.
   *
   * \param[in] node_label Node label to check.
   * \return If the node is terminal (has any missing children).
   */
  bool is_terminal(const std::string &node_label) {
    std::string child_label = ' ' + node_label;
    for (auto c : this->valid_character_chars) {
      child_label[0] = c;

      if (this->is_excluded(child_label)) {
        return true;
      }
    }
    return false;
  }

  /**! \brief Checks if removed_index is the only child that is missing.
   * \details
   * The node_index has to correspond to the parent of removed_index.
   * Iterates through the valid characters of the reverse children of
   * node_index, and checks if any child other than removed_index is excluded.
   *
   * \param[in] node_index Index of node to check.
   * \param[in] removed_index Index of child to node, recently excluded.
   * \return True if removed index is the only missing child of node_index
   */
  bool became_terminal(const std::string &node_label,
                       const std::string &removed_label) {
    for (auto c : this->valid_character_chars) {
      std::string child_label = c + node_label;

      if (this->is_excluded(child_label) && child_label != removed_label) {
        return false;
      }
    }

    return is_terminal(node_label);
  }

  std::tuple<bool, hashmap_value<alphabet_t>&> is_included_get_state(const std::string &label) {
    auto iter = this->counts.find(label);
    if (iter != this->counts.end()) {
      bool is_included = iter->second.is_included;
      return {is_included, iter->second};
    }

    hashmap_value<alphabet_t> v{};
    return {false, v};
  }

  bool is_included(const std::string &label) {
    return this->counts.find(label) != this->counts.end() &&
           this->counts[label].is_included;
  }

  bool is_excluded(const std::string &label) {
    return this->counts.find(label) == this->counts.end() ||
           !this->counts[label].is_included;
  }

  /**! \brief Append string for a node to the output stream.
   * \details
   * The format is Node: (index / 2) [ <reverse_children_count> ] count [
   * <children_count>] [ <reverse child indices> ]
   *
   * \param[in] label The node label to operate on
   * \param tree_string The output stream to write to.
   */
  void append_node_string(
      const std::string &label,
      robin_hood::unordered_map<std::string, size_t> &iteration_order_indices,
      std::ostringstream &tree_string) {
    // old pst-classifier compatibility:
    // tree_string << "Node: " << iteration_order_indices[label];

    tree_string << "Node: ";
    if (label.empty()) {
      tree_string << " # ";
    } else {
      tree_string << " " << label << " ";
    }

    append_reverse_child_counts(label, tree_string);

    tree_string << " " << this->get_count(label) << " ";

    append_child_counts(label, tree_string);

    //    tree_string << " ";

    append_reverse_children(label, iteration_order_indices, tree_string);

    tree_string << std::endl;
  }

  /**! \brief Append the count of forward children to the output stream.
   * \details
   * The format is [ x y z ... ]
   *
   * \param[in] node_index The node index to operate on.
   * \param tree_string The output stream to write to.
   */
  void append_child_counts(const std::string &label,
                           std::ostringstream &tree_string) {
    auto child_counts = this->get_prob_estimated_child_counts(label);

    std::vector<size_t> output(seqan3::alphabet_size<alphabet_t>, -1);
    for (auto char_rank : this->valid_characters) {
      output[char_rank] = child_counts[char_rank];
    };

    tree_string << "[ ";

    for (auto count : output) {
      if (count == -1) {
        continue;
      }
      tree_string << count << " ";
    }

    tree_string << "]";
  }

  /**! \brief Append the count of reverse children to the output stream.
   * \details
   * The format is [ x y z ... ]
   *
   * \param[in] node_label The node label to operate on.
   * \param tree_string The output stream to write to.
   */
  void append_reverse_child_counts(const std::string &node_label,
                                   std::ostringstream &tree_string) {
    std::vector<size_t> output(seqan3::alphabet_size<alphabet_t>, -1);

    for (auto char_rank : this->valid_characters) {
      alphabet_t c = seqan3::assign_rank_to(char_rank, alphabet_t{});

      std::string child_label = c.to_char() + node_label;

      auto count = this->get_count(child_label);

      output[char_rank] = count;
    };

    tree_string << "[ ";
    for (auto count : output) {
      if (count == -1) {
        continue;
      }
      tree_string << count << " ";
    }
    tree_string << "]";
  }

  /**! \brief Append the index of reverse children to the output stream.
   * \details
   * For excluded nodes, a -1 is appended, for included nodes, the
   * internal node index / 2 is used (which to be fair is arbitrary).
   * The format is [ x y z ... ]
   *
   * \param[in] node_index The node index to operate on.
   * \param tree_string The output stream to write to.
   */
  void append_reverse_children(
      const std::string &node_label,
      robin_hood::unordered_map<std::string, size_t> &iteration_order_indices,
      std::ostringstream &tree_string) {
    std::vector<size_t> output(seqan3::alphabet_size<alphabet_t>, -2);

    for (auto char_rank : this->valid_characters) {
      alphabet_t c = seqan3::assign_rank_to(char_rank, alphabet_t{});

      std::string child_label = c.to_char() + node_label;

      if (child_label.empty() || this->is_excluded(child_label)) {
        output[char_rank] = max_size;
      } else {
        // old pst-classifier compatibility:
        // output[char_rank] = iteration_order_indices[child_label];
        output[char_rank] = char_rank;
      }
    };

    tree_string << "[ ";
    for (auto count : output) {
      if (count == -2) {
        continue;
      } else if (count == max_size) {
        tree_string << "-1 ";
      } else {
        tree_string << count << " ";
      }
    }
    tree_string << "]";
  }

  /**! \brief Counts the nodes in the tree.
   *
   * \return The number of nodes in the tree.
   */
  size_t nodes_in_tree() {
    std::atomic<size_t> n_terminal_nodes = 0;

    for (auto &[node, v] : this->counts) {
      bool included = v.is_included;
      if (included) {
        n_terminal_nodes += 1;
      }
    }

    return n_terminal_nodes;
  }

  /**! \brief Determine if a node should be included in the tree.
   * \details
   * A node is included if the following three conditions are true:
   * - The length of the label is shorter than the max_depth of the tree.
   * - The label occurs at least as often as the specified min frequency.
   * - Every character in the label is included in the valid_characters.
   *
   * \param label_start start index of label in sequence.
   * \param label_end end index of label in sequence.
   * \param count number of times the label occurs in the sequence.
   * \return true if the node should be included.
   */
  bool include_node(size_t sequence_index, size_t lcp, size_t edge_lcp,
                    size_t count) {
    auto label_start = sequence_index - lcp;
    auto label_end = sequence_index + edge_lcp;
    auto label_length = label_end - label_start;

    // All characters before sequence_index have already been checked.
    return label_length <= this->max_depth && count >= this->freq &&
           label_valid(sequence_index, label_end);
  }

  /**! \brief Checks if a label is valid.
   *
   * \param[in] label_start starting index of the label in the sequence.
   * \param[in] label_end ending index of the label in the sequence.
   * \return true if the label is valid, false otherwise.
   */
  bool label_valid(size_t label_start, size_t label_end) {
    for (auto i = label_start; i < label_end; i++) {
      auto character_rank = this->get_character_rank(i);
      if (this->valid_characters.find(character_rank) ==
          this->valid_characters.end()) {
        return false;
      }
    }

    return true;
  }

  /** \copydoc lst::LazySuffixTree.skip_node()
   */
  bool skip_node(const std::string &node_label) {
    return this->is_excluded(node_label);
  }

  /**! \brief Iterates through all possible suffix tree children of the node
   * with label.
   *
   * \param[in] label string label of node to iterate children of.
   * \param f function to call on every child.
   */
  void iterate_children(const std::string &label,
                        const std::function<void(const std::string &)> &f) {
    for (auto char_rank : this->valid_characters) {
      alphabet_t c = seqan3::assign_rank_to(char_rank, alphabet_t{});

      std::string child_label = label + c.to_char();
      f(child_label);
    }
  }

  /**! \brief Iterates through all possible prefix tree children of the node
   * with label.
   *
   * \param[in] label string label of node to iterate children of.
   * \param f function to call on every child.
   */
  void iterate_pst_children(const std::string &label,
                            const std::function<void(std::string)> &f) {
    for (auto char_rank : this->valid_characters) {
      alphabet_t c = seqan3::assign_rank_to(char_rank, alphabet_t{});

      std::string child_label{c.to_char() + label};

      if (is_included(child_label)) {
        f(std::move(child_label));
      }
    }
  }

  /** Parses a line in the tree format.
   *
   * Adds the node label, counts and probabilities corresponding to the line.
   *
   * \param line Line to parse.
   * \param characters The valid characters, has to agree with input.
   */
  void parse_line(const std::string &line,
                  std::vector<seqan3::gapped<alphabet_t>> &characters) {
    if (line.substr(0, 5) == "Node:") {
      parse_node(line, characters);
    } else if (line.substr(0, 5) == "Name:") {
      this->id = line.substr(6);
    }
  }

  void parse_node(const std::string &line,
                  std::vector<seqan3::gapped<alphabet_t>> &characters) {
    std::stringstream line_stream{line};

    std::string node_label{"-1"};
    size_t node_count;

    std::string word, prev_word;

    bool found_label = false;
    bool found_count = false;
    bool found_probs = false;
    size_t prob_index = 0;

    while (line_stream >> word) {
      if (word == "[" && !found_label) {
        node_label = prev_word;
        if (node_label == "#") {
          node_label = "";
        }
        found_label = true;

      } else if (prev_word == "]" && !found_count) {
        node_count = std::stoi(word);
        this->counts[node_label] = {node_count, {}, true};
        found_count = true;

      } else if (found_count && !found_probs) {
        if (word[0] == ']') {
          found_probs = true;
        } else if (word != "[") {
          double child_count = std::stof(word);
          auto c = characters[prob_index];
          auto char_rank = seqan3::to_rank(c);
          this->counts[node_label].next_symbol_probabilities[char_rank] = child_count;
          prob_index++;
        }
      }
      prev_word = std::move(word);
    }
    convert_counts_to_probabilities(node_label);
  }

  /** Utility for parsing, converts the raw counts initially stored in the
   * probabilities map to probabilities.
   *
   * @param node_label Label to convert counts for.
   */
  void convert_counts_to_probabilities(std::string &node_label) {
    double child_sum =
        std::accumulate(this->counts[node_label].next_symbol_probabilities.begin(),
                        this->counts[node_label].next_symbol_probabilities.end(),
                        this->pseudo_count_amount * 4);

    for (size_t i = 0; i < seqan3::alphabet_size<alphabet_t>; i++) {
      this->counts[node_label].next_symbol_probabilities[i] =
          (double(this->counts[node_label].next_symbol_probabilities[i]) +
           this->pseudo_count_amount) /
          child_sum;
    }
  }

  size_t get_count(const std::string &label) {
    auto search = this->counts.find(label);

    if (search == this->counts.end()) {
      return 0;
    } else {
      return search->second.count;
    }
  }
};
} // namespace pst
