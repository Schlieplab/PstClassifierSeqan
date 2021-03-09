#include "../src/probabilistic_suffix_tree_map.hpp"
#include <unordered_map>

std::unordered_map<std::string, size_t>
get_label_count_map(pst::ProbabilisticSuffixTree<seqan3::dna5> &tree) {
  std::unordered_map<std::string, size_t> map{};

  static std::mutex labels_mutex{};

  tree.breadth_first_iteration(
      [&](size_t node_index, size_t lcp, size_t edge_lcp, size_t node_count,
          lst::details::alphabet_array<size_t, seqan3::dna5> &child_counts)
          -> bool {
        std::lock_guard labels_lock{labels_mutex};
        auto label = tree.node_label(node_index, lcp, edge_lcp);

        map[label] = tree.get_counts(node_index);

        return true;
      });

  return map;
}