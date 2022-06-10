#pragma once

#include "seqan3/alphabet/nucleotide/dna5.hpp"
#include "seqan3/alphabet/views/all.hpp"
#include "seqan3/utility/views/convert.hpp"

#include "indicators/block_progress_bar.hpp"
#include <highfive/H5File.hpp>

#include "pst/probabilistic_suffix_tree_map.hpp"

using tree_t = pst::ProbabilisticSuffixTreeMap<seqan3::dna5>;

std::vector<tree_t> get_trees(HighFive::File &file,
                              const double pseudo_count_amount) {
  const std::string DATASET_NAME("signatures");

  HighFive::DataSet dataset = file.getDataSet(DATASET_NAME);

  std::vector<std::string> result_string_list;
  dataset.read(result_string_list);

  indicators::BlockProgressBar bar{
      indicators::option::BarWidth{50},
      indicators::option::PrefixText{"Parsing trees"},
      indicators::option::ShowElapsedTime{true},
      indicators::option::ShowRemainingTime{true}};

  std::vector<tree_t> trees{};
  ulong n_trees = dataset.getDimensions()[0];
  for (ulong i = 0; i < n_trees; i++) {
    dataset.select(HighFive::ElementSet({i})).read(result_string_list);
    //    std::cout << result_string_list[0] << std::endl;
    trees.push_back(tree_t{result_string_list[0], pseudo_count_amount});
    bar.set_progress((double(i) / double(n_trees)) * 100.0);
  }

  bar.set_progress(100.0);
  bar.mark_as_completed();
  return trees;
}

std::vector<tree_t>
get_trees_from_directory(const std::filesystem::path &directory,
                         const double pseudo_count_amount) {
  std::vector<tree_t> trees{};
  for (auto const &dir_entry : std::filesystem::directory_iterator{directory}) {
    trees.push_back(tree_t{dir_entry.path(), pseudo_count_amount});
  }

  return trees;
}

std::vector<tree_t> get_trees(std::filesystem::path &path,
                              const double pseudo_count_amount) {
  if (std::filesystem::is_directory(path)) {
    return get_trees_from_directory(path, pseudo_count_amount);
  } else if (path.extension() == ".h5" || path.extension() == ".hdf5") {
    HighFive::File file{path, HighFive::File::ReadOnly};
    return get_trees(file, pseudo_count_amount);
  } else if (path.extension() == ".tree" || path.extension() == ".bintree") {
    pst::ProbabilisticSuffixTreeMap<seqan3::dna5> tree{path,
                                                       pseudo_count_amount};

    return std::vector<tree_t>{std::move(tree)};
  } else {
    std::stringstream message;
    message << "Error: " << path
            << "has invalid extension, and can't be parsed.";

    throw std::invalid_argument(message.str());
  }
}

std::string to_string(std::vector<seqan3::dna5> &sequence) {
  std::string seq =
      sequence | seqan3::views::to_char | seqan3::views::to<std::string>;
  return seq;
}