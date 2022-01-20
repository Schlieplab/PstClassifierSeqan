#include <algorithm>
#include <cstdlib>
#include <functional>
#include <iostream>
#include <optional>
#include <string>
#include <thread>

#include <Eigen/Dense>

#include <highfive/H5File.hpp>

#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/argument_parser/argument_parser.hpp>

#include <seqan3/std/filesystem>

#include <indicators/cursor_control.hpp>
#include <indicators/dynamic_progress.hpp>
#include <indicators/progress_bar.hpp>

#include "distances/cv.hpp"
#include "distances/d2.hpp"
#include "distances/d2star.hpp"
#include "distances/dvstar.hpp"
#include "distances/parallelize.hpp"
#include "probabilistic_suffix_tree_map.hpp"

using tree_t = pst::ProbabilisticSuffixTreeMap<seqan3::dna5>;
using matrix_t = Eigen::MatrixXd;

struct input_arguments {
  std::string distance_name{"dvstar"};
  size_t order{6};
  size_t background_order{0};
  std::filesystem::path filepath{""};
  std::filesystem::path filepath_to{""};
  std::filesystem::path scores{""};
  double pseudo_count_amount{1.0};
};

input_arguments parse_cli_arguments(int argc, char *argv[]) {
  input_arguments arguments{};

  seqan3::argument_parser parser{"pst-distance-calculation", argc, argv,
                                 seqan3::update_notifications::off};
  parser.info.short_description =
      "Calculate distances between PSTs saved in a hdf5 file.";

  parser.add_option(arguments.filepath, 'p', "path",
                    "Path to hdf5 file where PSTs are stored.");

  parser.add_option(arguments.filepath_to, 'y', "path-to",
                    "Path to hdf5 file where PSTs are stored.");

  parser.add_option(arguments.scores, 's', "scores-path",
                    "Path to hdf5 file where scores will be stored.");

  parser.add_option(arguments.pseudo_count_amount, 'm', "pseudo-count-amount",
                    "Size of pseudo count for probability estimation. See e.g. "
                    "https://en.wikipedia.org/wiki/Additive_smoothing .");

  parser.add_option(arguments.distance_name, 'n', "distance-name",
                    "Name of distance function.  Must be one of 'd2', "
                    "'d2star', 'dvstar', 'cv' and 'cv-estimation'");
  parser.add_option(arguments.order, 'o', "order",
                    "Length of contexts in some distance calculations.");
  parser.add_option(arguments.background_order, 'b', "background-order",
                    "Length of background in some distance calculations.");

  try {
    parser.parse();
  } catch (seqan3::argument_parser_error const &ext) {
    std::cout << "[PARSER ERROR] " << ext.what() << '\n';
    return arguments;
  }

  return arguments;
}

void calculate_vector_slice(
    size_t start_index, size_t stop_index, matrix_t &distances,
    std::vector<Eigen::VectorXd> &vectors,
    const std::function<float(tree_t &, tree_t &)> &fun) {
  for (size_t i = start_index; i < stop_index; i++) {
    for (size_t j = 0; j < vectors.size(); j++) {
      distances(i, j) =
          pst::distances::details::cosine_dissimilarity(vectors[i], vectors[j]);
    }
  }
}

void calculate_slice(size_t start_index, size_t stop_index, matrix_t &distances,
                     std::vector<tree_t> &trees, std::vector<tree_t> &trees_to,
                     const std::function<float(tree_t &, tree_t &)> &fun) {

  for (size_t i = start_index; i < stop_index; i++) {
    for (size_t j = 0; j < trees_to.size(); j++) {
      distances(i, j) = fun(trees[i], trees_to[j]);
    }
  }
}

void calculate_slice_with_progress(
    size_t start_index, size_t stop_index, matrix_t &distances,
    std::vector<tree_t> &trees, std::vector<tree_t> &trees_to,
    const std::function<float(tree_t &, tree_t &)> &fun,
    indicators::DynamicProgress<indicators::ProgressBar> &bars, int bars_i) {

  for (size_t i = start_index; i < stop_index; i++) {
    for (size_t j = 0; j < trees_to.size(); j++) {
      distances(i, j) = fun(trees[i], trees_to[j]);
    }
    float progress = float(i - start_index) / float(stop_index - start_index);
    bars[bars_i].set_progress(progress * 100);
  }

  bars[bars_i].mark_as_completed();
}

void calculate_slice(size_t start_index, size_t stop_index, matrix_t &distances,
                     std::vector<tree_t> &trees,
                     const std::function<float(tree_t &, tree_t &)> &fun) {
  calculate_slice(start_index, stop_index, distances, trees, trees, fun);
}

std::tuple<std::function<float(tree_t &, tree_t &)>, std::string>
parse_distance_function(input_arguments &arguments) {
  if (arguments.distance_name == "cv") {
    auto fun = [&](auto &left, auto &right) {
      return pst::distances::cv<seqan3::dna5>(left, right,
                                              arguments.background_order);
    };
    return {fun, "cv-" + std::to_string(arguments.background_order)};
  } else if (arguments.distance_name == "cv-estimation") {
    auto fun = [&](auto &left, auto &right) {
      return pst::distances::cv_estimation<seqan3::dna5>(
          left, right, arguments.order, arguments.background_order);
    };
    return {fun, "cv-estimation-" + std::to_string(arguments.order) + "-" +
                     std::to_string(arguments.background_order)};
  } else if (arguments.distance_name == "dvstar") {
    auto fun = [&](auto &left, auto &right) {
      return pst::distances::dvstar<seqan3::dna5>(left, right,
                                                  arguments.background_order);
    };
    return {fun, "dvstar-" + std::to_string(arguments.background_order)};
  } else if (arguments.distance_name == "d2star") {
    auto fun = [&](auto &left, auto &right) {
      return pst::distances::d2star<seqan3::dna5>(left, right,
                                                  arguments.background_order);
    };
    return {fun, "d2star-" + std::to_string(arguments.background_order)};
  } else if (arguments.distance_name == "d2") {
    auto fun = [&](auto &left, auto &right) {
      return pst::distances::d2<seqan3::dna5>(left, right);
    };
    return {fun, "d2"};
  }

  throw std::invalid_argument("Invalid distance function name.");
}

std::vector<Eigen::VectorXd> get_composition_vectors(std::vector<tree_t> &trees,
                                                     size_t order,
                                                     size_t background_order) {
  std::cout << "getting contexts..." << std::endl;
  auto contexts = pst::distances::details::get_all_contexts<seqan3::dna5>(
      order, trees[0].valid_characters);

  std::vector<Eigen::VectorXd> vectors{};
  std::cout << "calculating cvs..." << std::endl;
  for (size_t i = 0; i < trees.size(); i++) {
    vectors.push_back(pst::distances::details::composition_vector<seqan3::dna5>(
        trees[i], contexts, background_order));
  }

  return vectors;
}

matrix_t calculate_distances(
    std::vector<tree_t> &trees, std::vector<tree_t> &trees_to,
    input_arguments &arguments,
    const std::function<float(tree_t &, tree_t &)> &distance_fun) {
  std::cout << "calculating distances... " << trees.size() << " x "
            << trees_to.size() << std::endl;

  matrix_t distances(trees.size(), trees_to.size());

  if (arguments.distance_name == "cv-estimation") {
    std::vector<Eigen::VectorXd> vectors{};
    vectors = get_composition_vectors(trees, arguments.order,
                                      arguments.background_order);
    auto fun = [&](size_t start_index, size_t stop_index) {
      calculate_vector_slice(start_index, stop_index, std::ref(distances),
                             std::ref(vectors), distance_fun);
    };
    pst::parallelize::parallelize(trees.size(), fun);

  } else {
    auto fun = [&](size_t start_index, size_t stop_index,
                   indicators::DynamicProgress<indicators::ProgressBar> &bars,
                   int bar_i) {
      calculate_slice_with_progress(
          start_index, stop_index, std::ref(distances), std::ref(trees),
          std::ref(trees_to), distance_fun, bars, bar_i);
    };
    pst::parallelize::parallelize_with_progress(trees.size(), fun);
  }

  return distances;
}

std::vector<tree_t> get_trees(HighFive::File &file,
                              const double pseudo_count_amount) {
  const std::string DATASET_NAME("signatures");

  HighFive::DataSet dataset = file.getDataSet(DATASET_NAME);

  std::vector<std::string> result_string_list;
  dataset.read(result_string_list);

  std::vector<tree_t> trees{};
  std::cout << "parsing trees..." << std::endl;
  std::transform(result_string_list.begin(), result_string_list.end(),
                 std::back_inserter(trees),
                 [&pseudo_count_amount](std::string &tree) -> tree_t {
                   return tree_t{tree, pseudo_count_amount};
                 });

  return trees;
}

int main(int argc, char *argv[]) {
  input_arguments arguments = parse_cli_arguments(argc, argv);
  const auto [distance_fun, distance_name_with_args] =
      parse_distance_function(arguments);

  std::vector<tree_t> trees;
  std::vector<tree_t> trees_to;

  if (arguments.filepath.extension() == ".h5" ||
      arguments.filepath.extension() == ".hdf5") {
    HighFive::File file{arguments.filepath, HighFive::File::ReadOnly};
    trees = get_trees(file, arguments.pseudo_count_amount);
  } else if (arguments.filepath.extension() == ".tree" ||
             arguments.filepath.extension() == ".bintree") {
    pst::ProbabilisticSuffixTreeMap<seqan3::dna5> tree{
        arguments.filepath, arguments.pseudo_count_amount};

    trees = std::vector<tree_t>{std::move(tree)};
  }

  if (arguments.filepath_to.empty()) {
    trees_to = trees;
  } else if (arguments.filepath_to.extension() == ".h5" ||
             arguments.filepath_to.extension() == ".hdf5") {
    HighFive::File file{arguments.filepath_to, HighFive::File::ReadOnly};
    trees_to = get_trees(file, arguments.pseudo_count_amount);
  } else if (arguments.filepath_to.extension() == ".tree" ||
             arguments.filepath_to.extension() == ".bintree") {
    pst::ProbabilisticSuffixTreeMap<seqan3::dna5> tree{
        arguments.filepath_to, arguments.pseudo_count_amount};

    trees_to = std::vector<tree_t>{std::move(tree)};
  }

  matrix_t distances =
      calculate_distances(trees, trees_to, arguments, distance_fun);
  std::cout << "Done with distances" << std::endl;

  std::cout << "done with distances" << std::endl;

  if (arguments.scores.empty()) {
    for (int i = 0; i < trees.size(); i++) {
      for (int j = 0; j < trees_to.size(); j++) {
        std::cout << distances(i, j) << " ";
      }
      std::cout << std::endl;
    }
  } else if (arguments.scores.extension() == ".h5" ||
             arguments.scores.extension() == ".hdf5") {
    HighFive::File file{arguments.scores,
                        HighFive::File::ReadWrite | HighFive::File::Create};

    if (!file.exist("distances")) {
      file.createGroup("distances");
    }
    auto distance_group = file.getGroup("distances");

    if (!distance_group.exist(distance_name_with_args)) {
      std::vector<size_t> dims{trees.size(), trees_to.size()};
      distance_group.createDataSet<float>(distance_name_with_args,
                                          HighFive::DataSpace(dims));
    }

    auto distance_data_set = distance_group.getDataSet(distance_name_with_args);
    distance_data_set.write(distances);
  }

  return EXIT_SUCCESS;
}
