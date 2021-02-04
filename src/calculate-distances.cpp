#include <algorithm>
#include <cstdlib>
#include <filesystem>
#include <functional>
#include <iostream>
#include <string>
#include <thread>

#include <Eigen/Dense>

#include <highfive/H5File.hpp>

#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/argument_parser/argument_parser.hpp>

#include "distances/cv.hpp"
#include "distances/parallelize.hpp"
#include "probabilistic_suffix_tree_map.hpp"

using tree_t = pst::ProbabilisticSuffixTreeMap<seqan3::dna5>;
using matrix_t = Eigen::MatrixXd;

struct input_arguments {
  std::string distance_name{"cv"};
  size_t order{6};
  size_t background_order{2};
  std::filesystem::path filepath{""};
};

input_arguments parse_cli_arguments(int argc, char *argv[]) {
  input_arguments arguments{};

  seqan3::argument_parser parser{"pst-distance-calculation", argc, argv, false};
  parser.info.short_description =
      "Calculate distances between PSTs saved in a hdf5 file.";

  parser.add_option(arguments.filepath, 'p', "path",
                    "Path to hdf5 file where PSTs are stored.");
  parser.add_option(
      arguments.distance_name, 'n', "distance-name",
      "Name of distance function.  Must be one of 'cv' and 'cv-estimation'");
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
                     std::vector<tree_t> &trees,
                     const std::function<float(tree_t &, tree_t &)> &fun) {
  for (size_t i = start_index; i < stop_index; i++) {
    for (size_t j = 0; j < trees.size(); j++) {
      distances(i, j) = fun(trees[i], trees[j]);
    }
  }
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
    vectors.push_back(pst::distances::composition_vector<seqan3::dna5>(
        trees[i], contexts, background_order));
  }

  return vectors;
}

matrix_t calculate_distances(
    std::vector<tree_t> &trees, input_arguments &arguments,
    const std::function<float(tree_t &, tree_t &)> &distance_fun) {
  std::cout << "calculating distances... " << trees.size() << " x "
            << trees.size() << std::endl;

  matrix_t distances(trees.size(), trees.size());

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
    auto fun = [&](size_t start_index, size_t stop_index) {
      calculate_slice(start_index, stop_index, std::ref(distances),
                      std::ref(trees), distance_fun);
    };
    pst::parallelize::parallelize(trees.size(), fun);
  }

  return distances;
}

std::vector<tree_t> get_trees(HighFive::File &file) {
  const std::string DATASET_NAME("signatures");

  HighFive::DataSet dataset = file.getDataSet(DATASET_NAME);

  std::vector<std::string> result_string_list;
  dataset.read(result_string_list);

  std::vector<tree_t> trees{};
  std::cout << "parsing trees..." << std::endl;
  std::transform(result_string_list.begin(), result_string_list.end(),
                 std::back_inserter(trees),
                 [](std::string &tree) -> tree_t { return tree_t{tree}; });

  return trees;
}

int main(int argc, char *argv[]) {
  input_arguments arguments = parse_cli_arguments(argc, argv);
  const auto [distance_fun, distance_name_with_args] =
      parse_distance_function(arguments);

  HighFive::File file{arguments.filepath, HighFive::File::ReadWrite};

  auto trees = get_trees(file);
  matrix_t distances = calculate_distances(trees, arguments, distance_fun);

  std::cout << "done with distances" << std::endl;
  if (!file.exist("distances")) {
    file.createGroup("distances");
  }
  auto distance_group = file.getGroup("distances");

  if (!distance_group.exist(distance_name_with_args)) {
    std::vector<size_t> dims{trees.size(), trees.size()};
    distance_group.createDataSet<float>(distance_name_with_args,
                                        HighFive::DataSpace(dims));
  }

  auto distance_data_set = distance_group.getDataSet(distance_name_with_args);
  distance_data_set.write(distances);

  return EXIT_SUCCESS;
}
