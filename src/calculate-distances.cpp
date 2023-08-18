#include <algorithm>
#include <cstdlib>
#include <functional>
#include <iostream>
#include <optional>
#include <string>
#include <thread>

#include <Eigen/Dense>

#include <highfive/H5Easy.hpp>

#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/argument_parser/argument_parser.hpp>

#include <seqan3/std/filesystem>

#include <indicators/cursor_control.hpp>
#include <indicators/dynamic_progress.hpp>
#include <indicators/progress_bar.hpp>

#include "pst/distances/cv.hpp"
#include "pst/distances/d2.hpp"
#include "pst/distances/d2star.hpp"
#include "pst/distances/dvstar.hpp"
#include "pst/distances/kl_divergence.hpp"
#include "pst/distances/other_distances.hpp"
#include "pst/distances/parallelize.hpp"
#include "pst/probabilistic_suffix_tree_map.hpp"

#include "io_utils.hpp"

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

  parser.add_option(
      arguments.distance_name, 'n', "distance-name",
      "Name of distance function.  Must be one of 'd2', "
      "'d2star', 'dvstar', 'nearest-dvstar', 'penalized-dvstar', 'kl', "
      "'kl-both', 'nll', 'nll-background', 'cv' and 'cv-estimation'");
  parser.add_option(arguments.order, 'o', "order",
                    "Length of contexts/sequence in distances 'cv-estimation', "
                    "'nll', and 'kl'.");
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
    indicators::DynamicProgress<indicators::BlockProgressBar> &bars) {
  std::string text = "Distances from " + std::to_string(start_index) + " to " +
                     std::to_string(stop_index);

  indicators::BlockProgressBar bar{indicators::option::BarWidth{50},
                                   indicators::option::PrefixText{text},
                                   indicators::option::ShowElapsedTime{true},
                                   indicators::option::ShowRemainingTime{true}};
  //  auto bars_i = bars.push_back(bar);

  double progress_per_tick =
      1 / (float(stop_index - start_index) * trees_to.size());
  double progress = 0;
  for (size_t i = start_index; i < stop_index; i++) {
    for (size_t j = 0; j < trees_to.size(); j++) {
      distances(i, j) = fun(trees[i], trees_to[j]);
      progress += progress_per_tick;

      if ((int(progress * 1000) % 10) == 0) {
        //          bars[bars_i].set_progress(progress * 100);
      }
    }
  }
  //  bars[bars_i].mark_as_completed();
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
  } else if (arguments.distance_name == "penalized-dvstar") {
    auto fun = [&](auto &left, auto &right) {
      return pst::distances::penalized_dvstar<seqan3::dna5>(
          left, right, arguments.background_order);
    };
    return {fun,
            "penalized-dvstar-" + std::to_string(arguments.background_order)};
  } else if (arguments.distance_name == "nearest-dvstar") {
    auto fun = [&](auto &left, auto &right) {
      return pst::distances::nearest_dvstar<seqan3::dna5>(
          left, right, arguments.background_order);
    };
    return {fun,
            "nearest-dvstar-" + std::to_string(arguments.background_order)};
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
  } else if (arguments.distance_name == "kl") {
    auto fun = [&](auto &left, auto &right) {
      return pst::distances::symmetric_kl_divergence<seqan3::dna5>(
          left, right, arguments.order);
    };
    return {fun, "kl-" + std::to_string(arguments.order)};
  } else if (arguments.distance_name == "kl-both") {
    auto fun = [&](auto &left, auto &right) {
      return pst::distances::symmetric_kl_divergence_both<seqan3::dna5>(left,
                                                                        right);
    };
    return {fun, "kl-both"};
  } else if (arguments.distance_name == "kl-both-background") {
    auto fun = [&](auto &left, auto &right) {
      return pst::distances::symmetric_kl_divergence_both_background<
          seqan3::dna5>(left, right);
    };
    return {fun, "kl-both-background"};
  } else if (arguments.distance_name == "jensenshannon") {
    auto fun = [&](auto &left, auto &right) {
      return pst::distances::jensenshannon<seqan3::dna5>(left, right);
    };
    return {fun, arguments.distance_name};
  } else if (arguments.distance_name == "canberra") {
    auto fun = [&](auto &left, auto &right) {
      return pst::distances::canberra<seqan3::dna5>(left, right);
    };
    return {fun, arguments.distance_name};
  } else if (arguments.distance_name == "correlation") {
    auto fun = [&](auto &left, auto &right) {
      return pst::distances::correlation<seqan3::dna5>(left, right);
    };
    return {fun, arguments.distance_name};
  } else if (arguments.distance_name == "jaccard") {
    auto fun = [&](auto &left, auto &right) {
      return pst::distances::jaccard_estimation<seqan3::dna5>(left, right);
    };
    return {fun, arguments.distance_name};
  } else if (arguments.distance_name == "jaccard-prob") {
    auto fun = [&](auto &left, auto &right) {
      return pst::distances::jaccard_estimation_prob<seqan3::dna5>(left, right);
    };
    return {fun, arguments.distance_name};
  } else if (arguments.distance_name == "nll") {
    auto fun = [&](auto &left, auto &right) {
      return pst::distances::negative_log_likelihood<seqan3::dna5>(
          left, right, arguments.order);
    };
    return {fun, "nll-" + std::to_string(arguments.order)};
  } else if (arguments.distance_name == "nll-background") {
    auto fun = [&](auto &left, auto &right) {
      return pst::distances::negative_log_likelihood_background<seqan3::dna5>(
          left, right, arguments.order, arguments.background_order);
    };
    return {fun, "nll-" + std::to_string(arguments.order) + "-background-" +
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
    auto fun =
        [&](size_t start_index, size_t stop_index,
            indicators::DynamicProgress<indicators::BlockProgressBar> &bars) {
          calculate_slice_with_progress(start_index, stop_index,
                                        std::ref(distances), std::ref(trees),
                                        std::ref(trees_to), distance_fun, bars);
        };
    pst::parallelize::parallelize_with_progress(trees.size(), fun);
  }

  return distances;
}

int main(int argc, char *argv[]) {
  input_arguments arguments = parse_cli_arguments(argc, argv);
  const auto [distance_fun, distance_name_with_args] =
      parse_distance_function(arguments);

  std::vector<tree_t> trees =
      get_trees(arguments.filepath, arguments.pseudo_count_amount);
  std::vector<tree_t> trees_to;

  if (arguments.filepath_to.empty()) {
    trees_to = trees;
  } else {
    trees_to = get_trees(arguments.filepath_to, arguments.pseudo_count_amount);
  }

  matrix_t distances =
      calculate_distances(trees, trees_to, arguments, distance_fun);

  std::vector<std::string> ids{};
  std::transform(trees.begin(), trees.end(), std::back_inserter(ids),
                 [](tree_t &tree) -> std::string { return tree.id; });

  std::vector<std::string> ids_to{};
  std::transform(trees_to.begin(), trees_to.end(), std::back_inserter(ids_to),
                 [](tree_t &tree) -> std::string { return tree.id; });

  if (arguments.scores.empty()) {
    for (int i = 0; i < trees.size(); i++) {
      for (int j = 0; j < trees_to.size(); j++) {
        std::cout << distances(i, j) << " ";
      }
      std::cout << std::endl;
    }
  } else if (arguments.scores.extension() == ".h5" ||
             arguments.scores.extension() == ".hdf5") {
    std::cout << "Writing to file..." << std::endl;
    H5Easy::File file{arguments.scores, HighFive::File::OpenOrCreate};

    H5Easy::dump(file, "/distances/" + distance_name_with_args, distances,
                 H5Easy::DumpMode::Overwrite);
    H5Easy::dump(file, "/ids", ids, H5Easy::DumpMode::Overwrite);
    H5Easy::dump(file, "/ids_to", ids_to, H5Easy::DumpMode::Overwrite);

    std::cout << "Wrote distances to: " << arguments.scores.string()
              << std::endl;
  }

  return EXIT_SUCCESS;
}
