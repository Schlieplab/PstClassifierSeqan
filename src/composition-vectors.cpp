#include <Eigen/Dense>
#include <algorithm>
#include <functional>
#include <highfive/H5File.hpp>
#include <iostream>
#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/argument_parser/argument_parser.hpp>
#include <string>
#include <thread>

#include <seqan3/std/filesystem>

#include "distances/composition_vectors.hpp"
#include "distances/cv.hpp"
#include "probabilistic_suffix_tree_map.hpp"

using tree_t = pst::ProbabilisticSuffixTreeMap<seqan3::dna5>;

struct input_arguments {
  std::string type{"cv"};
  size_t order{6};
  size_t background_order{0};
  std::filesystem::path filepath{""};
};

input_arguments parse_cli_arguments(int argc, char *argv[]) {
  input_arguments arguments{};

  seqan3::argument_parser parser{"vlmc-to-composition-vectors", argc, argv,
                                 false};
  parser.info.short_description = "Turns a VLMC into a composition vector.";

  parser.add_option(arguments.filepath, 'p', "path",
                    "Path to hdf5 file where PSTs are stored.");
  parser.add_option(
      arguments.type, 't', "type",
      "Type of composition vector.  Must be one of 'cv' and 'cv-estimation'");
  parser.add_option(arguments.order, 'o', "order",
                    "Length of contexts in types.");
  parser.add_option(arguments.background_order, 'b', "background-order",
                    "Length of background.");

  try {
    parser.parse();
  } catch (seqan3::argument_parser_error const &ext) {
    std::cout << "[PARSER ERROR] " << ext.what() << '\n';
    return arguments;
  }

  return arguments;
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

std::vector<std::vector<double>>
get_composition_vectors(std::vector<tree_t> &trees, size_t order,
                        size_t background_order) {
  std::cout << "getting contexts..." << std::endl;
  auto contexts = pst::distances::details::get_all_contexts<seqan3::dna5>(
      order, trees[0].valid_characters);

  std::vector<std::vector<double>> vectors{};
  std::cout << "calculating cvs..." << std::endl;

  std::transform(
      trees.begin(), trees.end(), std::back_inserter(vectors),
      [&](tree_t &tree) -> std::vector<double> {
        Eigen::VectorXd vector =
            pst::distances::composition_vector_state_probability_scaled<
                seqan3::dna5>(tree, contexts, background_order);

        std::vector<double> vec(vector.size());
        for (size_t i = 0; i < vector.size(); i++) {
          vec[i] = vector(i);
        }
        return vec;
      });

  return vectors;
}

int main(int argc, char *argv[]) {
  input_arguments arguments = parse_cli_arguments(argc, argv);

  HighFive::File file{arguments.filepath, HighFive::File::ReadWrite};

  auto trees = get_trees(file);
  auto vectors = get_composition_vectors(trees, arguments.order,
                                         arguments.background_order);

  file.unlink("composition-vectors");
  if (!file.exist("composition-vectors")) {
    std::vector<size_t> dims{vectors.size(),
                             static_cast<size_t>(vectors[0].size())};
    std::cout << dims[0] << " " << dims[1] << std::endl;
    file.createDataSet<double>("composition-vectors",
                               HighFive::DataSpace(dims));
  }
  auto vector_dataset = file.getDataSet("composition-vectors");

  std::cout << "write..." << std::endl;
  vector_dataset.write(vectors);

  return EXIT_SUCCESS;
}
