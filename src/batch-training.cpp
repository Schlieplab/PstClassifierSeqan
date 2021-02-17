#include "distances/cv.hpp"
#include "kl_tree_map.hpp"
#include "probabilistic_suffix_tree_map.hpp"
#include <Eigen/Dense>
#include <algorithm>
#include <filesystem>
#include <functional>
#include <highfive/H5File.hpp>
#include <iostream>
#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/argument_parser/argument_parser.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/io/sequence_file/input.hpp>
#include <string>
#include <thread>

using tree_t = pst::ProbabilisticSuffixTreeMap<seqan3::dna5>;
using matrix_t = Eigen::MatrixXd;

struct input_arguments {
  std::filesystem::path fasta_path{};
  size_t max_depth{15};
  size_t min_count{100};
  float threshold{1.2};
  size_t number_of_parameters{192};
  std::string pruning_method{"cutoff"};
  std::string estimator{"KL"};
  bool multi_core{true};
  int parallel_depth{1};

  std::filesystem::path h5_path{""};
};

struct my_traits : seqan3::sequence_file_input_default_traits_dna {
  template <typename alph>
  using sequence_container =
      std::vector<alph>; // must be defined as a template!
};

input_arguments parse_cli_arguments(int argc, char *argv[]) {
  input_arguments arguments{};

  seqan3::argument_parser parser{"Pst-classifier-batch", argc, argv, false};
  parser.info.short_description = "Train on all sequences in a fasta file.";

  parser.add_option(arguments.h5_path, 'o', "path",
                    "Output path to hdf5 file where PSTs are to be stored.");
  parser.add_positional_option(arguments.fasta_path, "path to fasta file.");

  parser.add_option(arguments.max_depth, 'd', "max-depth",
                    "Max depth of the built probabilistic suffix tree.  "
                    "Corresponds to the max length of each branch/context.");

  parser.add_option(
      arguments.min_count, 'c', "min-count",
      "Minimum number of time each node/context has to appear in the "
      "string to be included in the probabilistic suffix tree.");

  parser.add_option(arguments.threshold, 'k', "threshold",
                    "Threshold for the pruning stage of the algorithm.  "
                    "Smaller value gives larger tree.");

  parser.add_option(arguments.number_of_parameters, 'n', "number-of-parameters",
                    "For the 'parameters' puning-method, the number of "
                    "parameters to prune the tree until.");

  parser.add_option(arguments.estimator, 'e', "estimator",
                    "estimator used to determine which states should be "
                    "pruned. Only current option is 'KL'.");

  parser.add_option(arguments.pruning_method, 'p', "pruning-method",
                    "pruning method to use. Either 'cutoff' for pruning until "
                    "the threshold is reached, or 'parameters' to prune until "
                    "a certain number of parameters have been reached.");

  // Multiprocessing
  parser.add_flag(arguments.multi_core, 'm', "multi-core",
                  "Enable Multi-core utilisation.");

  parser.add_option(
      arguments.parallel_depth, 's', "parallel-depth",
      "If multi-core, will spawn a new thread per node up until the "
      "parallel-depth."
      "Higher value increase the number of threads spawned. Default 1");

  try {
    parser.parse();
  } catch (seqan3::argument_parser_error const &ext) {
    std::cout << "[PARSER ERROR] " << ext.what() << '\n';
    return arguments;
  }

  return arguments;
}

std::string train(lst::details::sequence_t<seqan3::dna5> sequence,
                  std::string id, size_t max_depth, size_t min_count,
                  float threshold, size_t number_of_parameters,
                  std::string pruning_method, std::string estimator,
                  bool multi_core, int parallel_depth) {

  if (estimator == "KL") {
    pst::KullbackLieblerTreeMap<seqan3::dna5> pst{id,
                                                  sequence,
                                                  max_depth,
                                                  min_count,
                                                  threshold,
                                                  number_of_parameters,
                                                  pruning_method,
                                                  multi_core,
                                                  parallel_depth};
    pst.construct_tree();
    return pst.to_tree();
  } else {
    return "";
  }
}

int main(int argc, char *argv[]) {
  input_arguments arguments = parse_cli_arguments(argc, argv);

  std::vector<std::string> trees{};
  seqan3::sequence_file_input<my_traits> file_in{arguments.fasta_path};

  std::cout << "Training" << std::endl;
  for (auto &[seq, id, qual] : file_in) {
    std::cout << id << std::endl;
    std::string tree = train(
        lst::details::sequence_t<seqan3::dna5>{std::move(seq)}, id,
        arguments.max_depth, arguments.min_count, arguments.threshold,
        arguments.number_of_parameters, arguments.pruning_method,
        arguments.estimator, arguments.multi_core, arguments.parallel_depth);
    trees.push_back(tree);
  }

  HighFive::File file{arguments.h5_path,
                      HighFive::File::OpenOrCreate | HighFive::File::ReadWrite};

  if (!file.exist("signatures")) {
    file.createDataSet<std::string>("signatures",
                                    HighFive::DataSpace::From(trees));
  }

  auto signatures_dataset = file.getDataSet("signatures");
  signatures_dataset.write(trees);

  std::cout << "done training" << std::endl;

  return EXIT_SUCCESS;
}
