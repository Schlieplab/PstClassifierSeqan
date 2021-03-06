#include <chrono>
#include <cstdlib>
#include <string>
#include <vector>

#include <seqan3/alphabet/all.hpp>
#include <seqan3/argument_parser/argument_parser.hpp>
#include <seqan3/argument_parser/exceptions.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/io/sequence_file/input.hpp>

#include <seqan3/std/filesystem>

#include "kl_tree.hpp"
#include "kl_tree_map.hpp"
#include "probabilistic_suffix_tree_map.hpp"

using seqan3::operator""_dna5;

struct input_arguments {
  size_t max_depth{15};
  size_t min_count{100};
  float threshold{3.9075};
  size_t number_of_parameters{192};
  std::string algorithm_method{"hashmap"};
  std::string pruning_method{"cutoff"};
  std::string estimator{"KL"};
  lst::details::sequence_t<seqan3::dna5> sequence{};
  std::string id{};
  bool multi_core{false};
  int parallel_depth{1};
  std::filesystem::path out_path{""};
};

struct my_traits : seqan3::sequence_file_input_default_traits_dna {
  template <typename alph>
  using sequence_container =
      std::vector<alph>; // must be defined as a template!
};

input_arguments parse_cli_arguments(int argc, char *argv[]) {
  std::filesystem::path filename{};

  input_arguments arguments{};

  seqan3::argument_parser parser{"Pst-Classifier", argc, argv, false};
  parser.info.short_description = "Build PST/VLMC on the given fasta file.";

  parser.add_positional_option(filename, "path to fasta file.");

  parser.add_option(
      arguments.out_path, 'o', "out-path",
      "Name of output file.  The suffix '.tree' will be added if not present. "
      "Will write to stdout if not specified.");

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
                    "parameters to prune the tree until.  If there are "
                    "multiple nodes with the same estimator value, "
                    "which ones are returned is not defined.");

  parser.add_option(arguments.estimator, 'e', "estimator",
                    "estimator used to determine which states should be "
                    "pruned. Only current option is 'KL'.");

  parser.add_option(arguments.pruning_method, 'p', "pruning-method",
                    "pruning method to use. Either 'cutoff' for pruning until "
                    "the threshold is reached, or 'parameters' to prune until "
                    "a certain number of parameters have been reached.");

  parser.add_option(arguments.algorithm_method, 'a', "algorithm-method",
                    "Algorithm to use. Either 'hashmap', which stores the "
                    "k-mers in a hashmap, or 'tree' which will store "
                    "the k-mers in a tree and requires suffix links.");

  // Multiprocessing
  parser.add_flag(arguments.multi_core, 'm', "multi-core",
                  "Enable Multi-core utilisation.");

  parser.add_option(
      arguments.parallel_depth, 's', "parallel-depth",
      "If multi-core, will spawn a new thread per node up until the "
      "parallel-depth. "
      "Higher value increase the number of threads spawned. Default 1");

  try {
    parser.parse();
  } catch (seqan3::argument_parser_error const &ext) {
    std::cout << "[PARSER ERROR] " << ext.what() << '\n';
    return arguments;
  }

  seqan3::sequence_file_input<my_traits> file_in{filename};

  for (auto &[seq, id, qual] : file_in) {
    std::move(seq.begin(), seq.end(), std::back_inserter(arguments.sequence));
    arguments.sequence.push_back('N'_dna5);

    arguments.id += id;
    arguments.id += "|";
  }
  arguments.id.pop_back();
  arguments.sequence.pop_back();

  return arguments;
}

std::string train(lst::details::sequence_t<seqan3::dna5> sequence,
                  std::string id, size_t max_depth, size_t min_count,
                  float threshold, size_t number_of_parameters,
                  std::string pruning_method, std::string algorithm,
                  std::string estimator, bool multi_core, int parallel_depth) {

  if (estimator == "KL" && algorithm == "hashmap") {
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
  } else if (estimator == "KL") {
    pst::KullbackLieblerTree<seqan3::dna5> pst{id,
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

  std::string tree = train(
      arguments.sequence, arguments.id, arguments.max_depth,
      arguments.min_count, arguments.threshold, arguments.number_of_parameters,
      arguments.pruning_method, arguments.algorithm_method, arguments.estimator,
      arguments.multi_core, arguments.parallel_depth);

  if (arguments.out_path == "") {
    std::cout << tree << std::endl;
  } else {
    if (!arguments.out_path.has_extension()) {
      arguments.out_path.replace_extension(".tree");
    }

    std::ofstream ofs{arguments.out_path, std::ios::out};
    ofs << tree << std::endl;
    ofs.close();
  }

  return EXIT_SUCCESS;
}
