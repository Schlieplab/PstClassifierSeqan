#include <string>
#include <vector>

#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/argument_parser/all.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/io/sequence_file/input.hpp>
#include <seqan3/range/container/bitcompressed_vector.hpp>
#include <seqan3/range/view/char_to.hpp>

#include "probabilistic_suffix_tree.hpp"

struct input_arguments {
  size_t max_depth{15};
  size_t min_count{100};
  float threshold{1.2};
  size_t number_of_parameters{192};
  std::string pruning_method{"cutoff"};
  std::string estimator{"PS"};
  std::vector<seqan3::bitcompressed_vector<seqan3::dna5>> sequences{};
  std::vector<std::string> ids{};
};

struct my_traits : seqan3::sequence_file_input_default_traits_dna {
  template <typename alph>
  using sequence_container =
      seqan3::bitcompressed_vector<alph>; // must be defined as a template!
};

input_arguments parse_cli_arguments(int argc, char *argv[]) {
  std::string filename{};

  input_arguments arguments{};

  seqan3::argument_parser parser{"Count kmers in the fasta file.", argc, argv,
                                 false};

  parser.add_positional_option(filename, "path to fasta file.");

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
                    "pruned. Either 'KL' or 'PS'.");

  parser.add_option(arguments.pruning_method, 'p', "pruning-method",
                    "pruning method to use. Either 'cutoff' for pruning until "
                    "the threshold is reached, or 'parameters' to prune until "
                    "a certain number of parameters have been reached.");

  try {
    parser.parse();
  } catch (seqan3::parser_invalid_argument const &ext) {
    seqan3::debug_stream << "[PARSER ERROR] " << ext.what() << '\n';
    return arguments;
  }
  seqan3::debug_stream << "The text was: " << filename << "\n";

  seqan3::sequence_file_input<my_traits> file_in{filename};

  for (auto &[seq, id, qual] : file_in) {
    arguments.sequences.push_back(seq);
    arguments.ids.push_back(id);
  }

  return arguments;
}

std::string train(seqan3::bitcompressed_vector<seqan3::dna5> sequence,
                  std::string id, size_t max_depth, size_t min_count,
                  float threshold, size_t number_of_parameters,
                  std::string pruning_method, std::string estimator) {
  pst::ProbabilisticSuffixTree<seqan3::dna5> pst{id,
                                                 sequence,
                                                 max_depth,
                                                 min_count,
                                                 threshold,
                                                 number_of_parameters,
                                                 pruning_method,
                                                 estimator};

  return pst.to_tree();
}

int main(int argc, char *argv[]) {
  input_arguments arguments = parse_cli_arguments(argc, argv);

  seqan3::debug_stream << "Building index" << std::endl;

  std::string tree = train(arguments.sequences[0], arguments.ids[0],
                           arguments.max_depth, arguments.min_count,
                           arguments.threshold, arguments.number_of_parameters,
                           arguments.pruning_method, arguments.estimator);
  std::cout << tree << std::endl;

  return 0;
}