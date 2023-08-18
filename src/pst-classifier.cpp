#include <chrono>
#include <cstdlib>
#include <string>
#include <vector>

#include <seqan3/alphabet/all.hpp>
#include <seqan3/argument_parser/all.hpp>
#include <seqan3/argument_parser/exceptions.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/io/sequence_file/input.hpp>

#include <filesystem>

#include <vlmc_from_kmers/build_vlmc.hpp>

#include "pst/kl_tree.hpp"
#include "pst/kl_tree_map.hpp"
#include "pst/probabilistic_suffix_tree_map.hpp"

using seqan3::operator""_dna5;

struct input_arguments {
  size_t max_depth{15};
  size_t min_count{100};
  float threshold{3.9075};
  size_t number_of_parameters{192};
  std::string algorithm_method{"hashmap"};
  std::string pruning_method{"cutoff"};
  std::string estimator{"KL"};
  std::filesystem::path fasta_path{};
  bool multi_core{false};
  int parallel_depth{1};
  std::filesystem::path out_path{""};
  double pseudo_count_amount{1.0};
  std::filesystem::path tmp_path{"./tmp"};
  std::string in_or_out_of_core{"internal"};
};

struct my_traits : seqan3::sequence_file_input_default_traits_dna {
  template <typename alph>
  using sequence_container =
      std::vector<alph>; // must be defined as a template!
};

input_arguments parse_cli_arguments(int argc, char *argv[]) {
  input_arguments arguments{};

  seqan3::argument_parser parser{"Pst-Classifier", argc, argv,
                                 seqan3::update_notifications::off};
  parser.info.short_description = "Build PST/VLMC on the given fasta file.";

  parser.add_positional_option(arguments.fasta_path, "path to fasta file.");

  parser.add_option(arguments.out_path, 'o', "out-path",
                    "Name of output file.  The suffix '.tree' or '.bintree' "
                    "will be added if not present. ");

  parser.add_option(arguments.tmp_path, 't', "tmp-path",
                    "Path to temporary directory.  Required for the 'kmers' "
                    "method. Defaults to './tmp'.");

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
                    "Algorithm for construction of the VLMC."
                    "Either 'kmers' which builds the VLMC from kmers, "
                    "'hashmap', which stores the k-mers in a hashmap, "
                    "or 'tree' which will store "
                    "the k-mers in a lazy suffix tree.");

  parser.add_option(arguments.in_or_out_of_core, 'i', "in-or-out-of-core",
                    "For the 'kmers' method, specifies if the algorithm should "
                    "RAM ('internal') or disk ('external').");

  parser.add_option(arguments.pseudo_count_amount, 'j', "pseudo-count-amount",
                    "Size of pseudo count for probability estimation, does "
                    "not work for algorithm-method=='tree'. See e.g. "
                    "https://en.wikipedia.org/wiki/Additive_smoothing .");

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

  return arguments;
}

std::string train(const std::filesystem::path &fasta_path, size_t max_depth,
                  size_t min_count, float threshold,
                  size_t number_of_parameters, std::string pruning_method,
                  std::string algorithm, std::string estimator, bool multi_core,
                  int parallel_depth, const double pseudo_count_amount) {

  lst::details::sequence_t<seqan3::dna5> sequence{};
  std::string id{};

  seqan3::sequence_file_input<my_traits> file_in{fasta_path};

  for (auto &[seq, id_, qual] : file_in) {
    std::move(seq.begin(), seq.end(), std::back_inserter(sequence));
    sequence.push_back('N'_dna5);

    id += id_;
    id += "|";
  }
  id.pop_back();
  sequence.pop_back();

  if (estimator == "KL" && algorithm == "hashmap") {
    pst::KullbackLieblerTreeMap<seqan3::dna5> pst{id,
                                                  sequence,
                                                  max_depth,
                                                  min_count,
                                                  threshold,
                                                  number_of_parameters,
                                                  pruning_method,
                                                  multi_core,
                                                  parallel_depth,
                                                  pseudo_count_amount};
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

  if (arguments.algorithm_method == "kmers") {

    vlmc::Core core;
    if (arguments.in_or_out_of_core == "internal") {
      core = vlmc::Core::in;
    } else if (arguments.in_or_out_of_core == "external") {
      core = vlmc::Core::out;
    } else {
      std::cerr << "--in-or-out-of-core parameter given invalid value.  Only "
                   "'internal' or 'external' is allowed."
                << std::endl;
      return EXIT_FAILURE;
    }
    std::filesystem::path out_path = arguments.out_path;
    if (out_path == "") {
      out_path = std::filesystem::path(".tmp.bintree");
    } else if (!out_path.has_extension()) {
      out_path.replace_extension(".bintree");
    } else if (out_path.extension() != ".bintree") {
      std::cerr << "out path has invalid extension, should be '.bintree'."
                << std::endl;
      return EXIT_FAILURE;
    }

    vlmc::build_vlmc(arguments.fasta_path, arguments.max_depth,
                     arguments.min_count, arguments.threshold, out_path,
                     arguments.tmp_path, core, arguments.pseudo_count_amount);

    if (arguments.out_path == "") {
      vlmc::dump_path(out_path, std::filesystem::path{});
    }
  } else {

    std::string tree =
        train(arguments.fasta_path, arguments.max_depth, arguments.min_count,
              arguments.threshold, arguments.number_of_parameters,
              arguments.pruning_method, arguments.algorithm_method,
              arguments.estimator, arguments.multi_core,
              arguments.parallel_depth, arguments.pseudo_count_amount);

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
  }

  return EXIT_SUCCESS;
}
