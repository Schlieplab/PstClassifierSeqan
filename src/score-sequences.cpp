#include <algorithm>
#include <functional>
#include <iostream>
#include <string>
#include <thread>

#include <indicators/cursor_control.hpp>
#include <indicators/dynamic_progress.hpp>
#include <indicators/progress_bar.hpp>

#include <highfive/H5File.hpp>

#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/argument_parser/argument_parser.hpp>
#include <seqan3/io/sequence_file/input.hpp>

#include <filesystem>

#include "pst/distances/composition_vectors.hpp"
#include "pst/distances/negative_log_likelihood.hpp"
#include "pst/distances/parallelize.hpp"
#include "pst/distances/score.hpp"
#include "pst/probabilistic_suffix_tree_map.hpp"

#include "io_utils.hpp"

using tree_t = pst::ProbabilisticSuffixTreeMap<seqan3::dna5>;

struct input_arguments {
  size_t background_order{0};
  bool score_both_sequence_directions{false};
  bool background_adjusted{false};
  std::filesystem::path filepath{};
  std::filesystem::path outpath{};
  std::filesystem::path sequence_list{};
  double pseudo_count_amount{1.0};
};

input_arguments parse_cli_arguments(int argc, char *argv[]) {
  input_arguments arguments{};

  seqan3::argument_parser parser{"score-sequences", argc, argv,
                                 seqan3::update_notifications::off};
  parser.info.short_description = "Calculates the negative log-likelihood of a "
                                  "set of sequences for a set of signatures.";

  parser.add_option(arguments.filepath, 'p', "path",
                    "Path to hdf5 or .tree file where PSTs are stored.");
  parser.add_option(arguments.outpath, 'o', "out-path",
                    "Path to hdf5 file where scores will be stored.");
  parser.add_option(arguments.sequence_list, 's', "sequence-list",
                    "Path to a file with paths to sequences, or a fasta file "
                    "which will be scored.");
  parser.add_option(arguments.background_order, 'b', "background-order",
                    "Length of background for log-likelihood.");
  parser.add_flag(arguments.score_both_sequence_directions, 'r',
                  "score-forward-and-reverse",
                  "Flag to signify that the scoring function should score both "
                  "the forward and the reverse complement of the sequence.");
  parser.add_flag(arguments.background_adjusted, 'a', "background-adjusted",
                  "Flag to signify background-adjustment.");

  parser.add_option(arguments.pseudo_count_amount, 'm', "pseudo-count-amount",
                    "Size of pseudo count for probability estimation. See e.g. "
                    "https://en.wikipedia.org/wiki/Additive_smoothing .");

  try {
    parser.parse();
  } catch (seqan3::argument_parser_error const &ext) {
    std::cout << "[PARSER ERROR] " << ext.what() << '\n';
    return arguments;
  }

  return arguments;
}

std::function<double(tree_t &, std::vector<seqan3::dna5> &)>
parse_scoring_function(size_t n_sequences, bool score_both_sequence_directions,
                       bool background_adjusted, int background_order) {
  pst::distances::details::scoring::score_signature<seqan3::dna5>
      transition_score_fun =
          pst::distances::details::scoring::log_transition_prob<seqan3::dna5>;

  if (background_adjusted) {
    transition_score_fun = pst::distances::details::scoring::
        specialise_background_log_transition_prob<seqan3::dna5>(
            background_order);
  }

  std::function<double(tree_t &, std::vector<seqan3::dna5> &)> score_fun;

  if (n_sequences == 1) {
    if (score_both_sequence_directions) {
      score_fun = [transition_score_fun](
                      tree_t &tree,
                      std::vector<seqan3::dna5> &sequence) -> double {
        return pst::distances::negative_log_likelihood_symmetric_p<
            seqan3::dna5>(tree, sequence, transition_score_fun);
      };
    } else {
      score_fun = [transition_score_fun](
                      tree_t &tree,
                      std::vector<seqan3::dna5> &sequence) -> double {
        return pst::distances::negative_log_likelihood_p<seqan3::dna5>(
            tree, sequence, transition_score_fun);
      };
    }
  } else {
    if (score_both_sequence_directions) {
      score_fun = [transition_score_fun](
                      tree_t &tree,
                      std::vector<seqan3::dna5> &sequence) -> double {
        return pst::distances::negative_log_likelihood_symmetric_<seqan3::dna5>(
            tree, sequence, transition_score_fun);
      };
    } else {
      score_fun = [transition_score_fun](
                      tree_t &tree,
                      std::vector<seqan3::dna5> &sequence) -> double {
        return pst::distances::negative_log_likelihood<seqan3::dna5>(
            tree, sequence, transition_score_fun);
      };
    }
  }

  return score_fun;
}

std::tuple<std::vector<std::vector<double>>, std::vector<std::string>>
score_sequences_paths(std::vector<tree_t> &trees,
                      std::vector<std::filesystem::path> &sequence_list,
                      size_t background_order,
                      bool score_both_sequence_directions,
                      bool background_adjusted) {
  std::vector<std::vector<double>> all_scores{};
  std::vector<std::string> ids{};

  for (auto &path : sequence_list) {
    std::vector<seqan3::dna5_vector> sequences{};

    if (path.extension() == ".fastq") {
      seqan3::sequence_file_input file_in{std::ifstream{path},
                                          seqan3::format_fastq{}};
      for (auto &[seq, id, qual] : file_in) {
        ids.push_back(std::move(id));
        sequences.push_back(std::move(seq));
      }
    } else {
      seqan3::sequence_file_input file_in{std::ifstream{path},
                                          seqan3::format_fasta{}};
      for (auto &[seq, id, qual] : file_in) {
        ids.push_back(std::move(id));
        sequences.push_back(std::move(seq));
      }
    }

    std::vector<std::vector<double>> scores(sequences.size(),
                                            std::vector<double>(trees.size()));

    std::function<double(tree_t &, std::vector<seqan3::dna5> &)> score_fun =
        parse_scoring_function(sequences.size(), score_both_sequence_directions,
                               background_adjusted, background_order);

    indicators::DynamicProgress<indicators::ProgressBar> bars{};
    bars.set_option(indicators::option::HideBarWhenComplete{true});

    auto fun = [&](size_t start_index, size_t stop_index) {
      pst::score_sequences_slice_with_progress(
          start_index, stop_index, std::ref(scores), std::ref(trees),
          std::ref(sequences), score_fun, std::ref(bars));
    };
    pst::parallelize::parallelize(sequences.size(), fun);

    std::move(scores.begin(), scores.end(), std::back_inserter(all_scores));
  }

  return {all_scores, ids};
}

int main(int argc, char *argv[]) {
  input_arguments arguments = parse_cli_arguments(argc, argv);

  std::vector<tree_t> trees =
      get_trees(arguments.filepath, arguments.pseudo_count_amount);

  if (!std::filesystem::exists(arguments.filepath)) {
    std::cerr << "Error: " << arguments.filepath << " is not a file.";
    return EXIT_FAILURE;
  }
  if (!std::filesystem::exists(arguments.sequence_list)) {
    std::cerr << "Error: " << arguments.sequence_list << " is not a file.";
    return EXIT_FAILURE;
  }

  std::vector<std::filesystem::path> sequence_list{};

  if (arguments.sequence_list.extension() == ".fasta" ||
      arguments.sequence_list.extension() == ".fastq" ||
      arguments.sequence_list.extension() == ".fna" ||
      arguments.sequence_list.extension() == ".fa" ||
      arguments.sequence_list.extension() == ".gz") {
    sequence_list.push_back(arguments.sequence_list);
  } else if (arguments.sequence_list.extension() == ".txt") {
    std::ifstream infile(arguments.sequence_list);
    std::string path;
    while (std::getline(infile, path)) {
      sequence_list.emplace_back(path);
    }
  } else {
    std::cerr << "Error: " << arguments.sequence_list
              << " has invalid extension, and can't be parsed.";
    return EXIT_FAILURE;
  }

  auto [scores, ids] = score_sequences_paths(
      trees, sequence_list, arguments.background_order,
      arguments.score_both_sequence_directions, arguments.background_adjusted);

  if (arguments.outpath.empty()) {
    for (int i = 0; i < trees.size(); i++) {
      for (int j = 0; j < sequence_list.size(); j++) {
        std::cout << scores[i][j] << " ";
      }
      std::cout << std::endl;
    }
  } else {
    HighFive::File out_file{arguments.outpath,
                            HighFive::File::ReadWrite | HighFive::File::Create};

    if (!out_file.exist("scores")) {
      std::vector<size_t> dims{scores.size(),
                               static_cast<size_t>(scores[0].size())};
      out_file.createDataSet<double>("scores", HighFive::DataSpace(dims));
    }
    auto scores_dataset = out_file.getDataSet("scores");
    scores_dataset.write(scores);

    if (!out_file.exist("ids")) {
      out_file.createDataSet<std::string>("ids",
                                          HighFive::DataSpace::From(ids));
    }
    auto ids_dataset = out_file.getDataSet("ids");
    ids_dataset.write(ids);
  }

  return EXIT_SUCCESS;
}
