#include <algorithm>
#include <functional>
#include <iostream>
#include <string>
#include <thread>

#include <highfive/H5File.hpp>

#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/argument_parser/argument_parser.hpp>
#include <seqan3/io/sequence_file/input.hpp>

#include <seqan3/std/filesystem>

#include "distances/composition_vectors.hpp"
#include "distances/negative_log_likelihood.hpp"
#include "distances/parallelize.hpp"
#include "distances/score.hpp"
#include "probabilistic_suffix_tree_map.hpp"

using tree_t = pst::ProbabilisticSuffixTreeMap<seqan3::dna5>;

struct input_arguments {
  size_t background_order{0};
  bool score_both_sequence_directions{false};
  std::filesystem::path filepath{};
  std::filesystem::path outpath{};
  std::filesystem::path sequence_list{};
};

input_arguments parse_cli_arguments(int argc, char *argv[]) {
  input_arguments arguments{};

  seqan3::argument_parser parser{"score-sequences", argc, argv, false};
  parser.info.short_description = "Calculates the negative log-likelihood of a "
                                  "set of sequences for a set of signatures.";

  parser.add_option(arguments.filepath, 'p', "path",
                    "Path to hdf5 or .tree file where PSTs are stored.");
  parser.add_option(arguments.outpath, 'o', "out-path",
                    "Path to hdf5 file where scores will be stored.");
  parser.add_option(arguments.sequence_list, 's', "sequence-list",
                    "Path to a file with paths to sequences or a fasta file "
                    "to be scored.");
  parser.add_option(arguments.background_order, 'b', "background-order",
                    "Length of background for log-likelihood.");
  parser.add_flag(arguments.score_both_sequence_directions, 'r',
                  "score-forward-and-reverse",
                  "Flag to signify that the scoring function should score both "
                  "the forward and the reverse complement of the sequence.");

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
  std::transform(result_string_list.begin(), result_string_list.end(),
                 std::back_inserter(trees),
                 [](std::string &tree) -> tree_t { return tree_t{tree}; });

  return trees;
}

std::tuple<std::vector<std::vector<double>>, std::vector<std::string>>
score_sequences_paths(std::vector<tree_t> &trees,
                      std::vector<std::string> &sequence_list,
                      size_t background_order,
                      bool score_both_sequence_directions) {
  std::vector<std::vector<double>> all_scores{};
  std::vector<std::string> ids{};

  size_t sequence_idx = 0;

  for (auto &path : sequence_list) {
    seqan3::sequence_file_input file_in{std::ifstream{path},
                                        seqan3::format_fasta{}};
    std::vector<seqan3::dna5_vector> sequences{};
    for (auto &[seq, id, qual] : file_in) {
      ids.push_back(std::move(id));
      sequences.push_back(std::move(seq));
    }
    std::vector<std::vector<double>> scores(sequences.size(),
                                            std::vector<double>(trees.size()));

    auto fun = [&](size_t start_index, size_t stop_index) {
      std::function<double(tree_t &, std::vector<seqan3::dna5> &)> score_fun;

      if (sequences.size() == 1) {
        if (score_both_sequence_directions) {
          score_fun =
              pst::distances::negative_log_likelihood_symmetric_p<seqan3::dna5>;
        } else {
          score_fun = pst::distances::negative_log_likelihood_p<seqan3::dna5>;
        }
      } else {
        if (score_both_sequence_directions) {
          score_fun =
              pst::distances::negative_log_likelihood_symmetric<seqan3::dna5>;
        } else {
          score_fun = [&](tree_t &tree,
                          std::vector<seqan3::dna5> &sequence) -> double {
            return pst::distances::negative_log_likelihood<seqan3::dna5>(
                tree, sequence);
          };
        }
      }
      pst::score_sequences_slice(start_index, stop_index, std::ref(scores),
                                 std::ref(trees), std::ref(sequences),
                                 score_fun);
    };

    pst::parallelize::parallelize(sequences.size(), fun);

    std::move(scores.begin(), scores.end(), std::back_inserter(all_scores));
  }

  return {all_scores, ids};
}

int main(int argc, char *argv[]) {
  input_arguments arguments = parse_cli_arguments(argc, argv);

  std::vector<tree_t> trees;

  if (arguments.filepath.extension() == ".h5" ||
      arguments.filepath.extension() == ".hdf5") {

    HighFive::File file{arguments.filepath, HighFive::File::ReadOnly};
    trees = get_trees(file);
  } else if (arguments.filepath.extension() == ".tree") {
    pst::ProbabilisticSuffixTreeMap<seqan3::dna5> tree{arguments.filepath};

    trees = std::vector<tree_t>{std::move(tree)};
  }

  std::vector<std::string> sequence_list{};

  if (arguments.sequence_list.extension() == ".fasta" ||
      arguments.sequence_list.extension() == ".fna" ||
      arguments.sequence_list.extension() == ".fa") {
    sequence_list.push_back(arguments.sequence_list);
  } else if (arguments.sequence_list.extension() == ".txt") {
    std::ifstream infile(arguments.sequence_list);
    std::string path;
    while (std::getline(infile, path)) {
      sequence_list.push_back(path);
    }
  }

  auto [scores, ids] =
      score_sequences_paths(trees, sequence_list, arguments.background_order,
                            arguments.score_both_sequence_directions);

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
    out_file.createDataSet<std::string>("ids", HighFive::DataSpace::From(ids));
  }
  auto ids_dataset = out_file.getDataSet("ids");
  ids_dataset.write(ids);

  return EXIT_SUCCESS;
}
