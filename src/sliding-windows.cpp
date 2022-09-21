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

#include <seqan3/std/filesystem>

#include "pst/distances/composition_vectors.hpp"
#include "pst/distances/negative_log_likelihood.hpp"
#include "pst/distances/parallelize.hpp"
#include "pst/distances/score.hpp"
#include "pst/probabilistic_suffix_tree_map.hpp"

#include "io_utils.hpp"
#include "pst/distances/sliding-windows.hpp"

using tree_t = pst::ProbabilisticSuffixTreeMap<seqan3::dna5>;

struct input_arguments {
  size_t background_order{0};
  bool score_both_sequence_directions{false};
  bool background_adjusted{false};
  std::filesystem::path filepath{};
  std::filesystem::path outpath{};
  std::filesystem::path sequence_list{};
  double pseudo_count_amount{1.0};
  int window_size{300};
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

  parser.add_option(arguments.window_size, 'i', "window-size",
                    "Size of windows.");
  try {
    parser.parse();
  } catch (seqan3::argument_parser_error const &ext) {
    std::cout << "[PARSER ERROR] " << ext.what() << '\n';
    return arguments;
  }

  return arguments;
}

std::function<std::vector<double>(tree_t &, std::vector<seqan3::dna5> &)>
parse_sliding_window_function(size_t n_sequences, bool background_adjusted,
                              int background_order, int window_size) {
  pst::distances::details::scoring::score_signature<seqan3::dna5>
      transition_score_fun =
          pst::distances::details::scoring::log_transition_prob<seqan3::dna5>;

  std::function<std::vector<double>(tree_t &, std::vector<seqan3::dna5> &)>
      sliding_window_fun;

  if (background_adjusted) {
    sliding_window_fun =
        [window_size, background_order](
            tree_t &tree,
            std::vector<seqan3::dna5> &sequence) -> std::vector<double> {
      auto seq = to_string(sequence);
      std::string reverse_seq(seq);
      std::reverse(reverse_seq.begin(), reverse_seq.end());

      auto windows = pst::distances::sliding_windows_background<seqan3::dna5>(
          tree, seq, window_size, background_order);
      auto reverse_windows =
          pst::distances::sliding_windows_background<seqan3::dna5>(
              tree, reverse_seq, window_size, background_order);

      for (int i = 0; i < windows.size(); i++) {
        windows[i] = std::min(windows[i],
                              reverse_windows[reverse_windows.size() - 1 - i]);
      }

      return windows;
    };
  } else {
    sliding_window_fun =
        [window_size](
            tree_t &tree,
            std::vector<seqan3::dna5> &sequence) -> std::vector<double> {
      auto seq = to_string(sequence);
      std::string reverse_seq(seq);
      std::reverse(reverse_seq.begin(), reverse_seq.end());
      auto windows =
          pst::distances::sliding_windows<seqan3::dna5>(tree, seq, window_size);
      auto reverse_windows = pst::distances::sliding_windows<seqan3::dna5>(
          tree, reverse_seq, window_size);

      for (int i = 0; i < windows.size(); i++) {
        std::cout << windows[i] << " " << reverse_windows[reverse_windows.size() - 1 - i] << std::endl;
        windows[i] = std::min(windows[i],
                              reverse_windows[reverse_windows.size() - 1 - i]);
      }
      return windows;
    };
  }

  return sliding_window_fun;
}

std::tuple<std::vector<std::vector<std::vector<double>>>,
           std::vector<std::string>>
score_sliding_windows(std::vector<tree_t> &trees,
                      std::vector<std::filesystem::path> &sequence_list,
                      size_t background_order, int window_size,
                      bool background_adjusted) {
  std::vector<std::vector<std::vector<double>>> all_scores{};
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

    std::vector<std::vector<std::vector<double>>> scores(
        sequences.size(), std::vector<std::vector<double>>(trees.size()));

    std::function<std::vector<double>(tree_t &, std::vector<seqan3::dna5> &)>
        sliding_window_function =
            parse_sliding_window_function(sequences.size(), background_adjusted,
                                          background_order, window_size);

    auto fun = [&](size_t start_index, size_t stop_index) {
      for (size_t i = 0; i < trees.size(); i++) {
        for (size_t j = start_index; j < stop_index; j++) {
          scores[j][i] = sliding_window_function(trees[i], sequences[j]);
        }
      }
    };
    pst::parallelize::parallelize(sequences.size(), fun);

    std::move(scores.begin(), scores.end(), std::back_inserter(all_scores));
  }

  return {all_scores, ids};
}

int main(int argc, char *argv[]) {
  input_arguments arguments = parse_cli_arguments(argc, argv);

  std::vector<tree_t> trees;

  if (arguments.outpath.empty()) {
    std::cerr << "Error: " << arguments.outpath << " is not a file.";
    return EXIT_FAILURE;
  }
  if (!std::filesystem::exists(arguments.filepath)) {
    std::cerr << "Error: " << arguments.filepath << " is not a file.";
    return EXIT_FAILURE;
  }
  if (!std::filesystem::exists(arguments.sequence_list)) {
    std::cerr << "Error: " << arguments.sequence_list << " is not a file.";
    return EXIT_FAILURE;
  }

  if (arguments.filepath.extension() == ".h5" ||
      arguments.filepath.extension() == ".hdf5") {

    HighFive::File file{arguments.filepath, HighFive::File::ReadOnly};
    trees = get_trees(file, arguments.pseudo_count_amount);
  } else if (arguments.filepath.extension() == ".tree" ||
             arguments.filepath.extension() == ".bintree") {
    pst::ProbabilisticSuffixTreeMap<seqan3::dna5> tree{
        arguments.filepath, arguments.pseudo_count_amount};

    trees = std::vector<tree_t>{std::move(tree)};
  } else {
    std::cerr << "Error: " << arguments.filepath
              << " has invalid extension, and can't be parsed.";
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

  auto [scores, ids] = score_sliding_windows(
      trees, sequence_list, arguments.background_order, arguments.window_size,
      arguments.background_adjusted);

  HighFive::File out_file{arguments.outpath, HighFive::File::OpenOrCreate};

  if (!out_file.exist("windows")) {
    std::vector<size_t> dims{scores.size(),
                             static_cast<size_t>(scores[0].size())};
    out_file.createGroup("windows");
  }
  auto windows_group = out_file.getGroup("windows");
  for (int i = 0; i < ids.size(); i++) {
    auto id_ = ids[i];
    if (!windows_group.exist(id_)) {

      windows_group.createGroup(id_);
    }
    auto id_group = windows_group.getGroup(id_);
    for (int j = 0; j < trees.size(); j++) {
      auto tree_name = trees[j].id;
      if (!id_group.exist(tree_name)) {
        std::vector<size_t> dims{scores.size(),
                                 static_cast<size_t>(scores[0].size())};
        id_group.createGroup(tree_name);
      }
      auto tree_group = id_group.getGroup(tree_name);
      if (!tree_group.exist("window_scores")) {
        tree_group.createDataSet<double>(
            "window_scores", HighFive::DataSpace::From(scores[i][j]));
      }
      auto scores_dataset = tree_group.getDataSet("window_scores");
      scores_dataset.write(scores[i][j]);
    }
  }

  return EXIT_SUCCESS;
}
