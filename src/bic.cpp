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

#include "distances/negative_log_likelihood.hpp"
#include "kl_tree.hpp"
#include "kl_tree_map.hpp"
#include "probabilistic_suffix_tree_map.hpp"

using seqan3::operator""_dna5;

struct input_arguments {
  lst::details::sequence_t<seqan3::dna5> sequence{};
  std::string id{};
};

struct my_traits : seqan3::sequence_file_input_default_traits_dna {
  template <typename alph>
  using sequence_container =
      std::vector<alph>; // must be defined as a template!
};

input_arguments parse_cli_arguments(int argc, char *argv[]) {
  std::filesystem::path filename{};

  input_arguments arguments{};

  seqan3::argument_parser parser{"BIC", argc, argv,
                                 seqan3::update_notifications::off};
  parser.info.short_description =
      "Calculate BIC for a sequence under the variable-length Markov chain.";

  parser.add_positional_option(filename, "path to fasta file.");

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

size_t percentile_frequency(pst::KullbackLieblerTreeMap<seqan3::dna5> &tree,
                            float percentile) {
  std::vector<size_t> counts(tree.counts.size());
  int idx = 0;
  for (auto &[_c, v] : tree.counts) {
    auto &[count, _p, included] = v;
    if (included) {
      continue;
    }
    counts[idx] = count;
    idx++;
  }

  std::sort(counts.begin(), counts.end());

  size_t percentile_edge = counts.size() * percentile;
  return counts[percentile_edge];
}

std::tuple<double, size_t, double, size_t>
bic_score(std::string &id, lst::details::sequence_t<seqan3::dna5> &sequence,
          size_t min_count, size_t max_depth, float threshold) {
  pst::KullbackLieblerTreeMap<seqan3::dna5> tree{
      id, sequence, max_depth, min_count, threshold, true, 2};
  tree.construct_tree();
  double log_likelihood = pst::distances::log_likelihood(tree, sequence);

  double n_parameters = 3 * tree.count_terminal_nodes();

  double bic_value =
      n_parameters * std::log(sequence.size()) - 2 * log_likelihood;

  auto fifth_percentile_frequency = percentile_frequency(tree, 0.05);

  std::cout << n_parameters << " " << bic_value << " " << min_count << " "
            << max_depth << " " << threshold << " " << log_likelihood << " "
            << fifth_percentile_frequency << std::endl;

  return {bic_value, n_parameters, log_likelihood, fifth_percentile_frequency};
}

std::tuple<int, int, double, int, float>
bic(std::string &id, lst::details::sequence_t<seqan3::dna5> &sequence) {
  std::vector<int> max_depths{3,  4,  5,  6,  7,  8,  9,  10, 11,
                              12, 13, 14, 15, 16, 17, 18, 19, 20};
  std::vector<int> min_counts{2,   3,   4,   5,    6,    7,    8,    9,   10,
                              20,  30,  40,  50,   60,   70,   80,   90,  100,
                              125, 150, 175, 200,  250,  300,  400,  500, 600,
                              700, 800, 900, 1000, 2000, 3000, 4000, 5000};
  std::vector<float> thresholds{3.9075};

  std::vector<std::tuple<int, int, double, size_t, float, float, size_t>>
      results{};

  // This search can be made much quicker by modifying the
  // Probabilistic suffix tree map class to allow iterative support pruning,
  // i.e., computing the tree for maximal parameters and then removing nodes.
  for (auto min_count : min_counts) {
    for (auto max_depth : max_depths) {
      for (auto threshold : thresholds) {
        // Don't search parameters than previously have shown to never be "the
        // best" to save computation
        if (sequence.size() > 387424359 && min_count < 10) {
          continue;
        } else if (sequence.size() > 3099706404 && min_count < 100) {
          continue;
        } else if (sequence.size() > 6699723695 && min_count < 500) {
          continue;
        }

        if (sequence.size() < 3099706404 && min_count > 15) {
          continue;
        }

        auto [score, n_parameters, log_likelihood, fifth_percentile_frequency] =
            bic_score(id, sequence, min_count, max_depth, threshold);

        results.emplace_back(min_count, max_depth, score, n_parameters,
                             threshold, log_likelihood,
                             fifth_percentile_frequency);
      }
    }
  }

  std::stringstream ss;
  ss << "bic-" << id.substr(0, 20) << ".csv";
  std::string s = ss.str();

  std::filesystem::path path{s};

  bool add_header = !std::filesystem::exists(path);

  std::ofstream ofs(path, std::ios::app);

  if (add_header) {
    ofs << "n_params,bic_score,min_count,max_depth,threshold,log_likelihood,"
           "fifth_percentile_frequency"
        << std::endl;
  }

  auto [best_min_count, best_max_depth, best_bic, best_n_params, best_threshold,
        best_log, _] = results[0];
  for (auto [min_count, max_depth, score, n_params, threshold, log_likelihood,
             fifth_percentile_frequency] : results) {
    ofs << n_params << "," << score << "," << min_count << "," << max_depth
        << "," << threshold << "," << log_likelihood << ","
        << fifth_percentile_frequency << std::endl;

    if (score < best_bic) {
      best_bic = score;
      best_max_depth = max_depth;
      best_min_count = min_count;
      best_n_params = n_params;
      best_threshold = threshold;
      best_log = log_likelihood;
    }
  }
  ofs.close();

  return {best_min_count, best_max_depth, best_bic, best_n_params, best_log};
}

int main(int argc, char *argv[]) {
  input_arguments arguments = parse_cli_arguments(argc, argv);

  auto [min_count, max_depth, score, n_params, log_likelihood] =
      bic(arguments.id, arguments.sequence);
  std::cout << "min count: " << min_count << " max depth: " << max_depth
            << " bic: " << score << " parameters: " << n_params
            << " Log likelihood: " << log_likelihood
            << " log sequence: " << std::log(arguments.sequence.size())
            << std::endl;

  return EXIT_SUCCESS;
}
