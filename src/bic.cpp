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

  seqan3::argument_parser parser{"BIC", argc, argv, false};
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

std::tuple<float, size_t>
bic_score(std::string &id, lst::details::sequence_t<seqan3::dna5> &sequence,
          size_t min_count, size_t max_depth, float threshold) {
  pst::KullbackLieblerTreeMap<seqan3::dna5> pst{
      id, sequence, max_depth, min_count, threshold, true, 2};
  pst.construct_tree();

  double likelihood = -pst::distances::negative_log_likelihood(pst, sequence) *
                      double(sequence.size());

  float n_parameters = pst.count_terminal_nodes();

  double bic_value = n_parameters * std::log(sequence.size()) - 2 * likelihood;
  std::cout << n_parameters << " " << bic_value << " " << min_count << " "
            << max_depth << " " << threshold << std::endl;

  return {bic_value, n_parameters};
}

std::tuple<int, int, float, int>
bic(std::string &id, lst::details::sequence_t<seqan3::dna5> &sequence) {
  std::vector<int> max_depths{9, 12, 15};
  std::vector<int> min_counts{2, 3, 4, 5, 6, 7, 8, 9, 10, 20, 50, 100, 1000};
  std::vector<float> thresholds{3.9075};

  std::vector<std::tuple<int, int, float, size_t, float>> results{};

  for (auto min_count : min_counts) {
    for (auto max_depth : max_depths) {
      for (auto threshold : thresholds) {
        auto [score, n_parameters] =
            bic_score(id, sequence, min_count, max_depth, threshold);
        results.emplace_back(min_count, max_depth, score, n_parameters,
                             threshold);
      }
    }
  }

  std::stringstream ss;
  ss << "bic-" << id << ".csv";
  std::string s = ss.str();

  std::filesystem::path path{s};
  std::ofstream ofs(path);

  ofs << "n_params,bic_score,min_count,max_depth,threshold" << std::endl;

  auto [best_min_count, best_max_depth, best_bic, best_n_params,
        best_threshold] = results[0];
  for (auto [min_count, max_depth, score, n_params, threshold] : results) {
    if (score < best_bic) {
      best_bic = score;
      best_max_depth = max_depth;
      best_min_count = min_count;
      best_n_params = n_params;
      best_threshold = threshold;
    }
    ofs << n_params << "," << score << "," << min_count << "," << max_depth
        << "," << threshold << std::endl;
  }
  ofs.close();

  return {best_min_count, best_max_depth, best_bic, best_n_params};
}

int main(int argc, char *argv[]) {
  input_arguments arguments = parse_cli_arguments(argc, argv);

  auto [min_count, max_depth, score, n_params] =
      bic(arguments.id, arguments.sequence);
  std::cout << "min count: " << min_count << " max depth: " << max_depth
            << " bic: " << score << " parameters: " << n_params << std::endl;

  return EXIT_SUCCESS;
}
