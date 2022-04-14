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

#include <vlmc_from_kmers/build_vlmc.hpp>
#include <vlmc_from_kmers/kmer_container.hpp>

#include "pst/distances/negative_log_likelihood.hpp"
#include "pst/kl_tree.hpp"
#include "pst/kl_tree_map.hpp"
#include "pst/probabilistic_suffix_tree_map.hpp"

using seqan3::operator""_dna5;

struct input_arguments {
  lst::details::sequence_t<seqan3::dna5> sequence{};
  std::string id{};
  double threshold{3.9075};
  int min_max_depth{3};
  int max_max_depth{9};
  int min_min_count{2};
  int max_min_count{5000};
  std::filesystem::path filename{};
  std::filesystem::path out_path{};
  std::filesystem::path tmp_path{"test"};
};

struct my_traits : seqan3::sequence_file_input_default_traits_dna {
  template <typename alph>
  using sequence_container =
      std::vector<alph>; // must be defined as a template!
};

input_arguments parse_cli_arguments(int argc, char *argv[]) {
  input_arguments arguments{};

  seqan3::argument_parser parser{"BIC", argc, argv,
                                 seqan3::update_notifications::off};
  parser.info.short_description =
      "Calculate BIC for a sequence under the variable-length Markov chain.";

  parser.add_positional_option(arguments.filename, "path to fasta file.");
  parser.add_positional_option(
      arguments.out_path,
      "path to output file, has to have '.bintree' as extension.");
  parser.add_positional_option(arguments.tmp_path, "path to tmp path.");

  parser.add_option(arguments.threshold, 'k', "threshold",
                    "Threshold for the pruning stage of the algorithm.");

  parser.add_option(arguments.min_min_count, 'c', "min-min-count",
                    "Minimum min count to test.");
  parser.add_option(arguments.max_min_count, 'b', "max-min-count",
                    "Maximum min count to test.");

  parser.add_option(arguments.min_max_depth, 'd', "min-max-depth",
                    "Minimum max depth to test.");
  parser.add_option(arguments.max_max_depth, 'e', "max-max-depth",
                    "Maximum max depth to test.");

  try {
    parser.parse();
  } catch (seqan3::argument_parser_error const &ext) {
    std::cout << "[PARSER ERROR] " << ext.what() << '\n';
    return arguments;
  }

  seqan3::sequence_file_input<my_traits> file_in{arguments.filename};

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

struct Result {
  int min_count;
  int max_depth;
  double threshold;
  double bic_score;
  double aic_score;
  double aicc_score;
  double log_likelihood;
  uint64 n_parameters;
  size_t observed_min_count;
  size_t observed_max_depth;

  void output_result() const {
    std::cout << "{"
              << "\"Min count\": " << this->min_count << ","
              << "\"Max depth\": " << this->max_depth << ","
              << "\"Threshold\": " << this->threshold << ","
              << "\"BIC score\": " << this->bic_score << ","
              << "\"AIC score\": " << this->aic_score << ","
              << "\"AICc score\": " << this->aicc_score << ","
              << "\"Log likelihood\": " << this->log_likelihood << ","
              << "\"N parameters\": " << this->n_parameters << ","
              << "\"Observed min count\": " << this->observed_min_count << ","
              << "\"Observed max depth\": " << this->observed_max_depth << "}"
              << std::endl;
  }

  void write_csv_line(std::ofstream &ofs) const {
    ofs << n_parameters << ",";
    ofs << bic_score << ",";
    ofs << aic_score << ",";
    ofs << aicc_score << ",";
    ofs << min_count << ",";
    ofs << max_depth << ",";
    ofs << threshold << ",";
    ofs << log_likelihood << ",";
    ofs << observed_min_count << ",";
    ofs << observed_max_depth << std::endl;
  }

  static void write_header(std::ofstream &ofs) {
    ofs << "n_params,";
    ofs << "bic_score,";
    ofs << "aic_score,";
    ofs << "aicc_score,";
    ofs << "min_count,";
    ofs << "max_depth,";
    ofs << "threshold,";
    ofs << "log_likelihood,";
    ofs << "observed_min_count,";
    ofs << "observed_max_depth";
    ofs << std::endl;
  }
};

std::tuple<size_t, size_t>
get_min_count_max_depth(pst::ProbabilisticSuffixTreeMap<seqan3::dna5> &tree) {
  // This would have been cleaner as two loops, but when tree.counts is large
  // that would take quite some time.

  size_t min_count = std::numeric_limits<int>::max();
  size_t max_depth = 0;

  for (auto &[c, v] : tree.counts) {
    auto &[count, _p, included] = v;
    if (included) {
      //      std::cout << c <<  " " << c.size() << " " << count << std::endl;
      min_count = std::min(count, min_count);
      max_depth = std::max(c.size(), max_depth);
    }
  }

  return {min_count, max_depth};
}

Result bic_score(pst::KullbackLieblerTreeMap<seqan3::dna5> &tree,
                 lst::details::sequence_t<seqan3::dna5> &sequence,
                 int min_count, int max_depth, double threshold) {
  double log_likelihood = pst::distances::log_likelihood(tree, sequence);

  double n_parameters = 3 * tree.count_terminal_nodes();

  double bic_value =
      n_parameters * std::log(sequence.size()) - 2 * log_likelihood;

  double aic_value = 2 * n_parameters - 2 * log_likelihood;

  double aicc_value =
      aic_value + (2 * std::pow(n_parameters, 2.0) + 2 * n_parameters) /
                      (sequence.size() - n_parameters - 1);

  auto [obs_min_count, obs_max_depth] = get_min_count_max_depth(tree);

  Result res{min_count,     max_depth,    threshold,      bic_value,
             aic_value,     aicc_value,   log_likelihood, n_parameters,
             obs_min_count, obs_max_depth};

  return res;
}

std::vector<int> determine_min_counts(const int min_min_count,
                                      const int max_min_count) {
  std::vector<int> min_counts{};

  int current_min_count = min_min_count;

  while (current_min_count < max_min_count) {
    min_counts.push_back(current_min_count);

    if (current_min_count < 10) {
      current_min_count += 1;
    } else if (current_min_count < 100) {
      current_min_count += 10;
    } else if (current_min_count < 1000) {
      current_min_count += 100;
    } else {
      current_min_count += 1000;
    }
  }

  min_counts.push_back(max_min_count);

  return min_counts;
}

double g() {
  static double kl = 1.2;
  kl += 0.1;
  return kl;
}

std::tuple<int, int, double> find_best_parameters(
    const std::filesystem::path &fasta_path, const int min_max_depth,
    const int max_max_depth, const int min_min_count, const int max_min_count,
    const double start_threshold, const std::filesystem::path &tree_path,
    const std::filesystem::path &tmp_path, const vlmc::Core &in_or_out_of_core,
    lst::details::sequence_t<seqan3::dna5> &sequence) {

  std::vector<Result> results{};

  Result best_aic_result{0,   0,   0.0, 0.0, std::numeric_limits<double>::max(),
                         0.0, 0.0, 0,   0};
  Result best_aicc_result{
      0, 0, 0.0, 0.0, 0.0, std::numeric_limits<double>::max(), 0.0, 0, 0};
  Result best_bic_result{
      0, 0, 0.0, std::numeric_limits<double>::max(), 0.0, 0.0, 0.0, 0, 0};

  std::vector<double> thresholds(500);
  std::generate(thresholds.begin(), thresholds.end(), g);
  std::sort(thresholds.begin(), thresholds.end());

  std::vector<int> min_counts =
      determine_min_counts(min_min_count, max_min_count);

  std::cout << "Building vlmc" << std::endl;
  vlmc::build_vlmc(fasta_path, max_max_depth, min_min_count, thresholds[0],
                   tree_path, tmp_path, in_or_out_of_core);

  std::cout << "Finding min BIC" << std::endl;
  for (auto threshold : thresholds) {
    for (auto min_count : min_counts) {
      pst::KullbackLieblerTreeMap<seqan3::dna5> tree{tree_path};
      tree.multi_core = false;
      for (int max_depth = max_max_depth; max_depth >= min_max_depth;
           max_depth--) {
        // Prune for parameter settings
        tree.reprune_support(min_count, max_depth);
        tree.reprune_similarity(threshold);

        auto res = bic_score(tree, sequence, min_count, max_depth, threshold);

        if (res.bic_score < best_bic_result.bic_score) {
          best_bic_result = res;
        }

        if (res.aic_score < best_aic_result.aic_score) {
          best_aic_result = res;
        }

        if (res.aicc_score < best_aicc_result.aicc_score) {
          best_aicc_result = res;
        }

        results.push_back(res);

        res.output_result();
      }
    }
  }

  std::stringstream ss;
  ss << "bic-" << fasta_path.stem().string() << ".csv";
  std::string s = ss.str();

  std::filesystem::path path{tree_path.parent_path() / s};

  bool add_header = !std::filesystem::exists(path);

  std::ofstream ofs(path, std::ios::app);

  if (add_header) {
    Result::write_header(ofs);
  }

  for (auto result : results) {
    result.write_csv_line(ofs);
  }
  ofs.close();

  std::cout << "---------- Best BIC ----------" << std::endl;
  best_bic_result.output_result();
  std::cout << "---------- Best AIC ----------" << std::endl;
  best_aic_result.output_result();
  std::cout << "---------- Best AICc ----------" << std::endl;
  best_aicc_result.output_result();

  return {best_bic_result.max_depth, best_bic_result.min_count,
          best_bic_result.threshold};
}

int main(int argc, char *argv[]) {
  input_arguments arguments = parse_cli_arguments(argc, argv);

  vlmc::configure_stxxl(arguments.tmp_path);

  find_best_parameters(arguments.filename, arguments.min_max_depth,
                       arguments.max_max_depth, arguments.min_min_count,
                       arguments.max_min_count, arguments.threshold,
                       arguments.out_path, arguments.tmp_path, vlmc::Core::out,
                       arguments.sequence);

  return EXIT_SUCCESS;
}
