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

#include <build_vlmc.hpp>
#include <kmer_container.hpp>

#include "distances/negative_log_likelihood.hpp"
#include "kl_tree.hpp"
#include "kl_tree_map.hpp"
#include "probabilistic_suffix_tree_map.hpp"

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
  std::filesystem::path tmp_path{};
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
  size_t fifth_percentile;

  void output_result() const {
    std::cout << "Min count: " << this->min_count
              << ". Max depth: " << this->max_depth
              << ". Threshold: " << this->threshold
              << ". BIC score: " << this->bic_score
              << ". AIC score: " << this->aic_score
              << ". AICc score: " << this->aicc_score
              << ". Log likelihood: " << this->log_likelihood
              << ". N parameters: " << this->n_parameters
              << ". Fifth percentile: " << this->fifth_percentile << std::endl;
  }
};

size_t percentile_frequency(pst::ProbabilisticSuffixTreeMap<seqan3::dna5> &tree,
                            float percentile) {
  std::vector<size_t> counts{};
  for (auto &[_c, v] : tree.counts) {
    auto &[count, _p, included] = v;
    if (!included) {
      continue;
    }
    counts.push_back(count);
  }

  std::sort(counts.begin(), counts.end());

  size_t percentile_edge = counts.size() * percentile;
  return counts[percentile_edge];
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

  auto fifth_percentile_frequency = percentile_frequency(tree, 0.05);

  Result res{min_count,      max_depth,    threshold,
             bic_value,      aic_value,    aicc_value,
             log_likelihood, n_parameters, fifth_percentile_frequency};

  return res;
}

Result bic(std::string &id, lst::details::sequence_t<seqan3::dna5> &sequence,
           double threshold, const std::filesystem::path &out_path) {
  std::vector<int> max_depths{3,  4,  5,  6,  7,  8,  9,  10, 11,
                              12, 13, 14, 15, 16, 17, 18, 19, 20};
  std::vector<int> min_counts{1,   2,   3,   4,   5,    6,    7,    8,    9,
                              10,  20,  30,  40,  50,   60,   70,   80,   90,
                              100, 125, 150, 175, 200,  250,  300,  400,  500,
                              600, 700, 800, 900, 1000, 2000, 3000, 4000, 5000};
  std::vector<float> thresholds{1.2, 3.9075};

  std::vector<Result> results{};
  Result best_result{0, 0, 0.0, std::numeric_limits<double>::max(), 0.0, 0, 0};

  // This search can be made much quicker by modifying the
  // Probabilistic suffix tree map class to allow iterative support pruning,
  // i.e., computing the tree for maximal parameters and then removing nodes.
  for (auto min_count : min_counts) {
    for (auto max_depth : max_depths) {
      for (auto threshold : thresholds) {
        // Don't search parameters than previously have shown to never be "the
        // best" to save computation
        //        if (sequence.size() > 387424359 && min_count < 10) {
        //          continue;
        //        } else if (sequence.size() > 3099706404 && min_count < 100)
        //        {
        //          continue;
        //        } else if (sequence.size() > 6699723695 && min_count < 500)
        //        {
        //          continue;
        //        }
        //
        //        if (sequence.size() < 3099706404 && min_count > 15) {
        //          continue;
        //        }
        pst::KullbackLieblerTreeMap<seqan3::dna5> tree{
            id, sequence, max_depth, min_count, threshold, true, 2};
        tree.construct_tree();

        auto res = bic_score(tree, sequence, min_count, max_depth, threshold);

        if (res.bic_score < best_result.bic_score) {
          best_result = res;
        }
        results.push_back(res);
      }
    }
  }

  std::stringstream ss;
  ss << "bic-" << id.substr(0, 20) << ".csv";
  std::string s = ss.str();

  std::filesystem::path path{out_path / s};

  bool add_header = !std::filesystem::exists(path);

  std::ofstream ofs(path, std::ios::app);

  if (add_header) {
    ofs << "n_params,bic_score,min_count,max_depth,threshold,log_likelihood,"
           "fifth_percentile_frequency"
        << std::endl;
  }

  for (auto [min_count, max_depth, threshold, score, aic_score, aicc_score,
             log_likelihood, n_params, fifth_percentile] : results) {
    ofs << n_params << "," << score << "," << min_count << "," << max_depth
        << "," << threshold << "," << log_likelihood << "," << fifth_percentile
        << std::endl;
  }
  ofs.close();

  std::cout << "---------- Best bic ----------" << std::endl;
  best_result.output_result();

  return best_result;
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

  std::vector<double> thresholds{start_threshold};
  std::sort(thresholds.begin(), thresholds.end());

  std::vector<int> min_counts =
      determine_min_counts(min_min_count, max_min_count);

  auto kmc_db_path = vlmc::run_kmc(fasta_path, max_max_depth + 1, tmp_path,
                                   in_or_out_of_core, 2);

  vlmc::build_vlmc_from_kmc_db(fasta_path, max_max_depth, min_min_count,
                               thresholds[0], tree_path, tmp_path,
                               in_or_out_of_core, kmc_db_path);
  vlmc::remove_kmc_files(kmc_db_path);

  for (auto threshold : thresholds) {
    for (auto min_count : min_counts) {
      pst::KullbackLieblerTreeMap<seqan3::dna5> tree{tree_path};
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
  ss << "bic-" << fasta_path.stem() << ".csv";
  std::string s = ss.str();

  std::filesystem::path path{tree_path.parent_path() / s};

  bool add_header = !std::filesystem::exists(path);

  std::ofstream ofs(path, std::ios::app);

  if (add_header) {
    ofs << "n_params,bic_score,aic_score,aicc_score,min_count,max_depth,"
           "threshold,log_likelihood,"
           "fifth_percentile_frequency"
        << std::endl;
  }

  for (auto [min_count, max_depth, threshold, bic_score, aic_score, aicc_score,
             log_likelihood, n_params, fifth_percentile] : results) {
    ofs << n_params << "," << bic_score << "," << aic_score << "," << aicc_score
        << "," << min_count << "," << max_depth << "," << threshold << ","
        << log_likelihood << "," << fifth_percentile << std::endl;
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

  //  bic(arguments.id, arguments.sequence, arguments.threshold,
  //      arguments.out_path.parent_path());

  //  auto hashmap_done = std::chrono::steady_clock::now();

  //  std::chrono::duration<double> kmc_seconds = kmc_done - kmc_start;
  //  std::cout << "KMC time: " << kmc_seconds.count() << "s\n";
  //  std::chrono::duration<double> hashmap_seconds = hashmap_done - kmc_done;
  //  std::cout << "Hashmap time: " << hashmap_seconds.count() << "s\n";

  return EXIT_SUCCESS;
}
