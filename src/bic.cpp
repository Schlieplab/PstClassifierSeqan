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
  double log_likelihood;
  uint64 n_parameters;
  double fifth_percentile;

  void output_result() const {
    std::cout << "Min count: " << this->min_count
              << ". Max depth: " << this->max_depth
              << ". Threshold: " << this->threshold
              << ". Bic score: " << this->bic_score
              << ". Log likelihood: " << this->log_likelihood
              << ". N parameters: " << this->n_parameters
              << ". Fifth percentile: " << this->fifth_percentile << std::endl;
  }
};

size_t percentile_frequency(pst::ProbabilisticSuffixTreeMap<seqan3::dna5> &tree,
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

Result bic_score(std::string &id,
                 lst::details::sequence_t<seqan3::dna5> &sequence,
                 size_t min_count, size_t max_depth, float threshold) {
  pst::KullbackLieblerTreeMap<seqan3::dna5> tree{
      id, sequence, max_depth, min_count, threshold, true, 2};
  tree.construct_tree();
  double log_likelihood = pst::distances::log_likelihood(tree, sequence);

  double n_parameters = 3 * tree.count_terminal_nodes();

  double bic_value =
      n_parameters * std::log(sequence.size()) - 2 * log_likelihood;

  auto fifth_percentile_frequency = percentile_frequency(tree, 0.05);

  Result res{min_count,
             max_depth,
             threshold,
             bic_value,
             log_likelihood,
             n_parameters,
             fifth_percentile_frequency};

  res.output_result();

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

        auto res = bic_score(id, sequence, min_count, max_depth, threshold);

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

  for (auto [min_count, max_depth, threshold, score, log_likelihood, n_params,
             fifth_percentile] : results) {
    ofs << n_params << "," << score << "," << min_count << "," << max_depth
        << "," << threshold << "," << log_likelihood << "," << fifth_percentile
        << std::endl;
  }
  ofs.close();

  std::cout << "---------- Best bic ----------" << std::endl;
  best_result.output_result();

  return best_result;
}

std::tuple<int, int, double> find_best_parameters_bic(
    const std::filesystem::path &fasta_path, const int max_max_depth,
    const int min_min_count, const std::filesystem::path &tree_path,
    const std::filesystem::path &tmp_path, const vlmc::Core &in_or_out_of_core,
    lst::details::sequence_t<seqan3::dna5> &sequence) {

  std::vector<Result> results{};

  Result best_result{0, 0, 0.0, std::numeric_limits<double>::max(), 0.0, 0, 0};

  std::vector<double> thresholds{1.2, 3.9075};
  std::vector<int> min_counts{2,   3,   4,   5,    6,    7,    8,    9,   10,
                              20,  30,  40,  50,   60,   70,   80,   90,  100,
                              125, 150, 175, 200,  250,  300,  400,  500, 600,
                              700, 800, 900, 1000, 2000, 3000, 4000, 5000};

  for (int max_depth = 3; max_depth < max_max_depth; max_depth++) {
    const int kmer_size = max_depth + 1;

    auto kmc_db_path =
        run_kmc(fasta_path, kmer_size, tmp_path, in_or_out_of_core, 2);

    for (auto min_count : min_counts) {
      if (min_count < min_min_count) {
        continue;
      }
      if (max_depth > 12 && min_count < 5) {
        continue;
      }
      for (auto threshold : thresholds) {
        vlmc::build_vlmc_from_kmc_db(fasta_path, max_depth, min_count,
                                     threshold, tree_path, tmp_path,
                                     in_or_out_of_core, kmc_db_path);

        auto parse_to_tree_start = std::chrono::steady_clock::now();

        pst::ProbabilisticSuffixTreeMap<seqan3::dna5> tree{tree_path};

        auto parse_to_tree_end = std::chrono::steady_clock::now();

        std::chrono::duration<double> to_tree_seconds =
            parse_to_tree_end - parse_to_tree_start;
        std::cout << "Parse to tree time: " << to_tree_seconds.count() << "s\n";

        auto nll_start = std::chrono::steady_clock::now();

        double log_likelihood = pst::distances::log_likelihood(tree, sequence);

        auto nll_done = std::chrono::steady_clock::now();

        std::chrono::duration<double> nll_seconds = nll_done - nll_start;
        std::cout << "NLL time: " << nll_seconds.count() << "s\n";

        uint64 n_parameters = 3 * tree.count_terminal_nodes();

        double bic_score =
            n_parameters * std::log(sequence.size()) - 2 * log_likelihood;

        auto fifth_percentile_frequency = percentile_frequency(tree, 0.05);

        Result res{min_count,
                   max_depth,
                   threshold,
                   bic_score,
                   log_likelihood,
                   n_parameters,
                   fifth_percentile_frequency};

        if (bic_score < best_result.bic_score) {
          best_result = res;
        }

        results.push_back(res);

        res.output_result();
      }
    }

    vlmc::remove_kmc_files(kmc_db_path);
  }

  std::stringstream ss;
  ss << "bic-" << fasta_path.stem() << ".csv";
  std::string s = ss.str();

  std::filesystem::path path{tree_path.parent_path() / s};

  bool add_header = !std::filesystem::exists(path);

  std::ofstream ofs(path, std::ios::app);

  if (add_header) {
    ofs << "n_params,bic_score,min_count,max_depth,threshold,log_likelihood,"
           "fifth_percentile_frequency"
        << std::endl;
  }

  for (auto [min_count, max_depth, threshold, score, log_likelihood, n_params,
             fifth_percentile] : results) {
    ofs << n_params << "," << score << "," << min_count << "," << max_depth
        << "," << threshold << "," << log_likelihood << "," << fifth_percentile
        << std::endl;
  }
  ofs.close();

  std::cout << "---------- Best bic ----------" << std::endl;
  best_result.output_result();

  return {best_result.max_depth, best_result.min_count, best_result.threshold};
}

int main(int argc, char *argv[]) {
  input_arguments arguments = parse_cli_arguments(argc, argv);

  vlmc::configure_stxxl(arguments.tmp_path);

  find_best_parameters_bic(arguments.filename, 20, 100, arguments.out_path,
                           arguments.tmp_path, vlmc::Core::out,
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
