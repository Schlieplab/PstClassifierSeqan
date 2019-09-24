#include <benchmark/benchmark.h>

#include <string>
#include <vector>

#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/io/sequence_file/input.hpp>
#include <seqan3/range/container/bitcompressed_vector.hpp>

#include "../src/probabilistic_suffix_tree.hpp"

struct my_traits : seqan3::sequence_file_input_default_traits_dna {
  template <typename alph>
  using sequence_container =
      seqan3::bitcompressed_vector<alph>; // must be defined as a template!
};

static void HumanHerpesvirus5(benchmark::State &state) {
  std::string filename{"../NC_006273.2.fa"};
  seqan3::sequence_file_input<my_traits> file_in{filename};
  std::vector<seqan3::bitcompressed_vector<seqan3::dna5>> sequences{};
  std::vector<std::string> ids{};

  for (auto &[seq, id, qual] : file_in) {
    sequences.push_back(seq);
    ids.push_back(id);
  }

  size_t max_depth = 15;
  size_t min_count = 100;
  float threshold = 1.2;

  for (auto _ : state) {
    pst::ProbabilisticSuffixTree<seqan3::dna5> pst{
        ids[0], sequences[0], max_depth, min_count, threshold};
  }
}

BENCHMARK(HumanHerpesvirus5);

static void SaccharomycesCerevisiae(benchmark::State &state) {
  std::string filename{"../CM010781.1.fa"};
  seqan3::sequence_file_input<my_traits> file_in{filename};
  std::vector<seqan3::bitcompressed_vector<seqan3::dna5>> sequences{};
  std::vector<std::string> ids{};

  for (auto &[seq, id, qual] : file_in) {
    sequences.push_back(seq);
    ids.push_back(id);
  }

  size_t max_depth = 15;
  size_t min_count = 100;
  float threshold = 1.2;

  for (auto _ : state) {
    pst::ProbabilisticSuffixTree<seqan3::dna5> pst{
        ids[0], sequences[0], max_depth, min_count, threshold};
  }
}

BENCHMARK(SaccharomycesCerevisiae);

static void EColi(benchmark::State &state) {
  std::string filename{"../CP007136.1.fa"};
  seqan3::sequence_file_input<my_traits> file_in{filename};
  std::vector<seqan3::bitcompressed_vector<seqan3::dna5>> sequences{};
  std::vector<std::string> ids{};

  for (auto &[seq, id, qual] : file_in) {
    sequences.push_back(seq);
    ids.push_back(id);
  }

  size_t max_depth = 15;
  size_t min_count = 100;
  float threshold = 1.2;

  for (auto _ : state) {
    pst::ProbabilisticSuffixTree<seqan3::dna5> pst{
        ids[0], sequences[0], max_depth, min_count, threshold};
  }
}

BENCHMARK(EColi);

static void HumanChr17(benchmark::State &state) {
  std::string filename{"../NC_000017.11.fa"};
  seqan3::sequence_file_input<my_traits> file_in{filename};
  std::vector<seqan3::bitcompressed_vector<seqan3::dna5>> sequences{};
  std::vector<std::string> ids{};

  for (auto &[seq, id, qual] : file_in) {
    sequences.push_back(seq);
    ids.push_back(id);
  }

  size_t max_depth = 15;
  size_t min_count = 100;
  float threshold = 1.2;

  for (auto _ : state) {
    pst::ProbabilisticSuffixTree<seqan3::dna5> pst{
        ids[0], sequences[0], max_depth, min_count, threshold};
  }
}

BENCHMARK(HumanChr17);

static void HumanChr1(benchmark::State &state) {
  std::string filename{"../NC_000001.11.fa"};
  seqan3::sequence_file_input<my_traits> file_in{filename};
  std::vector<seqan3::bitcompressed_vector<seqan3::dna5>> sequences{};
  std::vector<std::string> ids{};

  for (auto &[seq, id, qual] : file_in) {
    sequences.push_back(seq);
    ids.push_back(id);
  }

  size_t max_depth = 15;
  size_t min_count = 100;
  float threshold = 1.2;

  for (auto _ : state) {
    pst::ProbabilisticSuffixTree<seqan3::dna5> pst{
        ids[0], sequences[0], max_depth, min_count, threshold};
  }
}

BENCHMARK(HumanChr1);

BENCHMARK_MAIN();