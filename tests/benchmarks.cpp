#include <benchmark/benchmark.h>

#include <string>
#include <vector>

#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/io/sequence_file/input.hpp>
#include <seqan3/range/container/bitcompressed_vector.hpp>

#include "../src/kl_tree.hpp"
#include "../src/probabilistic_suffix_tree.hpp"

struct my_traits : seqan3::sequence_file_input_default_traits_dna {
  template <typename alph>
  using sequence_container =
      std::vector<alph>; // must be defined as a template!
};

void test_benchmark(benchmark::State &state, std::string filename,
                    bool parallel, int parallel_depth) {
  pst::time_measurement = false;
  seqan3::sequence_file_input<my_traits> file_in{filename};
  std::vector<seqan3::bitcompressed_vector<seqan3::dna5>> sequences{};
  std::vector<std::string> ids{};

  for (auto &[seq, id, qual] : file_in) {
    sequences.push_back(
        seqan3::bitcompressed_vector<seqan3::dna5>{std::move(seq)});
    ids.push_back(id);
  }

  size_t max_depth = 15;
  size_t min_count = 100;
  float threshold = 1.2;

  for (auto _ : state) {
    pst::KullbackLieblerTree<seqan3::dna5> pst{
        ids[0],    sequences[0], max_depth,     min_count,
        threshold, parallel,     parallel_depth};
    pst.construct_tree();
  }
}

static void HumanHerpesvirus5(benchmark::State &state) {
  std::string filename{"../fasta/NC_006273.2.fa"};
  test_benchmark(state, filename, false, 0);
}

BENCHMARK(HumanHerpesvirus5);

static void SaccharomycesCerevisiae(benchmark::State &state) {
  std::string filename{"../fasta/CM010781.1.fa"};
  test_benchmark(state, filename, false, 0);
}

BENCHMARK(SaccharomycesCerevisiae);

static void EColi(benchmark::State &state) {
  std::string filename{"../fasta/CP007136.1.fa"};
  test_benchmark(state, filename, false, 0);
}

BENCHMARK(EColi);

static void HumanChr17(benchmark::State &state) {
  std::string filename{"../fasta/NC_000017.11.fa"};
  test_benchmark(state, filename, false, 0);
}

// BENCHMARK(HumanChr17);

static void HumanChr1(benchmark::State &state) {
  std::string filename{"../fasta/NC_000001.11.fa"};
  test_benchmark(state, filename, false, 0);
}

// BENCHMARK(HumanChr1);

static void HumanHerpesvirus5Parallel(benchmark::State &state) {
  std::string filename{"../fasta/NC_006273.2.fa"};
  test_benchmark(state, filename, true, 1);
}

BENCHMARK(HumanHerpesvirus5Parallel);

static void SaccharomycesCerevisiaeParallel(benchmark::State &state) {
  std::string filename{"../fasta/CM010781.1.fa"};
  test_benchmark(state, filename, true, 1);
}

BENCHMARK(SaccharomycesCerevisiaeParallel);

static void EColiParallel(benchmark::State &state) {
  std::string filename{"../fasta/CP007136.1.fa"};
  test_benchmark(state, filename, true, 1);
}

BENCHMARK(EColiParallel);

static void HumanChr17Parallel(benchmark::State &state) {
  std::string filename{"../fasta/NC_000017.11.fa"};
  test_benchmark(state, filename, true, 2);
}

// BENCHMARK(HumanChr17Parallel);

static void HumanChr1Parallel(benchmark::State &state) {
  std::string filename{"../fasta/NC_000001.11.fa"};
  test_benchmark(state, filename, true, 2);
}

// BENCHMARK(HumanChr1Parallel);

static void HumanHerpesvirus5Parallel2(benchmark::State &state) {
  std::string filename{"../fasta/NC_006273.2.fa"};
  test_benchmark(state, filename, true, 2);
}

BENCHMARK(HumanHerpesvirus5Parallel2);

static void SaccharomycesCerevisiaeParallel2(benchmark::State &state) {
  std::string filename{"../fasta/CM010781.1.fa"};
  test_benchmark(state, filename, true, 2);
}

BENCHMARK(SaccharomycesCerevisiaeParallel2);

static void EColiParallel2(benchmark::State &state) {
  std::string filename{"../fasta/CP007136.1.fa"};
  test_benchmark(state, filename, true, 2);
}

BENCHMARK(EColiParallel2);

static void HumanChr17Parallel2(benchmark::State &state) {
  std::string filename{"../fasta/NC_000017.11.fa"};
  test_benchmark(state, filename, true, 2);
}

// BENCHMARK(HumanChr17Parallel2);

static void HumanChr1Parallel2(benchmark::State &state) {
  std::string filename{"../fasta/NC_000001.11.fa"};
  test_benchmark(state, filename, true, 2);
}

// BENCHMARK(HumanChr1Parallel2);

BENCHMARK_MAIN();