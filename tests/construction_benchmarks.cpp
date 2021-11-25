#include <benchmark/benchmark.h>

#include <string>
#include <vector>

#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/io/sequence_file/input.hpp>

#include "../src/kl_tree_map.hpp"
#include "../src/probabilistic_suffix_tree_map.hpp"

struct my_traits : seqan3::sequence_file_input_default_traits_dna {
  template <typename alph>
  using sequence_container =
      std::vector<alph>; // must be defined as a template!
};


void test_benchmark_map(benchmark::State &state, std::string filename,
                        bool parallel, int parallel_depth) {
  seqan3::sequence_file_input<my_traits> file_in{filename};
  std::vector<lst::details::sequence_t<seqan3::dna5>> sequences{};
  std::vector<std::string> ids{};

  for (auto &[seq, id, qual] : file_in) {
    sequences.emplace_back(std::move(seq));
    ids.push_back(id);
  }

  size_t max_depth = 15;
  size_t min_count = 100;
  float threshold = 1.2;

  for (auto _ : state) {
    pst::KullbackLieblerTreeMap<seqan3::dna5> pst{
        ids[0],    sequences[0], max_depth,     min_count,
        threshold, parallel,     parallel_depth};
    pst.construct_tree();
  }
}

static void HumanHerpesvirus5Map(benchmark::State &state) {
  std::string filename{"../fasta/NC_006273.2.fa"};
  test_benchmark_map(state, filename, false, 0);
}

BENCHMARK(HumanHerpesvirus5Map);

static void SaccharomycesCerevisiaeMap(benchmark::State &state) {
  std::string filename{"../fasta/CM010781.1.fa"};
  test_benchmark_map(state, filename, false, 0);
}

BENCHMARK(SaccharomycesCerevisiaeMap);

static void EColiMap(benchmark::State &state) {
  std::string filename{"../fasta/CP007136.1.fa"};
  test_benchmark_map(state, filename, false, 0);
}

BENCHMARK(EColiMap);

static void HumanChr17Map(benchmark::State &state) {
  std::string filename{"../fasta/NC_000017.11.fa"};
  test_benchmark_map(state, filename, false, 0);
}

BENCHMARK(HumanChr17Map);

static void HumanChr1Map(benchmark::State &state) {
  std::string filename{"../fasta/NC_000001.11.fa"};
  test_benchmark_map(state, filename, false, 0);
}

// BENCHMARK(HumanChr1Map);

static void HumanHerpesvirus5ParallelMap(benchmark::State &state) {
  std::string filename{"../fasta/NC_006273.2.fa"};
  test_benchmark_map(state, filename, true, 1);
}

BENCHMARK(HumanHerpesvirus5ParallelMap);

static void SaccharomycesCerevisiaeParallelMap(benchmark::State &state) {
  std::string filename{"../fasta/CM010781.1.fa"};
  test_benchmark_map(state, filename, true, 1);
}

BENCHMARK(SaccharomycesCerevisiaeParallelMap);

static void EColiParallelMap(benchmark::State &state) {
  std::string filename{"../fasta/CP007136.1.fa"};
  test_benchmark_map(state, filename, true, 1);
}

BENCHMARK(EColiParallelMap);

static void HumanChr17ParallelMap(benchmark::State &state) {
  std::string filename{"../fasta/NC_000017.11.fa"};
  test_benchmark_map(state, filename, true, 2);
}

BENCHMARK(HumanChr17ParallelMap);

static void HumanChr1ParallelMap(benchmark::State &state) {
  std::string filename{"../fasta/NC_000001.11.fa"};
  test_benchmark_map(state, filename, true, 2);
}

// BENCHMARK(HumanChr1ParallelMap);

static void HumanHerpesvirus5Parallel2Map(benchmark::State &state) {
  std::string filename{"../fasta/NC_006273.2.fa"};
  test_benchmark_map(state, filename, true, 2);
}

BENCHMARK(HumanHerpesvirus5Parallel2Map);

static void SaccharomycesCerevisiaeParallel2Map(benchmark::State &state) {
  std::string filename{"../fasta/CM010781.1.fa"};
  test_benchmark_map(state, filename, true, 2);
}

BENCHMARK(SaccharomycesCerevisiaeParallel2Map);

static void EColiParallel2Map(benchmark::State &state) {
  std::string filename{"../fasta/CP007136.1.fa"};
  test_benchmark_map(state, filename, true, 2);
}

BENCHMARK(EColiParallel2Map);

static void HumanChr17Parallel2Map(benchmark::State &state) {
  std::string filename{"../fasta/NC_000017.11.fa"};
  test_benchmark_map(state, filename, true, 2);
}

BENCHMARK(HumanChr17Parallel2Map);

static void HumanChr1Parallel2Map(benchmark::State &state) {
  std::string filename{"../fasta/NC_000001.11.fa"};
  test_benchmark_map(state, filename, true, 2);
}

// BENCHMARK(HumanChr1Parallel2Map);

BENCHMARK_MAIN();