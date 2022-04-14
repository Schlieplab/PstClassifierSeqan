#include <benchmark/benchmark.h>

#include "pst/distances/cv.hpp"
#include "pst/distances/d2.hpp"
#include "pst/distances/dvstar.hpp"
#include "pst/distances/kl_divergence.hpp"
#include "pst/distances/negative_log_likelihood.hpp"
#include "pst/kl_tree.hpp"
#include "pst/kl_tree_map.hpp"

#include "random_sequence.hpp"

class D2StarBenchmarks : public benchmark::Fixture {
public:
  void SetUp(const ::benchmark::State &state) override {
    sakai = pst::KullbackLieblerTreeMap<seqan3::dna5>{sakai_path};
    ed1a = pst::KullbackLieblerTreeMap<seqan3::dna5>{ed1a_path};
  }
  std::filesystem::path sakai_path{"./../trees/Sakai.bintree"};
  std::filesystem::path ed1a_path{"./../trees/ED1a.bintree"};
  pst::KullbackLieblerTreeMap<seqan3::dna5> sakai;
  pst::KullbackLieblerTreeMap<seqan3::dna5> ed1a;
};

static void CV(benchmark::State &state) {
  std::filesystem::path first_path{"./../trees/CM008035.1.tree"};
  std::filesystem::path second_path{"./../trees/NC_009067.tree"};
  pst::KullbackLieblerTreeMap<seqan3::dna5> first{first_path};
  pst::KullbackLieblerTreeMap<seqan3::dna5> second{second_path};

  for (auto _ : state) {
    benchmark::DoNotOptimize(pst::distances::cv(first, second));
  }
}
BENCHMARK(CV);

static void D2(benchmark::State &state) {
  std::filesystem::path first_path{"./../trees/CM008035.1.tree"};
  std::filesystem::path second_path{"./../trees/NC_009067.tree"};
  pst::KullbackLieblerTreeMap<seqan3::dna5> first{first_path};
  pst::KullbackLieblerTreeMap<seqan3::dna5> second{second_path};

  for (auto _ : state) {
    benchmark::DoNotOptimize(pst::distances::d2(first, second));
  }
}
BENCHMARK(D2);

BENCHMARK_F(D2StarBenchmarks, Dvstar)
(benchmark::State &state) {
  for (auto _ : state) {
    benchmark::DoNotOptimize(pst::distances::dvstar(sakai, ed1a));
  }
}

BENCHMARK_DEFINE_F(D2StarBenchmarks, DvstarCore)
(benchmark::State &state) {
  auto background_order = state.range(0);
  for (auto _ : state) {
    benchmark::DoNotOptimize(
        pst::distances::details::dvstar::core_dvstar<seqan3::dna5>(
            sakai, ed1a, background_order));
  }
}
BENCHMARK_REGISTER_F(D2StarBenchmarks, DvstarCore)->DenseRange(0, 3, 1);

BENCHMARK_F(D2StarBenchmarks, DvstarIterateVectors)
(benchmark::State &state) {
  size_t number_of_entries =
      sakai.counts.size() * sakai.valid_characters.size();
  Eigen::VectorXd sakai_vector(number_of_entries);
  Eigen::VectorXd ed1a_vector(number_of_entries);

  for (auto _ : state) {
    Eigen::Index i = 0;
    pst::distances::details::iterate_included_in_both<seqan3::dna5>(
        sakai, ed1a, [&](auto &context, auto &sakai_v, auto &ed1a_v) {
          const auto background_context =
              pst::distances::details::get_background_context(context, 0);
          auto sakai_background_v = sakai.counts[background_context];
          auto ed1a_background_v = ed1a.counts[background_context];

          for (auto &char_rank : sakai.valid_characters) {
            sakai_vector(i) = pst::distances::details::dvstar::get_component(
                sakai, sakai_v, context, sakai_background_v, char_rank);

            ed1a_vector(i) = pst::distances::details::dvstar::get_component(
                ed1a, ed1a_v, context, ed1a_background_v, char_rank);
            i++;
          }
        });
  }
}

BENCHMARK_F(D2StarBenchmarks, DvstarDotProduct)
(benchmark::State &state) {

  for (auto _ : state) {
    size_t number_of_entries =
        sakai.counts.size() * sakai.valid_characters.size();
    Eigen::VectorXd sakai_vector(number_of_entries);
    Eigen::VectorXd ed1a_vector(number_of_entries);
    Eigen::Index i = 0;
    pst::distances::details::iterate_included_in_both<seqan3::dna5>(
        sakai, ed1a, [&](auto &context, auto &sakai_v, auto &ed1a_v) {
          const auto background_context =
              pst::distances::details::get_background_context(context, 0);
          auto sakai_background_v = sakai.counts[background_context];
          auto ed1a_background_v = ed1a.counts[background_context];

          for (auto &char_rank : sakai.valid_characters) {
            sakai_vector(i) = pst::distances::details::dvstar::get_component(
                sakai, sakai_v, context, sakai_background_v, char_rank);

            ed1a_vector(i) = pst::distances::details::dvstar::get_component(
                ed1a, ed1a_v, context, ed1a_background_v, char_rank);
            i++;
          }
        });
    benchmark::DoNotOptimize(sakai_vector.dot(ed1a_vector));
  }
}

BENCHMARK_F(D2StarBenchmarks, DvstarDirectDotProduct)
(benchmark::State &state) {
  for (auto _ : state) {
    double dot_product = 0.0;

    double sakai_norm = 0.0;
    double ed1a_norm = 0.0;

    pst::distances::details::iterate_included_in_both<seqan3::dna5>(
        sakai, ed1a, [&](auto &context, auto &sakai_v, auto &ed1a_v) {
          const auto background_context =
              pst::distances::details::get_background_context(context, 0);

          auto sakai_background_v = sakai.counts[background_context];
          auto ed1a_background_v = ed1a.counts[background_context];

          for (auto &char_rank : sakai.valid_characters) {
            double sakai_component_value =
                pst::distances::details::dvstar::get_component(
                    sakai, sakai_v, context, sakai_background_v, char_rank);

            double ed1a_component_value =
                pst::distances::details::dvstar::get_component(
                    ed1a, ed1a_v, context, ed1a_background_v, char_rank);

            dot_product += sakai_component_value * ed1a_component_value;
            sakai_norm += std::pow(sakai_component_value, 2.0);
            ed1a_norm += std::pow(ed1a_component_value, 2.0);
          }
        });

    sakai_norm = std::sqrt(sakai_norm);
    ed1a_norm = std::sqrt(ed1a_norm);
  }
}

BENCHMARK_F(D2StarBenchmarks, DvstarIterateInclude)
(benchmark::State &state) {
  for (auto _ : state) {
    size_t i = 0;
    pst::distances::details::iterate_included_in_both<seqan3::dna5>(
        sakai, ed1a, [&](auto &context, auto &left_v, auto &right_v) {
          for (auto &char_rank : sakai.valid_characters) {
            benchmark::DoNotOptimize(i++);
          }
        });
  }
}

BENCHMARK_F(D2StarBenchmarks, DvstarInclude)
(benchmark::State &state) {
  std::string context{"ACTAGA"};
  auto left_v = sakai.counts[context];
  for (auto _ : state) {
    size_t i = 0;
    pst::distances::details::is_included_in_both<seqan3::dna5>(sakai, ed1a, context, left_v);
  }
}

BENCHMARK_F(D2StarBenchmarks, DvstarIterateOneLookup)
(benchmark::State &state) {
  for (auto _ : state) {
    size_t i = 0;
    pst::distances::details::iterate_included_in_both<seqan3::dna5>(
        sakai, ed1a, [&](auto &context, auto &left_v, auto &right_v) {
          for (auto &char_rank : sakai.valid_characters) {
            sakai.counts[context].next_symbol_probabilities[char_rank];
            ed1a.counts[context].next_symbol_probabilities[char_rank];
            i++;
          }
        });
  }
}

BENCHMARK_F(D2StarBenchmarks, DvstarHead)
(benchmark::State &state) {
  Eigen::VectorXd vector(100000);

  for (auto _ : state) {
    benchmark::DoNotOptimize(vector.head(10000));
  }
}

BENCHMARK_DEFINE_F(D2StarBenchmarks, DvstarComponent)
(benchmark::State &state) {
  auto background_order = state.range(0);

  std::string context{"ACTAGA"};
  auto right_v = sakai.counts[context];
  bool right_included = right_v.is_included;

  const auto background_context =
      pst::distances::details::get_background_context(context,
                                                      background_order);

  for (auto _ : state) {
    benchmark::DoNotOptimize(pst::distances::details::dvstar::get_component(
        sakai, right_v, context, background_context, 0));
  }
}
BENCHMARK_REGISTER_F(D2StarBenchmarks, DvstarComponent)->DenseRange(0, 3, 1);

BENCHMARK_DEFINE_F(D2StarBenchmarks, DvstarBackgroundState)
(benchmark::State &state) {
  auto shared_contexts =
      pst::distances::details::get_shared_contexts(sakai, ed1a);
  auto background_order = state.range(0);

  auto vector = pst::distances::details::adjusted_transition_frequency_vector(
      sakai, shared_contexts, background_order);

  for (auto _ : state) {
    for (auto &context : shared_contexts) {
      benchmark::DoNotOptimize(pst::distances::details::get_background_context(
          context, background_order));
    }
  }
}
BENCHMARK_REGISTER_F(D2StarBenchmarks, DvstarBackgroundState)
    ->DenseRange(0, 3, 1);

BENCHMARK_F(D2StarBenchmarks, DvstarSingleIsIncluded)
(benchmark::State &state) {
  auto shared_contexts =
      pst::distances::details::get_shared_contexts(sakai, ed1a);

  for (auto _ : state) {
    benchmark::DoNotOptimize(sakai.is_included(shared_contexts[0]));
  }
}

static void SharedContexts(benchmark::State &state) {
  std::filesystem::path first_path{"./../trees/CM008035.1.tree"};
  std::filesystem::path second_path{"./../trees/NC_009067.tree"};
  pst::KullbackLieblerTreeMap<seqan3::dna5> first{first_path};
  pst::KullbackLieblerTreeMap<seqan3::dna5> second{second_path};

  for (auto _ : state) {
    benchmark::DoNotOptimize(
        pst::distances::details::get_shared_contexts(first, second));
  }
}
BENCHMARK(SharedContexts);

static void CompositionVector(benchmark::State &state) {
  std::filesystem::path first_path{"./../trees/CM008035.1.tree"};
  std::filesystem::path second_path{"./../trees/NC_009067.tree"};
  pst::KullbackLieblerTreeMap<seqan3::dna5> first{first_path};
  pst::KullbackLieblerTreeMap<seqan3::dna5> second{second_path};

  auto contexts = pst::distances::details::get_shared_contexts(first, second);

  for (auto _ : state) {
    benchmark::DoNotOptimize(
        pst::distances::details::composition_vector(first, contexts, 2));
  }
}
BENCHMARK(CompositionVector);

static void GetTerminalNodes(benchmark::State &state) {
  std::filesystem::path first_path{"./../trees/CM008035.1.tree"};
  pst::KullbackLieblerTreeMap<seqan3::dna5> first{first_path};

  for (auto _ : state) {
    benchmark::DoNotOptimize(first.get_terminal_nodes());
  }
}
BENCHMARK(GetTerminalNodes);

static void CountTerminalNodes(benchmark::State &state) {
  std::filesystem::path first_path{"./../trees/CM008035.1.tree"};
  pst::KullbackLieblerTreeMap<seqan3::dna5> first{first_path};

  for (auto _ : state) {
    benchmark::DoNotOptimize(first.count_terminal_nodes());
  }
}
BENCHMARK(CountTerminalNodes);

static void GetTransitionProbability(benchmark::State &state) {
  std::filesystem::path first_path{"./../trees/CM008035.1.tree"};
  pst::KullbackLieblerTreeMap<seqan3::dna5> first{first_path};

  std::string context{"GTGTGT"};
  size_t char_rank = 4;
  for (auto _ : state) {
    benchmark::DoNotOptimize(
        first.get_transition_probability(context, char_rank));
  }
}
BENCHMARK(GetTransitionProbability);

static void GetClosestState(benchmark::State &state) {
  std::filesystem::path first_path{"./../trees/CM008035.1.tree"};
  pst::KullbackLieblerTreeMap<seqan3::dna5> first{first_path};

  std::string context{"ACGGTGTGT"};
  for (auto _ : state) {
    benchmark::DoNotOptimize(first.get_closest_state(context));
  }
}
BENCHMARK(GetClosestState);

static void CVEstimation(benchmark::State &state) {
  std::filesystem::path first_path{"./../trees/CM008035.1.tree"};
  std::filesystem::path second_path{"./../trees/NC_009067.tree"};
  pst::KullbackLieblerTreeMap<seqan3::dna5> first{first_path};
  pst::KullbackLieblerTreeMap<seqan3::dna5> second{second_path};

  for (auto _ : state) {
    benchmark::DoNotOptimize(
        pst::distances::cv_estimation(first, second, 6, 2));
  }
}
BENCHMARK(CVEstimation);

static void GetAllContexts(benchmark::State &state) {
  std::filesystem::path first_path{"./../trees/CM008035.1.tree"};
  pst::KullbackLieblerTreeMap<seqan3::dna5> first{first_path};

  for (auto _ : state) {
    benchmark::DoNotOptimize(
        pst::distances::details::get_all_contexts<seqan3::dna5>(
            6, first.valid_characters));
  }
}
BENCHMARK(GetAllContexts);

static void NegativeLogLikelihood(benchmark::State &state) {
  using seqan3::operator""_dna5;
  std::filesystem::path first_path{"./../trees/CM008035.1.tree"};
  pst::KullbackLieblerTreeMap<seqan3::dna5> first{first_path};
  std::vector<seqan3::dna5> long_sequence{
      "CGCGCGGCGCGTCGTCGAGCTCGCCGCCCAGGGCGCGCCCCTCGGCGCGATCCGCGCCGCCCTCAACGACGTGACCCCCGCCGCCGACAAGGGCGAAGCCTTCGTCGGACGCCCCGACTGGCTGGGCGAGCTGTGGACCGCGCGCCGCACGGACCGGCCCCTGATCGACTCGATCACCAAAAAGGCCCTGCCCCGTGCGACCAAGGTCAAGGGCTGGCGCTGGAAGAAGCGCCCCGAGGTCGCCGACTACACGGGCAACAAAACCGAAATCTCCTCCAACGAGATCGAGACCGAGCCCGTCGAGGCCGCCGTCAAACGCATCGCCGCAGGGTGGGACACCGACCGCATTTTCGTCGACCTCGGCGACGGCGACATGATCGAGAGCCTGTGGGAGGGCGCCCGCGAGGACTACGCGATCAAGACCGAGGCCGCCGTCACCACCGGCCTCAAGACCGCAGCGACGAAGCTCACCGGCGCGCCCACCGAACTGGACAAGGCCCTCGTGGTCCTCGGCTCCAAGGCCGCCGCCATCGGCTCCCGCCTGAACTTCGTCGCCTTCGGCGCCGACGTGTGGAGCAAGTTCACCGCGCTGACCCGCGACCAGGTGCCGTGGTGGATCACCAACGGCGACCGCCTCAACCTCTCGACCGCGACCGGCGAGGTCAACGGCCTGCGCCTGTTCGTGGACCCGACCCTCGCGGCGGGCGACATCCTCGCGGGCGACACCCGGTCCGCGACGTTCTACGAGGAGCCGACCCCCATCAGGGTCAACGCCATCGACCTGCCCAAGGGCGGCGTGGACCTCGGCCTGTTCGGGTACCACGCCCTTCTCGTGAACGACCCGAACTCGTTGTTCATCATCACGGCGGGCTGACCCCATGACCCCCGACGACCTCGCCACGCGGGCCGCCGCGTGGGCGAAGCTCCCGGGCGGCGTGGACGACGCCATGAGGGCGTGCGCAGCCGCAGTGCACGCCCTCGTGGCCGCCCTGCCCGTCACGCAGGGCCGCCCCGCCTGGCGCGAGGACACGGCCCTCGGAGCGGTCATGCTCACCGCCCGCCTGCACCGCCGCCGCAACAGCCCGGCGGGCATCGAGTCCCTGACCGAGATGGGCGCGACCTACGTGAGCCGCTACGACAGCGACATCGCGCGCCTGCTGCGCATCGACGCCTTCGTCGGGCCCGTCGCCATCTGAGGGGGGCCACGAGATGAACCCGCTCTACGCGGCCGCCCAGGATGTGGCCGACATGCTCGCCGCGGCCGGAGTCCACACGGTCACCGACCCCAGGGACATCGAGCCGCCGTGCGCGTGGGTCAGCCCCAGCCGCATCGCCTACCCCACGCTCGCTGGCCGCCCCCGCACCGTCGAGTGGGAGGTGTACCTCATCGCACCCGACAGCGGCGCGCCCCTCTTCCCCCTCGGCGACCTCATCGACCGGGCCGCCACCGTCTTCCCCGGCATCGAGGCCCGCACCCTCGGCCTGACAATCCCCAACCTCAGCCCGGACCCCCTGCCCGCGATCACGTTCACCATCGAAACAGAAACGGACTAAACCCATGGCAGTGAAAACCCTCACCCTCGGCCCCGGCAAACTCAGCTTCGGCGCCCCCGAGTCCCTGACCCACGCCGCCGCCCAGGTCACCAAATGCGCCGTCAAGCCCACCGCGAAGCAGGGCGACTCCGTGGCCGTCCTGTCCGGCGACCGCGTGCCCGGCGACCGCACCGAAGCCGCGACCCTGGAGTTCACGATCTACCAGGACTTCGGCGAGGCCGAATCCTTCGTCGAATGGACCTGGGCCAACGCAGGGAAGGAACTCCCCTTCGAGTTCATCCCCGCCGACAAGCACGACAAGGCCGTGCGCGGCCGCGTCACGATCGAGCGGTCCGACATCGGCGGCGAGGTCGGCGTCAAGGTCACCGCCGACCTGGAGTTCACCTGCACCACCATGCCCACCATCGAGCCCAAAACCAAGATCGGGCACTGAGGTGGCCGACTACTCCGGCGTCAAGATCGACGGCGCGCGCCGCCTCCGCTCGACCCTGCGCAAAGCGGGCGCGGACATGCGCGACATGCGGGAGGTGAACCGCGTCGTCGCCGGCATCGTCGTCGGCGCGGCCACCGCCCGCGTCCCCCGACGCACCGGGGCCCTGGCCGCCACCGTGCGCGCAGGGGCCACCCAGGCCGCAGCCATCGGCCGCGCCGGGAACAACCGCCGCACCGGCGTCCCCTACGCCAACCCCATCCACTGGGGATGGCACCGCCACCGAATCCGCCCTAACCCGTTCCTCAGCCTCGCCGCCCAGGACACCGAACCCCAGTGGTTCGGCGTCTACGCCGACCGCATCGAACGCCTCATCAACAGCATCGAAGGAGCCTGACCCATGTCGAGCATCAAAGCCATCAACGTCGAGGTAGTCACCTCCGCCGTGACCGGCGACCTCGCCGCCGTCACCGTCCGCACCGACAACCGCGACCGCATCGCCTGGGACCTCGCGAGGGGCCGCAACAAATGGCCCCAAGCACAGGAGGCCCCCAGCCTGTGGGCCACCCACATCGCCTACACCGCGCTCCGCCGCACCGGCGAAGTCAGCTGCTCGTTCGAGGAGTTCTCCGAGGCAACTGTGAGCGCCGAACCCGAGGTCATCGACGTGGACCCTACCCGGACGGCGACCGCCGGGGCCTGATCGTCGCCCTGGCCCTCGCCACCCGCATCCCCATGAGCGAGTGGGAGACCCGCCCCGACGAGGACATCGCCACCGCACTGCAACTGCTAGAAGAGAGGAGGAGCTGACTTGGCGTCGAAAACCGCCATCCTGAGCGTCCGCGTCGTCTCCGACGTGAAGGACGCCACCAAGGGACTGGACGACGTGGCCGACAAGACCGGCCGCCTGGAGGACGGCCTCAAACGGGCCGCCGCCCCCGCCGGGATCGCCGTCGCCGCCCTCGCCGGGATCGGCAAGGCCGCCACCGACTCCGCCAGCGAGTTGCAGCAGAGCGCGGGCGCCGTCGAATCCGTATTCGGCGGGCACGCCGCCGCCGTCCAGGACGCCGCCAAGACCGCCGCCTCCAGCGTCGGCCTGGCAGCAAGCGAGTACCAGAACATGAGCGCGGTCCTGGGCGCCCAGCTCAAGAACATGGGCACCCCCATGGAGGACCTGGCCGGATCGACCCAGAACCTCATAGGCCTGGGCTCCGACCTCGCCGCCACCTTCGGGGGAACCACCGCCGACGCCGTGAGCGCCATCTCAGCCCTCCTCCGGGGCGAGCGCGACCCCATCGAGCGCTACGGCGTCTCGATCAAACAGTCGGATATCAACGCGCGTCTGGCCGCCGAGGGCATGGACAAGCTGGAAGGCGCGGCCAAGACCCAGGCCGAAGCCCAGGCCGCCCTCGCCCTGCTCACCGAGCAGACCGCATCTGCGCAAGGCCAGTTCGCGCGCGAGACCGACACGATGGCCGGGAGCCAGCAGATCGCCGCCGCCCAGTTCGAGAACGCAAAAGCCGCCCTCGGGGAGAAGCTGCTGCCCGTCGTCACGCAGTTCATGGAGGCCATGAGCGGGGCGGCTCAATGGGTCGCCCAGAACAGCGATGCGCTGCTCGTCCTCGGCGGCGTCGTCGGAACCATCGCGGGCGTGATCCTCGCCGCCAACGCCGCCATGGGCGTGTGGACCGCAGTCCAGACGACCGCCAGAGTCGCGACGGCCGCCTGGACCGGCGTCCAGGCCGCGTTCAACGCGGTCATGGCCCTGAACCCGATCACACTGGTGGTCATCGCCATCGGGGCCCTGGTCGCCGCCGTCGTCGTCGCCTACAACAAGTCCGAGGCGTTCCGCAACGCTGTATCCGCGCTGTGGGACGCCATCAAAGCGGGGGGCGGCTGGATAGTCGATCACGTCATCAAACCCATCGGAGACGCTTTCAACGCCGTTGTTGATGCCGTGAAATCCGTCTACGACTGGGTGAAGAACCTGTTCAGCGGATTCCAGCTCCCAGGCTGGCTATCGAGCGTACTCAGCTGGTTCGGGCTCGAAGCCCCCTCCGGCCCCGAATCCGGGGCCATCCTCGCCGCCACCGGGACAACGGACGCCCCCCTCGCGCGCCTCGCATCGTGGGCGCTCGCCCCCCGCACCGGCTCCAGCCCCACCCCCGCCGGCAGCGTCGTGAACATCACTGTGAACGGCGCCCTGGACCCTGACGCGGTGGCCCGCCAGATCGGCCGCATCCTGTCGCGCCGCGACCTCATCAACGGCACCGAGCAGATCGTGGGGGCGACCCTATGAGCGTCAGCGCAAGCCTCCGCGTCGCCGCAGGCGGCCTCGGCGGTGTCATCAACGCCGCCGCCGACAAGTACCCGACGACGGTGACCGTCCTCGACGACCTCACCGTCACGTGGGGCCGCGATAGCGTGGTCTCACACCCGGACCCGTCGTCAATGACCGCGACAATCGCCCTCGTGGACACTGTGCCCGACTGGCTGCGCGTCGGCGCGCTGGCGACCGTCAACGCGGTCGCCCGGACCGAGGAGTCCCAGCGGTCCTATATGCGGCTCCTGCCCTGGCGCGCCATCGAACCGGGCACCGGGTGGCGCCAGCAGGTCACCCCCGATCCGCCCGGCGCGTGGGTCGGGAGCCTCCCGGTGTTCGCCGCCGCCGGAGCAGACTCCGGGATCGGATGGTTCATAGCGCCAGGGGTGCAGCCCCCCTCCGATACCCCCGAGACCACCCAGTGGGCCGCGAACGCAAAGACGACAGCGGGCAAGCCAGTCACGTTCACAATCACCGTCCCCGAGCTCACCGGCGCGACCGTGCGCGCCTTCCCGCTCACTTACCGGCGCCCCGGCGGGCTCTACACCCGCGCCCCCGGCATGGCGATCGAGCTCTCCCCGGAGAAGTACACGCCGGGCACTGTCGAGTACTCGGGCACCTGGACCCCCGAGGCGACGGGCCTTTACGTCGGCGCCTATCTGCATATCCAGCTGCACAAGGCCCCGGCCTGGACGACGATCCCCCGCGAGCGGACCTGGCGGGCCGCACCCGGCACGTGGGCGGACGCCGGCGGCCGCGCAACCGTCACTGACGTTCACATCGCCGGGACCTCCGGCCACGTCGCTGAGCACGCCGTCGAGGTCTTCACCGGCCGGGTCCAGTCCCTCCGCGTCGAATGGTCCGAGCGCCTGTCCCGACCCATCGCGCGGATCACCGCAGTGGACAAGCTCGCCGATCTGAATGGTACCTACATCGGCGACACGCCGTGGGGCGAGGAATCCTGGAAGCTGCGCGCAGAACGCATCCTGAAACAGGCCCTCGGACCAGCAGACACACTGGAGGGCGAACCCGGCAATTGGCTGGGGACGATCCGCCCCCGAGACGTGGACCACCGCAGCGCCGGTGAGCTCATGAGGAACACACTCGCCTCGTGCGCCGCCACCGCATTCCCCGTCAGCTGGCGCAAATGGCGGGTAATCCCATTCATCTACAAGGGGAGCGATCAATCAATCACGATCCCCGGGCGCGCCATCCGCCGCGACGGCGTACAGGTCAGCACTGACGAATCCGCGAACATCTCGACAATTCAAGCCACGTATTTCGATGTGACCTACGACGGGAAAACCGGGAGGGTGAAAGACGTTATAGAGCGCACGACTACGCGGAAGAATACACCGGCAAATGAGGGACCCCCCAGGTCCATCAAAATGAAAACCGAACTATCCCGCAGTAACGAAGCGAGCGAACTCACCCGGATCATCGGGAAATACGTGAACGTGAACCAGTGGATCATCAGCGCCCTATCAGTGAAGCACGACCGGATCAGCGAGGACGCCCTTGTGCGCCTACTGTCCGCCACCGAGCGCATCGCCCAGCAGGTCGTCCTCACCGGCCTCCCACGATGGTTCCCAGCAGCGACGATGCGCGGCATCGTCATCGGCGGGTCCCTCACCATGCACCGCGGCCACTGGACCCCCACCCTCCGCATCGCAAACACACCCGACTAGAAAGAGTACCCATGCCATCGACAACCCCCCGGGGTCTTCCCTACGCAATCCCCACCGACGCCCAAGCCGCATTCCCCGACGCCGTGTCAAAACCCATCGCCGAATGGATCGAGGCGAACCTCCCGGTCATGCAAGCCGGGACCATCGCCTACCCCGCCCTCGGCTCCCAAGACCAGACGGGAGAGTACACAGTCACGTTTCCCAAACCCTTCCCAGTCACGCCCCGGATATTCATGCAAGCCGATAACCAGCGCCTCACAATCGCCGTATGGAATATCAGCCGCACCGGGTTCAAATGGATGGCCCGCAACAACAGCAACGGCAATTCGTCCTCTGGAGCGGCCTCGTGGTTCGCCGTTAGCGGCGCCACCGGACAGTAACGAAAGGAAACAAAGGAAATGACCACAGCGGTCGACGTGTTCACCGCCCGCCTCGCCTGGATGATGACCCAAGCCGACGGCGGCTACTCCCAGCCCAACCGCCTCGACGTGCGCCGCACGCGCGGCGTGTGGGACCCCGGCTTCCAGTTCGAGGGGGACTGCTCCTCCTGCGTCCTGGAGGCCGCCCACCAGGCGGGCCTGCCCACAGGCTCTGCGTCCTACACGGGCGACATGCGCGCGGGTCTGGAGGCCGTGGGATGGGCCGTCATCCCCTACGCCGCGACCGGCGGGGACCTTGACAACCTCGCCGACGGCGACGTGCTCCTATCCGAGGCCGCGAGCGGCGGCGTCGGCCATACCGGCGGCCTCATCCCCGGCGGCCTCGTCGCCGAGGCGTGGATCGACGGTCACGGAGACATCATGGGCTCCGCAGGCGGGGACGGGCCCGGCGACGACACCGGCGGGGAAACCCGCGCAGTGCCGTTCTATTCCCACCCCTACACAGTGCGGGGGCTCTGGACGCACGTCCTGCGCCCCCCAGCCCTCGACGCCGCAGACTCGCCCGCCGAACCCACCCCCACAACGAAAGGAATCCCCAATATGTTCGGAATCACCTACACGGCAAACGCCTTCGGCGGTATCACCGCCTACGTCCTCATCCACGAGTCCGCCGGTGCCGACGCCCTTGACCGCGTTCAGGCCCAGGTGTACAACAGCGTCCTTCCCAACGGCTTCACCGAGGTCCCTGAGCACCACGCCGAAATGCTCATCCGCGAGTCGTGGGTGCGCCACAACCGCATCGCCAACGCCGTCGCCGCGACCACTCGCGTAGACATCAACGAGGCCACCGCCCGCGTTCTCGCCGCCGTCAAGGAAGGAGCTGCCAAGTGAACGCCATCACGTCCCAGACCCCCGATGATCCCACGCCGCAGCCCATCTCCTGGCTTACACCCGCAGTGCGGCGTTACATCTACAACGTCACTATCGCCGCCCTCGGCGTCGCTCTCGTCTACGGCGTCGTCGATGGCCAGCACGCCGCCGCCTGGGAGGCCCTGGCCCTCGCTGTCGTGGGCCTCGCCCGCGCGCACGTCCCCGGAGACCCCCAATGAGCGACGCCTCGGCGGCCGTCGAGGTCATCGCCGCCATAGGCGGCCTCGGCGGCCTCGGGGCCGCGCTCTCAGCCGTCGCCTCCCTCATGGAGGCCCGCAGGGTCCGCGCCAGTATCCCCGCCGCCGCCGACCGCACCGAGGAGGCCATCGACGCCCTGCGCTCCGATGTCCGCGCCATCGACCGCCGCATCGGACACGAGCTCGGCGAAATCCGCCGAGCCGCCGACCGGGAACACGCCGACTATGACGCACGTCTCAGACGATTGGAGGGGTCATGAGTTGACACTCCCACTCAGGTGGTTATAAACTATAACCATCAGGAAGCCATAAGGGGCAAGCCTGAAACCTGAAGGGAGCGCGAAAATGCGCAAGTCAATCGAACTCACCTGGACCCCCGAAACCCGCGTTTGGGGCACAAGTGGGAACACCAGCGTGGCCGTCGGCACCGGAACCCTGGACGGCCGCCGCCTCGCCGTCTACGCCTTCCCCCAGTCCGATCACTGGTCGTTCTGGTCGCAAATCGAGCGTCCCGGCGGTGGTTCGACCTCGATAGAAATCAGGTCCACATTGCCCGCAGGCACCGTTCCATCGGTCCTCGGCCCGAACGGCGCCATCCGCGAAACGACGACAATCGAACTCTGACCCAACGCCCCGGCGCACACGCGCCGGGGCACTGTCACGGAAGGGCAAAGGATGAACGGCATCGAACTGCGCGCCCGCCGCGAAGCACTCGGGCTCTCACAAACCAAGTTCGCGAAAATGTGCGAAACGACTCAAGTGACCGTTTCCCGCTGGGAGAACGGCACCCGCGAACCAAGGAACGACATTGCAATACACCTGCTGATGGCAAATATCGAAGACGCCGCTATCGACCTCATCGAGGACCTGCTAGAGCTCGCCGAAGACGAAGAACTCCTAACAGCAACGCCCGACCTCCAGCTCACGGTGTACAACGACGAAGCCCGCTACGCCGCCGGGGAGCCCGTCTGGTCAAAACGCCTCCCCATGGAGACGCACCGCGTGTGCGCCGCCCGAGCCGCCGCCCTCCTCGGCGCCGAGGACGGCACACACGTCACACTCATCGAGGGCTGAGCGCCCTCACGCGATAGCCTCAACGACCTCGCGCACATCATCGTCGGGAACGAGGATGTAGCGAAGCGTCGTCGAGGGTGAAGCGTGCCCGAGCGCGCGCTGCACCGCGACGAGGTTCCTCGTGCGCGCGAACCCCGTCGAGGCGAAGGCGTGCCGGAGGGCGTGCATCGTCACGCCCTCCGGCAGCGCCCGGCCGACCAGCTTCCCCACCCACGCCGGGGAGAGGTGTCCATGGTCCGCCCCCGGGAAGAAGAACCCCGGGTCATGATCGAGCAGCTCATCGGCGAGAGAATGCGGGAGGGGGATCACGCGGGTTTTCCCGCCCTTTCCGTGGACCACGAGGGACCAGCCCGCCAGATCGCGCACAAGATCGCGCGTGTGCGCGCGGGCGACCTCACCCCGCCGCATCCCCAACTCCGCCGCCATGCGCACCATGAGACGCACTCGTGGATCCGTCGCCCGGCGCCCCACCGCGATGGCGCCCGGCGTCGCCGGCCTCGGCGCAGGGTCCGACTGTCTCACCGACGGCACCGGCGGCGCTACCTCGATGTATCCGACCCCCTGGGCCCACCGGTAGAACTGGTCGACGCTTTGATGCGCGCTCCGACGCGTATCCCGCGCCCAATCATGCGCCCCGGACCACTCGATCACCGTGAGCGGCCCAACCTCCCACGGGCCCGCCCGCAGATCGCGGGCGAACCTACTCACCCACTCGATCCGCAGTCGGATAGTCTCGGCCCGCCGGCCGGCCGCCGCCAGCGCTGTAGTCCACTCGCCTATAGGACCGGCCCATCCGGCGGGTACCGGTCGCGGTTTCATGCTCATGATGATTACCCTGCATCTATTCGGCCCCAGTGTCGCGTCATCACGCCGCCGGGACCCACGATGTAGGATCGAGCCGTGGGCCACCAATCCGTAGGTTGGGGGTTCGAGTCCCCCTGGGCCTACTCGCCCGCCCCAGCCCCGCCGGCTGGGGCGGTTGCCGCACTGCGCCGGTACCGACGCGCTGCCGACGGATGGAGGTCGCCATGTCCGACGACGCGGACAGGGCCCACGGCGCCCTGGCGGGACTCGCGCTGGGGGACGCCCTTGGCATGCCGACCCAGGCGATGACCGCCGATCAGATCAGGCTGACCTACGGGTGGGTGGACGCCCTGGTGCCAGCCGACGCCTCGCAGCCCTACGCGCCCGGCATGCCCGCCGGCAGCGTCACGGACGACACGGAGCAGGCGCTGCTCGTCGCCGGCCTGCTGGTATCGGGCGGGGGCGGCATCGACCCCCACGCCTTCTCCCGCGCCCTGCTGGACTGGGAGGACTCGATGGCGGCCCGCGGTTCCCTCGACCTGCTGGGCCCATCGACGAAGGCCGCCCTGGAACGCGTGCGGGCCGGGGAGGACCCCCTCCGCGTGGGCGGCGCGGGCACCACCAACGGCTCAGCGATGCGAGTCGCGCCCGTCGGGATCGCCTCCTCCACCCGGGATCCGCGTTTCGCCGACACCGTGTGGGAGTCGTGCCGCGTCACCCACGCCACCGAACAGGGCTTCCACGCCGCCGCGCTCGTGGCGGCAGCGGTCTCCCTCGGCATCGACGGAGCAGGGGCGGACAGCCCTTCGGACTCCGCCCGCGCCTCCTTGGAACGCGCCCTGGCCCTCGTGGAGGCACTCGGGCGCCGGGGGGCGCGGACGCCCCAGCCGGACGTGTGCGAGCGGACCCGCTACGCGCTGCGGTTCGCGCGCGCCCGCGACCCCGCCCCCGGTACTGCCGACGACGACCGGGCATTCGCCGGGGCACTGCGGGCACGCGTCGGCGCCTCCGTGGAGGCCGCCCAGTCCGTCCCCGCAGCATTCGCCATCGCCTGGCGCTACGCCGCCGATCCGTGGCGGGGCCTGTGCGTCGCCGCCAACCTCGGCGGTGACACCGACACGATCGGCGCTATCACCGGCGCCGTGCTCGGCGCCGCCCTGGGGGCCCGGTGCTGGCCCGCCCAGGAGCTGGAACGAGTGGAGGCCGTCTCCGGGCTGCGGCTGCGCGAGACCGCCGACGGTTTGCTCCGCCTGCGCGCCCACGGATCCCGACTGCCCGCCCACGGGGAGCCGGTCGCAGCACCGCAGGAGGGCAGGGTCGTCCTGCTCGGGCAGGTGGTCGTCGACCTCGCACTGCTGGCGCCGCGCGTGCCCGCTCCCGGCGGCGACGTGTTCGCAGAGGACGCGGGCATGCACGCGGGCGGGGGCTTCAACGTGCTGGCCGCTGCGCGCCGGATGGGAGCGGAGGCAGTGAGCCTGTCCGGCGTCGGGGACGGCGGATTCGCCTCGATCATCACCGCTGCGTTGGAGCGCATCGGCGCCTCCTGCGAGGGACCGCGCGTCGCGGGAACGGACTCGGGGTACTGCGTGGCCATCACGGACGGCGACGGCGAGCGCACCTTCGTCTCGACCAGGGGCGCGGAGGCCCGCCTGCCGCGCGGGTCGTGGTCCGCCCACGCGGCCCGCTTGCGCAGCGGGGACGTGGTGCACGTGGACGGCTACGCGCTGGCCCATCCGGCCAACACCGCAGCGCTGCGGGAGTTCCTCTCGGCGCACCTGCCCGCAGGGCTCCGCGCGATCGTCGACGTGTCGCCCGTCGTCGGCGATGTGGACCTCGACGACCTGCTTGCCCTGCGGGCCCTGGCCCCCCTGTGGTCCATGAACGAGCGCGAGGCGGGGATCCTCGCGGGCCGCCTCGCGCGGGCGTCCGCCGCTCCCCCGCACGGAGGCGCTCCCCCGGGGGAGGCGACACCACCGGCCGGAGCGGCCCCCGGGA"_dna5};

  for (auto _ : state) {
    benchmark::DoNotOptimize(
        pst::distances::negative_log_likelihood_symmetric<seqan3::dna5>(
            first, long_sequence));
  }
}

BENCHMARK(NegativeLogLikelihood);

static void HashMapScore(benchmark::State &state) {
  using seqan3::operator""_dna5;
  std::vector<seqan3::dna5> sequence = random_sequence(50000);

  pst::KullbackLieblerTreeMap<seqan3::dna5> tree{"HashMap", sequence, 7, 2,
                                                 3.9075,    true,     2};

  tree.construct_tree();

  for (auto _ : state) {
    benchmark::DoNotOptimize(
        pst::distances::negative_log_likelihood<seqan3::dna5>(tree, sequence));
  }
}

BENCHMARK(HashMapScore);

//static void TreeScore(benchmark::State &state) {
//  using seqan3::operator""_dna5;
//  std::vector<seqan3::dna5> sequence = random_sequence(50000);
//
//  pst::KullbackLieblerTree<seqan3::dna5> tree{"Tree", sequence, 7, 2,
//                                              3.9075, false,    2};
//  tree.construct_tree();
//
//  for (auto _ : state) {
//    benchmark::DoNotOptimize(
//        pst::distances::negative_log_likelihood<seqan3::dna5>(tree, sequence));
//  }
//}
//
//BENCHMARK(TreeScore);

static void KL_AllContexts(benchmark::State &state) {
  std::filesystem::path first_path{"./../trees/CM008035.1.tree"};
  pst::KullbackLieblerTreeMap<seqan3::dna5> first{first_path};

  for (auto _ : state) {
    benchmark::DoNotOptimize(
        pst::distances::details::get_all_contexts<seqan3::dna5>(
            6, first.valid_characters));
  }
}
BENCHMARK(KL_AllContexts);

static void KL_LikelihoodAllContexts(benchmark::State &state) {
  std::filesystem::path first_path{"./../trees/CM008035.1.tree"};
  pst::KullbackLieblerTreeMap<seqan3::dna5> first{first_path};

  auto all_contexts = pst::distances::details::get_all_contexts<seqan3::dna5>(
      6, first.valid_characters);

  for (auto _ : state) {
    for (auto &context : all_contexts) {
      benchmark::DoNotOptimize(
          pst::distances::details::likelihood_context(first, context));
    }
  }
}
BENCHMARK(KL_LikelihoodAllContexts);

static void KL_LikelihoodOneContext(benchmark::State &state) {
  std::filesystem::path first_path{"./../trees/CM008035.1.tree"};
  pst::KullbackLieblerTreeMap<seqan3::dna5> first{first_path};

  auto all_contexts = pst::distances::details::get_all_contexts<seqan3::dna5>(
      6, first.valid_characters);

  for (auto _ : state) {
    benchmark::DoNotOptimize(
        pst::distances::details::likelihood_context(first, all_contexts[0]));
  }
}
BENCHMARK(KL_LikelihoodOneContext);

static void KL_(benchmark::State &state) {
  std::filesystem::path first_path{"./../trees/CM008035.1.tree"};
  std::filesystem::path second_path{"./../trees/NC_009067.tree"};
  pst::KullbackLieblerTreeMap<seqan3::dna5> first{first_path};
  pst::KullbackLieblerTreeMap<seqan3::dna5> second{second_path};

  for (auto _ : state) {
    benchmark::DoNotOptimize(
        pst::distances::kl_divergence<seqan3::dna5>(first, second, 6));
  }
}
BENCHMARK(KL_);

BENCHMARK_MAIN();