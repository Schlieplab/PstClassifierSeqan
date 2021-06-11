#include <cstring>
#include <string>
#include <vector>

#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/alphabet/views/char_to.hpp>
#include <seqan3/std/ranges>
#include <seqan3/utility/views/to.hpp>

#include "kl_tree_map.hpp"
#include "probabilistic_suffix_tree_map.hpp"

lst::details::sequence_t<seqan3::dna5> to_dna(const char *sequence_) {
  std::string seq{sequence_};
  auto sequence_as_dna =
      seq | seqan3::views::char_to<seqan3::dna5> |
      seqan3::views::to<lst::details::sequence_t<seqan3::dna5>>;

  return sequence_as_dna;
}

std::string train_kl_cpp(const std::string &id_, const std::string &sequence_,
                         const size_t max_depth, const size_t min_count,
                         const float threshold, const bool multi_core,
                         const int parallel_depth) {
  auto sequence_as_dna =
      sequence_ | seqan3::views::char_to<seqan3::dna5> |
      seqan3::views::to<lst::details::sequence_t<seqan3::dna5>>;

  pst::KullbackLieblerTreeMap<seqan3::dna5> pst{
      id_,       sequence_as_dna, max_depth,     min_count,
      threshold, multi_core,      parallel_depth};
  pst.construct_tree();

  return pst.to_tree();
}

extern "C" {
static const char *train_kl(const char *id_, const char *sequence_,
                            size_t max_depth, size_t min_count, float threshold,
                            bool multi_core, int parallel_depth) {

  auto sequence = to_dna(sequence_);
  std::string id{id_};

  pst::KullbackLieblerTreeMap<seqan3::dna5> pst{
      id,        sequence,   max_depth,     min_count,
      threshold, multi_core, parallel_depth};
  pst.construct_tree();
  return strdup(pst.to_tree().c_str());
}

static const char *train_kl_parameters(const char *id_, const char *sequence_,
                                       size_t max_depth, size_t min_count,
                                       size_t n_parameters, bool multi_core,
                                       int parallel_depth) {

  std::string id = std::string(id_);
  auto sequence = to_dna(sequence_);

  pst::KullbackLieblerTreeMap<seqan3::dna5> pst{
      id,           sequence,   max_depth,     min_count,
      n_parameters, multi_core, parallel_depth};

  pst.construct_tree();
  return strdup(pst.to_tree().c_str());
}
}
