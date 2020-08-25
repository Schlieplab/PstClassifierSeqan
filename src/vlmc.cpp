#include <cstring>
#include <string>
#include <vector>

#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/range/views/all.hpp>
#include <seqan3/std/ranges>

#include "kl_tree.hpp"
#include "probabilistic_suffix_tree.hpp"
#include "ps_tree.hpp"

lst::details::sequence_t<seqan3::dna5> to_dna(const char *sequence_) {
  std::string seq = std::string{sequence_};
  lst::details::sequence_t<seqan3::dna5> seq_ =
      seq | seqan3::views::char_to<seqan3::dna5> |
      seqan3::views::to<lst::details::sequence_t<seqan3::dna5>>;
  return seq_;
}

extern "C" {
static const char *train_kl(const char *id_, const char *sequence_,
                            size_t max_depth, size_t min_count, float threshold,
                            bool multi_core, int parallel_depth) {

  auto sequence = to_dna(sequence_);
  std::string id{id_};

  pst::KullbackLieblerTree<seqan3::dna5> pst{
      id,        sequence,   max_depth,     min_count,
      threshold, multi_core, parallel_depth};
  pst.construct_tree();
  return strdup(pst.to_tree().c_str());
}

static const char *train_ps(const char *id_, const char *sequence_,
                            size_t max_depth, size_t min_count, bool multi_core,
                            int parallel_depth) {

  std::string id = std::string(id_);
  auto sequence = to_dna(sequence_);

  pst::PeresShieldsTree<seqan3::dna5> pst{
      id, sequence, max_depth, min_count, multi_core, parallel_depth};
  pst.construct_tree();
  return strdup(pst.to_tree().c_str());
}

static const char *train_ps_parameters(const char *id_, const char *sequence_,
                                       size_t max_depth, size_t min_count,
                                       size_t n_parameters, bool multi_core,
                                       int parallel_depth) {

  std::string id = std::string(id_);
  auto sequence = to_dna(sequence_);

  pst::PeresShieldsTree<seqan3::dna5> pst{
      id,           sequence,     max_depth,  min_count,
      n_parameters, "parameters", multi_core, parallel_depth};

  pst.construct_tree();
  return strdup(pst.to_tree().c_str());
}

static const char *train_kl_parameters(const char *id_, const char *sequence_,
                                       size_t max_depth, size_t min_count,
                                       size_t n_parameters, bool multi_core,
                                       int parallel_depth) {

  std::string id = std::string(id_);
  auto sequence = to_dna(sequence_);

  pst::KullbackLieblerTree<seqan3::dna5> pst{
      id,           sequence,   max_depth,     min_count,
      n_parameters, multi_core, parallel_depth};

  pst.construct_tree();
  return strdup(pst.to_tree().c_str());
}
}
