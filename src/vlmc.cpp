#include <cstring>
#include <vector>

#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/range/container/bitcompressed_vector.hpp>
#include <seqan3/range/views/to.hpp>

#include "kl_tree.hpp"
#include "probabilistic_suffix_tree.hpp"
#include "ps_tree.hpp"

extern "C" {
static const char *train_kl(const char *id_, const char *sequence_,
                            size_t max_depth, size_t min_count, float threshold,
                            bool multi_core, int parallel_depth) {

  std::string id{id_};
  std::string seq = std::string(sequence_);
  seqan3::bitcompressed_vector<seqan3::dna5> sequence =
      seq | seqan3::views::char_to<seqan3::dna5> |
      seqan3::views::to<seqan3::bitcompressed_vector<seqan3::dna5>>;

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
  std::string seq = std::string(sequence_);
  seqan3::bitcompressed_vector<seqan3::dna5> sequence =
      seq | seqan3::views::char_to<seqan3::dna5> |
      seqan3::views::to<seqan3::bitcompressed_vector<seqan3::dna5>>;

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
  std::string seq = std::string(sequence_);
  seqan3::bitcompressed_vector<seqan3::dna5> sequence =
      seq | seqan3::views::char_to<seqan3::dna5> |
      seqan3::views::to<seqan3::bitcompressed_vector<seqan3::dna5>>;

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
  std::string seq = std::string(sequence_);
  seqan3::bitcompressed_vector<seqan3::dna5> sequence =
      seq | seqan3::views::char_to<seqan3::dna5> |
      seqan3::views::to<seqan3::bitcompressed_vector<seqan3::dna5>>;

  pst::KullbackLieblerTree<seqan3::dna5> pst{
      id,           sequence,   max_depth,     min_count,
      n_parameters, multi_core, parallel_depth};

  pst.construct_tree();
  return strdup(pst.to_tree().c_str());
}
}
