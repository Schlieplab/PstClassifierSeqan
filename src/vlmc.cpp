#include <cstring>
#include <vector>

#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/range/view/char_to.hpp>

#include "probabilistic_suffix_tree.hpp"


extern "C" {
const char *train_kl(const char *id_, const char *sequence_, size_t max_depth,
                  size_t min_count, float threshold) {

  std::string id = std::string(id_);
  std::string seq = std::string(sequence_);
  seqan3::dna5_vector sequence = seq | seqan3::view::char_to<seqan3::dna5>;

  pst::ProbabilisticSuffixTree<seqan3::dna5> pst{id, sequence, max_depth,
                                                 min_count, threshold};

  return strdup(pst.to_tree().c_str());
}

const char *train_ps(const char *id_, const char *sequence_, size_t max_depth,
                  size_t min_count) {

  std::string id = std::string(id_);
  std::string seq = std::string(sequence_);
  seqan3::dna5_vector sequence = seq | seqan3::view::char_to<seqan3::dna5>;

  pst::ProbabilisticSuffixTree<seqan3::dna5> pst{id, sequence, max_depth,
                                                 min_count};

  return strdup(pst.to_tree().c_str());
}

const char *train_ps_parameters(const char *id_, const char *sequence_, size_t max_depth, size_t min_count, size_t n_parameters) {

  std::string id = std::string(id_);
  std::string seq = std::string(sequence_);
  seqan3::dna5_vector sequence = seq | seqan3::view::char_to<seqan3::dna5>;

  pst::ProbabilisticSuffixTree<seqan3::dna5> pst{id, sequence, max_depth,
                                                 min_count, 0.0, n_parameters, "parameters", "PS"};

  return strdup(pst.to_tree().c_str());
}

const char *train_kl_parameters(const char *id_, const char *sequence_, size_t max_depth, size_t min_count, size_t n_parameters) {

  std::string id = std::string(id_);
  std::string seq = std::string(sequence_);
  seqan3::dna5_vector sequence = seq | seqan3::view::char_to<seqan3::dna5>;

  pst::ProbabilisticSuffixTree<seqan3::dna5> pst{id, sequence, max_depth,
                                                 min_count, 0.0, n_parameters, "parameters", "KL"};

  return strdup(pst.to_tree().c_str());
}
}