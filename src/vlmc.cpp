#include <cstring>
#include <vector>

#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/range/container/bitcompressed_vector.hpp>
#include <seqan3/range/views/to.hpp>

#include "kl_tree.hpp"
#include "probabilistic_suffix_tree.hpp"
#include "ps_tree.hpp"
#include "pst-classifier.cpp"

extern "C" {
  const char *train_(const char *id_, const char *sequence_,
                    size_t max_depth, size_t min_count,
                    float threshold, size_t number_of_parameters,
                    const char *pruning_method_, const char *estimator_,
                    bool multi_core, size_t split_depth) {


    std::string id = std::string(id_);
    std::string seq = std::string(sequence_);
    std::string pruning_method = std::string(pruning_method_);
    std::string estimator = std::string(estimator_);
    seqan3::bitcompressed_vector<seqan3::dna5> sequence =
        seq | seqan3::views::char_to<seqan3::dna5> |
        seqan3::views::to<seqan3::bitcompressed_vector<seqan3::dna5>>;

    /*
    pst::KullbackLieblerTree<seqan3::dna5> pst{id, sequence, max_depth, min_count,
                                               threshold};
    pst.construct_tree();
     */
    std::string tree = train(sequence, id, max_depth, min_count,
                             threshold, number_of_parameters,
                             pruning_method, estimator,
                             multi_core, split_depth);

    return strdup(tree.c_str());
  }
}
