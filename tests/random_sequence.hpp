#include <random>
#include <vector>

#include <seqan3/alphabet/nucleotide/dna5.hpp>
using seqan3::operator""_dna5;

std::random_device rd;
std::mt19937 gen = std::mt19937{rd()};

std::vector<seqan3::dna5> random_sequence(size_t length) {
  std::vector<seqan3::dna5> local_sequence{"ACGT"_dna5};
  std::uniform_int_distribution<> distrib(0, 3);
  for (size_t i = 0; i < length; i++) {
    auto rand = distrib(gen);
    auto c_dna4 = seqan3::assign_rank_to(rand, seqan3::dna4{});
    auto c_dna5 = seqan3::assign_char_to(c_dna4.to_char(), seqan3::dna5{});
    local_sequence.push_back(c_dna5);
  }

  return local_sequence;
}