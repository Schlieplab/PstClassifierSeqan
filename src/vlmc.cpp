#include <string>
#include <vector>

#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/argument_parser/all.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/io/sequence_file/input.hpp>

#include "find_contexts.cpp"
#include "probabilistic_suffix_tree.hpp"

using namespace vlmc;

struct input_arguments {
  size_t max_depth{15};
  size_t min_count{10};
  std::vector<seqan3::dna5_vector> sequences{};
};

input_arguments parse_cli_arguments(int argc, char *argv[]) {
  std::string filename{};

  input_arguments arguments{};

  seqan3::argument_parser parser{"Count kmers in the fasta file.", argc, argv};

  parser.add_positional_option(filename, "path to fasta file.");

  parser.add_option(arguments.max_depth, 'd', "max-depth",
                    "max depth of the VLMC tree.");
  parser.add_option(arguments.min_count, 'c', "min-count",
                    "minimum number of time each context has to appear in the "
                    "string to be considered.");

  try {
    parser.parse();
  } catch (seqan3::parser_invalid_argument const &ext) {
    seqan3::debug_stream << "[PARSER ERROR] " << ext.what() << '\n';
    return arguments;
  }
  seqan3::debug_stream << "The text was: " << filename << "\n";

  seqan3::sequence_file_input file_in{filename};

  for (auto &[seq, id, qual] : file_in) {
    arguments.sequences.push_back(seq);
  }

  return arguments;
}

int main(int argc, char *argv[]) {
  input_arguments arguments = parse_cli_arguments(argc, argv);

  seqan3::debug_stream << "Building index" << std::endl;

  // seqan3::fm_index index{sequences[0]};
  // seqan3::bi_fm_index bi_index{arguments.sequences[0]};
  pst::ProbabilisticSuffixTree<seqan3::dna5> pst{
      arguments.sequences[0], arguments.max_depth, arguments.min_count, 1.2};
  pst.print();

  return 0;
}
