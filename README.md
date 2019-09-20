# VLMC construction using a Lazy Suffix Tree

Implements the VLMC(variable length Markov chain) algorithm described by [Schulz et al.](https://doi.org/10.1007/978-3-540-87361-7_26). In short, the VLMC of a sequence _S_ is constructed by building a probabilistic suffix tree (PST) in two phases: support pruning and similarity pruning. This PST has a direct correspondence to the VLMC.

In support pruning, we select only those branches which are at most _d_ deep (context is at most _d_ long), and occur at least _c_ times in the sequence. This is implemented in practice using a Lazy Suffix Tree (using the WOTD algorithm described by [Giegerich et al.](https://doi.org/10.1002/spe.535)). The suffix tree is extended with implicit nodes and suffix links. See [search/lazy_suffix_tree.hpp](src/search/lazy_suffix_tree.hpp) and [search/lazy_suffix_tree/](src/search/lazy_suffix_tree/) for the implementation details.

The similarity pruning proceeds to compute the forward probabilities (probability of a character _a_ occurring after a context _c_), and then prunes the tree bottom-up. The pruning uses either the [KL or PS](https://doi.org/10.2202/1544-6115.1214) estimator for calculation of which contexts should be kept in the tree.

The output is as a `.tree` file which for, each node, contains the forward and reverse counts of each child as well as the index of the (PST) children of that node.

## Usage

Example:

```cpp
#include "probabilistic_suffix_tree.hpp"
#include <seqan3/alphabet/nucleotide/dna5.hpp>

int main() {
    using seqan3::operator""_dna5;
    seqan3::dna5_vector sequence = "GATTACA"_dna5;

    pst::ProbabilisticSuffixTree<seqan3::dna5> vlmc{"Test", sequence, 2, 3};

    std::cout << vlmc.to_tree() << std::endl;

    return 0;
}
```

This will train a PST on `sequence` with a depth of `2`, and that occur at least `3` times.  The tree will then be pruned using the default parameters.

We also provide a cli, which can be used as follows:

```cpp
./pst-classifier CP007136.1.fa --min-count 100 --max-depth 15 --estimator "KL" --theshold 1.2
```

Where the first argument is a path to a fasta file.  For details on the arguments, run `./pst-classifier --help`.

## Build/Install

Building the executable and shared library should be fairly straight-forward (not sure about windows though).  Requires `cmake`, `make`, and a `c++17` (or later) compatible c++ compiler (e.g. `gcc`).

Create and go to a build directory:

```shell script
mkdir build
cd build
```

Configure the `cmake`-project and build.

```shell script
cmake -DCMAKE_BUILD_TYPE=Release ../src
make
```

This should yield two files in the `build` directory: `pst-classifier` which provides a command line interface to the PST/VLMC training and `libvlmc.so` which is a shared library, which can be used in a different project.  See e.g. [minimal_python_driver.py](minimal_python_driver.py).

## Python

See [minimal_python_driver.py](minimal_python_driver.py) for a minimal python interface. Requires the shared library to be built as specified in Build/Install.
