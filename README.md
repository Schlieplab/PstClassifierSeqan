# VLMC construction using a Lazy Suffix Tree

Implements the VLMC(variable length Markov chain) algorithm described by [Schulz et al.](https://doi.org/10.1007/978-3-540-87361-7_26). In short, the VLMC of a sequence _S_ is constructed by building a probabilistic suffix tree (PST) in two phases: support pruning and similarity pruning. This PST has a direct correspondence to the VLMC.

In support pruning, we select only those branches which are at most _d_ deep (_k_-mer is at most _d_ long), and occur at least _c_ times in the sequence. This is implemented in practice using a Lazy Suffix Tree (using the WOTD algorithm described by [Giegerich et al.](https://doi.org/10.1002/spe.535)). The suffix tree is extended with implicit nodes and suffix links. See [search/lazy_suffix_tree.hpp](src/search/lazy_suffix_tree.hpp) and [search/lazy_suffix_tree/](src/search/lazy_suffix_tree/) for the implementation details.

The similarity pruning proceeds to compute the forward probabilities (probability of a character _a_ occurring after a context _c_), and then prunes the tree bottom-up. The pruning uses the Kullback-Leibler estimator for calculation of which contexts should be kept in the tree.

The output is as a `.tree` file which for, each node, contains the forward and reverse counts of each child as well as the index of the (PST) children of that node.

## Usage

Example:

```cpp
#include "probabilistic_suffix_tree.hpp"
#include "distances/negative_log_likelihood.hpp"
#include <seqan3/alphabet/nucleotide/dna5.hpp>

int main() {
    using seqan3::operator""_dna5;
    seqan3::dna5_vector sequence = "GATTACA"_dna5;

    // Parameters: name, sequence, min_count, max_depth, Kullback-leibler threshold, parallel, parallel-depth
    pst::KullbackLieblerTreeMap<seqan3::dna5> vlmc{"Test", sequence, 2, 3, 3.9075, true, 2};

    std::cout << vlmc.to_tree() << std::endl;

    auto score = pst::distances::negative_log_likelihood(vlmc, sequence);

    std::cout << score << std::endl;

    return 0;
}
```


This will train a PST on `sequence` with a depth of `2`, and that occur at least `3` times.

We also provide a cli, which can be used as follows:

```cpp
./pst-classifier CP007136.1.fa --min-count 100 --max-depth 15 --estimator "KL" --theshold 3.9075 --multi-core --parallel-depth 2
```

Where the first argument is a path to a fasta file.  For details on the arguments, run `./pst-classifier --help`.

## Build/Install

Building the executable and shared library should be fairly straight-forward.  Requires `cmake`, `make`, `hdf5`, and a `c++17` (or later) compatible c++ compiler (e.g. `gcc`).

Update/download the SeqAn3 library (you can skip `--init` if you're only updating), and other dependencies:

```shell script
git submodule update --init --recursive
```

Create and go to a build directory:

```shell script
mkdir build
cd build
```

Configure the `cmake`-project and build.

```shell script
cmake -DCMAKE_BUILD_TYPE=Release ..
make
```

This should yield three files in the `build/src` directory:

* `pst-classifier` which provides a command line interface to the PST/VLMC training.
* `batch-training` which trains a VLMC for every sequence in a multi-fasta file, and saves the output to a h5 file.
* `score-sequences` which takes the signatures from the output of `batch-training` and computes the negative log-likelihood of as set of sequences, and saves the result to a (possibly) different h5 file.

## Python

We include a wrapper for the VLMC training and scoring in a package called `libvlmc`.
It can be built from the python-package directory, and requires `numpy`, `cython`, and `scikit-build`.

After ensuring the c++ package builds, the python package can be installed through:

```shell script
python -m pip install -r dev-requirements.txt
python setup.py bdist_wheel
python install --user dist/libvlmc-0.2-{Python_Version}-{Python_Version}-linux_x86_64.whl
```
