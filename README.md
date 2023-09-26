# VLMC construction using a Lazy Suffix Tree

Implements two VLMC (variable length Markov chain) construction algorithms, based on the ideas
by [Schulz et al.](https://doi.org/10.1007/978-3-540-87361-7_26). In short, the VLMC of a sequence _S_ is constructed by
building a probabilistic suffix tree (PST) in two phases: support pruning and similarity pruning. This PST has a direct
correspondence to the VLMC.

In support pruning, we select only those branches which are at most _d_ deep (_k_-mer is at most _d_ long), and occur at
least _c_ times in the sequence. This is implemented in practice using a Lazy Suffix Tree (using the WOTD algorithm
described by [Giegerich et al.](https://doi.org/10.1002/spe.535)). The suffix tree is extended with implicit nodes and
suffix links. See [search/lazy_suffix_tree.hpp](include/pst/search/lazy_suffix_tree.hpp)
and [search/lazy_suffix_tree/](include/pst/search/lazy_suffix_tree/) for the implementation details.

The similarity pruning proceeds to compute the forward probabilities (probability of a character _a_ occurring after a
context _c_), and then prunes the tree bottom-up. The pruning uses the Kullback-Leibler estimator for calculation of
which contexts should be kept in the tree.

The output is as a `.tree` file which for, each node, contains the forward and reverse counts of each child as well as
the index of the (PST) children of that node.

The algorithms are highly parallel, and at time of implementation, much faster than other current VLMC construction
algorithms, depending on parameter selection.

## Citation

Published in [BMC bioinformatics](https://doi.org/10.1186/s12859-021-04387-y).

## Usage

Pre-complied binaries for linux are available on GitHub under the [releases](https://github.com/Schlieplab/PstClassifierSeqan/releases).

The CLI can be used as follows. To train a VLMC on the file `CP007136.1.fa`, and compute the length-normalised negative log-likelihood:

```sh
./pst-classifier CP007136.1.fa --min-count 100 --max-depth 15 --threshold 3.9075 --multi-core --parallel-depth 2 > CP007136.1.tree

./pst-score-sequences -p CP007136.1.tree -s CP007136.1.fa
```

Where the first argument is a path to a fasta file. For details on the arguments, run e.g. `./pst-classifier --help`.

Equivalently, if `CP007136.1.fa` was a multi-fasta file,

```sh
./pst-batch-training CP007136.1.fa --min-count 100 --max-depth 15  --threshold 3.9075 --multi-core --parallel-depth 2 -o trees.h5

./pst-score-sequences -p trees.h5 -s CP007136.1.fa
```

To use the VLMC in a c++ project, include the corresponding header files (in include/), make sure SeqAn3 is installed. Example usage:

```cpp
#include "pst/probabilistic_suffix_tree.hpp"
#include "pst/distances/negative_log_likelihood.hpp"
#include <seqan3/alphabet/nucleotide/dna5.hpp>

int main() {
    using seqan3::operator""_dna5;
    seqan3::dna5_vector sequence = "GATTACA"_dna5;

    // Parameters: name, sequence, min_count, max_depth, Kullback-Leibler threshold, parallel, parallel-depth
    pst::KullbackLieblerTreeMap<seqan3::dna5> vlmc{"Test", sequence, 2, 3, 3.9075, true, 2};

    std::cout << vlmc.to_tree() << std::endl;

    auto score = pst::distances::negative_log_likelihood(vlmc, sequence);

    std::cout << score << std::endl;

    return 0;
}
```

This will train a PST on `sequence` with a depth of `3`, and include `k`-mers that occur at least `2` times, and then compute the length-normalised negative log-likelihood.

## Installation

Since we rely heavily on the SeqAn3 library, the main building dependencies are similar to theirs. Specifically, a
modern c++ compiler (full c++20 support, or GCC >= 7) is necessary.
See [the SeqAn3 quick-setup](https://docs.seqan.de/seqan/3-master-user/setup.html)
for the SeqAn3 specific dependencies and build instructions.

We also provide [singularity](https://sylabs.io/singularity/) and [docker](https://www.docker.com/) files, for easy
reproducibility and usage. The [`build_artifacts.sh`](build_artifacts.sh) script uses docker to build the binaries on debian 18.04 and outputs them in the `artifacts` folder.

### Linux

Dependencies (debian packages in parentheses):

- GCC >= 7 (g++) (__I've recently only been able to compile with GCC12__)
- CMake (`cmake`)
- HDF5 (`libhdf5-dev`)
- Eigen (`libeigen3-dev`)
- tbb (`libtbb-dev`)

Fetch the source code via git, or download the latest release.

Update/download the git submodules containing the `SeqAn3`, `eigen3`, `robin-hood-hashing`, `HighFive` libraries
(you can skip `--init` if you're only updating):

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
make pst-classifier pst-batch-training pst-score-sequences
```

This should yield three files in the `build/src` directory:

- `pst-classifier` which provides a command line interface to the PST/VLMC training.
- `pst-batch-training` which trains a VLMC for every sequence in a multi-fasta file, and saves the output to a hdf5
  file.
- `pst-score-sequences` which takes the signatures from the output of `pst-batch-training` or `pst-classifier` and
  computes the negative log-likelihood of as set of sequences, and saves the result to a hdf5 file.

Or for tests, build all everything:

```shell script
cmake -DCMAKE_BUILD_TYPE=Release ..
make
./tests/probabilistic_suffix_tree_map_tests
```

### macOS

Dependencies (homebrew packages in parentheses):

- GCC >= 7 (gcc) (__I've recently only been able to compile with GCC12__)
- CMake (cmake)
- HDF5 (hdf5)
- Eigen (eigen)
- tbb (tbb)

Other than the different packages, the instructions for Linux should work, but you may need to provide the paths to the installed gcc and g++ when you configure the cmake project, e.g.:
```shell
CC=/opt/homebrew/bin/gcc-12 CXX=/opt/homebrew/bin/g++-12 cmake -DCMAKE_BUILD_TYPE=Release ..
```

### Windows

Windows is not supported by SeqAn3, but through
the [Windows Subsystem for Linux](https://docs.microsoft.com/en-us/windows/wsl/about), the instructions for Linux should
work.

## Python

We include a wrapper for the VLMC training and scoring in a package called `libvlmc`. It can be built from the
python-package directory, and requires `numpy`, `cython`, and `scikit-build`.

After ensuring the c++ package builds, the python package can be installed through:

```shell script
python -m pip install -r dev-requirements.txt
python setup.py bdist_wheel
pip install --user dist/libvlmc-0.2-{Python_Version}-{Python_Version}-linux_x86_64.whl
```

Running the shell script [`python-cross-compile.sh`](python-cross-compile.sh) will use docker to build a version of the python package that runs on debian.
