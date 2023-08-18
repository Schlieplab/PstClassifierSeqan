# distutils: language = c++

from libcpp cimport bool
from libcpp.vector cimport vector
from libcpp.string cimport string

cdef extern from "../../src/vlmc.cpp":
    const char *train_kl_parameters(const char *id_, const char *sequence_,
                                    size_t max_depth, size_t min_count,
                                    size_t n_parameters, bool multi_core, int parallel_depth)

    const char *train_kl(const char *id_, const char *sequence_, size_t max_depth,
                         size_t min_count, float threshold, bool multi_core, int parallel_depth)

    const string train_kl_cpp(const string id_, const string sequence_, size_t max_depth,
                                                                               size_t min_count, float threshold, bool multi_core, int parallel_depth)


cpdef str train(str name, str sequence, int max_depth, int min_count, int threshold, bool multi_core, int parallel_depth):
    """Trains a VLMC with the given sequence and parameters.

    Parameters
    ----------
    name : str
        Unique name given to the model.
    sequence : str
        Sequence to train VLMC on.
    max_depth : int
        The max depth of the VLMC.
    min_count : int
        The minimum amount of time each word has to appear in the sequence
        to be considered for inclusion as a node in the VLMC.
    threshold :  float
        Threshold value to prune to (for KL estimator).
        
    Examples
    ________
    >>> train("NC_009334.1", "AGAATTTGTCTTGCTCTATTCACCGTTACTTTTCTTCTTGCCCGTTCTCTTTCTTAGTAT", 15, 2, 1.2)
    'Name: NC_009334.1\nDate: Fri Jul 10 10:18:33 2020\nTree: PST\nAlphabet: DNA5\nNumber(nodes): 5\nNumber(parameters):
    15\nNode: 0 # [ 8 15 7 30 ] 60 [ 8 15 7 30 ][ -1 1 -1 2 ]\nNode: 1 C [ 2 3 2 8 ] 15 [ 1 3 2 9 ][ -1 3 -1 -1 ]\n
    Node: 2 T [ 3 9 4 14 ] 30 [ 4 8 3 14 ][ -1 -1 -1 4 ]\nNode: 3 CC [ 0 0 0 0 ] 3 [ 1 2 3 1 ][ -1 -1 -1 -1 ]\n
    Node: 4 TT [ 2 6 2 4 ] 14 [ 2 5 3 4 ][ -1 -1 -1 -1 ]\n'
    
    The first few rows of the output are metadata, such as the name, date, number of nodes and parameters.  The rows
    starting with 'Node: ' are the nodes of the probabilistic suffix tree representation of the variable-lenght
    Markov chain.  The format of each row is:
    'Node: {index} {k-mer} [{reverse counts}] [{forward counts}][{indices of children}'.
    The reverse counts are the counts of the k-mer `k` with an extra character `c` added to the left: `ck`.  The forward
    counts are the counts with the character added to the right: `kc`.  The reverse counts are sometimes missing when 
    they are infrequent, but we don't use them for anything so that's okay.
    
    Returns
    -------
    str
        The vlmc in a tree format.  First rows are metadata, and following rows are nodes in the tree.
    """
    cdef string tree = train_kl_cpp(name.encode(), sequence.encode(), max_depth, min_count, threshold, multi_core, parallel_depth)
    return tree.decode("utf-8")


cdef extern from "../../include/pst/distances/score.hpp" namespace "pst":
    vector[vector[double]] score_cpp(const vector[string]& tree_strings, const vector[string]& sequence_list)
    vector[vector[double]] score_background_cpp(const vector[string]& tree_strings, const vector[string]& sequence_list)

cdef extern from "../../include/pst/distances/sliding-windows.hpp" namespace "pst":
    vector[vector[double]] sliding_windows_cpp(string tree_string, string sequence, vector[int] window_sizes)
    vector[vector[double]] sliding_windows_background_cpp(string tree_string, string sequence, vector[int] window_sizes, int background_order)

cpdef vector[vector[double]] score_sequences_cython(list trees, list sequence_list):
    cdef vector[string] trees_ = [tree.encode() for tree in trees]
    cdef vector[string] arr = [s.encode() for s in sequence_list]

    return score_cpp(trees_, arr)

cpdef vector[vector[double]] score_sequences_background_cython(list trees, list sequence_list):
    cdef vector[string] trees_ = [tree.encode() for tree in trees]
    cdef vector[string] arr = [s.encode() for s in sequence_list]

    return score_background_cpp(trees_, arr)

cpdef vector[vector[double]] sliding_windows_cython(str tree, str sequence, list window_sizes):
    cdef string tree_ = tree.encode()
    cdef string sequence_ = sequence.encode()
    cdef vector[int] window_sizes_ = window_sizes

    return sliding_windows_cpp(tree_, sequence_, window_sizes_)

cpdef vector[vector[double]] sliding_windows_background_cython(str tree, str sequence, list window_sizes, int background_order):
    cdef string tree_ = tree.encode()
    cdef string sequence_ = sequence.encode()
    cdef vector[int] window_sizes_ = window_sizes
    cdef int background_order_ = background_order

    return sliding_windows_background_cpp(tree_, sequence_, window_sizes_, background_order_)


cpdef score_sequences(trees, sequence_list):
    import numpy as np

    return np.array(score_sequences_cython(trees, sequence_list))


cpdef score_sequences_background(trees, sequence_list):
    import numpy as np

    return np.array(score_sequences_background_cython(trees, sequence_list))


cpdef sliding_windows(tree, sequence, window_sizes):
    import numpy as np

    return np.array(sliding_windows_cython(tree, sequence, window_sizes))


cpdef sliding_windows_background(tree, sequence, window_sizes, background_order):
    import numpy as np

    return np.array(sliding_windows_background_cython(tree, sequence, window_sizes, background_order))


cdef extern from "../../include/pst/distances/d2.hpp" namespace "pst::distances":
    double d2_cpp(string left_string, string right_string)

cdef extern from "../../include/pst/distances/d2star.hpp" namespace "pst::distances":
    double d2star_cpp(string left_string, string right_string, int background_order)

cdef extern from "../../include/pst/distances/dvstar.hpp" namespace "pst::distances":
    double dvstar_cpp(string left_string, string right_string, int background_order)

cdef extern from "../../include/pst/distances/cv.hpp" namespace "pst::distances":
    double cv_cpp(string left_string, string right_string, int background_order)

cpdef double d2(str left_tree, str right_tree):
    cdef string left = left_tree.encode()
    cdef string right = right_tree.encode()

    return d2_cpp(left, right)

cpdef double d2star(str left_tree, str right_tree, int background_order):
    cdef string left = left_tree.encode()
    cdef string right = right_tree.encode()

    return d2star_cpp(left, right, background_order)


cpdef double dvstar(str left_tree, str right_tree, int background_order):
    cdef string left = left_tree.encode()
    cdef string right = right_tree.encode()

    return dvstar_cpp(left, right, background_order)

cpdef double cv(str left_tree, str right_tree, int background_order):
    cdef string left = left_tree.encode()
    cdef string right = right_tree.encode()

    return cv_cpp(left, right, background_order)
