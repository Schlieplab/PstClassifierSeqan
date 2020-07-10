# distutils: language = c++

cdef extern from "../../src/vlmc.cpp":
    const char *train_kl_parameters(const char *id_, const char *sequence_,
                                    size_t max_depth, size_t min_count,
                                    size_t n_parameters)

    const char *train_kl(const char *id_, const char *sequence_, size_t max_depth,
                         size_t min_count, float threshold)

cpdef str train(str name, str sequence, int max_depth, int min_count, int threshold):
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
    cdef bytes tree = train_kl(name.encode(), sequence.encode(), max_depth, min_count, threshold)
    return tree.decode("utf-8")
