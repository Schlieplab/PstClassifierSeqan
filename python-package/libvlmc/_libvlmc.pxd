from libcpp cimport bool
from libcpp.vector cimport vector
from libcpp.string cimport string

cdef const char *train_kl_parameters(const char *id_, const char *sequence_,
                                size_t max_depth, size_t min_count,
                                size_t n_parameters, bool multi_core, int parallel_depth)

cdef const char *train_kl(const char *id_, const char *sequence_, size_t max_depth,
                     size_t min_count, float threshold, bool multi_core, int parallel_depth)

cpdef str train(str name, str sequence, int max_depth, int min_count, int threshold, bool multi_core, int parallel_depth)

cdef vector[vector[double]] score_cpp(const vector[string]& tree_strings, const vector[string]& sequence_list)

cpdef vector[vector[double]] score_sequences_cython(list trees, list sequence_list)
cpdef score_sequences(trees, sequence_list)
cpdef score_sequences_background(trees, sequence_list)

cpdef vector[vector[double]] sliding_windows_cython(str tree, str sequence, list window_sizes)
cpdef vector[vector[double]] sliding_windows_background_cython(str tree, str sequence, list window_sizes, int background_order)

cpdef sliding_windows(tree, sequence, window_sizes)
cpdef sliding_windows_background(tree, sequence, window_sizes, background_order)

cpdef double d2(str left_tree, str right_tree)
cpdef double d2star(str left_tree, str right_tree, int background_order)
cpdef double cv(str left_tree, str right_tree, int background_order)
