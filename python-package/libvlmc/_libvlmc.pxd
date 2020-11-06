from libcpp cimport bool
from libcpp.vector cimport vector
from libcpp.string cimport string

cdef const char *train_kl_parameters(const char *id_, const char *sequence_,
                                size_t max_depth, size_t min_count,
                                size_t n_parameters)
cdef const char *train_kl(const char *id_, const char *sequence_, size_t max_depth,
                     size_t min_count, float threshold)

cpdef str train(str name, str sequence, int max_depth, int min_count, int threshold)

cdef vector[vector[double]] score_cpp(vector[string] tree_strings, vector[string] sequence_list)

cpdef vector[vector[double]] score_sequences_cython(list hdf5_path, list sequence_list)

cpdef score_sequences(hdf5_path, sequence_list)
