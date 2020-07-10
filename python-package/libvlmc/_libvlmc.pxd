cdef const char *train_kl_parameters(const char *id_, const char *sequence_,
                                size_t max_depth, size_t min_count,
                                size_t n_parameters)
cdef const char *train_kl(const char *id_, const char *sequence_, size_t max_depth,
                     size_t min_count, float threshold)

cpdef str train(str name, str sequence, int max_depth, int min_count, int threshold)