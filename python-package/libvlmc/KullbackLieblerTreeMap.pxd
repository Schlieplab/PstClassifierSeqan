# distutils: language = c++

from libcpp.string cimport string
from libcpp.vector cimport vector
from libcpp cimport bool


cdef extern from "../../include/pst/kl_tree_map.hpp" namespace "pst":
    cdef cppclass KullbackLieblerTreeMap5:
        KullbackLieblerTreeMap5() except +
        string get_closest_state(string label)
        float get_transition_probability(string label, char character)
        string to_tree()
        string generate_sequence(int length)


cdef extern from "../../include/pst/kl_tree_map.hpp" namespace "pst":
    KullbackLieblerTreeMap5 create_kl_map(string id, string sequence, int max_depth, int freq, float cutoff_value_, int number_of_parameters, string pruning_method, bool multi_core, int parallel_depth)

cdef extern from "../../include/pst/distances/score.hpp" namespace "pst":
    vector[vector[double]] score_sequences(vector[string] tree_strings, vector[string] sequence_list, int background_order)

cdef extern from "../../include/pst/distances/d2.hpp" namespace "pst::distances":
    double d2(KullbackLieblerTreeMap5 left, KullbackLieblerTreeMap5 right)

cdef extern from "../../include/pst/distances/d2star.hpp" namespace "pst::distances":
    double d2star(KullbackLieblerTreeMap5 left, KullbackLieblerTreeMap5 right, int background_order)

cdef extern from "../../include/pst/distances/dvstar.hpp" namespace "pst::distances":
    double dvstar(KullbackLieblerTreeMap5 left, KullbackLieblerTreeMap5 right, int background_order)

cdef extern from "../../include/pst/distances/cv.hpp" namespace "pst::distances":
    double cv(KullbackLieblerTreeMap5 left, KullbackLieblerTreeMap5 right, int background_order)
