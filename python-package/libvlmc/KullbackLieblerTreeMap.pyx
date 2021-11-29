# distutils: language = c++

from libcpp.string cimport string
from libcpp cimport bool


from .KullbackLieblerTreeMap cimport KullbackLieblerTreeMap5, create_kl_map, d2, d2star, dvstar, cv, score_sequences

def create_vlmc(id: str, sequence: str, max_depth: int, freq: int, threshold: float, multi_core: bool = True, parallel_depth: int = 2):
    return pyKullbackLieblerTreeMap(id.encode(), sequence.encode(), max_depth, freq, threshold, multi_core, parallel_depth)

cdef class pyKullbackLieblerTreeMap:
    cdef KullbackLieblerTreeMap5 c_vlmc

    def __cinit__(self, string id, string sequence, int max_depth, int freq, float cutoff_value, bool multi_core, int parallel_depth):
        self.c_vlmc = create_kl_map(id, sequence, max_depth, freq, cutoff_value, 1, b"cutoff", multi_core, parallel_depth)

    def get_closest_state(self, label: str) -> str:
        return self.c_vlmc.get_closest_state(label.encode()).decode()

    def to_tree(self) -> str:
        return self.c_vlmc.to_tree().decode()

    def get_transition_probability(self, label: str, character: str):
        return self.c_vlmc.get_transition_probability(label, character)


# def negative_log_likelihood(left: list[pyKullbackLieblerTreeMap], right: list[str], background_order: int = 0) -> float:
#     return score_sequences([l.c_vlmc for l in left], [r.encode() for r in right], background_order)

def distance_d2(left: pyKullbackLieblerTreeMap, right: pyKullbackLieblerTreeMap) -> float:
    return d2(left.c_vlmc, right.c_vlmc)

def distance_d2star(left: pyKullbackLieblerTreeMap, right: pyKullbackLieblerTreeMap, background_order: int = 0) -> float:
    return d2star(left.c_vlmc, right.c_vlmc, background_order)

def distance_dvstar(left: pyKullbackLieblerTreeMap, right: pyKullbackLieblerTreeMap, background_order: int = 0) -> float:
    return dvstar(left.c_vlmc, right.c_vlmc, background_order)

def distance_cv(left: pyKullbackLieblerTreeMap, right: pyKullbackLieblerTreeMap, background_order: int = 0) -> float:
    return cv(left.c_vlmc, right.c_vlmc, background_order)
