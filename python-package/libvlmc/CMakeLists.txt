add_cython_target(_libvlmc CXX)
add_library(_libvlmc MODULE ${_libvlmc})
python_extension_module(_libvlmc)

add_cython_target(KullbackLieblerTreeMap CXX)
add_library(KullbackLieblerTreeMap MODULE ${KullbackLieblerTreeMap})
python_extension_module(KullbackLieblerTreeMap)

target_link_libraries(_libvlmc seqan3::seqan3)
target_link_libraries(_libvlmc Threads::Threads)
target_link_libraries(_libvlmc CountVLMC::CountVLMC)

target_link_libraries(KullbackLieblerTreeMap seqan3::seqan3)
target_link_libraries(KullbackLieblerTreeMap Threads::Threads)
target_link_libraries(KullbackLieblerTreeMap CountVLMC::CountVLMC)

install(TARGETS _libvlmc LIBRARY DESTINATION libvlmc)
install(TARGETS KullbackLieblerTreeMap LIBRARY DESTINATION libvlmc)
