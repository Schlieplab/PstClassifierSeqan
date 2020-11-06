add_cython_target(_libvlmc CXX)
add_library(_libvlmc MODULE ${_libvlmc})
python_extension_module(_libvlmc)

target_link_libraries(_libvlmc seqan3::seqan3)
target_link_libraries(_libvlmc Threads::Threads)
target_link_libraries(_libvlmc HighFive)

install(TARGETS _libvlmc LIBRARY DESTINATION libvlmc)
