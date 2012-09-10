

file(GLOB src_files *.f90)

add_library(ptc ${src_files})

set_target_properties(ptc PROPERTIES VERSION ${MADX_VERSION} SOVERSION 1)

madx_install_targets(ptc)
