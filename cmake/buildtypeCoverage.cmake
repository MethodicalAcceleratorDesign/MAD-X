
# Build type for Coverage checks:
set(CMAKE_CXX_FLAGS_DEBUGCOVERAGE "${CMAKE_CXX_FLAGS_DEBUG} -Wall -fprofile-arcs -ftest-coverage")
set(CMAKE_C_FLAGS_DEBUGCOVERAGE "${CMAKE_C_FLAGS_DEBUG} -Wall -fprofile-arcs -ftest-coverage")
set(CMAKE_Fortran_FLAGS_DEBUGCOVERAGE "${CMAKE_Fortran_FLAGS_DEBUG} -Wall -fprofile-arcs -ftest-coverage")
set(LINK_FLAGS_DEBUGCOVERAGE "${LINK_FLAGS_DEBUG} -Wall -fprofile-arcs -ftest-coverage")
