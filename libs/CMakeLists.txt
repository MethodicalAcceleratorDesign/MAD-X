
cmake_minimum_required(VERSION 2.6)

# We want to run coverage analysis of
# the files compiled here..
include(buildtypeCoverage)

add_subdirectory(ptc)
if(USE_GC)
   add_subdirectory(gc/gc)
endif()

