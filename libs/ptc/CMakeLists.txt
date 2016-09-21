cmake_minimum_required(VERSION 2.6)

# Some policies regarding configure
# warnings we do not want..
if(COMMAND cmake_policy)
  if(${CMAKE_VERSION} VERSION_GREATER 2.9)
    cmake_policy(SET CMP0042 NEW)
  endif()
endif(COMMAND cmake_policy)

# We want to run coverage analysis of
# the files compiled here..
include(buildtypeCoverage)

add_subdirectory(src)
