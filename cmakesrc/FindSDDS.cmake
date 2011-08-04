# - Find SDDS Libraries
# This module finds SDDS libraries if they exist, and sets the following
# variables:
#
#  SDDS_FOUND         - SDDS successfully found, including all libraries it depends on
#  SDDS_LIBRARIES     - SDDS libraries
#  SDDS_SEARCH_DIRS   - Additional directories to search in (these are prepended to the list)

message(STATUS "Looking for SDDS libraries")
set(SDDS_FOUND TRUE)
# set(SDDS_LIBRARIES )
set(SDDS_INCLUDE_DIR "")
# libraries needed to be found:
set(_SDDS_LIBS SDDS1c  SDDS1 rpnlib mdbmth mdblib gsl)

# directories to search in:
set(_dirs ${SDDS_SEARCH_DIRS} 
          /usr/lib/ 
          /usr/lib64/ 
          /usr/local/lib/ 
          /usr/local/lib64)

foreach(_lib ${_SDDS_LIBS})
    find_library(SDDS_${_lib}_LIB NAMES ${_lib} 
      PATHS ${_dirs}
      NO_DEFAULT_PATH
    )
    if(NOT SDDS_${_lib}_LIB)
        message(STATUS "Did not find lib${_lib}")
        set(SDDS_FOUND FALSE)
    else()
        set(SDDS_LIBRARIES ${SDDS_LIBRARIES} ${SDDS_${_lib}_LIB})
    endif()
    mark_as_advanced(SDDS_${_lib}_LIB)
endforeach()

if(CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
    # we trust the linker to find this library..
    set(SDDS_LIBRARIES ${SDDS_LIBRARIES} imf)
endif()

if(SDDS_FOUND)
    message(STATUS "Looking for SDDS libraries - found.")
else()
    if(SDDS_FIND_REQUIRED)
        message(FATAL_ERROR "Looking for SDDS libraries - not found.")
    else()
        message(STATUS "Looking for SDDS libraries - not found.")
    endif()
endif()
    
