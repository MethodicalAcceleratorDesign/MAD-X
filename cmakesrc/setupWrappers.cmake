###
#
# This file finds the Python binary and adds the commands to produce the sources/headers
# 
###


#execute python wrapper scripts (you need to be dependent on one of the output files or else this command will never be ran):
# Unsure about dependencies.. Might be an overkill this one.
find_package(PythonInterp REQUIRED)

if(PYTHON_VERSION_MAJOR) # for older cmake this variable is not defined..
    if(${PYTHON_VERSION_MAJOR} EQUAL 3)
        set(PYSCRPT_END "_py3")
    endif()
endif()
ADD_CUSTOM_COMMAND(
  OUTPUT c_wrappers.c c_wrappers.h c_prototypes.h c_wrappers_prototypes.h
  DEPENDS ${csrcfiles} ${fsrcfiles}
  COMMAND ${PYTHON_EXECUTABLE} wrap_C_calls${PYSCRPT_END}.py -o ${CMAKE_CURRENT_BINARY_DIR} 
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  COMMENT "Creating C wrapper files"
  )
ADD_CUSTOM_COMMAND(
  OUTPUT fortran_wrappers.c fortran_wrappers.h fortran_prototypes.h fortran_wrappers_prototypes.h
  DEPENDS ${csrcfiles} ${fsrcfiles}
  COMMAND ${PYTHON_EXECUTABLE} wrap_fortran_calls${PYSCRPT_END}.py --outdir=${CMAKE_CURRENT_BINARY_DIR} 
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  COMMENT "Creating fortran wrapper files"
  )
# execute python wrap_fortran_calls.py
# COMMAND python wrap_C_calls.py
