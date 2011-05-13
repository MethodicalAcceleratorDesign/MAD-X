
# list of c source files
FILE (GLOB csrcfiles *.c)
# list of fortran source files
FILE (GLOB fsrcfiles *.f90 *.F90)

macro

# add source files according to NTPSA option...
if (MADX_NTPSA )
  message("NTPSA turned on")
	file(GLOB to_remove c_dabnew.f90)
	list(REMOVE_ITEM fsrcfiles ${to_remove})
else (MADX_NTPSA )
	file(GLOB to_remove c_dabnew_berz.f90 c_tpsa_interface.F90)
	list(REMOVE_ITEM fsrcfiles ${to_remove})
	file(GLOB to_remove tpsa.cpp)
	list(REMOVE_ITEM csrcfiles ${to_remove})
endif  (MADX_NTPSA )

#execute python wrapper scripts (you need to be dependent on one of the output files or else this command will never be ran):
# Unsure about dependencies.. Might be an overkill this one.
ADD_CUSTOM_COMMAND(
  OUTPUT c_wrappers.c c_wrappers.h c_prototypes.h c_wrappers_prototypes.h
  DEPENDS ${csrcfiles} ${fsrcfiles}
  COMMAND python wrap_C_calls.py -o ${CMAKE_CURRENT_BINARY_DIR} 
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  COMMENT "Creating C wrapper files"
  )
ADD_CUSTOM_COMMAND(
  OUTPUT fortran_wrappers.c fortran_wrappers.h fortran_prototypes.h fortran_wrappers_prototypes.h
  DEPENDS ${csrcfiles} ${fsrcfiles}
  COMMAND python wrap_fortran_calls.py --outdir=${CMAKE_CURRENT_BINARY_DIR} 
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  COMMENT "Creating fortran wrapper files"
  )
# execute python wrap_fortran_calls.py
# COMMAND python wrap_C_calls.py


# main source files... 
set(srcfiles ${csrcfiles} ${fsrcfiles} ${CMAKE_CURRENT_BINARY_DIR}/fortran_wrappers.c ${CMAKE_CURRENT_BINARY_DIR}/c_wrappers.c )

# header files...
file(GLOB headerfiles *.h)

