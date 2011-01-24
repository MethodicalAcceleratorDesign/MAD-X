cmake_minimum_required(VERSION 2.6)

PROJECT(madX C CXX Fortran) 
#add name of languages used in the project as well

if(COMMAND cmake_policy)
  cmake_policy(SET CMP0003 NEW)
#   cmake_policy(SET CMP0014 OLD) #ignore not having CMakeLists.txt in subfolders (lib/lib64)
endif(COMMAND cmake_policy)

# project version
SET( ${PROJECT_NAME}_MAJOR_VERSION 4 )
SET( ${PROJECT_NAME}_MINOR_VERSION 01 )
SET( ${PROJECT_NAME}_PATCH_LEVEL 40 )

# add source dir as a place for CMake modules (e.g. FindGSL.cmake)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmakesrc") 


# project options
if ( NOT APPLE)
   OPTION( BUILD_SHARED_LIBS "Set to OFF to build static libraries"  OFF )
else ( NOT APPLE)
   OPTION( BUILD_SHARED_LIBS "Set to OFF to build static libraries"  ON )
endif ( NOT APPLE)
# OPTION( INSTALL_DOC "Set to OFF to skip build/install Documentation" OFF )

#Mad-X specific options (arch. specific options can be added in similar manner):
SET ( MADX_OPTION_ONLINE 0 CACHE VALUE "Set to 1 to build with Online model" )
SET ( MADX_OPTION_NTPSA 1 CACHE VALUE "Set to 1 to build with NTPSA" )
SET ( MADX_OPTION_FORCE_32 0 CACHE VALUE "Set to 1 to force 32bit build" )
SET ( MADX_OPTION_FORCE_FORCE_FORTRAN_COMPILER 0 CACHE VALUE "Set to 1 to force a specific fortran compiler, otherwise system default will be used" )

# Build type (defines different sets of flags)
SET(CMAKE_BUILD_TYPE RELEASE)
if ( NOT madX_PATCH_LEVEL EQUAL 00 )
 message( "Building a development version" )
 SET (BINARY_POSTFIX "_dev")
endif ( NOT madX_PATCH_LEVEL EQUAL 00 )

#My own build type (as an example):
#SET(CMAKE_BUILD_TYPE madbuild)
# SET(CMAKE_C_FLAGS_MADBUILD "-O4")

# include(FortranCInterface) # does cmake have automatic f/c interface??

if ( MADX_OPTION_FORCE_FORCE_FORTRAN_COMPILER EQUAL 1 )
  set(CMAKE_Fortran_COMPILER "g95")
  message ("Forcing fortran compiler: " ${CMAKE_Fortran_COMPILER})
  message ( WARNING "This is not recommended, rather use e.g. '-DCMAKE_Fortran_COMPILER=g95' as option when running cmake")
endif (MADX_OPTION_FORCE_FORCE_FORTRAN_COMPILER EQUAL 1)


#fortran compiler stuff... extensive example
# FFLAGS depend on the compiler

# SET(CMAKE_FORTRAN_COMPILER mpif77)
if (CMAKE_Fortran_COMPILER MATCHES "gfortran")
  message( WARNING " g95 or ifort recommended, run cmake with e.g. '-DCMAKE_Fortran_COMPILER=g95' option")
  set (CMAKE_Fortran_FLAGS_RELEASE " -g -funroll-loops -fno-range-check -fno-f2c -O2 ")
  set (CMAKE_Fortran_FLAGS_DEBUG   " -fno-f2c -O0 -g ")
  if ( NOT BUILD_SHARED_LIBS )
    set (CMAKE_Fortran_LINK_FLAGS   "${CMAKE_Fortran_LINK_FLAGS} -static ")
  endif ( NOT BUILD_SHARED_LIBS )
elseif (CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
  set (CMAKE_Fortran_FLAGS_RELEASE " -funroll-loops -f77rtl -O2 -assume noold_unit_star -D_INTEL_IFORT_SET_RECL")
  set (CMAKE_Fortran_FLAGS_DEBUG   " -f77rtl -O0 -g -assume noold_unit_star -D_INTEL_IFORT_SET_RECL")
  set (CMAKE_Fortran_LINK_FLAGS   "${CMAKE_Fortran_LINK_FLAGS} -static ")
elseif (CMAKE_Fortran_COMPILER MATCHES "lf95")
  message( WARNING " This compiler is not yet confirmed working properly with CMake")
  message( WARNING " On a 64 bit system you need to use the toolchain-file (see README) to get anywhere with the 32bit compiler.")
  SET ( MADX_OPTION_FORCE_32 1 )
  set (CMAKE_Fortran_FLAGS_RELEASE " --o2 --tp ")
  set (CMAKE_SKIP_RPATH ON)
  set (CMAKE_Fortran_FLAGS_DEBUG   " --info --f95 --lst -V -g  --ap --trace --trap --verbose ")
  set(CMAKE_SHARED_LIBRARY_LINK_Fortran_FLAGS "") #suppress rdynamic which doesn't work for lf95...
  set (CMAKE_Fortran_LINK_FLAGS   "${CMAKE_Fortran_LINK_FLAGS} -static ")
elseif (CMAKE_Fortran_COMPILER MATCHES "nagfor")
  message( WARNING " Make sure you use the same gcc as nagfor is compiled with, or linking WILL fail.")
  set (CMAKE_SKIP_RPATH ON)
  set (CMAKE_Fortran_FLAGS_RELEASE " -gline -maxcontin=100 -ieee=full -D_NAG ")
  set (CMAKE_Fortran_FLAGS_DEBUG   " -gline -maxcontin=100 -ieee=full -D_NAG -C=all -nan ")
  set(CMAKE_SHARED_LIBRARY_LINK_Fortran_FLAGS "") #suppress rdynamic which isn't recognized by nagfor...
  set (CMAKE_Fortran_LINK_FLAGS   "${CMAKE_Fortran_LINK_FLAGS} -Bstatic ")
elseif (CMAKE_Fortran_COMPILER MATCHES "g77")
  message( WARNING " This compiler is not yet confirmed working for mad-x")
  message( WARNING " g95 or ifort recommended, run cmake with e.g. '-DCMAKE_Fortran_COMPILER=g95' option")
  set (CMAKE_Fortran_FLAGS_RELEASE " -funroll-loops -fno-f2c -O3 ")
  set (CMAKE_Fortran_FLAGS_DEBUG   " -fno-f2c -O0 -g ")
elseif (CMAKE_Fortran_COMPILER MATCHES "g95")
  set (CMAKE_Fortran_FLAGS_RELEASE " -funroll-loops -fno-second-underscore -fshort-circuit -O2 ")
  set (CMAKE_Fortran_FLAGS_DEBUG   " -fno-second-underscore -O0 -g -Wall -pedantic -ggdb3")  
else (CMAKE_Fortran_COMPILER MATCHES "gfortran")
  message ("Fortran compiler full path: " ${CMAKE_Fortran_COMPILER})
  message ("Fortran compiler: " ${Fortran_COMPILER_NAME})
  message (WARNING " No optimized Fortran compiler flags are known for this compiler...")
  set (CMAKE_Fortran_FLAGS_RELEASE " -funroll-loops -fno-range-check -O2")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g")
endif (CMAKE_Fortran_COMPILER MATCHES "gfortran")
#end fortran compiler stuff...


#adding flags to our compilers:

# General compile flags:
set (CMAKE_C_FLAGS_DEBUG   " ${CMAKE_C_FLAGS_DEBUG} -Wall -pedantic ")
set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g -funroll-loops -D_CATCH_MEM -D_WRAP_FORTRAN_CALLS -D_WRAP_C_CALLS -D_FULL -I. -I${CMAKE_CURRENT_SOURCE_DIR} ")
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -funroll-loops -D_CATCH_MEM -D_WRAP_FORTRAN_CALLS -D_WRAP_C_CALLS -D_FULL -I. -I${CMAKE_CURRENT_SOURCE_DIR} ") #needed for c++ linking
set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -D_CATCH_MEM -D_WRAP_FORTRAN_CALLS -D_WRAP_C_CALLS -D_FULL -I. -I${CMAKE_CURRENT_SOURCE_DIR} ") 

if ( MADX_OPTION_FORCE_32 EQUAL 1 )
  message("Forced 32bit build" ) 
  set (CPACK_DEBIAN_PACKAGE_ARCHITECTURE "i386")
  set (CPACK_RPM_PACKAGE_ARCHITECTURE "i386")
  set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -m32")
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -m32")
  if (CMAKE_Fortran_COMPILER MATCHES "lf95")
    set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Wa,--32")
  elseif (CMAKE_Fortran_COMPILER MATCHES "nagfor")
    set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}  -Wc,-m32 -Wl,-m32 -abi=32")
  else (CMAKE_Fortran_COMPILER MATCHES "lf95")
    set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} --chk aesux")
  else (CMAKE_Fortran_COMPILER MATCHES "gfortran")
    set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -m32")
  endif (CMAKE_Fortran_COMPILER MATCHES "lf95")
  set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -m32 ")
  if (CMAKE_SYSTEM_NAME STREQUAL "Linux")
    LINK_DIRECTORIES(${CMAKE_SOURCE_DIR}/lib)
  endif (CMAKE_SYSTEM_NAME STREQUAL "Linux")
else ( MADX_OPTION_FORCE_32 EQUAL 1 )
   INCLUDE(CheckTypeSize)
   CHECK_TYPE_SIZE (long LONG_SIZE) #4 if 32bit, 8 if 64bit...
  if (${LONG_SIZE} EQUAL 4)
  message("32 bit system found" )
  message("${CMAKE_SYSTEM_PROCESSOR}")
   set (CPACK_DEBIAN_PACKAGE_ARCHITECTURE "i386")
   set (CPACK_RPM_PACKAGE_ARCHITECTURE "i386")
   set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -m32")
   set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -m32")
  if (CMAKE_Fortran_COMPILER MATCHES "lf95")
    set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} --chk aesux")
  endif (CMAKE_Fortran_COMPILER MATCHES "lf95")
   if (CMAKE_SYSTEM_NAME STREQUAL "Linux")
     LINK_DIRECTORIES(${CMAKE_SOURCE_DIR}/lib)
   endif (CMAKE_SYSTEM_NAME STREQUAL "Linux")
 elseif (${LONG_SIZE} EQUAL 8)
  message("64 bit system found")
	 #set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Wa,--64")
	 set (CPACK_DEBIAN_PACKAGE_ARCHITECTURE "amd64")
   set (CPACK_RPM_PACKAGE_ARCHITECTURE "x86_64")
	 set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -m64 ")
	 set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -m64 ")
  if (CMAKE_Fortran_COMPILER MATCHES "lf95")
    set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} --chk aefs")
  endif (CMAKE_Fortran_COMPILER MATCHES "lf95")
	 if (CMAKE_SYSTEM_NAME STREQUAL "Linux")
	   LINK_DIRECTORIES(${CMAKE_SOURCE_DIR}/lib64)
	 endif (CMAKE_SYSTEM_NAME STREQUAL "Linux")
  endif (${LONG_SIZE} EQUAL 4)
endif ( MADX_OPTION_FORCE_32 EQUAL 1 )

if (MADX_OPTION_ONLINE EQUAL 1)
  message("Online Model turned on" )
  set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -D_ONLINE ")
  set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -D_ONLINE ")
endif (MADX_OPTION_ONLINE EQUAL 1)

# OSX specifics:
if (APPLE)
  set(CMAKE_LIBRARY_PATH /usr/lib/ /usr/X11/lib/ ${CMAKE_LIBRARY_PATH})
endif (APPLE)


# C stuff:
get_filename_component (C_COMPILER_NAME ${CMAKE_C_COMPILER} NAME)
if (C_COMPILER_NAME MATCHES "gcc" AND NOT CMAKE_Fortran_COMPILER MATCHES "gfortran")
# this will probably crash on windows...
  execute_process(COMMAND ${C_COMPILER_NAME} -print-search-dirs
                OUTPUT_VARIABLE gccsearchdirs)
  string(REGEX REPLACE ".*libraries: =(.*)\n"  "\\1" gcclibs "${gccsearchdirs}")
  # need to do this many times because lf95 segfaults on lists with :
  string(REPLACE "/:/"  "/ -L/" gcclibs "${gcclibs}")
  # adding these to the linking process which is handled by a non-gnu fortran compiler in your case
  LINK_DIRECTORIES(${gcclibs}) 
endif (C_COMPILER_NAME MATCHES "gcc" AND NOT CMAKE_Fortran_COMPILER MATCHES "gfortran")
# end C stuff



# list of c source files
set(csrcfiles madxp.c gxx11c.c matchptcknobs.c rplot.c tpsa.cpp )
# list of fortran source files
set(fsrcfiles madx_main.F90 Sqb_accel_ptc.f90 Sr_spin.f90 Sra_fitting.f90 madx_ptc_intstate.f90 madx_ptc_setcavs.f90 madx_ptc_knobs.f90 madx_ptc_module.f90 madx_ptc_distrib.f90 madx_ptc_eplacement.f90 madx_ptc_normal.f90 madx_ptc_script.f90 madx_ptc_trackcavs.f90 madx_ptc_track_run.f90 madx_ptc_twiss.f90 match.f90 matchjc.f90 matchlib.f90 matchsa.f90 orbf.f90 plot.f90 poisson.f90 ptc_export_xml.f90 resindex.f90 run_madx.f90 sodd.f90 Spb_fake_gino_sub.f90 St_pointers.f90 survey.f90 timest.f90 timex.f90 touschek.f90 trrun.f90 twiss.f90 user2_photon.f90 wrap.f90 fortran_flush.F90 a_scratch_size.f90 b_da_arrays_all.f90 d_lielib.f90 util.f90 dynap.f90 emit.f90 gxx11.f90 h_definition.f90 ibsdb.f90 i_tpsa.f90 j_tpsalie.f90 k_tpsalie_analysis.f90 l_complex_taylor.f90 m_real_polymorph.f90 n_complex_polymorph.f90 o_tree_element.f90 Sa_extend_poly.f90 Sb_sagan_pol_arbitrary.f90 Sc_euclidean.f90 Sd_frame.f90 Se_status.f90 Sf_def_all_kinds.f90 Sg_sagan_wiggler.f90 Sh_def_kind.f90 Si_def_element.f90 Sk_link_list.f90 Sl_family.f90 Sm_tracking.f90 Sma0_beam_beam_ptc.f90 Sma_multiparticle.f90 Sn_mad_like.f90 So_fitting.f90 Sp_keywords.f90 Sq_orbit_ptc.f90)


# add source files according to NTPSA option...
if (MADX_OPTION_NTPSA EQUAL 1)
  message("NTPSA turned on")
 set(fsrcfiles ${fsrcfiles} c_dabnew_berz.f90 c_tpsa_interface.F90)
 set(csrcfiles ${csrcfiles} tpsa.cpp)
else (MADX_OPTION_NTPSA EQUAL 1)
 set(fsrcfiles ${fsrcfiles} c_dabnew.f90)
endif  (MADX_OPTION_NTPSA EQUAL 1)

#execute python wrapper scripts (you need to be dependent on one of the output files or else this command will never be ran):
ADD_CUSTOM_COMMAND(
  OUTPUT fortran_wrappers.c fortran_wrappers.h fortran_prototypes.h fortran_wrappers_prototypes.h
  DEPENDS ${fsrcfiles}
  COMMAND cp ${CMAKE_CURRENT_SOURCE_DIR}/*.c . && cp ${CMAKE_CURRENT_SOURCE_DIR}/*.f90 . && python ${CMAKE_CURRENT_SOURCE_DIR}/wrap_fortran_calls.py && cp ${CMAKE_CURRENT_SOURCE_DIR}/*.h .
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  COMMENT "Creating fortran wrapper files"
  )
# Unsure about dependencies.. Might be an overkill this one.
ADD_CUSTOM_COMMAND(
  OUTPUT c_wrappers.c c_wrappers.h c_prototypes.h c_wrappers_prototypes.h
  DEPENDS ${csrcfiles} ${fsrcfiles} fortran_wrappers.c
  COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/wrap_C_calls.py
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  COMMENT "Creating C wrapper files"
  )
# execute python wrap_fortran_calls.py
# COMMAND python wrap_C_calls.py


# main source files... 
# TODO: please check this list at some point!
set(srcfiles ${csrcfiles} ${fsrcfiles} fortran_wrappers.c c_wrappers.c )


if(APPLE)
  SET(MACOSX_BUNDLE_STARTUP_COMMAND madx${BINARY_POSTFIX})
  SET(MACOSX_BUNDLE_ICON_FILE "${CMAKE_CURRENT_SOURCE_DIR}/cmakesrc/MadX.icns")
  SET(MACOSX_BUNDLE_LONG_VERSION_STRING "MadX ${BINARY_POSTFIX} version ${madX_MAJOR_VERSION}.${madX_MINOR_VERSION}.${madX_PATCH_LEVEL}")
  SET(MACOSX_BUNDLE_BUNDLE_NAME "MadX${BINARY_POSTFIX}")
  SET(MACOSX_BUNDLE_GUI_IDENTIFIER "MadX${BINARY_POSTFIX}")
  # add icns to the .app/Resources with these two commands:
  SET(srcfiles ${srcfiles} ${CMAKE_CURRENT_SOURCE_DIR}/cmakesrc/MadX.icns)
  SET_SOURCE_FILES_PROPERTIES(${CMAKE_CURRENT_SOURCE_DIR}/cmakesrc/MadX.icns PROPERTIES MACOSX_PACKAGE_LOCATION Resources)
endif(APPLE)

add_executable(madx${BINARY_POSTFIX} MACOSX_BUNDLE ${srcfiles})
SET_TARGET_PROPERTIES(madx${BINARY_POSTFIX} PROPERTIES LINKER_LANGUAGE Fortran)

#I turn off search for libraries in case you are on Linux,
#to make sure we make use of the lib/lib64 folders
if(NOT CMAKE_SYSTEM_NAME STREQUAL "Linux")
#necessary to search for X11 for OSX instead of directly including it
FIND_PACKAGE(X11)
IF(X11_FOUND)
   message("Found X11 libraries")
   INCLUDE_DIRECTORIES(${X11_INCLUDE_DIR})
   TARGET_LINK_LIBRARIES(madx${BINARY_POSTFIX} ${X11_X11_LIB})
ENDIF(X11_FOUND)
else (NOT CMAKE_SYSTEM_NAME STREQUAL "Linux")
  target_link_libraries(madx${BINARY_POSTFIX} X11)
endif (NOT CMAKE_SYSTEM_NAME STREQUAL "Linux")

# Online libraries:
if (MADX_OPTION_ONLINE EQUAL 1)
target_link_libraries(madx${BINARY_POSTFIX} SDDS1c)
target_link_libraries(madx${BINARY_POSTFIX} SDDS1)
target_link_libraries(madx${BINARY_POSTFIX} gsl) # SDDS depends on gsl..
target_link_libraries(madx${BINARY_POSTFIX} rpnlib)
target_link_libraries(madx${BINARY_POSTFIX} mdbmth)
target_link_libraries(madx${BINARY_POSTFIX} mdblib)
endif (MADX_OPTION_ONLINE EQUAL 1)
target_link_libraries(madx${BINARY_POSTFIX} z)
# new additions from Frank:
target_link_libraries(madx${BINARY_POSTFIX} pthread)
target_link_libraries(madx${BINARY_POSTFIX} c)
target_link_libraries(madx${BINARY_POSTFIX} gcc_eh)


# Installation:

#FILE(GLOB files "${CMAKE_CURRENT_SOURCE_DIR}/*.h" "${CMAKE_CURRENT_BINARY_DIR}/*.h" "${CMAKE_CURRENT_BINARY_DIR}/*.mod")
#INSTALL(FILES ${files} DESTINATION include/madx)


set(APPS "\${CMAKE_INSTALL_PREFIX}/bin/madx${BINARY_POSTFIX}")  # paths to executables
set(DIRS "")

if(APPLE)
  set(APPS "\${CMAKE_INSTALL_PREFIX}/madx${BINARY_POSTFIX}.app")  # paths to executables
  set(DIRS "")
endif(APPLE)
INSTALL(TARGETS madx${BINARY_POSTFIX}
  BUNDLE DESTINATION .
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
)
if(APPLE OR WIN32) # I don't think this is supposed to have a function on GNU/Linux systems?
  INSTALL(CODE " 
    include(BundleUtilities) 
    fixup_bundle(\"${APPS}\"   \"\"   \"${DIRS}\") 
    " COMPONENT Runtime) 
endif(APPLE OR WIN32)

# CPACK stuff
 # build a CPack driven installer package
 set (CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/License.txt")
 
 set (CPACK_PACKAGE_CONTACT "Frank Schmidt <Frank.Schmidt@cern.ch>")   
 
 #Debian specific:
 set (CPACK_DEBIAN_PACKAGE_MAINTAINER
      "Yngve Inntjore Levinsen <Yngve.Inntjore.Levinsen@cern.ch>") # if this is not set, CPACK_PACKAGE_CONTACT is used instead..
 set (CPACK_DEBIAN_PACKAGE_DESCRIPTION "MadX is a program for accelerator design and simulation")
 set (CPACK_DEBIAN_PACKAGE_SECTION "science")
 set (CPACK_DEBIAN_PACKAGE_PRIORITY "extra")
 set (CPACK_DEBIAN_PACKAGE_RECOMMENDS "gnuplot")
 set (CPACK_DEBIAN_PACKAGE_DEPENDS 
      "libc6 (>= 2.3.1-6), libgcc1 (>= 1:4.1), zlib1g, libx11-xcb1, libxcb1, libxau6") # some version dependencies kept as examples
 
 # RPM Specific:
 set (CPACK_RPM_PACKAGE_SUMMARY
      "MadX is a program for accelerator design and simulation")
 set (CPACK_RPM_PACKAGE_RELEASE 1)
 set (CPACK_RPM_PACKAGE_LICENSE "custom")
 set (CPACK_RPM_PACKAGE_GROUP "Development/Tools")
 set(CPACK_RPM_PACKAGE_REQUIRES "libgcc >= 4.1.0, libxau >= 1.0.5") # I don't know the names of the packages...
 
 # Version:
 set (CPACK_PACKAGE_VERSION_MAJOR ${madX_MAJOR_VERSION})
 set (CPACK_PACKAGE_VERSION_MINOR ${madX_MINOR_VERSION})
 set (CPACK_PACKAGE_VERSION_PATCH ${madX_PATCH_LEVEL})
 
 
 # so that we can build dragndrop on osx:
 set(CPACK_BINARY_DRAGNDROP ON)
 include (CPack)
# End CPACK stuff


# Testing:
# enable dashboard scripting

enable_testing()
set (CTEST_PROJECT_NAME "MadX")
ADD_TEST(testruns "${CMAKE_CURRENT_SOURCE_DIR}/cmakesrc/ctests/testruns.sh" "${BINARY_POSTFIX}")
set_tests_properties (testruns 
  PROPERTIES PASS_REGULAR_EXPRESSION "Number of warnings: 0")
ADD_TEST(testsample "${CMAKE_CURRENT_SOURCE_DIR}/cmakesrc/ctests/testsample.sh" "${BINARY_POSTFIX}" "${CMAKE_CURRENT_SOURCE_DIR}/cmakesrc/ctests")
set_tests_properties (testsample 
  PROPERTIES FAIL_REGULAR_EXPRESSION "Number of errors: [1-9]") # The no warning check fails on this test..

