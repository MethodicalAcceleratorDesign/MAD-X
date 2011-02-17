# Installation:

#FILE(GLOB files "${CMAKE_CURRENT_SOURCE_DIR}/*.h" "${CMAKE_CURRENT_BINARY_DIR}/*.h" "${CMAKE_CURRENT_BINARY_DIR}/*.mod")
#INSTALL(FILES ${files} DESTINATION include/madx)

SET (CMAKE_Fortran_MODULE_DIRECTORY
    ${PROJECT_BINARY_DIR}/bin/fortran CACHE PATH "Single Directory for all fortran modules."
)

set(APPS "\${CMAKE_INSTALL_PREFIX}/bin/madx")  # paths to executables
set(DIRS "")

if(APPLE)
  set(APPS "\${CMAKE_INSTALL_PREFIX}/madx.app")  # paths to executables
  set(DIRS "")
endif(APPLE)
INSTALL(TARGETS madx madxlib
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

# This installs the header files to <prefix>/include/madX
INSTALL (FILES ${headerfiles} 
	${CMAKE_CURRENT_BINARY_DIR}/c_wrappers.h 
	${CMAKE_CURRENT_BINARY_DIR}/c_wrappers_prototypes.h
	${CMAKE_CURRENT_BINARY_DIR}/c_prototypes.h
	${CMAKE_CURRENT_BINARY_DIR}/fortran_prototypes.h
	${CMAKE_CURRENT_BINARY_DIR}/fortran_wrappers.h 
 		DESTINATION "include/${PROJECT_NAME}")

# CPACK stuff
 # build a CPack driven installer package
 set (CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/License.txt")
 if ( PROJECT_PATCH_LEVEL EQUAL 00 )
 	set (CPACK_PACKAGE_NAME "${PROJECT_NAME}")
 else ( PROJECT_PATCH_LEVEL EQUAL 00 )
 	set (CPACK_PACKAGE_NAME "${PROJECT_NAME}-dev")
 endif ( PROJECT_PATCH_LEVEL EQUAL 00 )
 # Version:
 set (CPACK_PACKAGE_VERSION_MAJOR ${PROJECT_MAJOR_VERSION})
 set (CPACK_PACKAGE_VERSION_MINOR ${PROJECT_MINOR_VERSION})
 set (CPACK_PACKAGE_VERSION_PATCH ${PROJECT_PATCH_LEVEL})

 set (CPACK_PACKAGE_CONTACT "Frank Schmidt <Frank.Schmidt@cern.ch>")   
 set (CPACK_PACKAGE_DESCRIPTION_SUMMARY "MadX is a program for accelerator design and simulation") 
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
 set (CPACK_RPM_PACKAGE_RELEASE 1)
 set (CPACK_RPM_PACKAGE_LICENSE "custom")
 set (CPACK_RPM_PACKAGE_GROUP "Development/Tools")
 set(CPACK_RPM_PACKAGE_REQUIRES "libgcc >= 4.1.0, libxau >= 1.0.5") # I don't know the names of the packages...
 # so that we can build dragndrop on osx:
 set(CPACK_BINARY_DRAGNDROP ON)
 include (CPack)
# End CPACK stuff
