# Installation:

#FILE(GLOB files "${CMAKE_CURRENT_SOURCE_DIR}/*.h" "${CMAKE_CURRENT_BINARY_DIR}/*.h" "${CMAKE_CURRENT_BINARY_DIR}/*.mod")
#INSTALL(FILES ${files} DESTINATION include/madx)

set(CMAKE_Fortran_MODULE_DIRECTORY
    ${PROJECT_BINARY_DIR}/include/fortran/madX CACHE PATH "Single Directory for all fortran modules."
)

install(TARGETS madxbin madx madx_fpp madx_ptc_core madx_ptc
  BUNDLE DESTINATION . COMPONENT Runtime
  RUNTIME DESTINATION bin COMPONENT Runtime
  LIBRARY DESTINATION lib COMPONENT Libraries
  ARCHIVE DESTINATION lib COMPONENT Libraries
)

# This installs the header files to <prefix>/include/madX
# We only want this in the development version...
if(NOT ${MADX_PATCH_LEVEL} EQUAL 00)
INSTALL (FILES ${headerfiles} 
  DESTINATION "include/${PROJECT_NAME}"
  COMPONENT Files)
endif()

INSTALL(FILES ${CMAKE_SOURCE_DIR}/License.txt 
    DESTINATION "share/doc/${PROJECT_NAME}${PKG_POSTFIX}"
    COMPONENT Files)
if(APPLE AND MADX_BUNDLE)
  set(APPS "\${CMAKE_INSTALL_PREFIX}/madx${BINARY_POSTFIX}")  # paths to executables
  set(DIRS "")
  INSTALL(CODE " 
    include(BundleUtilities) 
    fixup_bundle(\"${APPS}.app\"   \"\"   \"${DIRS}\")
    execute_process(COMMAND mv \"${APPS}.app/Contents\" \"${APPS}\"  )
    execute_process(COMMAND rm -rf \"${APPS}.app\" \"\${CMAKE_INSTALL_PREFIX}/lib\" )
    " COMPONENT Runtime) 
endif()

# CPACK stuff
 # build a CPack driven installer package
 set (CPACK_RESOURCE_FILE_LICENSE "${CMAKE_SOURCE_DIR}/License.txt")
 set (CPACK_PACKAGE_NAME "${PROJECT_NAME}${PKG_POSTFIX}")
 # Version:
 set (CPACK_PACKAGE_VERSION_MAJOR ${MADX_MAJOR_VERSION})
 set (CPACK_PACKAGE_VERSION_MINOR ${MADX_MINOR_VERSION})
 set (CPACK_PACKAGE_VERSION_PATCH ${MADX_PATCH_LEVEL})

 set (CPACK_PACKAGE_CONTACT "Frank Schmidt <Frank.Schmidt@cern.ch>")   
 set (CPACK_PACKAGE_DESCRIPTION_SUMMARY "MadX is a program for accelerator design and simulation") 
 #Debian specific:
 set (CPACK_DEBIAN_PACKAGE_MAINTAINER
      "Yngve Inntjore Levinsen <Yngve.Inntjore.Levinsen@cern.ch>") # if this is not set, CPACK_PACKAGE_CONTACT is used instead..
 set (CPACK_DEBIAN_PACKAGE_DESCRIPTION "MadX is a program for accelerator design and simulation")
 set (CPACK_DEBIAN_PACKAGE_SECTION "science")
 set (CPACK_DEBIAN_PACKAGE_PRIORITY "extra")
 set (CPACK_DEBIAN_PACKAGE_RECOMMENDS "gnuplot")
 if(NOT MADX_STATIC)
 set (CPACK_DEBIAN_PACKAGE_DEPENDS 
      "libc6 (>= 2.3.1-6), libgcc1 (>= 1:4.1), zlib1g, libx11-xcb1, libxcb1, libxau6") # some version dependencies kept as examples
 endif()
 
 # RPM Specific:
 set (CPACK_RPM_PACKAGE_RELEASE 1)
 set (CPACK_RPM_PACKAGE_LICENSE "custom")
 set (CPACK_RPM_PACKAGE_GROUP "Development/Tools")
 if (NOT MADX_STATIC)
 # set(CPACK_RPM_PACKAGE_REQUIRES "libgcc >= 4.1.0, libxau >= 1.0.5") # I don't know the names of the packages...
 endif ( NOT MADX_STATIC )
 
 # so that we can build dragndrop on osx:
 set(CPACK_BINARY_DRAGNDROP ON)
 include (CPack)
# End CPACK stuff
