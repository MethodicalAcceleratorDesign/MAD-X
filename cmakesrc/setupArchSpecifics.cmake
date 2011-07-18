
# OSX specifics:
if(APPLE)
  set(CMAKE_LIBRARY_PATH /usr/lib/ /usr/X11/lib/ ${CMAKE_LIBRARY_PATH})
    set(CMAKE_OSX_DEPLOYMENT_TARGET 10.5)
endif(APPLE)

if ( MADX_FORCE_32 OR ${CMAKE_SIZEOF_VOID_P} EQUAL 4 )

  message("32 bit build" ) 
  
  set (CPACK_DEBIAN_PACKAGE_ARCHITECTURE "i386")
	
	if ( MADX_STATIC )
	 set (CPACK_RPM_PACKAGE_ARCHITECTURE "noarch")
	 # I think this should be correct but it is not tested yet...
	 #set (CPACK_DEBIAN_PACKAGE_ARCHITECTURE "all")
	else ( MADX_STATIC )
	 if ( NOT MADX_FORCE_32 )
  	set (CPACK_RPM_PACKAGE_ARCHITECTURE "i686")
	 else ( NOT MADX_FORCE_32 )
		 message(WARNING "Don't use CPACK to generate RPM, it will make no sense for these settings.")
	 endif ( NOT MADX_FORCE_32 )
	endif ( MADX_STATIC )
  
  set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -m32")
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -m32")
  if (CMAKE_Fortran_COMPILER MATCHES "lf95")
    set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Wa,--32 --chk aesux ")
  elseif (CMAKE_Fortran_COMPILER MATCHES "nagfor")
    set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}  -Wc,-m32 -Wl,-m32 -abi=32")
  elseif (CMAKE_Fortran_COMPILER MATCHES "ifort")
    set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fp-model precise -m32")
  elseif (CMAKE_Fortran_COMPILER MATCHES "gfortran")
    set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -m32")
  endif (CMAKE_Fortran_COMPILER MATCHES "lf95")
  set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -m32 ")
  
  if (CMAKE_SYSTEM_NAME STREQUAL "Linux" AND MADX_STATIC)
    LINK_DIRECTORIES(${CMAKE_SOURCE_DIR}/lib)
  endif (CMAKE_SYSTEM_NAME STREQUAL "Linux" AND MADX_STATIC)
  
elseif (${CMAKE_SIZEOF_VOID_P} EQUAL 8)

  message("64 bit build")
  
	 set (CPACK_DEBIAN_PACKAGE_ARCHITECTURE "amd64")
   set (CPACK_RPM_PACKAGE_ARCHITECTURE "x86_64")
   
  if (CMAKE_Fortran_COMPILER MATCHES "lf95")
    set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} --chk aefs")
  endif (CMAKE_Fortran_COMPILER MATCHES "lf95")
  
  if (CMAKE_SYSTEM_NAME STREQUAL "Linux" AND MADX_STATIC)
	   LINK_DIRECTORIES(${CMAKE_SOURCE_DIR}/lib64)
  endif (CMAKE_SYSTEM_NAME STREQUAL "Linux" AND MADX_STATIC)
 endif ( MADX_FORCE_32   OR ${CMAKE_SIZEOF_VOID_P} EQUAL 4 )


# OSX specifics:
if (APPLE)
  set(CMAKE_LIBRARY_PATH /usr/lib/ /usr/X11/lib/ ${CMAKE_LIBRARY_PATH})
        set(CMAKE_OSX_DEPLOYMENT_TARGET 10.5)
endif (APPLE)
