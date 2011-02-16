
if ( MADX_FORCE_32 OR ${CMAKE_SIZEOF_VOID_P} EQUAL 4 )

  message("32 bit build" ) 
  
  set (CPACK_DEBIAN_PACKAGE_ARCHITECTURE "i386")
  set (CPACK_RPM_PACKAGE_ARCHITECTURE "i686")
  
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
  
  if (CMAKE_SYSTEM_NAME STREQUAL "Linux")
    LINK_DIRECTORIES(${CMAKE_SOURCE_DIR}/lib)
  endif (CMAKE_SYSTEM_NAME STREQUAL "Linux")
  
elseif (${CMAKE_SIZEOF_VOID_P} EQUAL 8)

  message("64 bit build")
  
	 set (CPACK_DEBIAN_PACKAGE_ARCHITECTURE "amd64")
   set (CPACK_RPM_PACKAGE_ARCHITECTURE "x86_64")
   
  if (CMAKE_Fortran_COMPILER MATCHES "lf95")
    set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} --chk aefs")
  endif (CMAKE_Fortran_COMPILER MATCHES "lf95")
  
	 if (CMAKE_SYSTEM_NAME STREQUAL "Linux")
	   LINK_DIRECTORIES(${CMAKE_SOURCE_DIR}/lib64)
	 endif (CMAKE_SYSTEM_NAME STREQUAL "Linux")
 endif ( MADX_FORCE_32   OR ${CMAKE_SIZEOF_VOID_P} EQUAL 4 )
