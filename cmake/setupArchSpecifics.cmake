if ( IS32BIT )
    message(STATUS "32 bit build" ) 
    
    set (CPACK_DEBIAN_PACKAGE_ARCHITECTURE "i386")
            
    if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
        if ( MADX_STATIC )
            set (CPACK_RPM_PACKAGE_ARCHITECTURE "noarch")
            # I think this should be correct but it is not tested yet...
            #set (CPACK_DEBIAN_PACKAGE_ARCHITECTURE "all")
        else ()
            if ( NOT MADX_FORCE_32 )
                set (CPACK_RPM_PACKAGE_ARCHITECTURE "i686")
            else ()
                    message(WARNING "Don't use CPACK to generate RPM, it will make no sense for these settings.")
            endif ()
        endif ()
    endif()
    
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -m32")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -m32")

    if (CMAKE_SYSTEM_NAME STREQUAL "Linux" AND MADX_STATIC)
        LINK_DIRECTORIES(${CMAKE_SOURCE_DIR}/lib32)
    endif ()
  
else()

    message(STATUS "64 bit build")
    
    set (CPACK_DEBIAN_PACKAGE_ARCHITECTURE "amd64")
    set (CPACK_RPM_PACKAGE_ARCHITECTURE "x86_64")
    
    if (CMAKE_SYSTEM_NAME STREQUAL "Linux" AND MADX_STATIC)
            LINK_DIRECTORIES(${CMAKE_SOURCE_DIR}/lib64)
    endif ()
endif ()


# OSX specifics:
if (APPLE)
    set(CMAKE_LIBRARY_PATH /usr/lib/ /usr/X11/lib/ ${CMAKE_LIBRARY_PATH})
    # This does not work for 10.7 or above, so commenting out for now
    # I think it is just some cpack related option..
    # set(CMAKE_OSX_DEPLOYMENT_TARGET 10.5)
endif (APPLE)

if (WIN32)
   if(NOT CMAKE_C_COMPILER_ID MATCHES "GNU")
      add_definitions("-D_WIN32")
   endif()
elseif(APPLE)
   add_definitions("-D_DARWIN")
elseif(CMAKE_SYSTEM_NAME STREQUAL "Linux")
   add_definitions("-D_LINUX")
endif()
