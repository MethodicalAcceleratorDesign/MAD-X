if ( IS32BIT )
    message(STATUS "32 bit build" ) 
    
    set (CPACK_DEBIAN_PACKAGE_ARCHITECTURE "i386")
            
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
    set(CMAKE_OSX_DEPLOYMENT_TARGET 10.5)
endif (APPLE)

# Windows specifics:
if (WIN32)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -D_WIN32  -D_CATCH_MEM_W")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D_WIN32")
else()
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -D_CATCH_MEM")
endif()

