
find_package(PythonInterp)
if(${CMAKE_VERSION} VERSION_GREATER 2.8 AND PYTHONINTERP_FOUND)
    option(INCLUDE_PYMAD "Build the external pymad project from pymad.github.com" OFF)
    if(INCLUDE_PYMAD)
        include(ExternalProject)
        include(ExternalProject-pymad)
        # For pymad to work, we need to force building shared libraries:
        set(BUILD_SHARED_LIBS ON CACHE BOOL "Build shared libraries (.so/.dll)" FORCE)
    endif()
endif()