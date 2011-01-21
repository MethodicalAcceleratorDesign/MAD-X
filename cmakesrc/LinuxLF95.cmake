INCLUDE(CMakeForceCompiler)

# this one is important
SET(CMAKE_SYSTEM_NAME Linux)
# this is necessary..
set(CMAKE_SYSTEM_PROCESSOR i386)

#enable_language(Fortran)

# where is the target environment 
set(CMAKE_FIND_ROOT_PATH
        /opt/lib32
        /
       )

# stuff needed for CMake 2.6 (on lxplus)
if(CMAKE_VERSION MATCHES "2.6.")
set(CMAKE_C_COMPILER gcc)
set(CMAKE_CXX_COMPILER g++)
set(CMAKE_Fortran_COMPILER lf95)

set(CMAKE_C_COMPILER_ENV_VAR CC)
set(CMAKE_CXX_COMPILER_ENV_VAR CXX)
set(CMAKE_Fortran_COMPILER_ENV_VAR FC)
elseif(CMAKE_VERSION MATCHES "2.8.")
 CMAKE_FORCE_Fortran_COMPILER(lf95 Lahey)
endif(CMAKE_VERSION MATCHES "2.6.")

# specify the cross compiler
#set(CMAKE_Fortran_COMPILER_ENV_VAR FC)


# search for programs in the build host directories
SET(CMAKE_FIND_ROOT_PATH_MODE_PROGRAM NEVER)
# for libraries and headers in the target directories
SET(CMAKE_FIND_ROOT_PATH_MODE_LIBRARY ONLY)
SET(CMAKE_FIND_ROOT_PATH_MODE_INCLUDE ONLY)


SET(CMAKE_Fortran_FLAGS_INIT " -Wa,--32 ")
SET(CMAKE_SHARED_LIBRARY_LINK_Fortran_FLAGS "")
SET(CMAKE_SKIP_RPATH ON)
SET(CMAKE_C_FLAGS_INIT " -m32 ")
SET(CMAKE_CXX_FLAGS_INIT " -m32 ")
