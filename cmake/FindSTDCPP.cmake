
if(NOT (WIN32 AND CMAKE_Fortran_COMPILER_ID STREQUAL "Intel") AND NOT APPLE)
   # On Debian Stable, stdc++ is not linked by default...
   set(STDCPP_LIBS stdc++ gcc_eh)
endif()
