# This script builds pymad as an external project...

find_package(PythonInterp REQUIRED)
find_package(CYTHON REQUIRED)

set(PYMAD_BUILD_COMMAND ${PYTHON_EXECUTABLE} setup.py build)
set(PYMAD_INSTALL_COMMAND ${PYTHON_EXECUTABLE} setup.py install)

get_target_property(LIBLOCATION madx LOCATION)

ExternalProject_Add(
   pymad
   GIT_REPOSITORY git://github.com/pymad/pymad.git
#    DEPENDS madx
   CMAKE_ARGS
   CONFIGURE_COMMAND "pwd"
   BINARY_DIR pymad/src/
   SOURCE_DIR pymad
   BUILD_COMMAND ${PYMAD_BUILD_COMMAND}
   INSTALL_COMMAND ${PYMAD_INSTALL_COMMAND}
)
