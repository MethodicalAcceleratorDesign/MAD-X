# This script builds pymad as an external project...

find_package(PythonInterp REQUIRED)
find_package(CYTHON REQUIRED)

set(PYMAD_BUILD_COMMAND ${PYTHON_EXECUTABLE} setup.py build)

ExternalProject_Add(
   pymad
   GIT_REPOSITORY git://github.com/pymad/pymad.git
#    DEPENDS madx
   CMAKE_ARGS
   CONFIGURE_COMMAND "pwd"
   BINARY_DIR pymad/src/
   SOURCE_DIR pymad
   BUILD_COMMAND ${PYMAD_BUILD_COMMAND}
   INSTALL_COMMAND ""
)

ExternalProject_Get_Property(pymad binary_dir)

set(PYMAD_INSTALL_COMMAND ${PYTHON_EXECUTABLE} setup.py install --prefix=${CMAKE_INSTALL_PREFIX})

install(CODE "execute_process(COMMAND ${PYMAD_INSTALL_COMMAND} WORKING_DIRECTORY ${binary_dir})" COMPONENT Runtime)