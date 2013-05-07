# This script builds pymad as an external project...

# It does not work at the moment because it requires that Mad-X is already installed...

find_package(PythonInterp REQUIRED)
find_package(CYTHON REQUIRED)

set(PYMAD_BUILD_COMMAND   ${PYTHON_EXECUTABLE} setup.py build --madxdir=${CMAKE_INSTALL_PREFIX})
set(PYMAD_INSTALL_COMMAND ${PYTHON_EXECUTABLE} setup.py install --prefix=${CMAKE_INSTALL_PREFIX})


ExternalProject_Add(
   pymad
   GIT_REPOSITORY git://github.com/pymad/pymad.git
   CMAKE_ARGS
   UPDATE_COMMAND git pull
   CONFIGURE_COMMAND ""
   BINARY_DIR pymad/src/
   SOURCE_DIR pymad
   BUILD_COMMAND ${PYMAD_BUILD_COMMAND}
   INSTALL_COMMAND  ""
)

ExternalProject_Get_Property(pymad binary_dir)

#install(CODE "execute_process(COMMAND ${PYMAD_BUILD_COMMAND}   WORKING_DIRECTORY ${binary_dir})" COMPONENT Runtime)
install(CODE "execute_process(COMMAND ${PYMAD_INSTALL_COMMAND} WORKING_DIRECTORY ${binary_dir})" COMPONENT Runtime)

