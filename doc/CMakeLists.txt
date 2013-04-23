
# Installs usrguide to <prefix>/share/madX/usrguide

# creates a command madx-doc which opens documentation in
# default system browser

option(MADX_INSTALL_DOC "Install documentation on supported systems" ON)

if(UNIX AND MADX_INSTALL_DOC)
   install(DIRECTORY usrguide DESTINATION share/madX)

   # Write the file in binary directory
   set(_script "${CMAKE_CURRENT_BINARY_DIR}/madx-doc")
   # Get correct open command
   if(APPLE)
      set(_doc_open_cmd "open")
   else()
      set(_doc_open_cmd "xdg-open")
   endif()

   # Write the content to the file..
   file(WRITE ${_script}  "#!/bin/sh\n")
   file(APPEND ${_script} "if [ \"\$\#\" == \"1\" ];then\n")
   file(APPEND ${_script} " ${_doc_open_cmd} ${CMAKE_INSTALL_PREFIX}/share/madX/usrguide/$1/$1.html\n")
   file(APPEND ${_script} "else\n")
   file(APPEND ${_script} " ${_doc_open_cmd} ${CMAKE_INSTALL_PREFIX}/share/madX/usrguide/uguide.html\n")
   file(APPEND ${_script} "fi\n")

   install(PROGRAMS ${_script} DESTINATION bin)
endif()
