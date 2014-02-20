
macro(madx_install_targets)
    # Installs targets 
    # to the correct locations...
    install(TARGETS ${ARGN}
        BUNDLE DESTINATION . 
        COMPONENT Runtime
        RUNTIME DESTINATION bin 
        COMPONENT Runtime
        LIBRARY DESTINATION lib 
        COMPONENT Libraries
        ARCHIVE DESTINATION lib 
        COMPONENT Libraries
    )
endmacro()

macro(madx_install_headers)
    # This installs the header files to <prefix>/include/madX
        install(FILES ${ARGN} 
                DESTINATION "include/${PROJECT_NAME}"
                COMPONENT Files)
endmacro()

