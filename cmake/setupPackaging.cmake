

if(APPLE)
    set(MACOSX_BUNDLE_STARTUP_COMMAND madx${BINARY_POSTFIX})
    set(MACOSX_BUNDLE_ICON_FILE "${CMAKE_SOURCE_DIR}/cmake/MadX.icns")
    set(MACOSX_BUNDLE_LONG_VERSION_STRING "MadX  version ${madX_MAJOR_VERSION}.${madX_MINOR_VERSION}.${madX_PATCH_LEVEL}")
    set(MACOSX_BUNDLE_BUNDLE_NAME "MadX${BINARY_POSTFIX}")
    set(MACOSX_BUNDLE_GUI_IDENTIFIER "MadX${BINARY_POSTFIX}")
    # add icns to the .app/Resources with these two commands:
    set(srcfiles ${srcfiles} ${CMAKE_SOURCE_DIR}/cmake/MadX.icns)
    set_source_files_properties(${CMAKE_SOURCE_DIR}/cmake/MadX.icns PROPERTIES MACOSX_PACKAGE_LOCATION Resources)
endif()
