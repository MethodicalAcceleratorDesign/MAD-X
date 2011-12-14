# In case we are on a Unix system we install syntax file for Kate
# Only in case of devel version, so packages don't conflict...
if(NOT ${PROJECT_PATCH_LEVEL} EQUAL 00 AND UNIX)
  install(FILES kate/madx.xml
    DESTINATION "share/apps/katepart/syntax/"
    OPTIONAL)
  install(FILES vim/madx.vim
    DESTINATION "share/vim/vimfiles/syntax/"
    OPTIONAL)
endif()
