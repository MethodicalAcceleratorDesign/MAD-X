# In case we are on a Unix system we install syntax file for Kate
# Only in case of devel version, so packages don't conflict...
if(UNIX)
  if(NOT APPLE)
    install(FILES kate/madx.xml
      DESTINATION "share/apps/katepart/syntax/"
      OPTIONAL)
  endif()
  install(FILES emacs/madx.el
    DESTINATION "share/emacs/site-lisp/mad-x/"
    OPTIONAL)
  install(FILES vim/madx.vim
    DESTINATION "share/vim/vimfiles/syntax/"
    OPTIONAL)
endif()
