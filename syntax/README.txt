For some text editors someone in the community have written a highlighting file, making it easier to read Mad-X scripts.



Vim: VIM is a popular editor on any platform. To use this file, copy madx.vim to ~/.vim/syntax/madx.vim and then add the following line to ~/.vimrc:
au BufNewFile,BufRead *.madx,*.seq,*.str,*.mask setf madx


Kate: Kate is a text editor for KDE, popular on Linux platforms. To use this file, copy madx.xml to one of the following folders (whichever is existing on your particular system):
   ~/.kde/share/apps/katepart/syntax
   ~/.kde3/share/apps/katepart/syntax
   ~/.kde4/share/apps/katepart/syntax

Emacs : Emacs is a text editor based on lisp. Check the file madx.el to load the madx highlighting mode (madx-mode).
   By default, any file with extension .madx will be automatically highlighted,
     but any buffer can be highlighted by doing in emacs "META-x", typing "madx-mode", and RETURN.
     NOTE: META is the emacs meta-character, it seems that in linux META is the same as ALT.

