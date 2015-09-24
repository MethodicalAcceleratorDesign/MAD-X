subroutine call_gino(GINO_COMMAND)
  implicit none
  CHARACTER(*) GINO_COMMAND

  Write(6,*) "  "
  Write(6,*) GINO_COMMAND
  Write(6,*) " This is not available for this version  of MAD-X"
  Write(6,*) " The code for the Windows Operating System is available at"
  Write(6,*) " http://mad.web.cern.ch/mad/PTC_proper/mad-xp-gino/"
  Write(6,*) "  All questions and complaints should be directed to "
  Write(6,*) "  Etienne Forest at jceepf@hushmail.com  "
  Write(6,*) "  "
End subroutine call_gino

subroutine open_gino_graphics
  implicit none
End subroutine open_gino_graphics

subroutine close_gino_graphics
  implicit none
End subroutine close_gino_graphics