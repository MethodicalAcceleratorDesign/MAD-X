# |
# o---------------------------------------------------------------------o
# |
# | MAD makefile - ld linker settings
# |
# o---------------------------------------------------------------------o
# |
# | Methodical Accelerator Design
# |
# | Copyright (c) 2011+ CERN, mad@cern.ch
# |
# | For more information, see http://cern.ch/mad
# |
# o---------------------------------------------------------------------o
# |
# | $Id$
# |

ifeq ($(LD),)
$(error Wrong linker settings (not configured?), 'ld' as default linker is not supported)
endif

# end of makefile
