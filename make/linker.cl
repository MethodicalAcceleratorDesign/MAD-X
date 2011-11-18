# |
# o---------------------------------------------------------------------o
# |
# | MAD makefile - cl linker settings
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

#
# linker flags
#

LDFLAGS += /nologo /O3 /extlnk:.o

#
# options flags
#

ifeq ($(DEBUG),yes)
LDFLAGS += /Zi /Yd
endif

ifeq ($(PROFILE),yes)
LDFLAGS +=
endif

#
# command translator
#

LD_tr = $(strip $(subst $(SPACE)-o , /Fe,$1))

# end of makefile
