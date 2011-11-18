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

CL_LD1 := -o%
CL_LD2 := /Fe%

LD_tr = $(strip $(subst $(SPACE)/Fe , /Fe,$(call trans,$(CL_LD1),$(CL_LD2),$1)))

# end of makefile
