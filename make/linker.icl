# |
# o---------------------------------------------------------------------o
# |
# | MAD makefile - icl/icc linker settings
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

LDFLAGS = /nologo /O$(NOPT) /extlnk:.o
LDLIBS  =

#
# options flags
#

ifeq ($(DEBUG),yes)
LDFLAGS += /debug:all
endif

ifeq ($(PROFILE),yes)
LDFLAGS += /Qprof-use
endif

ifeq ($(STATIC),yes)
LDFLAGS += /static
endif

ifeq ($(SHARED),yes)
LDFLAGS += /shared
endif

ifeq ($(PLUGIN),yes)
LDFLAGS += /MD # TODO
endif

#
# command translator
#

LD_tr = $(strip $(subst $(SPACE)/O0 , /Od ,$(subst $(SPACE)-o , /Fe,$1)))

# end of makefile
