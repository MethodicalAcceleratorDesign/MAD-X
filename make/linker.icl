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

LDFLAGS = /nologo /O$(NOPT) /Qdiag-disable:10161
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
LDFLAGS += /MT
endif

ifeq ($(SHARED),yes)
LDFLAGS += /MD
endif

ifeq ($(PLUGIN),yes)
LDFLAGS += /MD # todo
endif

ifeq ($(OPENMP),yes)
LDFLAGS += /Qopenmp $(if $(call eq,$(STATIC),yes),/Qopenmp-link:static,)
endif

#
# command translator
#

LD_tr = $(strip $(subst $(SPACE)/O0 , /Od ,$(subst $(SPACE)-o , /Fe,$1)))

# end of makefile
