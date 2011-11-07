# |
# o---------------------------------------------------------------------o
# |
# | MAD makefile - f95 compiler settings
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
# makedep
#

FDEP := $(FC) -M

#
# linker
#

LD      := $(FC)
LDFLAGS += -m$(ARCH)
LDLIBS  += -lstdc++

#
# basic flags
#

FFLAGS := -m$(ARCH) -std=f95 -Wall -pedantic -pipe -O3

#
# options flags
#

ifeq ($(DEBUG),yes)
FFLAGS += -g
endif

ifeq ($(PROFILE),yes)
FFLAGS  += -pg
LDFLAGS += -pg
endif

#
# extra flags
#

# none

# end of makefile
