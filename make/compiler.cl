# |
# o---------------------------------------------------------------------o
# |
# | MAD makefile - cl compiler settings (Visual Studio C++)
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

$(error Not Yet Supported)

#
# makedep
#

CDEP := $(CC) -MM

#
# compiler
#

CFLAGS   := -m$(ARCH) -O3
CXXFLAGS := -m$(ARCH) -std=c++0x -O3 -Wall

# trig too many warnings (> 100+) for the moment!
# CFLAGS   := -m$(ARCH) -std=c99 -Wall -O3

#
# options flags
#

ifeq ($(DEBUG),yes)
CFLAGS   += -g -gdwarf-2
CXXFLAGS += -g -gdwarf-2
endif

ifeq ($(PROFILE),yes)
CFLAGS   += -p
CXXFLAGS += -p
endif

#
# extra flags
#
 
CFLAGS   += -mp1 -fp-model precise
CXXFLAGS += -mp1 -fp-model precise

# end of makefile
