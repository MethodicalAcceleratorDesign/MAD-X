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

# must use mcpp, http://mcpp.sourceforge.net/
CDEP := $(CC) -MM

#
# compiler
#

CFLAGS   := -m$(ARCH) -std=c99   -Wall -O3 -c
CXXFLAGS := -m$(ARCH) -std=c++0x -Wall -O3 -c

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
