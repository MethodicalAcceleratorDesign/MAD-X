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

$(error $(CC) is not supported)

#
# makedep
#

# must use mcpp, http://mcpp.sourceforge.net/
# TODO: not supported by cl!!!
# CDEP := $(CC) /nolog /c /Zs /showIncludes

#
# compiler
#

CFLAGS   := /Za /Wall /O2 /c
CXXFLAGS := /Za /Wall /O2 /c

#
# options flags
#

ifeq ($(DEBUG),yes)
CFLAGS   += /Zi /Yd
CXXFLAGS += /Zi /Yd
endif

ifeq ($(PROFILE),yes)
CFLAGS   +=
CXXFLAGS +=
endif

#
# extra flags
#
 
CFLAGS   += /nologo /fp:precise /EHc
CXXFLAGS += /nologo /fp:precise /EHc

# end of makefile
