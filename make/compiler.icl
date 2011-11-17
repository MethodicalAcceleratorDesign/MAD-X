# |
# o---------------------------------------------------------------------o
# |
# | MAD makefile - icl/icc compiler settings
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

CDEP := $(CC) /nologo /QMM

#
# compiler
#

CFLAGS   := /Qstd=c99   /O3 /c
CXXFLAGS := /Qstd=c++0x /O3 /c

#CFLAGS   := /Qstd=c99   /Wall /Wcheck /Wp64 /O3 /c
#CXXFLAGS := /Qstd=c++0x /Wall /Wcheck /Wp64 /O3 /c

#
# options flags
#

ifeq ($(DEBUG),yes)
CFLAGS   += /debug:all
CXXFLAGS += /debug:all
endif

ifeq ($(PROFILE),yes)
CFLAGS   += /Qprof-use
CXXFLAGS += /Qprof-use
endif

#
# extra flags
#
 
CFLAGS   += /nologo /Qprec /fp:precise /EHc /Qrestrict
CXXFLAGS += /nologo /Qprec /fp:precise /EHc /Qrestrict

# end of makefile
