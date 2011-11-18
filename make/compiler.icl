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

#CDEP   := $(CC) /nologo /QMM
#CXXDEP := $(CC) /nologo /QMM

#
# compiler
#

CPPFLAGS += -D_INTEL_CL
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

CPPFLAGS += /D_CRT_SECURE_NO_WARNINGS 
CFLAGS   += /nologo /Qprec /fp:source /EHc /Qrestrict
CXXFLAGS += /nologo /Qprec /fp:source /EHc /Qrestrict

#
# command translator
#

ICL_CC1 := -D%  -I%  -o%
ICL_CC2 := /D%  /I%  /Fo%

CC_tr  = $(strip $(subst $(SPACE)/Fo , /Fo,$(call trans,$(ICL_CC1),$(ICL_CC2),$1)))
CXX_tr = $(CC_tr)

# end of makefile
