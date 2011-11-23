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

#
# makedep
#

ifneq ($(and $(SED),$(GREP)),)
CDEP   := $(CC) /nologo /c /Zs /showIncludes
CXXDEP := $(CDEP)

# CDEP output translator
CDEP_tr = | $(GREP) -i -F "$(call f1bs,$(CURDIR))" \
          | $(SED)  -e "s/$(call f2bs,$(CURDIR)/)//gi" \
                    -e "s/Note: including file:/$<:/gi" \
                    -e "s/\.c:/\.o:/g"
CXXDEP_tr = $(CDEP_tr)
endif

#
# compiler
#

CPPFLAGS += -D_MCC
CFLAGS   := /O2 /c
CXXFLAGS := /O2 /c

# CFLAGS   := /Wall /O2 /c

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

CPPFLAGS += /D_CRT_SECURE_NO_WARNINGS /Dinline=__inline
CFLAGS   += /nologo /fp:precise /Zm1000 /EHsc
CXXFLAGS += /nologo /fp:precise /Zm1000 /EHsc

#
# command translator
#

CL_CC1 := -D%  -I%
CL_CC2 := /D%  /I%

CC_tr  = $(strip $(subst $(SPACE)-o , /Fo,$(call trans,$(CL_CC1),$(CL_CC2),$1)))
CXX_tr = $(CC_tr)

# end of makefile
