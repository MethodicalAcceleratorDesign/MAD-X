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

#####################
# CL specific
#

CL_CC1 := -D%  -I% /O3 /O0
CL_CC2 := /D%  /I% /O2 /Od

###############
# C language
#

ifeq ($(CCNAME),cl)

#
# makedep
#

ifneq ($(and $(SED),$(GREP)),)
CDEP = $(CC) /nologo /c /Zs /showIncludes
CDEP_tr = | $(GREP) -i -F "$(call f1bs,$(CURDIR))" \
          | $(SED)  -e "s/$(call f2bs,$(CURDIR)/)//gi" \
                    -e "s/Note: including file:/$<:/gi" \
                    -e "s/\.c:/\.o:/g"
endif

#
# compiler
#

CFLAGS = /O$(NOPT) /c /D_MCC # /Wall

#
# options flags
#

ifeq ($(DEBUG),yes)
CFLAGS += /Zi /Yd
endif

ifeq ($(PROFILE),yes)
CFLAGS +=
endif

#
# extra flags
#

CFLAGS += /nologo /fp:precise /Zm1000 /EHsc /D_CRT_SECURE_NO_WARNINGS /Dinline=__inline

#
# command translator
#

CC_tr = $(strip $(subst $(SPACE)-o , /Fo,$(call trans,$(CL_CC1),$(CL_CC2),$1)))

endif

###############
# C++ language
#

ifeq ($(CXXNAME),cl)

#
# makedep
#

ifneq ($(and $(SED),$(GREP)),)
CXXDEP = $(CC) /nologo /c /Zs /showIncludes
CXXDEP_tr = | $(GREP) -i -F "$(call f1bs,$(CURDIR))" \
            | $(SED)  -e "s/$(call f2bs,$(CURDIR)/)//gi" \
                      -e "s/Note: including file:/$<:/gi" \
                      -e "s/\.c:/\.o:/g"
endif

#
# compiler
#

CXXFLAGS = /O$(NOPT) /c /D_MCC # /Wall

#
# options flags
#

ifeq ($(DEBUG),yes)
CXXFLAGS += /Zi /Yd
endif

ifeq ($(PROFILE),yes)
CXXFLAGS +=
endif

#
# extra flags
#

CXXFLAGS += /nologo /fp:precise /Zm1000 /EHsc /D_CRT_SECURE_NO_WARNINGS /Dinline=__inline

#
# command translator
#

CXX_tr = $(strip $(subst $(SPACE)-o , /Fo,$(call trans,$(CL_CC1),$(CL_CC2),$1)))

endif

# end of makefile
