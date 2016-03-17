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

#####################
# ICC specific
#

#
# preprocessor flags
#

CPPFLAGS += -D_ICC -D_ICL

#
# command translator
#

ICL_CC1 := -D%  -I% /O0
ICL_CC2 := /D%  /I% /Od

###############
# C language
#

ifeq ($(CCNAME),icl)

#
# makedep
#

ifneq ($(SED),)
CDEP = $(CC) /nologo /Zs /QMM
CDEP_tr = | $(SED) -e "s/$(call f2bs,$(CURDIR)/)//gi" -e "s+\\+/+g" -e "s/\.obj:/\.o:/g"
else
$(warning cannot compute files dependencies)
endif

#
# compiler
#

CFLAGS = /Qstd=c99 /Wall /Wcheck /Wp64 /O$(NOPT) /c

#
# diagnostics
#

CFLAGS += /D_CRT_SECURE_NO_WARNINGS /Qdiag-disable:3280,2259,1572,981 # /Qdiag-enable:sc2

#
# options flags
#

ifeq ($(DEBUG),yes)
CFLAGS += /Z7 /debug:full /traceback /RTC1 /RTCc /Qcheck-pointers:rw /Qcheck-pointers-dangling:all /check:conversions,stack,uninit
endif

ifeq ($(PROFILE),yes)
CFLAGS += /Qprof-use
endif

ifeq ($(OPENMP),yes)
CFLAGS += /Qopenmp
endif

#
# extra flags
#

CFLAGS += /nologo /Qprec /fp:strict /EHc /Qrestrict

#
# command translator
#

CC_tr = $(strip $(subst $(SPACE)-o , /Fo,$(call trans,$(ICL_CC1),$(ICL_CC2),$1)))

endif

###############
# C++ language
#

ifeq ($(CXXNAME),icl)

#
# makedep
#

ifneq ($(SED),)
CXXDEP = $(CXX) /nologo /Zs /QMM
CXXDEP_tr = | $(SED) -e "s/$(call f2bs,$(CURDIR)/)//gi" -e "s+\\+/+g" -e "s/\.obj:/\.o:/g"
else
$(warning cannot compute files dependencies)
endif

#
# compiler
#

CXXFLAGS = /Qstd=c++0x /Wall /Wcheck /Wp64 /O$(NOPT) /c

#
# diagnostics
#

CXXFLAGS += /D_CRT_SECURE_NO_WARNINGS /Qdiag-disable:3280,2259,1572,981 # /Qdiag-enable:sc2

#
# options flags
#

ifeq ($(DEBUG),yes)
CXXFLAGS += /Z7 /debug:full /traceback /RTC1 /RTCc /Qcheck-pointers:rw /Qcheck-pointers-dangling:all /check:conversions,stack,uninit
endif

ifeq ($(PROFILE),yes)
CXXFLAGS += /Qprof-use
endif

ifeq ($(OPENMP),yes)
CXXFLAGS += /Qopenmp
endif

ifeq ($(OPENMP),yes)
CXXFLAGS += /Qopenmp
endif

#
# extra flags
#

CXXFLAGS += /nologo /Qprec /fp:strict /EHc /Qrestrict

#
# command translator
#

CXX_tr = $(strip $(subst $(SPACE)-o , /Fo,$(call trans,$(ICL_CC1),$(ICL_CC2),$1)))

endif

# end of makefile
