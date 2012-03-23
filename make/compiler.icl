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

ifneq ($(SED),)
CDEP   = $(CC)  /nologo /Zs /QMM $(addprefix /I,$(CC_DIR))
CXXDEP = $(CXX) /nologo /Zs /QMM $(addprefix /I,$(CXX_DIR))

# CDEP output translator
CDEP_tr = | $(SED) -e "s/$(call f2bs,$(CURDIR)/)//gi" \
                   -e "s/\.obj:/\.o:/g"
CXXDEP_tr = $(CDEP_tr)
endif

#
# compiler
#

CPPFLAGS += -D_ICC
CFLAGS   = /Qstd=c99   /Wall /Wcheck /Wp64 /O$(NOPT) /c
CXXFLAGS = /Qstd=c++0x /Wall /Wcheck /Wp64 /O$(NOPT) /c

#
# diagnostics
#

CFLAGS   += /Qdiag-disable:2259,1572,981 # /Qdiag-enable:sc2
CXXFLAGS += /Qdiag-disable:2259,1572,981 # /Qdiag-enable:sc2

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
CFLAGS   += /nologo /Qprec /fp:strict /EHc /Qrestrict $(addprefix /I,$(CC_DIR))
CXXFLAGS += /nologo /Qprec /fp:strict /EHc /Qrestrict $(addprefix /I,$(CXX_DIR))

#
# command translator
#

ICL_CC1 := -D%  -I% /O0
ICL_CC2 := /D%  /I% /Od

CC_tr  = $(strip $(subst $(SPACE)-o , /Fo,$(call trans,$(ICL_CC1),$(ICL_CC2),$1)))
CXX_tr = $(CC_tr)

# end of makefile
