# |
# o---------------------------------------------------------------------o
# |
# | MADX makefile
# |
# o---------------------------------------------------------------------o
# |
# | MADX as a library
# |
# | For more information, see http://cern.ch/mad
# |
# o---------------------------------------------------------------------o
# |
# | $Id$
# |

# Note: this makefile defines a "virtual project"
#   - is not supposed to be used directly
#   - assume that the madx environment is already setup
#   - add settings to the original madx environment (PROJECT, SHARED, DESTDIR)
#   - call/reload the madx Makefile again

# For makefile documentation, please read make/README
# For information and bug report, please contact mad@cern.ch

# debug
ifeq ($(DEBUG_MAKE),yes)
$(info MAKEFLAGS = $(MAKEFLAGS))
$(info MAKEOVERRIDES = $(MAKEOVERRIDES))
endif

###################
# Project settings

PROJECT := libmadx

SHARED  := yes

#################
# Build settings
#

# add the archiver
AR := ar

#################
# Recursive call
#

.PHONY: $(PROJECT)

export

$(PROJECT):
	@ $(MAKE) -C ../.. $(MAKEFLAGS) $(MAKEOVERRIDES) \
                     PROJECT=$(PROJECT) SHARED=$(SHARED) DESTDIR=libs/madx/$(DESTDIR)

# end of makefile
