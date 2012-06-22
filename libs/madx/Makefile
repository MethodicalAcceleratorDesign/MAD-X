# |
# o---------------------------------------------------------------------o
# |
# | MADX makefile
# |
# o---------------------------------------------------------------------o
# |
# | Polymorphic Tracking Code from Etienne Forest
# |
# | For more information, see http://cern.ch/mad
# |
# o---------------------------------------------------------------------o
# |
# | $Id$
# |

# For makefile documentation, please read make/README
# For information and bug report, please contact mad@cern.ch

# debug
ifeq ($(DEBUG_MAKE),yes)
$(warning MAKEFLAGS = $(MAKEFLAGS))
$(warning MAKEOVERRIDES = $(MAKEOVERRIDES))
endif

###################
# Project settings

PROJECT := libmadx

SHARED  := yes

#################
# Build settings
#

AR := ar

#################
# Recursive call
#

.PHONY: any

export

any:
	@ $(MAKE) -C ../.. $(MAKEFLAGS) $(MAKEOVERRIDES) \
                     PROJECT=$(PROJECT) SHARED=$(SHARED) DEST_DIR="libs/madx/"

# end of makefile
