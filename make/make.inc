# |
# o---------------------------------------------------------------------o
# |
# | MAD makefile - includes
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

# make utilities
include $(makedir)/make.lib

# make config
include $(makedir)/make.cfg

# project files
-include $(FILE_VER)
-include $(FILE_CPP)

$(call exists, $(makedir)/compiler.$(CCNAME), \
unable to include settings for compiler CC=$(CCNAME) (make [options] infoinc); \
check CC or CCNAME for correct compiler configuration file)
include $(makedir)/compiler.$(CCNAME)
-include $(FILE_C)

$(call exists, $(makedir)/compiler.$(CXXNAME), \
unable to include settings for compiler CXX=$(CXXNAME) (make [options] infoinc); \
check CXX or CXXNAME for correct compiler configuration file)
ifneq ($(CCNAME),$(CXXNAME))
include $(makedir)/compiler.$(CXXNAME)
endif
-include $(FILE_CXX)

$(call exists, $(makedir)/compiler.$(FCNAME), \
unable to include settings for compiler FC=$(FCNAME) (make [options] infoinc); \
check FC or FCNAME for correct compiler configuration file)
ifneq ($(CCNAME),$(FCNAME))
include $(makedir)/compiler.$(FCNAME)
endif
-include $(FILE_F90)

-include $(FILE_LIB)
-include $(FILE_SYS)

# linker specific settings
$(call exists, $(makedir)/linker.$(LDNAME), \
unable to include settings for linker LD=$(LDNAME) (make [options] infoinc); \
check LD or LDNAME for correct linker configuration file)
include $(makedir)/linker.$(LDNAME)

# system specific settings
-include $(makedir)/system.$(OSTYPE)

# user's extra stuff
-include $(FILE_USR)

# compilers, linker and depend rules
include $(makedir)/compiler.rules
include $(makedir)/linker.rules
include $(makedir)/depend.rules

# cleaning & debugging
include $(makedir)/clean.rules
include $(makedir)/infos.rules

# end of makefile
