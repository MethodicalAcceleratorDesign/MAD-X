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

-include $(makedir)/compiler.$(CCNAME)
-include $(FILE_C)

-include $(makedir)/compiler.$(CXXNAME)
-include $(FILE_CXX)

-include $(makedir)/compiler.$(FCNAME)
-include $(FILE_F90)

-include $(FILE_LIB)
-include $(FILE_SYS)

# linker specific settings
-include $(makedir)/linker.$(LDNAME)

# system specific settings
-include $(makedir)/system.$(OSTYPE)

# compilers, linker and depend rules
include $(makedir)/compiler.rules
include $(makedir)/linker.rules
include $(makedir)/depend.rules

# cleaning & debugging
include $(makedir)/clean.rules
include $(makedir)/infos.rules

# user's extra stuff
-include $(FILE_USR)

# end of makefile
